
import time
import logging
import numpy
import queue

import scipy.optimize

from sklearn import cluster, datasets, mixture
from sklearn.neighbors import kneighbors_graph
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt




def select_samples(times, reads, noise, n_samples, speedy=False):

    if (not speedy):
        # time differences between reads
        dt = times.reshape((-1, 1)) - times.reshape((1, -1))
        df = reads.reshape((-1, 1)) - reads.reshape((1, -1))
        d_noise = noise.reshape((-1, 1)) + noise.reshape((1, -1))
    else:
        # let's cut down the number of samples to use for pairwise fitting

        # use all samples
        n_select_x = numpy.arange(numpy.min([30, n_samples]))
        n_select_y = numpy.arange(numpy.min([30, n_samples]))
        pass
        if (n_samples > 30):
            # use first 50, plus every 3rd after wards
            n_select_x = numpy.append(n_select_x, numpy.arange(30, numpy.min([150, n_samples, 4])))
            n_select_x = numpy.append(n_select_x, numpy.arange(30, numpy.min([151, n_samples, 3])))
            pass
        if (n_samples > 150):
            #  use first 50, plus every 5th to 150, then fewer based on prime numbers to make sure we get
            #  unique pairings
            n_select_x = numpy.append(n_select_x, numpy.arange(150, n_samples, 7))
            n_select_x = numpy.append(n_select_x, numpy.arange(152, n_samples, 11))
        # always include the last sample
        n_select_x = numpy.append(n_select_x, [-1])
        n_select_y = numpy.append(n_select_y, [-1])

        dt = times[n_select_x].reshape((-1, 1)) - times[n_select_y].reshape((1, -1))
        df = reads[n_select_x].reshape((-1, 1)) - reads[n_select_y].reshape((1, -1))
        d_noise = noise[n_select_x].reshape((-1, 1)) + noise[n_select_y].reshape((1, -1))

    return dt,df,d_noise




def worker__fit_rauscher2007(
        shmem_cube_corrected, shmem_results, cube_shape, results_shape,
        jobqueue, read_times, workername=None, group_cutoff=None, recombine=False, **kwargs,
):
    logger = logging.getLogger(workername if workername is not None else "Rauscher2007worker")
    logger.debug("Starting worker")

    cube_linearized = numpy.ndarray(shape=cube_shape, dtype=numpy.float32,
                             buffer=shmem_cube_corrected.buf)
    cube_results = numpy.ndarray(shape=results_shape, dtype=numpy.float32,
                             buffer=shmem_results.buf)

    while (True):
        try:
            job = jobqueue.get()
            if (job is None):
                # this is the termination signal
                jobqueue.task_done()
                break
        except queue.Empty as e:
            logger.warning("job queue empty (%s)" % (e))
            break

        y = job
        logger.debug("Starting to fit pairwise slopes for row %d" % (y))

        t1 = time.time()
        # get correction from data

        for x in range(4, cube_results.shape[2]-4 ):

            reads = cube_linearized[:,y,x]
            noise = numpy.sqrt(reads) # TODO: THIS NEEDS FIXING
            times = numpy.array(read_times)

            if (recombine):
                reads, recombine_corrections = recombine_signals(times, reads)

            good_reads = numpy.isfinite(reads) & numpy.isfinite(noise) & numpy.isfinite(times) & (times >= 0)
            if (group_cutoff is not None and group_cutoff>0):
                good_reads[group_cutoff:] = False

            if (numpy.sum(good_reads) < 2):
                # not enough data to do anything with
                continue

            times = times[good_reads]
            reads = reads[good_reads]
            noise = noise[good_reads]
            n = numpy.sum(good_reads)

            rauscher_b = (n * numpy.sum(times*reads) - numpy.sum(times)*numpy.sum(reads)) / (n*numpy.sum(times**2) - numpy.sum(times)**2)
            rauscher_a = (numpy.sum(times**2)*numpy.sum(reads) - numpy.sum(times)*numpy.sum(times*reads)) / (n*numpy.sum(times**2) - numpy.sum(times)**2)

            max_t = numpy.max(times)
            n_samples = times.shape[0]

            # get the pair-wise sampling of all parameters (times, reads, and noise)
            weighted = rauscher_b
            _med = 0
            _sigma = 0
            n_useful_pairs = n

            cube_results[:,y,x] = [weighted, _med, _sigma, n_useful_pairs, max_t]

        t2 = time.time()
        logger.debug("Fitting rauscher2007 for row %d done after %.3f seconds" % (y, t2-t1))

        jobqueue.task_done()

    logger.debug("Shutting down")
    shmem_results.close()
    shmem_cube_corrected.close()



def __fit_linear_regression(p, times):
    return p[0] + p[1]*times
def __fit_linear_regression_error(p, times, reads, noise):
    fit = __fit_linear_regression(p, times)
    sigma = (reads - fit) / noise
    return sigma**2

def worker__fit_linear_regression(
        shmem_cube_corrected, shmem_results, cube_shape, results_shape,
        jobqueue, read_times, workername=None, group_cutoff=None,
        recombine=False, **kwargs,
):
    logger = logging.getLogger(workername if workername is not None else "LinearRegressionWorker")
    logger.debug("Starting worker")

    cube_linearized = numpy.ndarray(shape=cube_shape, dtype=numpy.float32,
                             buffer=shmem_cube_corrected.buf)
    cube_results = numpy.ndarray(shape=results_shape, dtype=numpy.float32,
                             buffer=shmem_results.buf)

    while (True):
        try:
            job = jobqueue.get()
            if (job is None):
                # this is the termination signal
                jobqueue.task_done()
                break
        except queue.Empty as e:
            logger.warning("job queue empty (%s)" % (e))
            break

        y = job
        logger.debug("Starting to fit linear regressions for row %d" % (y))

        t1 = time.time()
        # get correction from data

        for x in range(4, cube_results.shape[2]-4 ):

            reads = cube_linearized[:,y,x]
            noise = numpy.sqrt(reads) # TODO: THIS NEEDS FIXING
            times = numpy.array(read_times)

            good_reads = numpy.isfinite(reads) & numpy.isfinite(noise) & numpy.isfinite(times) & (times >= 0) & (noise>0)
            if (group_cutoff is not None and group_cutoff>0):
                good_reads[group_cutoff:] = False

            n = numpy.sum(good_reads)
            if (n < 2):
                # not enough data to do anything with
                continue

            if (n > 10):
                # we have enough data, skip the first N reads
                N = 3
                good_reads[:N] = False
                n -= N

            times = times[good_reads]
            reads = reads[good_reads]
            noise = noise[good_reads]

            if (recombine):
                try:
                    new_reads, recombine_corrections = recombine_signals(times, reads)
                    reads = new_reads
                except Exception as e:
                    numpy.savetxt("recombine_dump_%d_%d.dat" % (x,y), numpy.array([times, reads]).T)
                    logger.warning("Recombine failed for x/y = %d/%d: %s" % (x,y, str(e)))

            # slope = (n * numpy.sum(times*reads) - numpy.sum(times)*numpy.sum(reads)) / (n*numpy.sum(times**2) - numpy.sum(times)**2)
            # intercept = (numpy.sum(times**2)*numpy.sum(reads) - numpy.sum(times)*numpy.sum(times*reads)) / (n*numpy.sum(times**2) - numpy.sum(times)**2)
            slope, intercept = fit_line(times, reads)
            p0 = [intercept, slope]

            best_fit = [numpy.NaN, numpy.NaN]
            perr = [numpy.NaN, numpy.NaN]
            if (n >= 2):
                # Only attempt to fit a slope if we have at least 2 points
                try:
                    results = scipy.optimize.leastsq(
                        func=__fit_linear_regression_error,
                        x0=p0,
                        args=(times, reads, noise),
                        full_output=True,
                    )
                    best_fit, pcov, infodict, mesg, ier = results
                    perr = numpy.sqrt(numpy.diag(pcov))
                except:
                    pass

            max_t = numpy.max(times)
            n_samples = times.shape[0]

            # get the pair-wise sampling of all parameters (times, reads, and noise)
            weighted = best_fit[1]
            _med = slope
            _sigma = perr[1]
            n_useful_pairs = n #perr[0]

            cube_results[:,y,x] = [weighted, _med, _sigma, n_useful_pairs, max_t]

        t2 = time.time()
        logger.debug("Fitting rauscher2007 for row %d done after %.3f seconds" % (y, t2-t1))

        jobqueue.task_done()

    logger.debug("Shutting down")
    shmem_results.close()
    shmem_cube_corrected.close()




def fit_pairwise_slope_samples(dt, df, d_noise):

    useful_pairs = (dt > 0) & numpy.isfinite(df)
    n_useful_pairs = numpy.sum(useful_pairs)
    rates = (df / dt)[useful_pairs]
    noises = d_noise[useful_pairs]
    #     print(rates.shape)

    # if (y >= 1890 and y <= 1907 and x >= 1405 and x <= 1414):
    #     print("Dump %d %d" % (x, y))
    #     numpy.savetxt("dump_%d_%d.p1" % (x, y),
    #                   numpy.array([times, reads, noise]).T)
    #     numpy.savetxt("dump_%d_%d.p2" % (x, y),
    #                   numpy.array([rates, noises]).T)

    try:
        good = numpy.isfinite(rates) & numpy.isfinite(noises)
        n_good = numpy.sum(good)

        for it in range(3):
            _stats = numpy.nanpercentile(rates[good], [16, 50, 84])
            _med = _stats[1]
            _sigma = 0.5 * (_stats[2] - _stats[0])
            new_good = good & (rates > _med - 3 * _sigma) & (rates < _med + 3 * _sigma)
            n_good = numpy.sum(new_good)
            if (n_good < 10):
                # if after this iteration we are left we too few good results then
                # skip this iteration step and use all remaining data
                break
            else:
                good = new_good

        weights = 1. / noises
        weighted = numpy.sum((rates * weights)[good]) / numpy.sum(weights[good])

        # cube_results[:, y, x] = [weighted, _med, _sigma, n_useful_pairs, max_t]
    except Exception as e:
        # logger.debug("Encountered exception in pairfitting for x=%d, y=%d: %s" % (x, y, str(e)))
        weighted, _med, _sigma, n_useful_pairs = numpy.NaN, numpy.NaN, numpy.NaN, numpy.NaN

    return weighted, _med, _sigma, n_useful_pairs




def worker__fit_pairwise_slopes(
        shmem_cube_corrected, shmem_results, cube_shape, results_shape,
        jobqueue, read_times, speedy=False, workername=None, **kwargs,
):

    logger = logging.getLogger(workername if workername is not None else "PairSlopesWorker")
    logger.debug("Starting worker")

    cube_linearized = numpy.ndarray(shape=cube_shape, dtype=numpy.float32,
                             buffer=shmem_cube_corrected.buf)
    cube_results = numpy.ndarray(shape=results_shape, dtype=numpy.float32,
                             buffer=shmem_results.buf)

    while (True):
        try:
            job = jobqueue.get()
            if (job is None):
                # this is the termination signal
                jobqueue.task_done()
                break
        except queue.Empty as e:
            logger.warning("job queue empty (%s)" % (e))
            break

        y = job
        logger.debug("Starting to fit pairwise slopes for row %d" % (y))

        t1 = time.time()
        # get correction from data

        for x in range(4, cube_results.shape[2]-4 ):

            reads = cube_linearized[:,y,x]
            noise = numpy.sqrt(reads) # TODO: THIS NEEDS FIXING
            times = numpy.array(read_times)

            good_reads = numpy.isfinite(reads) & numpy.isfinite(noise) & numpy.isfinite(times) & (times >= 0)
            if (numpy.sum(good_reads) < 1):
                # not enough data to do anything with
                continue

            times = times[good_reads]
            reads = reads[good_reads]
            noise = noise[good_reads]

            max_t = numpy.max(times)
            n_samples = times.shape[0]

            # get the pair-wise sampling of all parameters (times, reads, and noise)
            dt,df,d_noise = select_samples(times, reads, noise, n_samples, speedy)

            #
            weighted, _med, _sigma, n_useful_pairs = fit_pairwise_slope_samples(dt,df,d_noise)
            cube_results[:,y,x] = [weighted, _med, _sigma, n_useful_pairs, max_t]

        t2 = time.time()
        logger.debug("Fitting pair-slopes for row %d done after %.3f seconds" % (y, t2-t1))

        jobqueue.task_done()

    logger.debug("Shutting down")
    shmem_results.close()
    shmem_cube_corrected.close()



def fit_line(t,r):
    """
    Fits a straight line to a set of x & y (or t and r) datapoints.

    :param t: x coordinate (read times)
    :param r: y coordinate (read fluxes)
    :return: slope, intercept

    """
    valid = numpy.isfinite(t) & numpy.isfinite(r)
    n = numpy.sum(valid)
    _t = t[valid]
    _r = r[valid]
    slope = (n * numpy.sum(_t*_r) - numpy.sum(_t)*numpy.sum(_r)) / (n*numpy.sum(_t**2) - numpy.sum(_t)**2)
    intercept = (numpy.sum(_t**2)*numpy.sum(_r) - numpy.sum(_t)*numpy.sum(_t*_r)) / (n*numpy.sum(_t**2) - numpy.sum(_t)**2)
    return slope,intercept

def noprint(*args):
    pass

def recombine_signals(times, reads, max_iterations=5, debug=False):
    """
    Helper function to recombine URG sequences with jumps due to cosmics or alternating "bias"-levels

    :param times: read time (seconds or read number)
    :param reads: reference pixel corrected intensity
    :return:
    """

    # reads = numpy.array(reads)
    # reads[700:] += 100

    myprint = noprint
    mysave = noprint
    if (debug):
        myprint = print
        mysave = numpy.savetxt

    # ensure only valid data
    bad_data = ~numpy.isfinite(times) | ~numpy.isfinite(reads)
    if (bad_data.any()):
        raise ValueError("Invalid data encountered during signal recombination")

    n_clusters = 10

    corrected_reads = numpy.array(reads)
    corrections = numpy.zeros_like(times, dtype=float)
    iteration = 0

    db_eps = 0.5
    # max_iterations = 3

    while (n_clusters > 1 and iteration < max_iterations):
        iteration += 1
        myprint("\n\n\n", "ITERATION", iteration)

        # fig, axs = plt.subplots(ncols=2, figsize=(12, 3))
        # ax = axs[0]
        # fig.suptitle("iteration %d" % (iteration))
        # ax.set_title("raw/working data")
        # ax.scatter(times, corrected_reads, s=1, c='black', label='raw')

        # prepare data for analysis
        raw_X = numpy.array([times, corrected_reads]).T
        scaled_X = StandardScaler().fit_transform(raw_X)
        mysave("scaled_X", scaled_X)
        myprint("scaled data", scaled_X.shape)

        # perform a clustering analysis on the raw sequence
        db = cluster.DBSCAN(min_samples=10, eps=db_eps)
        db.fit(scaled_X)
        myprint("DBSCAN cluster", db)

        labels = db.labels_
        myprint(labels)
        #  time.sleep(0.5)

        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        if (n_clusters_ < 1):
            # nothing to do, either all noise or not enough data
            break

        n_noise_ = list(labels).count(-1)
        myprint("Estimated number of clusters: %d" % n_clusters_)
        myprint("Estimated number of noise points: %d" % n_noise_)

        unique_labels = numpy.array(list(set(labels)))
        unique_labels = unique_labels[unique_labels >= 0]
        myprint("unique labels:", unique_labels)
        core_samples_mask = numpy.zeros_like(labels, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        core_samples_mask[:] = True
        label_count = numpy.array([numpy.sum(labels == label) for label in unique_labels])
        myprint("label counts", label_count)
        myprint("SAVING DATA DUMP", iteration)
        mysave("iteration_%d.all" % iteration, numpy.array([
            times, corrected_reads, scaled_X[:,0], scaled_X[:,1], labels, core_samples_mask,
        ]).T)

        # sort by number of reads in each sequence(
        sort_by_size = numpy.argsort(label_count)[::-1]
        myprint("label count, sorted by size", label_count[sort_by_size])
        labels_sorted_by_size = unique_labels[sort_by_size]
        sorted_size = label_count[sort_by_size]
        myprint("SORTED (by size) LABELS", labels_sorted_by_size)

        min_times = numpy.zeros_like((unique_labels))
        max_times = numpy.zeros_like((unique_labels))
        mean_times = numpy.zeros_like((unique_labels))
        for i_label, label in enumerate(labels_sorted_by_size):
            myprint("CHECKING idx=", i_label, "   // label=", label)
            this_sequence = (labels == label) #& core_samples_mask
            this_times = times[this_sequence]
            this_reads = corrected_reads[this_sequence]
            myprint("label=", label, "   nreads=", this_times.shape[0])
            min_times[i_label] = numpy.min(this_times)
            max_times[i_label] = numpy.max(this_times)
            mean_times[i_label] = numpy.nanmean(this_times)
            mean_read = numpy.nanmean(this_reads)
            # ax.scatter(mean_times[i_label], mean_read, c='white', s=300, alpha=0.5, zorder=90)
            # ax.annotate("%s" % (label), xy=(mean_times[i_label], mean_read),
            #             xytext=(mean_times[i_label], mean_read),
            #             xycoords='data', zorder=91, fontsize='large', fontweight='bold',
            #             verticalalignment='center', horizontalalignment='center')

        myprint("label // min-time // max-time // mean-time\n",
                numpy.array([unique_labels, min_times, max_times, mean_times]).T)
        myprint("MAIN SEQUENCE label", labels_sorted_by_size[0])

        main_label = labels_sorted_by_size[0]
        main_sequence = (labels == main_label) & core_samples_mask
        myprint(main_sequence)
        main_times = times[main_sequence]
        main_reads = corrected_reads[main_sequence]
        main_slope, main_intercept = fit_line(main_times, main_reads)
        myprint(main_slope, main_intercept)

        mysave("it_%d_main.txt" % iteration, numpy.array([main_times, main_reads]).T)
        myprint("main_slope/intercept:", main_slope, main_intercept)

        # ax.scatter(main_times, main_reads, s=40, alpha=0.1, color='red', label='main seq')
        # ax.plot(main_times, main_times * main_slope + main_intercept, color='red')

        #         corrected_data |= this_sequence

        if (unique_labels.shape[0] == 1):
            # we only have one sequence, so this is all we need to do
            break

        # check if we have any overlapping sequences
        overlapping = (mean_times > min_times[0]) & (mean_times < max_times[0])
        myprint("OVERLAPPING:", overlapping)
        if (numpy.sum(overlapping[1:]) > 0):
            # we have at least one overlapping sequence
            merge_label = labels_sorted_by_size[overlapping][1]
            myprint("merging label:", merge_label)
            to_merge = (labels == merge_label)
            to_merge_core = to_merge & core_samples_mask
            merge_times = times[to_merge_core]
            merge_reads = corrected_reads[to_merge_core]

            # estimate the mean offset
            merge_offsets = (merge_times * main_slope + main_intercept) - merge_reads
            myprint(numpy.mean(merge_offsets))
            myprint(merge_offsets)
            mean_merge_offset = numpy.mean(merge_offsets)
            corrections[to_merge] += mean_merge_offset
            corrected_reads[to_merge] += mean_merge_offset

            myprint("we have merged a sequence, restarting")
            max_iterations += 1
            continue

        myprint("No overlapping sequences found, continuing")
        #
        # now more overlapping sequences, so let's look for the closest one
        # note that index 0 is the main sequence itself
        #
        sequence_closeness = numpy.argsort(numpy.fabs(mean_times - mean_times[0]))
        myprint("    (delta)-time)", mean_times - mean_times[0])
        myprint("fabs(delta)-time)", numpy.fabs(mean_times - mean_times[0]))
        myprint("labels:", labels_sorted_by_size[sequence_closeness])
        myprint("closeness:", sequence_closeness)
        closest_sequence = sequence_closeness[1]
        myprint("closest:", closest_sequence, labels_sorted_by_size[closest_sequence])
        closest_label = labels_sorted_by_size[closest_sequence]

        # check out the other sequences
        next_sequence = (labels == closest_label)  #
        next_sequence_core = next_sequence & core_samples_mask
        next_times = times[next_sequence_core]
        next_reads = corrected_reads[next_sequence_core]
        # ax.scatter(next_times, next_reads, s=15, alpha=0.2, label='next')
        mysave("it_%d_next.txt" % iteration, numpy.array([next_times, next_reads]).T)
        next_slope, next_intercept = fit_line(next_times, next_reads)
        myprint("next_slope/intercept:", next_slope, next_intercept)
        # ax.plot(next_times, next_times * next_slope + next_intercept, label='next fit')
        # ax.plot(main_times, main_times * main_slope + main_intercept, label='main fit')

        # ax.legend()
        # now find connection point between main and next sequence
        if (numpy.min(next_times) > numpy.max(main_times)):
            # this sequence comes later
            midpoint = 0.5 * (numpy.min(next_times) + numpy.max(main_times))
        elif (numpy.max(next_times) < numpy.min(main_times)):
            # this sequence comes earlier
            midpoint = 0.5 * (numpy.max(next_times) + numpy.min(main_times))
        else:
            myprint("This shouldn't happen!")
            myprint("    main:", numpy.min(main_times), numpy.max(main_times))
            myprint("    next:", numpy.min(next_times), numpy.max(next_times))
            break
        myprint("mid-point", midpoint)

        connection_point_main = main_slope * midpoint + main_intercept
        connection_point_next = next_slope * midpoint + next_intercept
        # if (~numpy.isfinite(next_slope) or ~numpy.isfinite(next_intercept)):
        #
        # else:
        #     sequence_offset
        sequence_offset = connection_point_next - connection_point_main
        myprint(connection_point_main, connection_point_next, "-->", sequence_offset)

        corrected_reads[next_sequence] -= sequence_offset
        corrections[next_sequence] -= sequence_offset

        # axs[1].scatter(times, corrected_reads, s=2)
        # axs[1].set_title("corrected data after iteration %d" % (iteration))
        # ax.scatter(next_times, next_reads-sequence_offset, s=15, alpha=0.2)
        # break


    return corrected_reads, corrections


