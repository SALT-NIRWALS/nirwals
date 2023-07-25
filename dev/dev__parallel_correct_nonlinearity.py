#!/usr/bin/env python3

import os
import sys
import astropy.io.fits as pyfits
import numpy
import multiprocessing
import multiprocessing.shared_memory
import multiparlog as mplog
import logging
import queue


def worker_nonlinearity(shmem_raw, shmem_corrected, shmem_nonlinearity_factors,
                        data_shape, nonlin_shape,
                        workqueue,
                        workername="NonLinApply"):

    logger = logging.getLogger(workername)

    # make shared memory accessible here
    _raw = numpy.ndarray(
            shape=data_shape, dtype=numpy.float32,
            buffer=shmem_raw.buf,
        )
    _corrected = numpy.ndarray(
            shape=data_shape, dtype=numpy.float32,
            buffer=shmem_corrected.buf,
        )
    _nonlin = numpy.ndarray(
            shape=nonlin_shape, dtype=numpy.float32,
            buffer=shmem_nonlinearity_factors.buf,
        )

    # derive polynomial degree from number of factors
    poly_order = _nonlin.shape[0] - 1

    while (True):
        try:
            job = workqueue.get(timeout=5)
        except queue.Empty as e:
            job = None

        if (job is None):
            logger.info("Received shutdown command")
            workqueue.task_done()
            break

        y = job
        logger.debug("Working on line %d" % (y))

        linecube = _raw[:,y,:]
        linefactors = _nonlin[:,y,:]

        outbuf = numpy.zeros_like(linecube)
        for p in range(poly_order):
            # logger.info("poly-order %d (img^%d)" % (p, poly_order-p))
            outbuf += linefactors[p] * numpy.power(linecube, poly_order-p)

        _corrected[:,y,:] = outbuf[:,:]
        workqueue.task_done()

    shmem_raw.close()
    shmem_corrected.close()
    shmem_nonlinearity_factors.close()
    logger.debug("Shutting down")




if __name__ == "__main__":


    mplog.setup_logging(debug_filename="debug.log",
                        log_filename="run_analysis.log")
    mpl_logger = logging.getLogger('matplotlib')
    mpl_logger.setLevel(logging.WARNING)

    logger = logging.getLogger("RunAnalysis")

    raw_fn = sys.argv[1]
    nonlin_fn = sys.argv[2]
    output_fn = sys.argv[3]

    #
    # load input data cube and read data into shared memory
    #
    logger.info("Reading input datacube")
    raw_hdu = pyfits.open(raw_fn)
    raw_cube = raw_hdu[1].data
    print(raw_cube.shape)
    print(raw_cube.itemsize, raw_cube.size, raw_cube.itemsize * raw_cube.size / 2**30, "GB")

    logger.info("Allocating shared memory: raw cube")
    data_shape = raw_cube.shape
    shmem_raw = multiprocessing.shared_memory.SharedMemory(
        name='raw_datacube', create=True,
        size=(raw_cube.itemsize * raw_cube.size),
    )
    logger.info("Copying datacube to shared memory")
    _raw = numpy.ndarray(shape=data_shape, dtype=numpy.float32, buffer=shmem_raw.buf)
    _raw[:,:,:] = raw_cube[:,:,:]

    # allocate buffer for corrected data
    logger.info("Allocating shared memory: linearized cube")
    shmem_linearized = multiprocessing.shared_memory.SharedMemory(
        name='linearized_datacube', create=True,
        size=(raw_cube.itemsize * raw_cube.size),
    )
    _linearized = numpy.ndarray(shape=data_shape, dtype=numpy.float32, buffer=shmem_linearized.buf)
    _linearized[:,:,:] = 0.

    logger.info("Reading nonlin poly ")
    nonlin_hdu = pyfits.open(nonlin_fn)
    nonlin_cube = nonlin_hdu[0].data
    nonlin_shape = nonlin_cube.shape
    logger.info("Allocating shared memory: nonlinearity corrections")
    shmem_nonlin = multiprocessing.shared_memory.SharedMemory(
        name='nonlinearity_corrections', create=True,
        size=(nonlin_cube.itemsize * nonlin_cube.size),
    )
    logger.info("Copying datacube to shared memory")
    _nonlin = numpy.ndarray(shape=nonlin_shape, dtype=numpy.float32, buffer=shmem_nonlin.buf)
    _nonlin[:,:,:] = nonlin_cube[:,:,:]

    # Prepare jobs
    jobqueue = multiprocessing.JoinableQueue()
    for y in range(2048):
        jobqueue.put(y)

    # setup and start worker processes
    worker_processes = []
    n_workers = multiprocessing.cpu_count()
    for n in range(n_workers):
        p = multiprocessing.Process(
            target=worker_nonlinearity,
            kwargs=dict(shmem_raw=shmem_raw,
                        shmem_corrected=shmem_linearized,
                        shmem_nonlinearity_factors=shmem_nonlin,
                        data_shape=data_shape, nonlin_shape=nonlin_shape,
                        workqueue=jobqueue,
                        workername="NonLinApply-%03d" % (n+1)),
            daemon=True
        )
        jobqueue.put(None)
        p.start()
        worker_processes.append(p)

    # wait for work to be done
    logger.info("Working for all jobs to finish")
    jobqueue.join()
    for p in worker_processes:
        p.join()

    logger.info("Saving results to %s" % (output_fn))
    pyfits.PrimaryHDU(data=_linearized).writeto(output_fn, overwrite=True)

    logger.info("Freeing up shared memory")
    # release shared memory
    for shmem in [shmem_nonlin, shmem_raw, shmem_linearized]:
        shmem.close()
        shmem.unlink()

    logger.info("All done!")

