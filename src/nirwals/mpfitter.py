#!/usr/bin/env python3

import scipy
import numpy
import multiprocessing

import nirwals_reduce

import sys



def mpworker(rss, job_queue, result_queue, executable_function):
    print("Worker has started")
    while (True):
        job = job_queue.get()
        if (job is None):
            job_queue.task_done()
            print("Worker shutting down")
            break

        [x, y] = job
        if ((y % 10) == 0 and x == 0): print(job)

        result = executable_function(rss, x, y)
        result_queue.put((x, y, result))

        job_queue.task_done()


def fit_fct(rss, x, y):
    print(rss.filebase, x, y)
    return((x,y,[1,x,y]))





if __name__ == "__main__":

    # get some data to play with
    fn = sys.argv[1]
    rss = rss_reduce.RSS(fn, max_number_files=40, )
    rss.load_all_files()
    # rss.reduce(write_dumps=False)

    iy,ix = numpy.indices((rss.image_stack.shape[1], rss.image_stack.shape[2]))
    print(iy.shape)

    job_queue = multiprocessing.JoinableQueue()
    result_queue = multiprocessing.Queue()
    ixy = numpy.dstack([ix,iy]).reshape((-1,2))
    print(ixy.shape)
    print(ixy[:10])

    # prepare jobqueue
    print("preparing jobs")
    for _xy in ixy:
        job_queue.put(_xy)

    # prepare workers
    print("Creating workers")
    worker_processes = []
    for n in range(3):
        p = multiprocessing.Process(
            target=mpworker,
            kwargs=dict(rss=rss,
                        job_queue=job_queue,
                        result_queue=result_queue,
                        executable_function=fit_fct,
                        )
        )
        p.daemon = True
        p.start()
        worker_processes.append(p)

        job_queue.put(None)

    # wait for completion
    print("Waiting for completion")
    job_queue.join()

    print("all done!")