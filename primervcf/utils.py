# std library import
import multiprocessing
import subprocess
import signal
import os
import sys

def del_files(filename_l):
    for filename in filename_l:
        try:
            os.remove(filename)
        except OSError:
            pass


def parmap(f, X, nprocs=multiprocessing.cpu_count()):
    """
    a function to use muti map inside a function
    modified from stackoverflow, 3288595
    :param f:
    :param X:
    :param nprocs: core, if not given, use all core
    :return:
    """
    q_in = multiprocessing.Queue(1)
    q_out = multiprocessing.Queue()

    proc = [multiprocessing.Process(target=fun, args=(f, q_in, q_out))
            for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i, x)) for i, x in enumerate(X)]
    [q_in.put((None, None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]

    [p.join() for p in proc]

    return [x for i, x in sorted(res)]


def fun(f, q_in, q_out):
    """
    for parmap
    :param f:
    :param q_in:
    :param q_out:
    :return:
    """
    while True:
        i, x = q_in.get()
        if i is None:
            break
        q_out.put((i, f(x)))


def myexe(cmd, timeout=10):
    """
    a simple wrap of the shell
    mainly used to run the bwa mem mapping and samtool orders
    """
    def setupAlarm():
        signal.signal(signal.SIGALRM, alarmHandler)
        signal.alarm(timeout)

    def alarmHandler(signum, frame):
        sys.exit(1)

    proc=subprocess.Popen(cmd, shell=True, preexec_fn=setupAlarm,
                 stdout=subprocess.PIPE, stderr=subprocess.PIPE,cwd=os.getcwd())
    #print("Running:  ",  cmd)
    out, err=proc.communicate()
    #print("stderr out:  ",err, "\n","return code:  ",proc.returncode)
    return out