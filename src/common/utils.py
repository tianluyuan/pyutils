import glob
import sys
import os
import hashlib
import zlib
import subprocess
import shlex
import time
import subprocess
import multiprocessing as mp


def md5_checksum(filepath):
    """ Stolen with love from:
    http://joelverhagen.com/blog/2011/02/md5-hash-of-file-in-python/

    Reads in a file_path and returns the MD5 hash
    """
    with open(filepath, 'rb') as fh:
        m = hashlib.md5()
        while True:
            data = fh.read(8192)
            if not data:
                break
            m.update(data)
        return m.hexdigest()


def adler32_checksum(filepath):
    """ Modified from md5_checksum function.  This uses the adler32
    algorithm to generate a running checksum.

    This should be faster than md5 ideally.
    """
    with open(filepath, 'rb') as fh:
        curr_chksum = 0
        while True:
            data = fh.read(8192)
            if not data:
                break
            curr_chksum = zlib.adler32(data, curr_chksum)

        return curr_chksum


def find_and_replace(template_path, out_path, **kwargs):
    """ Use a template_path file to generated an out_path file
    by replacing items in the kwargs dict.

    kwargs: keys will be replaced by their values
    """
    print('Making:', out_path)
    with open(out_path, 'w') as out_file:
        # This will be the template by which we build the output
        with open(template_path, 'r') as template_file:
            for line in template_file:
                for key, value in kwargs.items():
                    line = line.replace(key, str(value))

                out_file.write(line)


def make_dirs_if_needed(dir_path):
    if dir_path != '' and not os.path.isdir(dir_path):
        print('making directory ', dir_path)
        os.makedirs(dir_path)


def make_lfcdirs_if_needed(lfn_dir):
    lfc_ls = subprocess.Popen(['lfc-ls', lfn_dir])
    lfc_ls.communicate()
    if lfc_ls.returncode == 1:
        print('making directory:', lfn_dir)
        subprocess.Popen(['lfc-mkdir', '-p', lfn_dir]).communicate()


def chunk(l, n):
    """ Returns list l chunked into lists of length n
    """
    if n < 1:
        n = 1

    return [l[i:i+n] for i in range(0, len(l), n)]


def get_n_jobs(user='tianlu', state='a'):
    """ returns the number of jobs running by 'user' with state 'state'
    """
    qstat = 'qstat -u {user} -s {state}'.format(user=user,
                                                state=state)

    proc_stat = subprocess.Popen(shlex.split(qstat),
                                 stdout=subprocess.PIPE)

    # grep user so that we ignore the first two lines outputted by qstat
    # when counting with 'wc -l'
    proc_grep = subprocess.Popen(['grep', user],
                                 stdin=proc_stat.stdout,
                                 stdout=subprocess.PIPE)

    proc_wc = subprocess.Popen(['wc', '-l'],
                               stdin=proc_grep.stdout,
                               stdout=subprocess.PIPE)

    return int(proc_wc.communicate()[0])


def queue_check(user, limit_queued=100, limit_running=float('inf')):
    """ Checks to see if the number of queued jobs is greater than
    limit_queued for user, and if there are suspended jobs.

    If so, wait 60 seconds and check again.
    """
    while True:
        n_running = get_n_jobs(user, state='r')
        n_queued = get_n_jobs(user, state='p')
        n_suspended = get_n_jobs(user, state='s')

        print('Currently there are...')
        print(n_running, 'jobs running')
        print(n_queued, 'jobs queued')
        print(n_suspended, 'jobs suspended')

        if (n_queued <= limit_queued and
            n_running <= limit_running and
            n_suspended == 0):
            break
        else:
            if n_suspended > 0:
                print('there are jobs suspended so submission will wait')
            elif n_queued > limit_queued:
                print(('the queue is greater than '+str(limit_queued)+''
                       ' so submission will wait'))
            elif n_running > limit_running:
                print(('the running jobs is greater than '+str(limit_running)+''
                       ' so submission will wait'))


            time.sleep(60)


def prompt_yes_no(query):
    """ Prompts the user for yes or no.  Tries to emulate ipython exit prompt.
    """
    bool_str = input(query+' ')
    if bool_str:
        resp = bool_str.capitalize()[0]
        if resp == 'Y':
            return True
        elif resp == 'N':
            return False
        else:
            return prompt_yes_no(query)
    else:
        return True



def almost_equal_relative_and_abs(a, b,
                                  max_diff=sys.float_info.epsilon,
                                  max_rel_diff=sys.float_info.epsilon):
    """ Compares two floating poitn numbers using a relative epsilon
    method and an absulute method if close to zero
    """
    abs_diff = abs(a-b)

    # case where numbers are near zero needs to be absolute
    if abs_diff < max_diff:
        return True

    # relative comparison
    if abs_diff < max(abs(a), abs(b)) * max_rel_diff:
        return True

    return False


def process_cmd(cmd, dry_run):
    print(cmd+'\n')

    if not dry_run:
        return subprocess.Popen(cmd, shell=True).communicate()

    print('=====================================================')


def process_cmd_star(cmd_tup):
    """ convert `f([1,2])` to `f(1,2)` call.  For multiprocessing Pool.
    http://stackoverflow.com/questions/5442910/python-multiprocessing-pool-map-for-multiple-arguments
    """
    return process_cmd(*cmd_tup)


def pool_wrapper(func, iterprocs, nprocs=None, chunksize=1):
    """ Wrapper for submitting multithreaded jobs to the pool

    need to check len(processes) because if len==0 the pool hangs
    http://bugs.python.org/issue6433
    """
    listprocs = list(iterprocs)
    if len(listprocs) > 0:
        # optimizes for the number of cpu threads if we don't pass argument to
        # Pool
        proc_pool = mp.Pool(nprocs)
        proc_pool.map(func, listprocs, chunksize)
        proc_pool.close()
        proc_pool.join()


def split_path(p):
    """ Most robust way to split paths, instead of using str.split('/')
    http://stackoverflow.com/questions/3167154/how-to-split-a-dos-path-into-its-components-in-python
    """
    a,b = os.path.split(os.path.normpath(p))
    return (split_path(a) if len(a) and len(b) else []) + [b]


def get_cwd():
    """ Since we're working in a mixed sl5/sl6 environment paths
    beginning with /amd/... can cause issues as sl6 uses autofs, which
    does not rely on the automounter for nfs shares.

    Try to use 'pawd' command first, on error revert to os.getcwd()
    """
    return (subprocess.getoutput('pawd')
            if not subprocess.getstatusoutput('pawd')[0] # 0 on success
            else os.getcwd())


def abspath(path):
    """ Removes malformed /amd/blah from full path
    """
    pawd = get_cwd()
    return os.path.normpath(os.path.join(pawd,
                                         os.path.relpath(
                                             os.path.expanduser(path))))


def glob_newest(match_str):
    """ Returns the newest file from list of files matching match_str.
    Useful for when we want a single file from a set of possible matches.
    """
    try:
        return max(glob.iglob(match_str), key=os.path.getctime)
    except ValueError:
        print('Error in utils.glob_newest: No matching files for', match_str)


def frange(x, y, jump):
    """Simple substitute for numpy.arange. numpy doesn't work on certain
    cluster machines
    """
    while x < y:
        yield x
        x += jump


def z_to_pdg(Z, A):
    """ Converts atomic Z to pdg
    """
    return 1000000000 + Z*10000 + A*10


def timefn(fn, *args, **kwargs):
    """ Runs fn and prints out the elapsed time

    returns result of fn call
    """
    if hasattr(fn, '__call__'):
        start = time.clock()
        ret = fn(*args, **kwargs)
        end = time.clock()
        print('Time elapsed for {0}: {1}'.format(fn.__name__,
                                             end - start))

        return ret
    else:
        print('Warning <utils.timefn>: non-callable passed')


def midpoints(l):
    """ Returns the midpoints between pairs in list l
    """
    import numpy as np
    narr = np.array(l)
    return (narr[:-1]+narr[1:])/2.
