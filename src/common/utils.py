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
    """
    Compute the MD5 checksum of a file.

    See http://joelverhagen.com/blog/2011/02/md5-hash-of-file-in-python/

    Parameters
    ----------
    filepath : str
        Path to the file for which the MD5 checksum is to be calculated.

    Returns
    -------
    str
        The MD5 checksum of the file.
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
    """
    Compute the Adler-32 checksum of a file.

    Modified from the md5_checksum function. This uses the Adler-32
    algorithm to generate a running checksum. This should be faster than MD5 ideally.

    Parameters
    ----------
    filepath : str
        Path to the file for which the Adler-32 checksum is to be calculated.

    Returns
    -------
    int
        The Adler-32 checksum of the file.
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
    """
    Replace placeholders in a template file and save the output to a new file.

    Parameters
    ----------
    template_path : str
        Path to the template file.
    out_path : str
        Path where the output file will be saved.
    **kwargs : dict
        Dictionary containing placeholder keys and their replacement values.

    Returns
    -------
    None
    """
    print('Making:', out_path)
    with open(out_path, 'w') as out_file:
        with open(template_path, 'r') as template_file:
            for line in template_file:
                for key, value in kwargs.items():
                    line = line.replace(key, str(value))
                out_file.write(line)


def make_dirs_if_needed(dir_path):
    """
    Create a directory if it does not already exist.

    Parameters
    ----------
    dir_path : str
        Path to the directory.

    Returns
    -------
    None
    """
    if dir_path != '' and not os.path.isdir(dir_path):
        print('making directory ', dir_path)
        os.makedirs(dir_path)


def make_lfcdirs_if_needed(lfn_dir):
    """
    Create a logical file catalog (LFC) directory if it does not already exist.

    Parameters
    ----------
    lfn_dir : str
        Path to the logical file catalog directory.

    Returns
    -------
    None
    """
    lfc_ls = subprocess.Popen(['lfc-ls', lfn_dir])
    lfc_ls.communicate()
    if lfc_ls.returncode == 1:
        print('making directory:', lfn_dir)
        subprocess.Popen(['lfc-mkdir', '-p', lfn_dir]).communicate()


def chunk(l, n):
    """
    Split a list into smaller chunks of a specified size.

    Parameters
    ----------
    l : list
        The list to be split.
    n : int
        Size of each chunk.

    Returns
    -------
    list of list
        List of chunks.
    """
    if n < 1:
        n = 1

    return [l[i:i+n] for i in range(0, len(l), n)]


def get_n_jobs(user='tianlu', state='a'):
    """
    Get the number of jobs running for a given user and state.

    Parameters
    ----------
    user : str, optional
        Username to check for jobs. Default is 'tianlu'.
    state : str, optional
        State of the jobs ('r' for running, 'p' for pending, etc.). Default is 'a'.

    Returns
    -------
    int
        The number of jobs.
    """
    qstat = f'qstat -u {user} -s {state}'

    proc_stat = subprocess.Popen(shlex.split(qstat), stdout=subprocess.PIPE)

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
    """
    Check if the number of queued or running jobs exceeds the specified limits.

    If the queue limit is exceeded or there are suspended jobs, the function waits 60 seconds
    and checks again.

    Parameters
    ----------
    user : str
        Username to check for jobs.
    limit_queued : int, optional
        Maximum number of queued jobs allowed. Default is 100.
    limit_running : int, optional
        Maximum number of running jobs allowed. Default is infinity.

    Returns
    -------
    None
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
    """
    Prompt the user with a yes/no question.

    Parameters
    ----------
    query : str
        The question to prompt the user.

    Returns
    -------
    bool
        True if the user responds with 'yes', False otherwise.
    """
    bool_str = input(f"{query} ")
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


def almost_equal_relative_and_abs(a, b, max_diff=sys.float_info.epsilon, max_rel_diff=sys.float_info.epsilon):
    """
    Compare two floating-point numbers using both relative and absolute tolerances.

    Parameters
    ----------
    a : float
        The first number.
    b : float
        The second number.
    max_diff : float, optional
        Maximum absolute difference. Default is machine epsilon.
    max_rel_diff : float, optional
        Maximum relative difference. Default is machine epsilon.

    Returns
    -------
    bool
        True if the numbers are considered nearly equal, False otherwise.
    """
    abs_diff = abs(a-b)

    # Case where numbers are near zero needs to be absolute comparison
    if abs_diff < max_diff:
        return True

    # Relative comparison
    if abs_diff < max(abs(a), abs(b)) * max_rel_diff:
        return True

    return False


def process_cmd(cmd, dry_run):
    """
    Execute a shell command.

    Parameters
    ----------
    cmd : str
        The command to execute.
    dry_run : bool
        If True, the command is not executed.

    Returns
    -------
    tuple
        The output of the command, if executed.
    """
    print(f"{cmd}\n")

    if not dry_run:
        return subprocess.Popen(cmd, shell=True).communicate()

    print('=====================================================')


def process_cmd_star(cmd_tup):
    """
    Convert `f([1,2])` to `f(1,2)` call for multiprocessing Pool.

    Reference:
    http://stackoverflow.com/questions/5442910/python-multiprocessing-pool-map-for-multiple-arguments

    Parameters
    ----------
    cmd_tup : tuple
        Command and arguments to execute.

    Returns
    -------
    tuple
        The output of the command, if executed.
    """
    return process_cmd(*cmd_tup)


def pool_wrapper(func, iterprocs, nprocs=None, chunksize=1):
    """
    Submit multithreaded jobs to a multiprocessing pool.

    Reference:
    http://bugs.python.org/issue6433

    Parameters
    ----------
    func : callable
        The function to execute.
    iterprocs : iterable
        Iterable of processes to run.
    nprocs : int, optional
        Number of processes in the pool. Default is None.
    chunksize : int, optional
        Number of tasks to assign to each worker. Default is 1.

    Returns
    -------
    None
    """
    listprocs = list(iterprocs)
    if len(listprocs) > 0:
        proc_pool = mp.Pool(nprocs)
        proc_pool.map(func, listprocs, chunksize)
        proc_pool.close()
        proc_pool.join()


def split_path(p):
    """
    Split a file path into its components.

    Parameters
    ----------
    p : str
        The file path to split.

    Returns
    -------
    list of str
        List of path components.

    Reference
    ---------
    http://stackoverflow.com/questions/3167154/how-to-split-a-dos-path-into-its-components-in-python
    """
    a, b = os.path.split(os.path.normpath(p))
    return (split_path(a) if len(a) and len(b) else []) + [b]


def get_cwd():
    """
    Get the current working directory, handling edge cases for automounted paths.

    Reference:
    SL5/SL6 mixed environment paths, where paths beginning with /amd/... may cause issues 
    due to SL6's autofs.

    Returns
    -------
    str
        The current working directory.
    """
    return (subprocess.getoutput('pawd')
            if not subprocess.getstatusoutput('pawd')[0]  # 0 on success
            else os.getcwd())


def glob_newest(match_str):
    """
    Get the newest file matching a pattern.

    Parameters
    ----------
    match_str : str
        The pattern to match.

    Returns
    -------
    str
        The newest file matching the pattern, or None if no files match.
    """
    try:
        return max(glob.iglob(match_str), key=os.path.getctime)
    except ValueError:
        print('Error in utils.glob_newest: No matching files for', match_str)


def frange(x, y, jump):
    """
    Generate a range of floating-point numbers.

    Parameters
    ----------
    x : float
        Start of the range.
    y : float
        End of the range.
    jump : float
        Step size.

    Yields
    ------
    float
        The next number in the range.
    """
    while x < y:
        yield x
        x += jump


def z_to_pdg(Z, A):
    """
    Convert atomic number and mass number to a PDG code.

    Parameters
    ----------
    Z : int
        Atomic number.
    A : int
        Mass number.

    Returns
    -------
    int
        The PDG code.
    """
    return 1000000000 + Z*10000 + A*10


def timefn(fn, *args, **kwargs):
    """
    Measure the execution time of a function.

    Parameters
    ----------
    fn : callable
        The function to measure.
    *args : tuple
        Positional arguments to pass to the function.
    **kwargs : dict
        Keyword arguments to pass to the function.

    Returns
    -------
    Any
        The return value of the function.

    Notes
    -----
    Prints the elapsed time for the function execution.
    """
    if hasattr(fn, '__call__'):
        start = time.clock()
        ret = fn(*args, **kwargs)
        end = time.clock()
        print(f'Time elapsed for {fn.__name__}: {end - start}')

        return ret
    else:
        print('Warning <utils.timefn>: non-callable passed')


def midpoints(l):
    """
    Compute the midpoints between pairs in a list.

    Parameters
    ----------
    l : list of float
        List of values.

    Returns
    -------
    numpy.ndarray
        Array of midpoints.
    """
    import numpy as np
    narr = np.array(l)
    return (narr[:-1]+narr[1:])/2.
