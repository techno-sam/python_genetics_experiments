from .genome_rspy import *

def _get_num_threads(num_threads: int = None) -> int:
    import os
    import multiprocessing
    if num_threads is not None:
        return num_threads
    if "PST_NUM_THREADS" in os.environ:
        return int(os.environ["PST_NUM_THREADS"])
    if "NUM_THREADS" in os.environ:
        return int(os.environ["NUM_THREADS"])
    if "MKL_NUM_THREADS" in os.environ:
        return int(os.environ["MKL_NUM_THREADS"])
    return multiprocessing.cpu_count()

def _interruptable_worker(send_end, func, args, kwargs):
    try:
        res = func(*args, **kwargs)
    except Exception as e:
        send_end.send((False, e))
    else: # If no exception
        send_end.send((True, res))

def _interruptable(func, *args, **kwargs):
    import multiprocessing as mp
    recv_end, send_end = mp.Pipe(False)
    p = mp.Process(target=_interruptable_worker, args=(send_end, func, args, kwargs))
    try:
        p.start()
        p.join()
    except Exception as e:
        # print("Terminating")
        p.terminate()
        # print("Terminating: joining")
        p.join()
        # print("Terminated")
        raise e from None

    success, res = recv_end.recv()
    if success:
        return res
    else:
        raise res

_parse_chromosomes = parse_chromosomes
def parse_chromosomes(fasta_path: str, storage_dir: str, num_threads: int = None) -> list[str]:
    return _interruptable(_parse_chromosomes, fasta_path, storage_dir, _get_num_threads(num_threads))
parse_chromosomes.__doc__ = _parse_chromosomes.__doc__

HeatDeathError.__doc__ = "This will never complete. No, seriously."

__doc__ = genome_rspy.__doc__
if hasattr(genome_rspy, "__all__"):
    __all__ = genome_rspy.__all__

del genome_rspy
