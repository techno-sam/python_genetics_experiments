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

def _interruptable_worker(send_end, result_converter, func, args, kwargs):
    try:
        res = func(*args, **kwargs)
    except Exception as e:
        send_end.send((False, e))
    else: # If no exception
        send_end.send((True, result_converter(res)))

def _interruptable(func, args=(), kwargs={}, result_converter = lambda x: x):
    import multiprocessing as mp
    recv_end, send_end = mp.Pipe(False)
    p = mp.Process(target=_interruptable_worker, args=(send_end, result_converter, func, args, kwargs))
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
def parse_chromosomes(fasta_path: str, storage_dir: str, index_size: int = 8, num_threads: int = None) -> list[tuple[str, Exception | None]]:
    return _interruptable(
            _parse_chromosomes,
            args=(fasta_path, storage_dir, index_size, _get_num_threads(num_threads))
    )
parse_chromosomes.__doc__ = _parse_chromosomes.__doc__

_HitRecord = HitRecord

class HitRecord:
    """Record of a read-chromosome match"""
    __slots__ = ("_read_name", "_chromosome_start", "_len")

    def __init__(self, read_name: str, chromosome_start: int, len_: int = 0):
        self._read_name: str = read_name
        self._chromosome_start: int = chromosome_start
        self._len: int = len_

    @property
    def read_name(self) -> str:
        """name of the read that was found"""
        return self._read_name

    @property
    def chromosome_start(self) -> int:
        """where in the chromosome the match starts (inclusive)"""
        return self._chromosome_start

    @property
    def chromosome_end(self) -> int:
        """where in the chromosome the match ends (exclusive)"""
        return self._chromosome_start + self._len

    def __repr__(self) -> str:
        return f"HitRecord('{self.read_name}', {self.chromosome_start}, {self._len})"

    def __str__(self) -> str:
        return f"HitRecord({self.read_name} {self.chromosome_start} - {self.chromosome_end})"

def _convert_hit_records(old: list[_HitRecord]) -> list[HitRecord]:
    return [HitRecord(o.read_name, o.chromosome_start) for o in old]

_search_chromosome = search_chromosome
def search_chromosome(reads_path: str, chromosome_name: str, storage_dir: str, kmer_length: int = 31, num_threads: int = None) -> list[HitRecord]:
    out = _interruptable(
            _search_chromosome,
            args=(reads_path, chromosome_name, storage_dir, kmer_length, _get_num_threads(num_threads)),
            result_converter=_convert_hit_records
    )
    for hr in out:
        hr._len = kmer_length
    return out
search_chromosome.__doc__ = _search_chromosome.__doc__

HeatDeathError.__doc__ = "This will never complete. No, seriously."
IndexSizeError.__doc__ = "Invalid size for chromosome index"
InvalidNucleotideCharacterError.__doc__ = "Invalid nucleotide character"

__doc__ = genome_rspy.__doc__
if hasattr(genome_rspy, "__all__"):
    __all__ = genome_rspy.__all__

del genome_rspy

__version__ = "0.1.0"
__author__ = "Sam Wagenaar (techno-sam)"
__copyright__ = "Copyright 2024, Sam Wagenaar"
__license__ = "GPL"
__maintainer__ = "Sam Wagenaar"
__status__ = "Development"
