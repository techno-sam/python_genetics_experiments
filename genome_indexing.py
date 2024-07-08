import random
import time
import sys
from utils import random_sequence, describe_bytes

random.seed(0)


def index(genome: str, chunk_length: int = 2) -> dict[str, list[int]]:
    index_dict = {}
    for i in range(0, len(genome) - chunk_length + 1):
        chunk = genome[i:i + chunk_length]
        if chunk not in index_dict:
            index_dict[chunk] = []

        index_dict[chunk].append(i)

    return index_dict


if __name__ == "__main__":
    n = 10_000#_000
    chunk_len = 2

    genome_ = random_sequence(n)

    start = time.perf_counter_ns()
    indexed = index(genome_, chunk_length=chunk_len)
    end = time.perf_counter_ns()

    print(f"Indexing random genome of length {n:,} with chunk length {chunk_len} took {(end - start)/1.0e6:.2f} ms")
    print(f"Index size: {describe_bytes(sys.getsizeof(indexed) + sum(sys.getsizeof(v) for v in indexed.values()))}")
    if n <= 10_000:
        for k, v in indexed.items():
            print(k, v)