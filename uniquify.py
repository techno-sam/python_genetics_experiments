"""
Python Genetics Experiments
Copyright (C) 2024  Sam Wagenaar

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import random
import time
from typing import Optional
from colorama import Fore
from utils import random_sequence

random.seed(0)


def highlight_string(v: str, ranges: list[tuple[int, int, Fore]]):
    """
    Highlight the given string with the given ranges.

    :param v: The string to be highlighted.
    :param ranges: The ranges to be highlighted. Each range should be a tuple of (start, end, color).
    :return: The highlighted string.
    """
    # must handle overlapping ranges (later ranges take higher priority, e.g. [(0, 4, Fore.Red), (2, 3, Fore.Green)] would result in a red highlight with a green splotch in the middle)
    color_stack = []

    color_ops: list[tuple[int, Optional[Fore]]] = []
    for start, end, color in ranges:
        color_ops.append((start, color))
        color_ops.append((end, None))

    color_ops_by_idx = {idx: [] for idx, _ in color_ops}
    for idx, color in color_ops:
        color_ops_by_idx[idx].append(color)

    out = ""
    last_color = None
    for i in range(len(v)):
        if i in color_ops_by_idx:
            for op in color_ops_by_idx[i]:
                if op:
                    color_stack.append(op)
                else:
                    color_stack.pop()

        if color_stack:
            color = color_stack[-1]
        else:
            color = None

        if color != last_color:
            if color:
                out += color
            else:
                out += Fore.RESET
            last_color = color

        out += v[i]

    return out + Fore.RESET


def unique_context(sequence: str, target_idx: int, original_context_length: int = 5, verbose: bool = False) -> str:
    # try to get context centered around the target index, but shift if near edges
    start_idx = max(0, target_idx - original_context_length // 2)
    end_idx = min(len(sequence), start_idx + original_context_length)
    start_idx = max(0, end_idx - original_context_length)

    original_context = sequence[start_idx:end_idx]
    if verbose:
        print(f"Original context: {original_context} (sequence[{start_idx}:{end_idx}])", end=" ")
        basic_highlight = [(start_idx, end_idx, Fore.BLUE), (target_idx, target_idx + 1, Fore.RED)]
        print(highlight_string(sequence, basic_highlight))
    else:
        basic_highlight = []

    found_ranges: list[tuple[int, int]] = []

    end = 0
    while True:
        start = sequence.find(original_context, end)
        if start == -1:
            break
        end = start + len(original_context)
        found_ranges.append((start, end))

    found_ranges.remove((start_idx, end_idx))

    if verbose:
        print(f"Found {len(found_ranges)} occurrences. {found_ranges}")
        print(highlight_string(sequence, [(start, end, Fore.GREEN) for start, end in found_ranges] + basic_highlight))

    # start extending the context
    new_start = start_idx
    new_end = end_idx
    while len(found_ranges) > 0:
        candidates: list[tuple[int, int]] = []
        if new_start > 0:
            candidates.append((-1, 0))
        if new_end < len(sequence):
            candidates.append((0, 1))

        best_candidate: tuple[int, int] = candidates[0]
        best_candidate_score = len(found_ranges)
        best_new_found_ranges: list[tuple[int, int]] = found_ranges
        for cs_ofs, ce_ofs in candidates:
            candidate_score = 0
            new_ranges: list[tuple[int, int]] = []
            for range_start, range_end in found_ranges:
                if range_start + cs_ofs >= 0 and range_end + ce_ofs <= len(sequence) and sequence[new_start+cs_ofs] == sequence[range_start+cs_ofs] and sequence[new_end+ce_ofs-1] == sequence[range_end+ce_ofs-1]:
                    candidate_score += 1
                    new_ranges.append((range_start+cs_ofs, range_end+ce_ofs))
            if candidate_score < best_candidate_score:
                best_candidate = (cs_ofs, ce_ofs)
                best_candidate_score = candidate_score
                best_new_found_ranges = new_ranges

        new_start += best_candidate[0]
        new_end += best_candidate[1]
        found_ranges = best_new_found_ranges
        if verbose:
            print(f"New context: {sequence[new_start:new_end]}", end=" ")
            print(highlight_string(
                sequence,
                [(start, end, Fore.GREEN) for start, end in found_ranges]
                + [(new_start, new_end, Fore.BLUE), (target_idx, target_idx + 1, Fore.RED)]
            ), end=" ")
            print("Score:", best_candidate_score)

    final_context = sequence[new_start:new_end]

    final_count = sequence.count(final_context)
    assert final_count == 1, f"Final count {final_count} is not 1!"
    if verbose:
        print(f"Final context: {final_context}[{new_start}:{new_end}]", end=" ")
        print(highlight_string(sequence, [(new_start, new_end, Fore.BLUE), (target_idx, target_idx + 1, Fore.RED)]))

    return final_context



if __name__ == "__main__":
    seq = random_sequence(1_000_000)
    original_print = print
    print(seq)

    time_start = time.time()
    uc = unique_context(seq, 5, original_context_length=30)
    time_end = time.time()
    print(f"Time taken: {time_end - time_start:.2f}s")
    print(f"Result: {uc}")