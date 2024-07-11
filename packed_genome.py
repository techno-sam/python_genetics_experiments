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

import math
import typing
import numpy as np
from typing import Optional

try:
    from tqdm import tqdm
except ImportError:
    tqdm = lambda x, **kwargs: x


class PackedSequence:
    _I = "ACGT"
    _LOOKUP = {n: i for i, n in enumerate(_I)}

    __slots__ = ["_packed", "_length"]

    def __init__(self, nucleotides: typing.Union[str, 'PackedSequence'], verbose: bool = False):
        if isinstance(nucleotides, PackedSequence):
            self._packed = bytearray(nucleotides._packed)
            self._length = nucleotides._length
        elif isinstance(nucleotides, str):
            self._packed = bytearray()
            if verbose:
                print("> packing...", end="")

            memoized_bytes: dict[str, int] = {}

            tqdm_ = tqdm(range(0, len(nucleotides), 4), desc="packing") if verbose else range(0, len(nucleotides), 4)
            for i in tqdm_:
                packed_byte = 0
                if i + 4 <= len(nucleotides):
                    key = nucleotides[i:i+4]
                    if key in memoized_bytes:
                        packed_byte = memoized_bytes[key]
                    else:
                        for j in range(4):
                            packed_byte <<= 2
                            packed_byte |= PackedSequence._LOOKUP[nucleotides[i + j]]
                        memoized_bytes[key] = packed_byte
                else:
                    for j in range(4):
                        if i + j < len(nucleotides):
                            packed_byte <<= 2
                            packed_byte |= PackedSequence._LOOKUP[nucleotides[i + j]]
                        else:
                            packed_byte <<= 2
                self._packed.append(packed_byte)

                #if verbose and i % 1000 == 0:
                #    print(f"\r> packing... {i}/{len(nucleotides)} ({i/len(nucleotides)*100:.2f}%)", end="")
            self._length = len(nucleotides)
            if verbose:
                print()
        else:
            raise TypeError(f"Unsupported type for 'nucleotides': '{type(nucleotides).__name__}'")

    def __len__(self):
        return self._length

    def __getitem__(self, idx: int) -> str:
        byte_idx = idx // 4
        shift = 6 - 2 * (idx % 4)
        return PackedSequence._I[(self._packed[byte_idx] >> shift) & 0b11]

    @staticmethod
    def _from_variants(variants: list[bytearray], start: int, length: int) -> 'PackedSequence':
        variant = variants[start % 4]
        packed = variant[start // 4:start // 4 + (length + 3) // 4]

        pg = PackedSequence("")
        pg._packed = bytearray(packed)
        pg._length = length
        return pg

    def subsections(self, chunk_length: int) -> typing.Iterable['PackedSequence']:
        # need 4 variants: this, shifted by 1, 2, 3
        variants = [np.array(self._packed)]
        for i in range(1, 4):
            variant = np.array(variants[-1], dtype=np.uint8)
            ## shift all the bytes, incorporating upper bits from next byte
            variant_copy = np.copy(variant)
            variant <<= 2
            variant &= 0xff
            variant[:-1] |= variant_copy[1:] >> 6
            variants.append(variant)

        variants = [bytearray(variant) for variant in variants]

        for i in range(0, len(self) - chunk_length + 1):
            yield PackedSequence._from_variants(variants, i, chunk_length)

        return StopIteration

    def __repr__(self):
        return f"PackedSequence({self._length} {str(self)})"

    def __str__(self):
        return "".join(self[i] for i in range(len(self)))

    def is_in(self, other: typing.Union[str, 'PackedSequence']) -> bool:
        """
        Check if this sequence is in the other sequence.
        """
        if isinstance(other, str):
            return str(self) in other
        elif isinstance(other, PackedSequence):
            return self in other
        else:
            raise TypeError(f"Unsupported operand type(s) for 'is_in': '{type(other).__name__}' and '{type(self).__name__}'")

    def find(self, sub: typing.Union[str, 'PackedSequence'], start: typing.Optional[int] = None, end: typing.Optional[int] = None) -> int:
        """
        S.find(sub[, start[, end]]) -> int

        Return the lowest index in S where subsequence sub is found,
        such that sub is contained within S[start:end].  Optional
        arguments start and end are interpreted as in slice notation.

        Return -1 on failure.
        """
        if isinstance(sub, str):
            return str(self).find(sub, start, end)
        elif isinstance(sub, PackedSequence):
            return str(self).find(str(sub), start, end)
        else:
            raise TypeError(f"Unsupported operand type(s) for 'find': '{type(sub).__name__}' and '{type(self).__name__}'")

    def __contains__(self, item: typing.Union[str, 'PackedSequence']) -> bool:
        if isinstance(item, str):
            return item in str(self)
        elif isinstance(item, PackedSequence):
            return str(item) in str(self)
        else:
            raise TypeError(f"Unsupported operand type(s) for 'in': '{type(item).__name__}' and '{type(self).__name__}'")


class PreVariedPackedSequence(PackedSequence):

    __slots__ = ["_variants"]

    def __init__(self, nucleotides: typing.Union[str, 'PackedSequence'], verbose: bool = False):
        super().__init__(nucleotides, verbose=verbose)

        if isinstance(nucleotides, PreVariedPackedSequence):
            self._variants = [bytearray(variant) for variant in nucleotides._variants]
        else:
            if verbose:
                print("> varying...", end="")
            self._variants: list[bytearray] = [self._packed]
            for i in range(1, 4):
                if verbose:
                    print(f"\r> varying... {i+1}/4", end="")
                variant = np.array(self._variants[-1], dtype=np.uint8)
                variant_copy = np.copy(variant)
                variant <<= 2
                variant &= 0xff
                variant[:-1] |= variant_copy[1:] >> 6
                self._variants.append(bytearray(variant))
            if verbose:
                print()

    def __repr__(self):
        return f"PreVariedPackedSequence({self._length} {str(self)})"

    def subsections(self, chunk_length: int) -> typing.Iterable['PackedSequence']:
        for i in range(0, len(self) - chunk_length + 1):
            yield PackedSequence._from_variants(self._variants, i, chunk_length)

        return StopIteration


class IndexedPackedSequence(PreVariedPackedSequence):

    __slots__ = ["_index", "_chunk_length"]

    def __init__(
            self,
            nucleotides: typing.Union[str, 'PackedSequence'],
            chunk_length: int = 2,
            no_index_last: Optional[int] = None,
            verbose: bool = False
    ):
        """Create a packed sequence with an index for faster searching.

        :param nucleotides: string of format 'ACGTGGATACA' or PackedSequence
        :param chunk_length: length of the chunks to index, longer chunks result in faster searching but require more memory
        :param no_index_last: if not None, the last 'no_index_last' nucleotides will not be indexed
        :param verbose: if True, print progress
        """
        if verbose:
            print("> packing and varying...")
        super().__init__(nucleotides, verbose=verbose)
        self._index: dict[bytes, list[int]] = {}
        self._chunk_length = chunk_length

        skip_last = chunk_length
        if no_index_last is not None:
            skip_last = max(skip_last, no_index_last)

        if verbose:
            print("> indexing...", end="")
        max_idx = len(self) - skip_last + 1
        tqdm_ = tqdm(range(0, max_idx), desc="indexing") if verbose else range(0, max_idx)
        for i in tqdm_:
            sample_from = self._variants[i % 4]
            chunk = bytearray()
            for j in range(math.ceil(chunk_length / 4)):
                chunk.append(sample_from[i//4 + j])
            # zero out the extra bits
            if chunk_length % 4:
                chunk[len(chunk)-1] &= (0xff << (8 - 2 * (chunk_length % 4))) & 0xff
            # chunk has to be converted to a hashable type
            chunk = bytes(chunk)
            if chunk not in self._index:
                self._index[chunk] = []
            self._index[chunk].append(i)

            # if verbose and i % 10_000 == 0:
            #     print(f"\r> indexing... {i}/{max_idx} ({i/max_idx*100:.2f}%)", end="")
        if verbose:
            print()

    def __repr__(self):
        return f"IndexedPackedSequence({self._length} {str(self)})"

    def debug_index(self):
        print("Index (Note: if chunk length is not a multiple of 4, there may be trailing 'filler' A's):")
        for k, v in self._index.items():
            # convert to bytearray first
            if isinstance(k, int):
                k = bytearray([k])

            chunk = "".join(PackedSequence._I[(k[j // 4] >> (6 - 2 * (j % 4)) & 0b11)] for j in range(len(k) * 4))
            print("\t"+chunk, v)

    def find(self, sub: typing.Union[str, 'PackedSequence'], start: typing.Optional[int] = None, end: typing.Optional[int] = None) -> int:
        """
        S.find(sub[, start[, end]]) -> int

        Return the lowest index in S where subsequence sub is found,
        such that sub is contained within S[start:end].  Optional
        arguments start and end are interpreted as in slice notation.

        Return -1 on failure.
        """
        if isinstance(sub, str):
            return self.find(PackedSequence(sub), start, end)
        elif isinstance(sub, PackedSequence):
            if len(sub) == 0:
                raise ValueError("Cannot search for an empty sequence, the return value is ambiguous.")
            if len(self) == 0:
                return -1

            if self._chunk_length > len(sub):
                print("Warning: chunk length is greater than the length of the item, falling back to slow implementation.")
                return str(self).find(str(sub), start, end)

            key = sub._packed[:math.ceil(self._chunk_length / 4)]
            if self._chunk_length % 4:
                key[len(key)-1] &= (0xff << (8 - 2 * (self._chunk_length % 4))) & 0xff
            key = bytes(key)

            if key not in self._index:
                return -1

            for i in self._index[key]:
                if start and i < start:
                    continue
                if end and i + len(sub) > end:
                    break

                comparison_section = self._variants[i % 4][i // 4:i // 4 + (len(sub) + 3) // 4]
                if len(sub) % 4:
                    comparison_section[len(comparison_section)-1] &= (0xff << (8 - 2 * (len(sub) % 4)) & 0xff)

                if sub._packed == comparison_section:
                    return i

            return -1
        else:
            raise TypeError(f"Unsupported operand type(s) for 'index_of': '{type(sub).__name__}' and '{type(self).__name__}'")

    def __contains__(self, item: typing.Union[str, 'PackedSequence']) -> bool:
        if isinstance(item, str):
            return PackedSequence(item) in self
        elif isinstance(item, PackedSequence):
            return self.find(item) != -1
        else:
            raise TypeError(f"Unsupported operand type(s) for 'in': '{type(item).__name__}' and '{type(self).__name__}'")


if __name__ == "__main__":
    import random
    from utils import random_sequence

    random.seed(0)

    DEBUG = False
    if DEBUG:
        packed_genome = PackedSequence("ACGTATTCCCGGATC")
        assert packed_genome._packed[0] == 0b00_01_10_11
        for i_ in range(len(packed_genome)):
            print(packed_genome[i_])
        for chunk_ in packed_genome.subsections(3):
            print(repr(chunk_))

        print("GTA" in packed_genome)
        #print(False in packed_genome)
        print(packed_genome.is_in("ACGTATTCCCGGATCACGTATTCCCGGATCACGTATTCCCGGATC"))

        indexed_genome = IndexedPackedSequence("ACGTATTCCCGGATCAAGTATTCCCGAATCACCTATTCCCGGATC", chunk_length=10)
        indexed_genome.debug_index()
        print("GTATTCC" in indexed_genome)
        print("GTATTCC" in PackedSequence(indexed_genome))


    # Simple example
    reference_sequence = random_sequence(500)

    # create an indexed genome, to accelerate searches
    reference_genome = IndexedPackedSequence(reference_sequence, chunk_length=4)
    reference_genome.debug_index()
    print(f"reference_genome = {reference_genome}")

    mutant_example = "GACTACTAGATT"
    nonexistent_example = "TTAGATCATCAG"
    print(f"'{mutant_example}' in reference_genome = {mutant_example in reference_genome}")
    print(f"Index of '{mutant_example}' in reference_genome = {reference_genome.find(mutant_example)}")
    print(f"'{nonexistent_example}' in reference_genome = {nonexistent_example in reference_genome}")
