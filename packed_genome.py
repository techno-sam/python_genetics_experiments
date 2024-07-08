import math
import typing
import random
from utils import random_sequence

random.seed(0)


class PackedGenome:
    _I = "ACGT"

    def __init__(self, nucleotides: typing.Union[str, 'PackedGenome']):
        if isinstance(nucleotides, PackedGenome):
            self._packed = bytearray(nucleotides._packed)
            self._length = nucleotides._length
        elif isinstance(nucleotides, str):
            self._packed = bytearray()
            for i in range(0, len(nucleotides), 4):
                packed_byte = 0
                for j in range(4):
                    if i + j < len(nucleotides):
                        packed_byte <<= 2
                        packed_byte |= PackedGenome._I.index(nucleotides[i + j])
                    else:
                        packed_byte <<= 2
                self._packed.append(packed_byte)
            self._length = len(nucleotides)
        else:
            raise TypeError(f"Unsupported type for 'nucleotides': '{type(nucleotides).__name__}'")

    def __len__(self):
        return self._length

    def __getitem__(self, idx: int) -> str:
        byte_idx = idx // 4
        shift = 6 - 2 * (idx % 4)
        return PackedGenome._I[(self._packed[byte_idx] >> shift) & 0b11]

    @staticmethod
    def _from_variants(variants: list['PackedGenome'], start: int, length: int) -> 'PackedGenome':
        variant = variants[start % 4]
        packed = variant._packed[start // 4:start // 4 + (length + 3) // 4]

        pg = PackedGenome("")
        pg._packed = bytearray(packed)
        pg._length = length
        return pg

    def subsections(self, chunk_length: int) -> typing.Iterable['PackedGenome']:
        # need 4 variants: this, shifted by 1, 2, 3
        variants = [self]
        for i in range(1, 4):
            variant = PackedGenome("")
            variant._packed = bytearray(variants[-1]._packed)
            # shift all the bytes, incorporating upper bits from next byte
            for j in range(len(variant._packed)):
                variant._packed[j] = (variant._packed[j] << 2) & 0xff
                if j + 1 < len(variant._packed):
                    variant._packed[j] |= (variant._packed[j + 1] >> 6)
            variant._length = variants[-1]._length - 1
            variants.append(variant)

        for i in range(0, len(self) - chunk_length + 1):
            yield PackedGenome._from_variants(variants, i, chunk_length)

        return StopIteration

    def __repr__(self):
        return f"PackedGenome({self._length} {str(self)})"

    def __str__(self):
        return "".join(self[i] for i in range(len(self)))

    def is_in(self, other: typing.Union[str, 'PackedGenome']) -> bool:
        """
        Check if this sequence is in the other sequence.
        """
        if isinstance(other, str):
            return str(self) in other
        elif isinstance(other, PackedGenome):
            return self in other
        else:
            raise TypeError(f"Unsupported operand type(s) for 'is_in': '{type(other).__name__}' and '{type(self).__name__}'")

    def find(self, sub: typing.Union[str, 'PackedGenome'], start: typing.Optional[int] = None, end: typing.Optional[int] = None) -> int:
        """
        S.find(sub[, start[, end]]) -> int

        Return the lowest index in S where subsequence sub is found,
        such that sub is contained within S[start:end].  Optional
        arguments start and end are interpreted as in slice notation.

        Return -1 on failure.
        """
        if isinstance(sub, str):
            return str(self).find(sub, start, end)
        elif isinstance(sub, PackedGenome):
            return str(self).find(str(sub), start, end)
        else:
            raise TypeError(f"Unsupported operand type(s) for 'find': '{type(sub).__name__}' and '{type(self).__name__}'")

    def __contains__(self, item: typing.Union[str, 'PackedGenome']) -> bool:
        if isinstance(item, str):
            return item in str(self)
        elif isinstance(item, PackedGenome):
            return str(item) in str(self)
        else:
            raise TypeError(f"Unsupported operand type(s) for 'in': '{type(item).__name__}' and '{type(self).__name__}'")


class _PreVariedPackedGenome(PackedGenome):
    def __init__(self, nucleotides: typing.Union[str, 'PackedGenome']):
        super().__init__(nucleotides)

        if isinstance(nucleotides, _PreVariedPackedGenome):
            self._variants = [bytearray(variant) for variant in nucleotides._variants]
        else:
            self._variants: list[bytearray] = [self._packed]
            for i in range(1, 4):
                variant = bytearray(self._variants[-1])
                for j in range(len(variant)):
                    variant[j] = (variant[j] << 2) & 0xff
                    if j + 1 < len(variant):
                        variant[j] |= (variant[j + 1] >> 6)
                self._variants.append(variant)

    def __repr__(self):
        return f"PreVariedPackedGenome({self._length} {str(self)})"


class IndexedPackedGenome(_PreVariedPackedGenome):
    def __init__(self, nucleotides: typing.Union[str, 'PackedGenome'], chunk_length: int = 2):
        super().__init__(nucleotides)
        self._index: dict[bytes, list[int]] = {}
        self._chunk_length = chunk_length

        for i in range(0, len(self) - chunk_length + 1):
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

    def __repr__(self):
        return f"IndexedPackedGenome({self._length} {str(self)})"

    def debug_index(self):
        print("Index (Note: if chunk length is not a multiple of 4, there may be trailing 'filler' A's):")
        for k, v in self._index.items():
            # convert to bytearray first
            if isinstance(k, int):
                k = bytearray([k])

            chunk = "".join(PackedGenome._I[(k[j//4] >> (6 - 2 * (j % 4)) & 0b11)] for j in range(len(k) * 4))
            print("\t"+chunk, v)

    def find(self, sub: typing.Union[str, 'PackedGenome'], start: typing.Optional[int] = None, end: typing.Optional[int] = None) -> int:
        """
        S.find(sub[, start[, end]]) -> int

        Return the lowest index in S where subsequence sub is found,
        such that sub is contained within S[start:end].  Optional
        arguments start and end are interpreted as in slice notation.

        Return -1 on failure.
        """
        if isinstance(sub, str):
            return self.find(PackedGenome(sub), start, end)
        elif isinstance(sub, PackedGenome):
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

    def __contains__(self, item: typing.Union[str, 'PackedGenome']) -> bool:
        if isinstance(item, str):
            return PackedGenome(item) in self
        elif isinstance(item, PackedGenome):
            return self.find(item) != -1
        else:
            raise TypeError(f"Unsupported operand type(s) for 'in': '{type(item).__name__}' and '{type(self).__name__}'")


if __name__ == "__main__":
    DEBUG = False
    if DEBUG:
        packed_genome = PackedGenome("ACGTATTCCCGGATC")
        assert packed_genome._packed[0] == 0b00_01_10_11
        for i_ in range(len(packed_genome)):
            print(packed_genome[i_])
        for chunk_ in packed_genome.subsections(3):
            print(repr(chunk_))

        print("GTA" in packed_genome)
        #print(False in packed_genome)
        print(packed_genome.is_in("ACGTATTCCCGGATCACGTATTCCCGGATCACGTATTCCCGGATC"))

        indexed_genome = IndexedPackedGenome("ACGTATTCCCGGATCAAGTATTCCCGAATCACCTATTCCCGGATC", chunk_length=10)
        indexed_genome.debug_index()
        print("GTATTCC" in indexed_genome)
        print("GTATTCC" in PackedGenome(indexed_genome))


    # Simple example
    reference_sequence = random_sequence(500)

    # create an indexed genome, to accelerate searches
    reference_genome = IndexedPackedGenome(reference_sequence, chunk_length=4)
    reference_genome.debug_index()
    print(f"reference_genome = {reference_genome}")

    mutant_example = "GACTACTAGATT"
    nonexistent_example = "TTAGATCATCAG"
    print(f"'{mutant_example}' in reference_genome = {mutant_example in reference_genome}")
    print(f"Index of '{mutant_example}' in reference_genome = {reference_genome.find(mutant_example)}")
    print(f"'{nonexistent_example}' in reference_genome = {nonexistent_example in reference_genome}")
