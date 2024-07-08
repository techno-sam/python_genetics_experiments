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

__all__ = ["describe_bytes", "random_sequence"]

import random

def describe_bytes(size: int) -> str:
    """
    Describe a size in bytes in human-readable format.

    :param size:    size in bytes
    """

    steps = ["bytes", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"]

    unit = steps.pop(0)

    while size >= 1000 and steps:
        size /= 1000
        unit = steps.pop(0)

    if unit == "bytes":
        return f"{size} {unit}"
    return f"{size:.2f} {unit}"


def random_sequence(length: int) -> str:
    return "".join(random.choices("ACGT", k=length))