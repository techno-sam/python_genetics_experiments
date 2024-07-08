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