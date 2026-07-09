from __future__ import annotations

import sys
from dataclasses import dataclass
from typing import Literal, TypedDict

if sys.version_info >= (3, 11):
    from typing import NotRequired
else:
    from typing_extensions import NotRequired

Company = Literal["23andme", "ancestrydna", "ftdna", "myheritage", "decodeme", "livingdnacsv", "all"]

RejectType = Literal[
    "invalid", "small", "incomplete", "grch36", "grch38",
    "notgrch37", "vcf", "zip", "ref", "imputation",
]


class SNPRecord(TypedDict):
    rsid: str
    chro: str
    pos: str
    call: str
    ref: str
    gav: str
    ba: NotRequired[bool]


@dataclass
class ProcessResult:
    status: str
    info: dict[str, int | float | str] | str
