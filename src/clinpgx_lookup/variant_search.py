"""Variant lookup against the PharmGKB API.

rsIDs (``rs...``) are looked up via the variant endpoint; everything else is
treated as a star allele and looked up via the haplotype endpoint.
"""

from typing import List

import requests
from loguru import logger
from pydantic import BaseModel

from .search_utils import calc_similarity

PGKB_BASE = "https://api.pharmgkb.org/v1/data"
DEFAULT_TIMEOUT = 10


class VariantSearchResult(BaseModel):
    raw_input: str
    id: str
    name: str
    url: str
    score: float
    source: str = "pharmgkb"


def _pgkb_query(endpoint: str, symbol: str) -> List[dict]:
    """Query a PharmGKB data endpoint by ``symbol`` and return the data list."""
    try:
        response = requests.get(
            f"{PGKB_BASE}/{endpoint}",
            params={"symbol": symbol.strip()},
            timeout=DEFAULT_TIMEOUT,
        )
    except requests.RequestException as exc:
        logger.warning(f"PharmGKB {endpoint} request failed: {exc}")
        return []
    if response.status_code != 200:
        return []
    return response.json().get("data") or []


def pgkb_star_allele_search(
    star_allele: str, top_k: int = 1
) -> List[VariantSearchResult]:
    data = _pgkb_query("haplotype", star_allele)
    return [
        VariantSearchResult(
            raw_input=star_allele,
            id=result["id"],
            name=result["symbol"],
            url=f"https://www.clinpgx.org/haplotype/{result['id']}",
            score=calc_similarity(star_allele, result["symbol"]),
        )
        for result in data[:top_k]
    ]


def pgkb_rsid_search(rsid: str, top_k: int = 1) -> List[VariantSearchResult]:
    data = _pgkb_query("variant", rsid)
    return [
        VariantSearchResult(
            raw_input=rsid,
            id=result["id"],
            name=result["symbol"],
            url=f"https://www.clinpgx.org/variant/{result['id']}",
            score=calc_similarity(rsid, result["symbol"]),
        )
        for result in data[:top_k]
    ]


class VariantLookup:
    """Look up variants (rsIDs and star alleles) in PharmGKB."""

    def rsid_lookup(
        self, rsid: str, threshold: float = 0.8, top_k: int = 1
    ) -> List[VariantSearchResult]:
        results = pgkb_rsid_search(rsid, top_k=top_k)
        results.sort(key=lambda result: result.score, reverse=True)
        return results[:top_k]

    def star_lookup(
        self, star_allele: str, threshold: float = 0.8, top_k: int = 1
    ) -> List[VariantSearchResult]:
        results = pgkb_star_allele_search(star_allele, top_k=top_k)
        results.sort(key=lambda result: result.score, reverse=True)
        return results[:top_k]

    def search(
        self, variant: str, threshold: float = 0.8, top_k: int = 1
    ) -> List[VariantSearchResult]:
        """Look up ``variant``, routing rsIDs and star alleles appropriately.

        ``threshold`` is accepted for API compatibility; matching is delegated
        to the PharmGKB symbol lookup.
        """
        variant = variant.strip()
        if not variant:
            return []
        if variant.lower().startswith("rs"):
            return self.rsid_lookup(variant, threshold=threshold, top_k=top_k)
        return self.star_lookup(variant, threshold=threshold, top_k=top_k)
