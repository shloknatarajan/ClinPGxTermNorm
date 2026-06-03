"""Drug/chemical lookup against the PharmGKB API, with RxNorm for fuzzy matching.

Flow:
1. Exact PharmGKB chemical lookup by name.
2. On a miss, fuzzy-match the input with RxNorm's ``approximateTerm`` endpoint,
   resolve the result to its ingredient (so brand names map to the ingredient
   PharmGKB indexes), and re-query PharmGKB by that name.
3. On a miss, return the RxNorm result itself as a fallback.
"""

from typing import List, Optional, Tuple

import requests
from loguru import logger
from pydantic import BaseModel

from .search_utils import calc_similarity

PGKB_BASE = "https://api.pharmgkb.org/v1/data"
RXNAV_BASE = "https://rxnav.nlm.nih.gov/REST"
DEFAULT_TIMEOUT = 10


class DrugSearchResult(BaseModel):
    raw_input: str
    id: str
    name: str
    url: str
    score: float
    source: str = "pharmgkb"


""" PharmGKB helpers """


def _pgkb_chemical_by_name(name: str) -> Optional[dict]:
    """Return the first PharmGKB chemical matching ``name`` exactly, or None."""
    try:
        response = requests.get(
            f"{PGKB_BASE}/chemical",
            params={"name": name, "view": "base"},
            timeout=DEFAULT_TIMEOUT,
        )
    except requests.RequestException as exc:
        logger.warning(f"PharmGKB chemical request failed: {exc}")
        return None
    if response.status_code != 200:
        return None
    data = response.json().get("data") or []
    return data[0] if data else None


""" RxNorm helpers """


def _rxnorm_approximate(term: str) -> Optional[dict]:
    """Return the top RxNorm ``approximateTerm`` candidate for ``term``, or None."""
    try:
        response = requests.get(
            f"{RXNAV_BASE}/approximateTerm.json",
            params={"term": term, "maxEntries": 1},
            timeout=DEFAULT_TIMEOUT,
        )
    except requests.RequestException as exc:
        logger.warning(f"RxNorm approximateTerm request failed: {exc}")
        return None
    if response.status_code != 200:
        return None
    candidates = response.json().get("approximateGroup", {}).get("candidate", [])
    return candidates[0] if candidates else None


def _rxnorm_ingredient(rxcui: str) -> Optional[Tuple[str, str]]:
    """Resolve an rxcui to its ingredient as ``(name, rxcui)``.

    Brand names resolve to their ingredient; ingredients resolve to themselves.
    """
    try:
        response = requests.get(
            f"{RXNAV_BASE}/rxcui/{rxcui}/related.json",
            params={"tty": "IN"},
            timeout=DEFAULT_TIMEOUT,
        )
    except requests.RequestException as exc:
        logger.warning(f"RxNorm related request failed: {exc}")
        return None
    if response.status_code != 200:
        return None
    for group in response.json().get("relatedGroup", {}).get("conceptGroup", []):
        if group.get("tty") == "IN":
            props = group.get("conceptProperties", [])
            if props:
                return props[0]["name"], props[0]["rxcui"]
    return None


class DrugLookup:
    """Look up drugs/chemicals in PharmGKB, falling back to RxNorm fuzzy search."""

    def _pgkb_result(self, raw_input: str, chemical: dict) -> DrugSearchResult:
        return DrugSearchResult(
            raw_input=raw_input,
            id=chemical["id"],
            name=chemical["name"],
            url=f"https://www.clinpgx.org/chemical/{chemical['id']}",
            score=calc_similarity(raw_input, chemical["name"]),
            source="pharmgkb",
        )

    def search(
        self, drug_name: str, threshold: float = 0.8, top_k: int = 1
    ) -> List[DrugSearchResult]:
        """Look up ``drug_name`` and return matching :class:`DrugSearchResult` records.

        ``threshold`` and ``top_k`` are accepted for API compatibility; lookups
        are resolved through exact PharmGKB matches and RxNorm fuzzy matching,
        so at most one best match is returned.
        """
        drug_name = drug_name.strip()
        if not drug_name:
            return []

        # 1. Exact PharmGKB name match.
        chemical = _pgkb_chemical_by_name(drug_name)
        if chemical:
            return [self._pgkb_result(drug_name, chemical)][:top_k]

        # 2. Fuzzy-match via RxNorm, then re-check PharmGKB.
        candidate = _rxnorm_approximate(drug_name)
        if not candidate:
            logger.warning(f"No RxNorm candidate found for '{drug_name}'")
            return []
        rxcui = candidate.get("rxcui")

        ingredient = _rxnorm_ingredient(rxcui) if rxcui else None
        candidate_names = []
        if ingredient:
            candidate_names.append(ingredient[0])
        if candidate.get("name"):
            candidate_names.append(candidate["name"])

        for name in candidate_names:
            chemical = _pgkb_chemical_by_name(name)
            if chemical:
                return [self._pgkb_result(drug_name, chemical)][:top_k]

        # 3. Fall back to the RxNorm result itself.
        rx_name = ingredient[0] if ingredient else candidate.get("name", "")
        if rxcui:
            logger.warning(
                f"'{drug_name}' resolved via RxNorm but not found in PharmGKB"
            )
            return [
                DrugSearchResult(
                    raw_input=drug_name,
                    id=f"RXN{rxcui}",
                    name=rx_name or "",
                    url=f"{RXNAV_BASE}/rxcui/{rxcui}",
                    score=calc_similarity(drug_name, rx_name) if rx_name else 0.0,
                    source="rxnorm",
                )
            ]
        return []
