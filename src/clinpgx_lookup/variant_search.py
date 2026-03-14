from pydantic import BaseModel
from typing import List, Optional
import requests
from clinpgx_lookup.search_utils import (
    calc_similarity,
    general_search,
    general_search_comma_list,
)
from clinpgx_lookup.data import get_tsv_path
import pandas as pd
import logging

logger = logging.getLogger(__name__)


class VariantSearchResult(BaseModel):
    raw_input: str
    id: str
    name: str
    url: str
    score: float
    source: str = "clinpgx"


def pgkb_star_allele_search(
    star_allele: str, threshold: float = 0.8, top_k: int = 1
) -> List[VariantSearchResult]:
    base_url = "https://api.pharmgkb.org/v1/data/haplotype?symbol="
    try:
        response = requests.get(base_url + star_allele, timeout=10)
    except requests.RequestException:
        logger.warning("PharmGKB API request failed for '%s'", star_allele)
        return []
    if response.status_code == 200:
        data = response.json()
        if data.get("data"):
            score = calc_similarity(star_allele, data["data"][0]["symbol"])
            return [
                VariantSearchResult(
                    raw_input=star_allele,
                    id=result["id"],
                    name=result["symbol"],
                    url=f"https://www.clinpgx.org/haplotype/{result['id']}",
                    score=score,
                    source="pharmgkb",
                )
                for result in data["data"][:top_k]
            ]
    return []


def pgkb_rsid_search(
    rsid: str, threshold: float = 0.8, top_k: int = 1
) -> List[VariantSearchResult]:
    base_url = "https://api.pharmgkb.org/v1/data/variant?symbol="
    try:
        response = requests.get(base_url + rsid.strip(), timeout=10)
    except requests.RequestException:
        logger.warning("PharmGKB API request failed for '%s'", rsid)
        return []
    if response.status_code == 200:
        data = response.json()
        if data.get("data"):
            score = calc_similarity(rsid, data["data"][0]["symbol"])
            return [
                VariantSearchResult(
                    raw_input=rsid,
                    id=result["id"],
                    name=result["symbol"],
                    url=f"https://www.clinpgx.org/variant/{result['id']}",
                    score=score,
                    source="pharmgkb",
                )
                for result in data["data"][:top_k]
            ]
    return []


def ncbi_rsid_search(rsid: str, top_k: int = 1) -> List[VariantSearchResult]:
    """Look up an rsID via NCBI dbSNP Variation Services API."""
    num = rsid.strip().lower().replace("rs", "")
    url = f"https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{num}"
    try:
        response = requests.get(url, timeout=10)
    except requests.RequestException:
        logger.warning("NCBI API request failed for '%s'", rsid)
        return []
    if response.status_code != 200:
        return []
    data = response.json()
    refsnp_id = data.get("refsnp_id", "")
    if not refsnp_id:
        return []
    canonical_rsid = f"rs{refsnp_id}"
    # Handle merged SNPs — follow the merge to the current rsID
    if "merged_snapshot_data" in data:
        merged_into = data["merged_snapshot_data"].get("merged_into", [])
        if merged_into:
            canonical_rsid = f"rs{merged_into[0]}"
    return [
        VariantSearchResult(
            raw_input=rsid,
            id=f"dbSNP:{canonical_rsid}",
            name=canonical_rsid,
            url=f"https://www.ncbi.nlm.nih.gov/snp/{canonical_rsid}",
            score=1.0 if canonical_rsid == rsid.strip().lower() else 0.9,
            source="dbsnp",
        )
    ][:top_k]


class VariantLookup:
    def __init__(self, data_path: Optional[str] = None) -> None:
        self.data_path = data_path or str(get_tsv_path("variants"))
        self._df: Optional[pd.DataFrame] = None

    @property
    def df(self) -> pd.DataFrame:
        if self._df is None:
            self._df = pd.read_csv(self.data_path, sep="\t")
        return self._df

    def _clinpgx_variant_search(
        self, variant: str, threshold: float = 0.8, top_k: int = 1
    ) -> List[VariantSearchResult]:
        """Search Variant Name and Synonyms columns."""
        results = general_search(
            self.df,
            variant,
            "Variant Name",
            "Variant ID",
            threshold=threshold,
            top_k=top_k,
        )
        results.extend(
            general_search_comma_list(
                self.df,
                variant,
                "Synonyms",
                "Variant ID",
                threshold=threshold,
                top_k=top_k,
            )
        )
        results.sort(key=lambda x: x["score"], reverse=True)
        if results:
            return [
                VariantSearchResult(
                    raw_input=variant,
                    id=result["Variant ID"],
                    name=result["Variant Name"],
                    url=f"https://www.clinpgx.org/variant/{result['Variant ID']}",
                    score=result["score"],
                )
                for result in results[:top_k]
            ]
        return []

    def star_lookup(
        self, star_allele: str, threshold: float = 0.8, top_k: int = 1
    ) -> List[VariantSearchResult]:
        """Search for star alleles (e.g., CYP2C19*2)."""
        results = pgkb_star_allele_search(star_allele, threshold=threshold, top_k=top_k)
        results.extend(
            self._clinpgx_variant_search(star_allele, threshold=threshold, top_k=top_k)
        )
        results.sort(key=lambda x: x.score, reverse=True)
        return results[:top_k] if results else []

    def rsid_lookup(
        self, rsid: str, threshold: float = 0.8, top_k: int = 1
    ) -> List[VariantSearchResult]:
        """Search for rsIDs (e.g., rs1000002).

        Fallback chain: local ClinPGx (exact) → PharmGKB API → NCBI dbSNP API.
        """
        # 1. Local ClinPGx data — exact match only
        results = self._clinpgx_variant_search(rsid, threshold=1.0, top_k=top_k)
        if not results:
            # 2. PharmGKB API
            results = pgkb_rsid_search(rsid, threshold=threshold, top_k=top_k)
        if not results:
            # 3. NCBI dbSNP — validates the rsID exists
            results = ncbi_rsid_search(rsid, top_k=top_k)
        return results[:top_k] if results else []

    def search(
        self, variant: str, threshold: float = 0.8, top_k: int = 1
    ) -> List[VariantSearchResult]:
        """Search for a variant. Auto-detects rsID vs star allele."""
        if variant.strip().startswith("rs"):
            return self.rsid_lookup(variant, threshold=threshold, top_k=top_k)
        else:
            return self.star_lookup(variant, threshold=threshold, top_k=top_k)
