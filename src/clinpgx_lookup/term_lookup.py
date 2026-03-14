from typing import List, Optional, Union
from clinpgx_lookup.drug_search import DrugLookup, DrugSearchResult
from clinpgx_lookup.variant_search import VariantLookup, VariantSearchResult
from clinpgx_lookup.gene_search import GeneLookup, GeneSearchResult
from clinpgx_lookup.phenotype_search import PhenotypeLookup, PhenotypeSearchResult
from clinpgx_lookup.chemical_search import ChemicalLookup, ChemicalSearchResult

SearchResult = Union[
    DrugSearchResult,
    VariantSearchResult,
    GeneSearchResult,
    PhenotypeSearchResult,
    ChemicalSearchResult,
]

ENTITY_TYPES = ("drug", "variant", "gene", "phenotype", "chemical")

# Lazily initialized singletons
_lookups: dict = {}


def _get_lookup(entity_type: str):
    if entity_type not in _lookups:
        if entity_type == "drug":
            _lookups[entity_type] = DrugLookup()
        elif entity_type == "variant":
            _lookups[entity_type] = VariantLookup()
        elif entity_type == "gene":
            _lookups[entity_type] = GeneLookup()
        elif entity_type == "phenotype":
            _lookups[entity_type] = PhenotypeLookup()
        elif entity_type == "chemical":
            _lookups[entity_type] = ChemicalLookup()
    return _lookups[entity_type]


def normalize(
    term: str,
    entity_type: Optional[str] = None,
    threshold: float = 0.8,
    top_k: int = 1,
) -> List[SearchResult]:
    """Normalize a clinical pharmacogenomics term to its best ClinPGx match.

    Args:
        term: The term to normalize (e.g., "aspirin", "rs1000002", "CYP2D6").
        entity_type: One of "drug", "variant", "gene", "phenotype", "chemical".
            If None, auto-detects based on the term format.
        threshold: Minimum similarity score (0-1) for fuzzy matching. Default 0.8.
        top_k: Number of top results to return. Default 1.

    Returns:
        List of search results, sorted by score descending.

    Examples:
        >>> from clinpgx_lookup import normalize
        >>> normalize("aspirin", entity_type="drug")
        >>> normalize("rs1000002")  # auto-detects as variant
        >>> normalize("CYP2D6", entity_type="gene")
    """
    if entity_type is not None:
        entity_type = entity_type.lower().strip()
        if entity_type not in ENTITY_TYPES:
            raise ValueError(
                f"Unknown entity_type '{entity_type}'. "
                f"Must be one of: {', '.join(ENTITY_TYPES)}"
            )
        lookup = _get_lookup(entity_type)
        return lookup.search(term, threshold=threshold, top_k=top_k)

    # Auto-detection
    stripped = term.strip()

    # rsIDs -> variant
    if stripped.lower().startswith("rs") and any(c.isdigit() for c in stripped):
        return _get_lookup("variant").search(stripped, threshold=threshold, top_k=top_k)

    # Star alleles (e.g., CYP2C19*2) -> variant
    if "*" in stripped:
        return _get_lookup("variant").search(stripped, threshold=threshold, top_k=top_k)

    # Default to drug (most common use case)
    return _get_lookup("drug").search(term, threshold=threshold, top_k=top_k)
