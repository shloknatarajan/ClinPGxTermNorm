from clinpgx_lookup import VariantLookup, VariantSearchResult


def test_variant_local_search():
    """Test local TSV search without hitting external APIs."""
    lookup = VariantLookup()
    results = lookup._clinpgx_variant_search("rs1045642", threshold=0.8, top_k=1)
    assert len(results) >= 1
    assert isinstance(results[0], VariantSearchResult)
    assert results[0].score >= 0.8
