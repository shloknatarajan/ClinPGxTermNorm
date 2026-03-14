from clinpgx_lookup import GeneLookup, GeneSearchResult


def test_gene_lookup_by_symbol():
    lookup = GeneLookup()
    results = lookup.search("CYP2D6", threshold=0.8, top_k=1)
    assert len(results) >= 1
    assert isinstance(results[0], GeneSearchResult)
    assert results[0].symbol == "CYP2D6"


def test_gene_lookup_no_match():
    lookup = GeneLookup()
    results = lookup.search("XYZNOTGENE999", threshold=0.95, top_k=1)
    assert results == []
