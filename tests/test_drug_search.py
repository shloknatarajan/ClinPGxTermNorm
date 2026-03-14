from clinpgx_lookup import DrugLookup, DrugSearchResult


def test_drug_lookup_aspirin():
    lookup = DrugLookup()
    results = lookup.search("aspirin", threshold=0.8, top_k=1)
    assert len(results) >= 1
    assert isinstance(results[0], DrugSearchResult)
    assert results[0].score >= 0.8
    assert results[0].name.lower() == "aspirin"


def test_drug_lookup_no_match():
    lookup = DrugLookup()
    results = lookup.clinpgx_lookup("xyznotadrug123", threshold=0.95, top_k=1)
    assert results == []


def test_drug_lookup_trade_name():
    lookup = DrugLookup()
    results = lookup._clinpgx_drug_alternatives_search(
        "Tylenol", threshold=0.8, top_k=1
    )
    assert len(results) >= 1
    assert isinstance(results[0], DrugSearchResult)
