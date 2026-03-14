import pytest
from clinpgx_lookup import normalize


def test_normalize_drug():
    results = normalize("aspirin", entity_type="drug")
    assert len(results) >= 1
    assert results[0].name.lower() == "aspirin"


def test_normalize_variant_autodetect():
    results = normalize("rs1045642")
    assert len(results) >= 1
    # Should auto-detect as variant
    assert "rs" in results[0].name.lower() or results[0].id.startswith("PA")


def test_normalize_gene():
    results = normalize("CYP2D6", entity_type="gene")
    assert len(results) >= 1
    assert results[0].score >= 0.8


def test_normalize_invalid_entity_type():
    with pytest.raises(ValueError, match="Unknown entity_type"):
        normalize("aspirin", entity_type="invalid")


def test_normalize_phenotype():
    results = normalize("diabetes", entity_type="phenotype", threshold=0.5)
    assert isinstance(results, list)
