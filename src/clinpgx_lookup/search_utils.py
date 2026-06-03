from difflib import SequenceMatcher


def calc_similarity(query: str, text: str) -> float:
    """Return a 0-1 fuzzy similarity ratio between two strings (case-insensitive)."""
    return SequenceMatcher(None, query.lower().strip(), text.lower().strip()).ratio()
