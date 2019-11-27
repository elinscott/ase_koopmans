import re

_pattern_test_data = []


def _define_pattern(pattern, example, *args):
    """Accepts a regular expression pattern, as well as an example
    string that the pattern should match. Returns the compiled
    pattern. Additionally, stores the compiled pattern and the
    example string in pattern_test_data for unit testing."""
    regex = re.compile(pattern, *args)
    _pattern_test_data.append((regex, example))
    return regex
