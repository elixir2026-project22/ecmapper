import pytest

from ecmapper.ec import is_complete, normalize_ec


@pytest.mark.parametrize(
    "raw, expected",
    [
        ("1.1.1.1", "1.1.1.1"),
        (" 1.14.11.13", "1.14.11.13"),
        ("EC:2.7.1.11", "2.7.1.11"),
        ("EC-2.3.1.122", "2.3.1.122"),
        ("ec 2.7.1.40", "2.7.1.40"),
        ("1.1.1.n1", "1.1.1.n1"),
        ("2.3.1.B43", "2.3.1.B43"),
        ("1.1.1.-", "1.1.1.-"),
        ("2.4.1", "2.4.1.-"),
        ("2.4", "2.4.-.-"),
        ("1.-.-.-", "1.-.-.-"),
        ("3.A.1.7.-", None),
        ("2.A.40.-.-", None),
        ("7.3.2.a", None),
        ("1.8.1.M1", None),
        ("6.3.4.16-", None),
        ("1.14.13..-", None),
        ("1.-.1.1", None),
        ("SPONTANEOUS", None),
        ("", None),
        ("1", None),
        ("R00200", None),
    ],
)
def test_normalize_ec(raw, expected):
    assert normalize_ec(raw) == expected


def test_is_complete():
    assert is_complete("1.1.1.1")
    assert is_complete("1.1.1.n1")
    assert not is_complete("1.1.1.-")
