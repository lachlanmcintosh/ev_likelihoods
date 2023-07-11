from typing import List

def get_ev_string(pre: int, mid: int, post: int) -> List[str]:
    ev_str = []
    if pre > 0:
        ev_str += ["A"] * pre
    if mid > -1:
        ev_str += ["G"]
    if mid > 0:
        ev_str += ["A"] * mid
    if post > -1:
        ev_str += ["G"]
    if post > 0:
        ev_str += ["A"] * post
    return ev_str


def get_ev_string_v2(pre: int, mid: int, post: int) -> List[str]:
    assert isinstance(pre, int) and pre >= -1, "Invalid pre value"
    assert isinstance(mid, int) and mid >= -1, "Invalid mid value"
    assert isinstance(post, int) and post >= -1, "Invalid post value"

    if pre == -1:
        assert mid == -1
        assert post == -1
    if mid == -1:
        assert post == -1

    ev_str = []
    ev_str += ["A"] * pre
    ev_str += ["G"] * (mid > -1)
    ev_str += ["A"] * mid
    ev_str += ["G"] * (post > -1)
    ev_str += ["A"] * post
    return ev_str


def test_consistency_between_implementations():
    test_cases = [
        (2, 0, 0),
        (0, 1, 0),
        (0, 0, 2),
        (2, 1, 2),
        (-1, -1, -1),
        (2, -1, -1),
        (0, -1, -1),
    ]

    for pre, mid, post in test_cases:
        result_v1 = get_ev_string(pre, mid, post)
        result_v2 = get_ev_string_v2(pre, mid, post)
        assert result_v1 == result_v2, f"Mismatch for inputs ({pre}, {mid}, {post}): {result_v1} vs {result_v2}"

test_consistency_between_implementations()

def test_get_ev_string():
    # valid inputs
    assert get_ev_string(2, 0, 0) == ["A", "A", "G", "G"]
    assert get_ev_string(0, 1, 0) == ["G", "A", "G"]
    assert get_ev_string(0, 0, 2) == ["G", "G", "A", "A"]
    assert get_ev_string(2, 1, 2) == ["A", "A", "G", "A", "G", "A", "A"]

    # pre is -1
    assert get_ev_string(-1, -1, -1) == []

    # mid is -1
    assert get_ev_string(2, -1, -1) == ["A", "A"]
    assert get_ev_string(0, -1, -1) == []
    assert get_ev_string(0, -1, -1) == []
    assert get_ev_string(-1, -1, -1) == []

test_get_ev_string()
