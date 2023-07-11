def path_code_to_pre_mid_post(path):
    """
    Takes a path string or a list of path strings containing 'G' characters as separators and returns a tuple of integers
    representing the number of elements before the first 'G', the number of elements between the first and second 'G',
    and the number of elements after the second 'G'.

    Args:
        path (str or list of str): The path string or list of path strings.

    Returns:
        tuple: A tuple containing three integers (pre, mid, post).

    Raises:
        TypeError: If the path is neither a string nor a list of strings.
        ValueError: If the path contains more than two 'G' characters.
    """
    pre, mid, post = 0, -1, -1

    if isinstance(path, str):
        bits = [int(x) for x in path.split("G") if x] + [-1, -1]
        pre, mid, post = bits[:3]
    elif isinstance(path, list) and all(isinstance(x, str) for x in path):
        n_g = path.count('G')
        if n_g >= 3:
            raise ValueError(f"Expected at most 2 'G' characters in the path, received {n_g}: {path}")

        if n_g == 0:
            pre = len(path)
        elif n_g == 1:
            idx = path.index('G')
            pre, mid = idx, len(path) - idx - 1
        else:  # n_g == 2
            idx1, idx2 = [i for i, x in enumerate(path) if x == 'G']
            pre, mid, post = idx1, idx2 - idx1 - 1, len(path) - idx2 - 1
    else:
        raise TypeError(f"Expected a string or a list of strings, received {type(path)}: {path}")

    return (pre, mid, post)

def path_code_to_path_length(path):
    (pre,mid,post) = path_code_to_pre_mid_post(path)
    return pre+mid+post+2

def path_code_num_anueploidy_epochs(path):
    (pre, mid, post) = path_code_to_pre_mid_post(path)
    return pre*(pre>0) + mid*(mid>0) + post*(post>0)


def test_path_code_to_pre_mid_post():
    def assert_result(path, expected_pre_mid_post):
        result = path_code_to_pre_mid_post(path)
        assert result == expected_pre_mid_post, f"Expected {expected_pre_mid_post}, got {result} for path {path}"

    # Test case 1
    assert_result("2G1G0", (2, 1, 0))

    # Test case 2
    assert_result("1G0G2", (1, 0, 2))

    # Test case 3
    assert_result("0G1G1", (0, 1, 1))

    # Test case 4
    assert_result("1G1", (1, 1, -1))

    # Test case 5 - List input with no G's
    assert_result(['A', 'A', 'A'], (3, -1, -1))

    # Test case 6 - List input with one G
    assert_result(['A', 'G', 'A', 'A'], (1, 2, -1))

    # Test case 7 - List input with two G's
    assert_result(['A', 'A', 'G', 'A', 'G', 'A'], (2, 1, 1))



test_path_code_to_pre_mid_post()
