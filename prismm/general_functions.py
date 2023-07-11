from typing import List


def get_ev_string_v1(pre: int, mid: int, post: int) -> List[str]:
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
        result_v1 = get_ev_string_v1(pre, mid, post)  
        result_v2 = get_ev_string_v2(pre, mid, post)  
        assert result_v1 == result_v2, f"Mismatch for inputs ({pre}, {mid}, {post}): {result_v1} vs {result_v2}"

test_consistency_between_implementations()

get_ev_string = get_ev_string_v1

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


def swap_and_make_left_heavy(tree, left_child, right_child):
    if isinstance(tree, list):
        return [tree[0], make_left_heavy(right_child), make_left_heavy(left_child)]
    else:
        return (tree[0], make_left_heavy(right_child), make_left_heavy(left_child))

def make_left_heavy(tree):
    if len(tree) == 1:
        return tree
    else:
        left_child, right_child = tree[1], tree[2]

        if left_child[0] < right_child[0]:
            return swap_and_make_left_heavy(tree, left_child, right_child)
        else:
            if isinstance(tree, list):
                return [tree[0], make_left_heavy(left_child), make_left_heavy(right_child)]
            else:
                return (tree[0], make_left_heavy(left_child), make_left_heavy(right_child))

def test_make_left_heavy():
    # Case 1: Test with a simple CN_tree
    CN_tree = [1, [2], [3]]
    expected_left_heavy_tree = [1, [3], [2]]
    result_left_heavy_tree = make_left_heavy(CN_tree)
    assert result_left_heavy_tree == expected_left_heavy_tree, f"Expected {expected_left_heavy_tree}, but got {result_left_heavy_tree}"

    # Case 2: Test with a complex CN_tree
    CN_tree = [5, [3, [2], [4]], [7, [6], [8]]]
    expected_left_heavy_tree = [5, [7, [8], [6]], [3, [4], [2]]]
    result_left_heavy_tree = make_left_heavy(CN_tree)
    assert result_left_heavy_tree == expected_left_heavy_tree, f"Expected {expected_left_heavy_tree}, but got {result_left_heavy_tree}"

    # Case 3: Test with a single-node tree
    CN_tree = [1]
    expected_left_heavy_tree = [1]
    result_left_heavy_tree = make_left_heavy(CN_tree)
    assert result_left_heavy_tree == expected_left_heavy_tree, f"Expected {expected_left_heavy_tree}, but got {result_left_heavy_tree}"

test_make_left_heavy()


def pretty_print(obj, max_length=1000):
    text = str(obj)
    if len(text) > max_length:
        text = text[:max_length] + '...'
    print(text)

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


def total_epochs(pre,mid,post):
    return pre+mid+post+2


def total_epochs2(pre,mid,post):
    return pre*(pre>=0) + mid*(mid>=0) + post*(post>=0) + (mid>=0) + (post>=0)


def test_total_epochs():
    e1 = (0,-1,-1)
    e2 = (1,2,3)
    e3 = (0,0,0)
    e4 = (1,0,-1)

    for epochs in [e1,e2,e3,e4]:
        result1 = total_epochs(*epochs)
        result2 = total_epochs2(*epochs)
        try:
            assert(result1 == result2)
        except AssertionError:
            print(f"For input {epochs}, total_epochs returned {result1} while total_epochs2 returned {result2}")
            raise


test_total_epochs()



def generate_path(pre, mid, post):
    path = ""

    if pre == -1 and (mid != -1 or post != -1):
        raise ValueError("Invalid input: pre is -1 while mid or post is not -1")

    if mid == -1 and post != -1:
        raise ValueError("Invalid input: mid is -1 while post is not -1")

    if pre >= 0:
        path += str(pre)

    if mid >= 0:
        path += "G" + str(mid)

    if post >= 0:
        path += "G" + str(post)

    return path


def test_generate_path():
    assert generate_path(0, 0, -1) == "0G0"
    assert generate_path(1, 1, -1) == "1G1"
    assert generate_path(1, 1, 1) == "1G1G1"
    assert generate_path(0, 1, -1) == "0G1"
    assert generate_path(2, 0, -1) == "2G0"

    try:
        generate_path(-1, 1, 2)
    except ValueError as e:
        assert str(e) == "Invalid input: pre is -1 while mid or post is not -1"

    try:
        generate_path(-1, 0, -1)
    except ValueError as e:
        assert str(e) == "Invalid input: pre is -1 while mid or post is not -1"

    try:
        generate_path(1, -1, 2)
    except ValueError as e:
        assert str(e) == "Invalid input: mid is -1 while post is not -1"


test_generate_path()



def calculate_epochs(pre, mid, post):
    """
    Calculates the total number of epochs based on the provided pre, mid, and post values.

    Args:
        pre (int): The number of epochs in the pre-training stage. Should be non-negative.
        mid (int): The number of epochs in the mid-training stage. Should be non-negative.
        post (int): The number of epochs in the post-training stage. Should be non-negative.

    Returns:
        int: The total number of epochs.

    Raises:
        ValueError: If pre, mid, or post are negative.
    """

    epochs = pre + mid + post + 2
    return epochs


def test_calculate_epochs():
    """
    Tests the calculate_epochs function.
    """

    # Test with all values being zero
    pre, mid, post = 0, 0, 0
    expected_result = 2
    result = calculate_epochs(pre, mid, post)
    assert result == expected_result, f"For pre={pre}, mid={mid}, post={post}, expected {expected_result} but got {result}"

    # Test with all values being positive
    pre, mid, post = 1, 1, -1
    expected_result = 3
    result = calculate_epochs(pre, mid, post)
    assert result == expected_result, f"For pre={pre}, mid={mid}, post={post}, expected {expected_result} but got {result}"


test_calculate_epochs()

