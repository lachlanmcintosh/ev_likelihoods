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
