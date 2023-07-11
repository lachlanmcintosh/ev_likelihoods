
def test_generate_trees():
    # Test case 1
    observed_CNs_1 = [3, 2]
    SNV_CNs_1 = [3, 2]
    expected_trees_1 = [(5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,)))]
    expected_trees_1 = [sort_tree(tree) for tree in expected_trees_1]
    result_trees_1 = generate_trees(observed_CNs_1, SNV_CNs_1)
    assert forests_are_equal(result_trees_1, expected_trees_1), f"Expected {expected_trees_1}, but got {result_trees_1}"

    # Test case 2
    observed_CNs_2 = [4, 1]
    SNV_CNs_2 = [2]
    expected_trees_2 = [(5, (4, (2, (1,), (1,)), (2, (1,), (1,))), (1,))]
    expected_trees_2 = [sort_tree(tree) for tree in expected_trees_2]
    result_trees_2 = generate_trees(observed_CNs_2, SNV_CNs_2)
    assert forests_are_equal(result_trees_2, expected_trees_2), f"Expected {expected_trees_2}, but got {result_trees_2}"

    observed_CNs_3 = [3, 0]
    SNV_CNs_3 = [3, 2, 1]
    result_trees_3 = generate_trees(observed_CNs_3, SNV_CNs_3)
    expected_trees_3 = [(3, (3, (2, (1,), (1,)), (1,)), (0,))]
    expected_trees_3 = [sort_tree(tree) for tree in expected_trees_3]
    assert forests_are_equal(result_trees_3, expected_trees_3), f"Expected {expected_trees_3}, but got {result_trees_3}"


def test_complete_trees():
    trees = [
        (2, (1,), (1,)),
        (3, (1,), (2,)),
        (4, (2,), (2,))
    ]
    
    completed_trees_result = complete_trees(trees)
    expected_result = [
        (2, (1,), (1,)),
        (3, (2, (1,), (1,)), (1,)),
        (4, (2, (1,), (1,)), (2, (1,), (1,)))
    ]
    
    assert completed_trees_result == expected_result, f"Expected {expected_result}, but got {completed_trees_result}"

    trees = [
        (6, (4,), (2,)), 
        (6, (4, (2,), (2,)), (2,))
    ]
    completed_trees_result = complete_trees(trees)
    expected_result = [
        (6, (4, (3, (2, (1,), (1,)), (1,)), (1,)), (2, (1,), (1,))), 
        (6, (4, (2, (1,), (1,)), (2, (1,), (1,))), (2, (1,), (1,)))
    ]

    assert completed_trees_result == expected_result, f"Expected {expected_result}, but got {completed_trees_result}"


def test_complete_tree():
    tree = (4,)
    completed_trees = complete_tree(tree)
    expected_trees = {
        (4, (3, (2, (1,), (1,)), (1,)), (1,)),
        (4, (2, (1,), (1,)), (2, (1,), (1,)))
    }
    assert completed_trees == expected_trees, f"Expected: {expected_trees}, Got: {completed_trees}"

    tree = (5,)
    completed_trees = complete_tree(tree)
    expected_trees = {
        (5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,))),
        (5, (4, (2, (1,), (1,)), (2, (1,), (1,))), (1,)),
        (5, (4, (3, (2, (1,), (1,)), (1,)), (1,)), (1,))
    }
    assert completed_trees == expected_trees, f"Expected: {expected_trees}, Got: {completed_trees}"

    tree = ((5, (3, (2,), (1,)), (2,)),)
    completed_trees = complete_tree(tree)
    expected_trees = {(5, (3, (2, (1,), (1,)), (1,)), (2, (1,), (1,)))}

    assert completed_trees == expected_trees, f"Expected: {expected_trees}, Got: {completed_trees}"


def test_insert_node():
    def sort_expected_output(trees):
        return [sort_tree(tree) for tree in trees]

    # Test Case 1: Inserting a node into an empty tree
    trees = [()]
    CN = 5
    result = insert_node(trees, CN)
    expected_output = [(5,)]
    assert result == sort_expected_output(expected_output), f"Expected {expected_output}, but got {result}"

    # Test Case 2: Inserting a node into a tree with only one node
    trees = [(10,)]
    CN = 3
    result = insert_node(trees, CN)
    expected_output = [(10, (7,), (3,))]
    assert result == sort_expected_output(expected_output), f"Expected {expected_output}, but got {result}"

    # Test Case 3: Inserting a node into a tree with multiple nodes
    trees = [(12, (5,), (7,))]
    CN = 2
    result = insert_node(trees, CN)
    expected_output = [
                        (12, (5, (3,), (2,)), (7,)),
                        (12, (5,), (7, (5,), (2,)))
                      ]
    assert result == sort_expected_output(expected_output), f"Expected {expected_output}, but got {result}"

    trees = [(7,)]
    CN = 10
    result = insert_node(trees, CN)
    expected_output = []
    assert result == sort_expected_output(expected_output), f"Expected {expected_output}, but got {result}"

    # Test Case 6: Inserting a node with a CN that's equal to the root node
    trees = [(5,)]
    CN = 5
    result = insert_node(trees, CN)
    expected_output = []
    assert result == sort_expected_output(expected_output), f"Expected {expected_output}, but got {result}"



def test_compare_trees_and_forests_are_equal():
    # Test cases for compare_trees
    tree1 = (5, (3,), (2, (1,), (1,)))
    tree2 = (5, (3,), (2, (1,), (1,)))
    tree3 = (5, (2, (1,), (1,)), (3,))
    tree4 = (5, (4,), (1, (1,), (0,)))
    tree5 = (5, (1, (1,), (0,)), (4,))

    assert compare_trees(tree1, tree2), "Expected tree1 and tree2 to be equal"
    assert compare_trees(tree1, tree3), "Expected tree1 and tree3 to be not equal"
    assert compare_trees(tree4, tree5), "Expected tree4 and tree5 to be not equal"

    # Test cases for forests_are_equal
    trees1 = [tree1, tree3]
    trees2 = [tree2, tree3]
    trees3 = [tree1, tree4]
    trees4 = [tree1, tree2, tree3, tree4]
    trees5 = [tree2, tree3, tree1, tree4]

    assert forests_are_equal(trees1, trees2), "Expected trees1 and trees2 to be equal"
    assert forests_are_equal(trees1, trees3), "Expected trees1 and trees3 to be not equal"
    assert forests_are_equal(trees4, trees5), "Expected trees4 and trees5 to be equal"


# Call the test function
test_compare_trees_and_forests_are_equal()

def test_sort_tree():
    # Test case 1: Empty tree
    tree = None
    assert sort_tree(tree) == None

    # Test case 2: Tree with a single value
    tree = (5,)
    assert sort_tree(tree) == (5,)

    # Test case 3: Tree with two values
    tree = (4, (2,))
    assert sort_tree(tree) == (4, (2,))

    # Test case 4: Tree with three values
    tree = (3, (5,), (2,))
    assert sort_tree(tree) == (3, (5,), (2,))

    # Test case 5: Tree with nested values
    tree = (7, (4, (8,), (6,)), (9,))
    assert sort_tree(tree) == (7, (9,), (4, (8,), (6,)))

test_sort_tree()