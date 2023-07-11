from utils import make_left_heavy

##### STEP 4; now we have SNV counts, make all possible trees that could explain those SNV counts for the given epoch structure “(pre,mid, post)”

##### To save redundancy and speed up computational time-complexity of discovering and sorting through all evolutionary 
##### trees that explain the SNV data we do not insert copy number 1’s when finding all tree structures too iterate over.
##### They can be considered implicit and inserted everywhere at the very end. 
##### Non unitary copy numbers are not necessarily 
##### The reason for this is that there is a bijection between tree structures that include copy number 1’s and 
##### have paired timeline arrays where it every copy number must have existed for non zero time except for copy number 1
##### and the tree structures that do not include copy number 1 as every leaf of every tree and force all nodes to have 
##### a non zero evolutionary history where SNVs were allowed to accumulate. The latter is easier to computationally discover.
##### The tree structures can be simply constructed in a tuple of tuples format like tree = (value, left subtree, right subtree). 

### some basic tree functions:
def sort_tree(tree):
    """
    Sorts a tree in ascending order.

    Args:
        tree (tuple or None): The tree to be sorted.

    Returns:
        tuple or None: The sorted tree.

    """
    if tree is None:
        return None

    if len(tree) == 1:
        return tree

    if len(tree) == 2:
        value, subtree = tree
        return (value, sort_tree(subtree))

    if len(tree) == 3:
        value, left_subtree, right_subtree = tree
        sorted_left = sort_tree(left_subtree)
        sorted_right = sort_tree(right_subtree)

        if sorted_left[0] > sorted_right[0]:
            return (value, sorted_left, sorted_right)
        else:
            return (value, sorted_right, sorted_left)




def compare_trees(tree1, tree2):
    tree1 = sort_tree(tree1)
    tree2 = sort_tree(tree2)
    len_tree1 = len(tree1)
    len_tree2 = len(tree2)

    if tree1 is None and tree2 is None:
        return True
    if tree1 is None or tree2 is None:
        return False

    # If the lengths of the trees are different, they are not equal
    if len_tree1 != len_tree2:
        return False

    # If the lengths of the trees are both 1, compare their values
    if len_tree1 == 1 and len_tree2 == 1:
        return tree1[0] == tree2[0]

    # If the lengths of the trees are both 2, compare their values
    if len_tree1 == 2 and len_tree2 == 2:
        if tree1[0] != tree2[0]:
            return False
        return compare_trees(tree1[1], tree2[1])

    # If the lengths of the trees are both 3, compare their values
    if len_tree1 == 3 and len_tree2 == 3:
        if tree1[0] != tree2[0]:
            return False
        return compare_trees(tree1[1], tree2[1]) and compare_trees(tree1[2], tree2[2])

    # If the code reaches this point, the trees are not equal
    return False

def forests_are_equal(trees1, trees2):
    if len(trees1) != len(trees2):
        return False
    for tree1 in trees1:
        found_match = False
        for tree2 in trees2:
            if compare_trees(tree1, tree2):
                found_match = True
                break
        if not found_match:
            return False
    return True


def insert_node(trees, CN):
    """
    Insert a node with the given copy number (CN) into every possible location in the input trees.

    :param trees: A list of binary trees.
    :param CN: The copy number to be inserted into the trees.
    :return: A list of new trees with the CN inserted once into each input tree.
    """

    # base case
    if trees == [] or trees == [()]:
        return [(CN,)]
    new_trees = []

    # otherwise
    for tree in trees:
        if len(tree) == 1 and CN < tree[0]:
            new_CNs = (CN, tree[0] - CN)
            new_trees.append(sort_tree((tree[0], (max(new_CNs),), (min(new_CNs),))))
        elif len(tree) == 3:
            for subtree in insert_node([tree[1]], CN):
                new_trees.append(sort_tree((tree[0], subtree, tree[2])))

            for subtree in insert_node([tree[2]], CN):
                new_trees.append(sort_tree((tree[0], tree[1], subtree)))

    return new_trees

def complete_tree(tree):
    """
    Given a tree structure as tuple, completes the tree by adding nodes and
    making it left-heavy. Returns a set of completed tree structures.

    :param tree: A tuple representing a tree structure with 1 to 3 elements.
    :return: A set of completed tree structures.
    :raises ValueError: If the input tree tuple has a length other than 1, 2, or 3.
    """
    def _complete_tree_three_elements(tree):
        left_trees = complete_tree(tree[1])
        right_trees = complete_tree(tree[2])

        completed_trees = set()
        for left_tree in left_trees:
            for right_tree in right_trees:
                completed_trees.add(make_left_heavy((tree[0], left_tree, right_tree)))
        return completed_trees


    def _complete_tree_two_elements(tree):
        left_trees = complete_tree(tree[1])
        completed_trees = set()

        for left_tree in left_trees:
            for i in range(1, tree[0] - left_tree[0]):
                right_trees = complete_tree((tree[0] - left_tree[0] - i,))
                for right_tree in right_trees:
                    completed_trees.add(make_left_heavy((tree[0], left_tree, right_tree)))
        return completed_trees


    def _complete_tree_one_element(tree):
        if tree[0] == 0 or tree[0] == 1:
            return {tree}

        completed_trees = set()
        for i in range(1, tree[0]):
            left_trees = complete_tree((i,))
            right_trees = complete_tree((tree[0] - i,))

            for left_tree in left_trees:
                for right_tree in right_trees:
                    completed_trees.add(make_left_heavy((tree[0], left_tree, right_tree)))
        return completed_trees

    if len(tree) == 3:
        return _complete_tree_three_elements(tree)
    elif len(tree) == 2:
        return _complete_tree_two_elements(tree)
    elif len(tree) == 1:
        return _complete_tree_one_element(tree)
    else:
        raise ValueError("Invalid tree structure. Expected tuple length between 1 and 3.")


def complete_trees(trees):
    """
    Given a list of trees, returns a new list with each tree completed.
    
    Args:
    trees (list): A list of trees, where each tree is a tuple.
    
    Returns:
    list: A list of completed trees.
    """
    result = set()

    for tree in trees:
        result = result.union(complete_tree(tree))

    return list(result)



def generate_trees(observed_CNs, SNV_CNs):
    """
    Generate all possible trees of every copy number from the multiplicity counts of the chromosomes and the SNVs.

    :param observed_CNs: A list of observed copy numbers.
    :param SNV_CNs: A list of copy numbers for SNVs.
    :return: A list of unique trees.
    """

    SNV_CNs.sort(reverse=True)
    observed_CNs.sort(reverse=True)

    # Initialize trees with the following tree for each chromosome
    trees = [(sum(observed_CNs), (max(observed_CNs),), (min(observed_CNs),))]

    for SNV_CN in SNV_CNs:
        if SNV_CN == 1:
            continue

        trees_with_new_node = insert_node(trees, SNV_CN)

        if not trees_with_new_node:
            assert SNV_CN in observed_CNs
            continue

        if SNV_CN in observed_CNs:
            trees += trees_with_new_node
        else:
            trees = trees_with_new_node

        while True:
            trees_with_node_inserted_again = insert_node(trees_with_new_node, SNV_CN)
            if not trees_with_node_inserted_again:
                break

            trees += trees_with_node_inserted_again
            trees_with_new_node = trees_with_node_inserted_again

    # Insert the "leaf nodes" into the tree - all of which are of CN 1
    trees = complete_trees(trees)

    trees = [sort_tree(tree) for tree in trees] 
    return trees


