from general_functions import *

def complete_tree(tree):
    if len(tree) == 3:
        left_trees = complete_tree(tree[1])
        right_trees = complete_tree(tree[2])

        completed_trees = set()
        for left_tree in left_trees:
            for right_tree in right_trees:
                completed_trees.add(make_left_heavy((tree[0], left_tree, right_tree)))
        return completed_trees

    elif len(tree) == 2:
        left_trees = complete_tree(tree[1])
        completed_trees = set()

        for left_tree in left_trees:
            for i in range(1, tree[0] - left_tree[0]):
                right_trees = complete_tree((tree[0] - left_tree[0] - i,))
                for right_tree in right_trees:
                    completed_trees.add(make_left_heavy((tree[0], left_tree, right_tree)))
        return completed_trees

    elif len(tree) == 1:
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

    else:
        raise ValueError("Invalid tree structure. Expected tuple length between 1 and 3.")

    return None

tree = (4,)
completed_trees = complete_tree(tree)
for t in completed_trees:
    print(t)


tree = (5,)
completed_trees = complete_tree(tree)
for t in completed_trees:
    print(t)
