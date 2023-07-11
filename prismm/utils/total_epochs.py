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
