def generate_path_from_pre_mid_post(pre, mid, post):
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


def test_generate_path_from_pre_mid_post():
    assert generate_path_from_pre_mid_post(0, 0, -1) == "0G0"
    assert genergenerate_path_from_pre_mid_postate_path(1, 1, -1) == "1G1"
    assert generate_path_from_pre_mid_post(1, 1, 1) == "1G1G1"
    assert generate_path_from_pre_mid_post(0, 1, -1) == "0G1"
    assert generate_path_from_pre_mid_post(2, 0, -1) == "2G0"

    try:
        generate_path_from_pre_mid_post(-1, 1, 2)
    except ValueError as e:
        assert str(e) == "Invalid input: pre is -1 while mid or post is not -1"

    try:
        generate_path_from_pre_mid_post(-1, 0, -1)
    except ValueError as e:
        assert str(e) == "Invalid input: pre is -1 while mid or post is not -1"

    try:
        generate_path_from_pre_mid_post(1, -1, 2)
    except ValueError as e:
        assert str(e) == "Invalid input: mid is -1 while post is not -1"


test_generate_path_from_pre_mid_post()
