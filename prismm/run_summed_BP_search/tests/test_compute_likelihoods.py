
def test_likelihoods_to_marginal_likelihoods():
    likelihoods = pd.DataFrame({
        "p_up": [0.1, 0.2, 0.3, 0.4],
        "p_down": [0.3, 0.2, 0.1, 0.0],
        "path": ["path1", "path2", "path3", "path4"],
        "likelihood": [0.25, 0.25, 0.25, 0.25]
    })

    max_number_of_solutions = 2
    default_paths = ["path1", "path4"]

    expected_result = pd.DataFrame({
        "p_up": [0.0, 0.0, 1.0, 4.0],
        "p_down": [1.0, 0.0, 3.0, 0.0],
        "path": ["path4", "path1", "path2", "path3"],
        "likelihood": [0.25, 0.25, 0.25, 0.25]
    })

    result = likelihoods_to_marginal_likelihoods(likelihoods, max_number_of_solutions, default_paths)
    pd.testing.assert_frame_equal(result, expected_result)

    # Additional test cases can be added here with different likelihoods, max_number_of_solutions, and default_paths values


def test_likelihoods_to_best_likelihood_by_path():
    likelihoods = pd.DataFrame({
        "p_up": [0.1, 0.2, 0.3, 0.4],
        "p_down": [0.3, 0.2, 0.1, 0.0],
        "path": ["path1", "path2", "path3", "path4"],
        "likelihood": [0.1, 0.2, 0.3, 0.4]
    })

    max_number_of_solutions = 2
    default_paths = ["path1", "path4"]

    expected_result = pd.DataFrame({
        "p_up": [0.4, 0.3, 0.1, 0.2],
        "p_down": [0.0, 0.1, 0.3, 0.2],
        "path": ["path4", "path3", "path1", "path2"],
        "likelihood": [0.4, 0.3, 0.1, 0.2]
    })

    result = likelihoods_to_best_likelihood_by_path(likelihoods, max_number_of_solutions, default_paths)
    pd.testing.assert_frame_equal(result, expected_result)

    # Additional test cases can be added here with different likelihoods, max_number_of_solutions, and default_paths values


def test_filter_rows_by_path_length():
    data = {
        "p_up": [33, 2, 97, 1, 1, 82, 0, 0, 0, 0],
        "p_down": [33, 2, 2, 1, 1, 13, 0, 0, 0, 0],
        "path": ["0G0", "1G0", "1", "0G1", "2G0", "0G0", "3G3", "3G4", "4G3", "5G0"],
        "likelihood": [0.989199, 0.001347, 0.001347, 0.000524, 0.000523, 0.000192, 0.000192, 0.000192, 0.000192, 0.000192]
    }

    df = pd.DataFrame(data)
    x = 2
    filtered_df = filter_rows_by_path_length(df, x)

    expected_filtered_data = {
        "p_up": [33, 2, 97, 1, 1, 82, 0],
        "p_down": [33, 2, 2, 1, 1, 13, 0],
        "path": ["0G0", "1G0", "1", "0G1", "2G0", "0G0", "5G0"],
        "likelihood": [0.989199, 0.001347, 0.001347, 0.000524, 0.000523, 0.000192, 0.000192]
    }

    expected_filtered_df = pd.DataFrame(expected_filtered_data)

    # Reset the index of both DataFrames
    filtered_df.reset_index(drop=True, inplace=True)
    expected_filtered_df.reset_index(drop=True, inplace=True)

    # Compare the filtered DataFrame with the expected output
    pd.testing.assert_frame_equal(filtered_df, expected_filtered_df)

#test_filter_rows_by_path_length()


def test_filter_out_p_up_p_down_equals_zero_row():
    data = {
        "p_up": [33, 2, 97, 1, 1, 82, 0, 0, 0, 0],
        "p_down": [33, 2, 2, 1, 1, 13, 0, 0, 0, 0],
        "path": ["0G0", "1G0", "1", "0G1", "2G0", "0G0", "3G3", "3G4", "4G3", "5G0"],
        "likelihood": [0.989199, 0.001347, 0.001347, 0.000524, 0.000523, 0.000192, 0.000192, 0.000192, 0.000192, 0.000192]
    }

    df = pd.DataFrame(data)
    filtered_df = filter_out_p_up_p_down_equals_zero_row(df)

    expected_filtered_data = {
        "p_up": [33, 2, 97, 1, 1, 82, 0],
        "p_down": [33, 2, 2, 1, 1, 13, 0],
        "path": ["0G0", "1G0", "1", "0G1", "2G0", "0G0", "5G0"],
        "likelihood": [0.989199, 0.001347, 0.001347, 0.000524, 0.000523, 0.000192, 0.000192]
    }

    expected_filtered_df = pd.DataFrame(expected_filtered_data)

    # Reset the index of both DataFrames
    filtered_df.reset_index(drop=True, inplace=True)
    expected_filtered_df.reset_index(drop=True, inplace=True)

    # Compare the filtered DataFrame with the expected output
    try:
        pd.testing.assert_frame_equal(filtered_df, expected_filtered_df)
    except AssertionError as e:
        print("Expected:")
        print(expected_filtered_df)
        print("Actual:")
        print(filtered_df)
        raise e