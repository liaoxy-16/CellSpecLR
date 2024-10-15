import numpy as np
import pandas as pd


def filter_CCC_matrix(matrix: np.ndarray, t: float, output_file: str):
    """
    Parameters:
    - matrix: Input 2D NumPy matrix
    - t: Threshold, where the element representing the top-t proportion by absolute value is set to 1
    - output_file: Path of the output CSV file

    Returns:
    - The processed matrix as a DataFrame
    """

    abs_matrix = np.abs(matrix)

    # Flatten the matrix and calculate the number of top-t elements based on absolute values
    flattened = abs_matrix.flatten()
    n_elements = len(flattened)
    top_t_percent_index = int(n_elements * t)

    # Sort the flattened matrix in descending order and get the indices of the top t% largest elements
    sorted_indices = np.argsort(flattened)[::-1]  # Sort in descending order
    top_indices = sorted_indices[:top_t_percent_index]

    # Create a new matrix with zeros
    new_matrix = np.zeros_like(abs_matrix)

    # Set the top t% largest elements to 1
    np.put(new_matrix, top_indices, 1)

    # Set the diagonal elements to 0
    np.fill_diagonal(new_matrix, 0)

    # Convert the matrix to a DataFrame and assign row and column indices
    row_index = np.arange(1, new_matrix.shape[0] + 1)
    col_index = np.arange(1, new_matrix.shape[1] + 1)
    new_matrix_df = pd.DataFrame(new_matrix, index=row_index, columns=col_index)

    # Save the processed matrix to the specified CSV file
    new_matrix_df.to_csv(output_file)

    return new_matrix_df


