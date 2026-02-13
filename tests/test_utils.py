import pytest
import numpy as np


np.random.seed(0)


def test_connectivity_matrix():

    from bionmf.utils import connectivity_matrix
    from collections import Counter
    from scipy.linalg import issymmetric


    X = np.random.randint(0, 10, size = 100)
    counts = Counter(X)
    conn = connectivity_matrix(X)

    # Assert matrix is symmetric
    assert(issymmetric(conn))
    assert(
        np.all(
            np.equal(
                # Row (or column) sum equals total number of occurrences of certain
                # cluster -> Counter(A)
                conn.sum(axis = 0),
                list(map(lambda x: counts[x], X))
            )
        )
    )


def test_sparse_rmse():

    from bionmf.utils import sparse_rmse
    import scipy.sparse as sp


    A = sp.random(100, 100, 0.5)
    B = sp.random(100, 100, 0.5)

    assert(sparse_rmse(A, A) == 0)
    assert(((A != B).nnz > 0) and sparse_rmse(A, B) > 0)
    assert(sparse_rmse(sp.coo_array([0]), sp.coo_array([1])) == 1)
