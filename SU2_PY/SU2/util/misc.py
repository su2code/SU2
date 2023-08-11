import numpy as np


def check_array(A, oned_as="row"):
    """ensures A is an array and at least of rank 2"""
    if not isinstance(A, np.ndarray):
        A = np.array(A)
    if np.rank(A) < 2:
        A = np.array(np.matrix(A))
        if oned_as == "row":
            pass
        elif oned_as == "col":
            A = A.T
        else:
            raise Exception("oned_as must be 'row' or 'col' ")

    return A
