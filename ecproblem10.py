def matrix_transpose(matrix):
    return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]

def matrix_multiply(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]

def matrix_identity(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]

def matrix_norm(matrix):
    return sum(matrix[i][j]**2 for i in range(len(matrix)) for j in range(len(matrix[0])))**0.5

def scalar_multiply(matrix, scalar):
    return [[matrix[i][j] * scalar for j in range(len(matrix[0]))] for i in range(len(matrix))]

def polar_decomposition(matrix):
    transpose = matrix_transpose(matrix)
    sym_part = matrix_multiply(transpose, matrix)
    U = [[1 if i == j else 0 for j in range(len(matrix))] for i in range(len(matrix))]
    for _ in range(10):  # Iterative refinement
        U = scalar_multiply(matrix_multiply(U, sym_part), 0.5)
    P = matrix_multiply(U, transpose)
    return U, P

def cholesky_decomposition(matrix):
    n = len(matrix)
    L = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            if i == j:
                L[i][j] = (matrix[i][i] - sum(L[i][k]**2 for k in range(j)))**0.5
            else:
                L[i][j] = (matrix[i][j] - sum(L[i][k] * L[j][k] for k in range(j))) / L[j][j]
    return L

def svd(matrix):
    transpose = matrix_transpose(matrix)
    AAT = matrix_multiply(matrix, transpose)
    ATA = matrix_multiply(transpose, matrix)

    U, eigen_vals_AAT = eigen_decomposition(AAT)
    V, eigen_vals_ATA = eigen_decomposition(ATA)

    singular_values = [eigen_val**0.5 for eigen_val in eigen_vals_AAT]

    return U, singular_values, matrix_transpose(V)


def eigen_decomposition(matrix):
    n = len(matrix)
    eigen_vals = [1] * n
    eigen_vecs = matrix_identity(n)
    for _ in range(100):
        for i in range(n):
            for j in range(n):
                eigen_vecs[i][j] = eigen_vals[j]
        eigen_vals = [sum(matrix[i][j] * eigen_vecs[j][i] for j in range(n)) for i in range(n)]
    return eigen_vecs, eigen_vals

def eigenvalues(matrix):
    coeffs = characteristic_polynomial(matrix)
    return abert_method(coeffs)

def abert_method(coeffs):
    n = len(coeffs) - 1
    roots = [complex(1, 1) for _ in range(n)]
    for _ in range(100):
        for i in range(n):
            f = sum(c * roots[i]**(n-j) for j, c in enumerate(coeffs))
            df = sum((n-j) * c * roots[i]**(n-j-1) for j, c in enumerate(coeffs[:-1]))
            roots[i] -= f / df
    return roots
def test_decompositions():
    matrix = [
        [2, 1],
        [1, 2]
    ]

    print("\n--- Polar Decomposition ---")
    U, P = polar_decomposition(matrix)
    print("U (Unitary Matrix):", U)
    print("P (Positive Semi-Definite Matrix):", P)

    print("\n--- Cholesky Decomposition ---")
    L = cholesky_decomposition(matrix)
    print("L (Lower Triangular Matrix):", L)

    print("\n--- Singular Value Decomposition ---")
    U, S, V = svd(matrix)
    print("U (Left Singular Vectors):", U)
    print("Singular Values:", S)
    print("V (Right Singular Vectors):", V)

test_decompositions()
