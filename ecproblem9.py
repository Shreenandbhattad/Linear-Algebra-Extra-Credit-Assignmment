def abert_method(coeffs):
    n = len(coeffs) - 1
    roots = [complex(1, 1) for _ in range(n)]  
    for _ in range(100):  
        for i in range(n):
            f = sum(c * roots[i]**(n-j) for j, c in enumerate(coeffs))  
            df = sum((n-j) * c * roots[i]**(n-j-1) for j, c in enumerate(coeffs[:-1]))  # Derivative
            roots[i] -= f / df  # Newton's method to update roots
    return roots

def characteristic_polynomial(matrix):
    n = len(matrix)
    coeffs = [1] + [0] * n
    for i in range(n):
        coeffs[i] = -sum(matrix[i][j] * matrix[j][i] for j in range(n))
    return coeffs

def minimal_polynomial(matrix):
    n = len(matrix)
    coeffs = [1] + [0] * n
    for i in range(n):
        coeffs[i] = -sum(matrix[i][j] * matrix[j][i] for j in range(n))
    return coeffs

def eigenvalues(matrix):
    coeffs = characteristic_polynomial(matrix)
    roots = abert_method(coeffs)
    return roots

def matrix_multiply(A, B):
    return [[sum(a * b for a, b in zip(A_row, B_col)) for B_col in zip(*B)] for A_row in A]

def is_similar(A, B):
    return eigenvalues(A) == eigenvalues(B)

def change_of_basis_matrix(A, B):
    if not is_similar(A, B):
        return None
    return matrix_multiply(B, inverse(A))

def eigen_bases(matrix, eigenvalue):
    n = len(matrix)
    eigen_vector = [0] * n
    for i in range(n):
        if matrix[i][i] == eigenvalue:
            eigen_vector[i] = 1
    return [eigen_vector]

def algebraic_multiplicity(matrix, eigenvalue):
    return sum(1 for eigenvalue in eigenvalues(matrix) if eigenvalue == eigenvalue)

def geometric_multiplicity(matrix, eigenvalue):
    return sum(1 for eigenvalue in eigenvalues(matrix) if eigenvalue == eigenvalue)

def inverse(matrix):
    n = len(matrix)
    augmented = [row + (1 if i == j else 0) for i, row in enumerate(matrix)]
    for i in range(n):
        if augmented[i][i] == 0:
            return None
        for j in range(i + 1, n):
            factor = augmented[j][i] / augmented[i][i]
            for k in range(2 * n):
                augmented[j][k] -= augmented[i][k] * factor
    for i in range(n - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            factor = augmented[j][i] / augmented[i][i]
            for k in range(2 * n):
                augmented[j][k] -= augmented[i][k] * factor
    for i in range(n):
        factor = augmented[i][i]
        for j in range(2 * n):
            augmented[i][j] /= factor
    return [row[n:] for row in augmented]

def test_eigen_methods():
    matrix = [
        [4, -1, 1],
        [2, 4, -2],
        [1, 1, 3]
    ]

    print("Eigenvalues:", eigenvalues(matrix)) 
    print("Characteristic Polynomial:", characteristic_polynomial(matrix))  
    print("Minimal Polynomial:", minimal_polynomial(matrix))  
    print("matrix similar to itself:", is_similar(matrix, matrix)) 
    print("Algebraic multiplicity of eigenvalue 1:", algebraic_multiplicity(matrix, 1)) 
    print("Geometric multiplicity of eigenvalue 1:", geometric_multiplicity(matrix, 1))  
    print("Eigen-basis for eigenvalue 1:", eigen_bases(matrix, 1))  

test_eigen_methods()
