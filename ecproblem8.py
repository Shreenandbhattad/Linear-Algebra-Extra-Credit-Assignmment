def dot_product(v1, v2):
    return sum(x * y for x, y in zip(v1, v2))

def is_orthogonal(v1, v2):
    return dot_product(v1, v2) == 0

def gram_schmidt(vectors):
    orthogonal_vectors = []
    for v in vectors:
        for u in orthogonal_vectors:
            proj = dot_product(v, u) / dot_product(u, u)
            v = [v_i - proj * u_i for v_i, u_i in zip(v, u)]
        orthogonal_vectors.append(v)
    return orthogonal_vectors

def qr_factorization(matrix):
    m, n = len(matrix), len(matrix[0])
    Q = []
    R = [[0] * n for _ in range(n)]
    for i in range(n):
        Q.append(matrix[i][:])
        for j in range(i):
            R[j][i] = dot_product(Q[j], Q[i])
            for k in range(m):
                Q[i][k] -= R[j][i] * Q[j][k]
        R[i][i] = (sum(x * x for x in Q[i])) ** 0.5
        for k in range(m):
            Q[i][k] /= R[i][i]
    return Q, R

def transpose(matrix):
    return list(zip(*matrix))

def matrix_multiply(A, B):
    return [
        [sum(a * b for a, b in zip(A_row, B_col)) for B_col in zip(*B)]
        for A_row in A
    ]

def moore_penrose_pseudoinverse(matrix):
    Q, R = qr_factorization(matrix)
    R_inv = invert_upper_triangle(R)
    return matrix_multiply(transpose(Q), R_inv)

def invert_upper_triangle(R):
    n = len(R)
    R_inv = [[0] * n for _ in range(n)]
    for i in range(n):
        R_inv[i][i] = 1 / R[i][i]
        for j in range(i - 1, -1, -1):
            R_inv[j][i] = -sum(R[j][k] * R_inv[k][i] for k in range(j + 1, n)) / R[j][j]
    return R_inv

def least_squares_solution(A, b):
    A_T = transpose(A)
    A_T_A = matrix_multiply(A_T, A)
    A_T_b = [sum(A_T[i][j] * b[j] for j in range(len(b))) for i in range(len(A_T))]
    R_inv = invert_upper_triangle(A_T_A)
    return [sum(R_inv[i][j] * A_T_b[j] for j in range(len(A_T_b))) for i in range(len(R_inv))]
