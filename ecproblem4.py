class Matrix:
    def __init__(self, n, m, field="real", rows=None):
        self.n = n
        self.m = m
        self.field = field.lower()
        self.rows = rows if rows else [[0] * m for _ in range(n)]

    def rref(self, show_operations=False):
        mat = [row[:] for row in self.rows]
        lead = 0
        row_count, column_count = len(mat), len(mat[0])
        for r in range(row_count):
            if lead >= column_count:
                return Matrix(self.n, self.m, self.field, mat)
            i = r
            while mat[i][lead] == 0:
                i += 1
                if i == row_count:
                    i = r
                    lead += 1
                    if column_count == lead:
                        return Matrix(self.n, self.m, self.field, mat)
            mat[i], mat[r] = mat[r], mat[i]
            lv = mat[r][lead]
            mat[r] = [mrx / lv for mrx in mat[r]]
            for i in range(row_count):
                if i != r:
                    lv = mat[i][lead]
                    mat[i] = [iv - lv * rv for rv, iv in zip(mat[r], mat[i])]
            lead += 1
        return Matrix(self.n, self.m, self.field, mat)

    def rank(self):
        rref_matrix = self.rref()
        return sum(any(cell != 0 for cell in row) for row in rref_matrix.rows)

    def size(self):
        return self.n, self.m

    def plu_decomposition(self):
        if self.n != self.m:
            raise ValueError("PLU decomposition is only valid for square matrices.")
        
        P = [[1 if i == j else 0 for j in range(self.n)] for i in range(self.n)]
        L = [[0 for _ in range(self.n)] for _ in range(self.n)]
        U = [row[:] for row in self.rows]  # Copy of the original matrix

        for i in range(self.n):
            max_row = max(range(i, self.n), key=lambda r: abs(U[r][i]))
            if i != max_row:
                U[i], U[max_row] = U[max_row], U[i]
                P[i], P[max_row] = P[max_row], P[i]
            L[i][i] = 1
            for j in range(i+1, self.n):
                if U[j][i] != 0:
                    L[j][i] = U[j][i] / U[i][i]
                    for k in range(i, self.n):
                        U[j][k] -= L[j][i] * U[i][k]

        return Matrix(self.n, self.n, self.field, P), Matrix(self.n, self.n, self.field, L), Matrix(self.n, self.n, self.field, U)

class LinearSystem:
    def __init__(self, A, b):
        self.A = A
        self.b = b
        if len(A) != len(b) or len(A[0]) != len(b):
            raise ValueError("Matrix A and vector b dimensions do not match.")

    def is_consistent(self):
        augmented_matrix = [self.A[i] + [self.b[i]] for i in range(len(self.A))]
        rref_matrix = Matrix(len(augmented_matrix), len(augmented_matrix[0]), "real", augmented_matrix).rref()
        for row in rref_matrix.rows:
            if all(x == 0 for x in row[:-1]) and row[-1] != 0:
                return False
        return True

    def gaussian_elimination(self):
        augmented_matrix = [self.A[i] + [self.b[i]] for i in range(len(self.A))]
        rref_matrix = Matrix(len(augmented_matrix), len(augmented_matrix[0]), "real", augmented_matrix).rref()
        return [row[-1] for row in rref_matrix.rows]

    def is_subspace(self, S1, S2):
        test_matrix = Matrix(len(S1), len(S1[0]), "real", S1)
        return test_matrix.rank() <= Matrix(len(S2), len(S2[0]), "real", S2).rank()

    def express_solution_set(self):
        augmented_matrix = [self.A[i] + [self.b[i]] for i in range(len(self.A))]
        rref_matrix = Matrix(len(augmented_matrix), len(augmented_matrix[0]), "real", augmented_matrix).rref()
        return rref_matrix

    def pl_decomposition_solution(self):
        P, L, U = Matrix(len(self.A), len(self.A[0]), "real", self.A).plu_decomposition()
        y = [0] * len(self.A)
        for i in range(len(self.A)):
            y[i] = self.b[i]
            for j in range(i):
                y[i] -= L.rows[i][j] * y[j]  
        x = [0] * len(self.A)
        for i in range(len(self.A)-1, -1, -1):
            x[i] = y[i]
            for j in range(i+1, len(self.A)):
                x[i] -= U.rows[i][j] * x[j]  
            x[i] /= U.rows[i][i] 
        return x

# Test cases
A = [[2, 1], [1, 3]]
b = [5, 8]
system = LinearSystem(A, b)

print("Is system consistent:", system.is_consistent()) 

print("Solution using Gaussian Elimination:", system.gaussian_elimination())  # Expected: [3, 2]

S1 = [[1, 0], [0, 1]]
S2 = [[1, 1], [1, 0]]
print("Is span of S1 a subspace of S2:", system.is_subspace(S1, S2))  # Expected: True

print("Solution set with free variables:")
solution_set = system.express_solution_set()
print(solution_set)

print("Solution using PLU decomposition:", system.pl_decomposition_solution())  # Expected: [3, 2]
