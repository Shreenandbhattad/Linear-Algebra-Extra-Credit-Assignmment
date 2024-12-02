class Matrix:
    def __init__(self, n, m, field="real", rows=None):
        self.n = n
        self.m = m
        self.field = field.lower()
        self.rows = rows if rows else [[0] * m for _ in range(n)]

    def size(self):
        return self.n, self.m

    def length_of_vector(self, vector):
        return (sum(x * x for x in vector)) ** 0.5

    def rank(self):
        rref_matrix = self.rref()
        return sum(any(cell != 0 for cell in row) for row in rref_matrix.rows)

    def nullity(self):
        return self.m - self.rank()

    def rref(self, show_operations=False):
        mat = [row[:] for row in self.rows]  # Copy rows to avoid altering original
        lead = 0
        row_count, column_count = len(mat), len(mat[0])

        for r in range(row_count):
            if lead >= column_count:
                return Matrix(self.n, self.m, self.field, mat)  # Return matrix with rows
            i = r
            while mat[i][lead] == 0:
                i += 1
                if i == row_count:
                    i = r
                    lead += 1
                    if column_count == lead:
                        return Matrix(self.n, self.m, self.field, mat)  # Return matrix with rows
            mat[i], mat[r] = mat[r], mat[i]
            lv = mat[r][lead]
            mat[r] = [mrx / lv for mrx in mat[r]]
            for i in range(row_count):
                if i != r:
                    lv = mat[i][lead]
                    mat[i] = [iv - lv * rv for rv, iv in zip(mat[r], mat[i])]
            lead += 1

            if show_operations:
                print(f"After row operation {r + 1}: {mat}")

        return Matrix(self.n, self.m, self.field, mat)  # Return matrix with rows

    def linearly_independent(self, vectors):
        if not vectors or len(vectors[0]) != len(vectors):
            raise ValueError("Vectors must form a square matrix for independence check.")
        test_matrix = Matrix(len(vectors), len(vectors[0]), self.field, vectors)
        return test_matrix.rank() == len(vectors)

    def dimension_of_span(self, vectors):
        test_matrix = Matrix(len(vectors[0]), len(vectors), self.field, list(zip(*vectors)))
        return test_matrix.rank()

    def basis_for_span(self, vectors):
        rref_matrix = self.rref()
        basis = []
        for row in rref_matrix.rows:
            if any(cell != 0 for cell in row):
                basis.append(row)
        return basis

    def rank_factorization(self):
        rref_matrix = self.rref()
        rank = self.rank()
        U = [[1 if i == j else 0 for j in range(self.m)] for i in range(rank)]
        return Matrix(rank, self.m, self.field, rref_matrix.rows[:rank]), Matrix(rank, self.m, self.field, U)

    def lu_decomposition(self):
        if self.n != self.m:
            raise ValueError("LU decomposition only for square matrices.")
        L = [[0] * self.n for _ in range(self.n)]
        U = [[0] * self.n for _ in range(self.n)]
        for i in range(self.n):
            for j in range(i, self.n):
                U[i][j] = self.rows[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
            for j in range(i, self.n):
                if i == j:
                    L[i][i] = 1
                else:
                    L[j][i] = (self.rows[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]
        return Matrix(self.n, self.n, self.field, L), Matrix(self.n, self.n, self.field, U)

    def plu_decomposition(self):
        if self.n != self.m:
            raise ValueError("PLU decomposition only for square matrices.")
        P = [[1 if i == j else 0 for j in range(self.n)] for i in range(self.n)]
        L, U = self.lu_decomposition()
        return Matrix(self.n, self.n, self.field, P), L, U
# Test Cases

#(a) 
m1 = Matrix(3, 3, "real", [[1, 2, 3], [4, 5, 6], [7, 8, 10]])
vector = [1, 2, 2]
print("Length of Vector:", m1.length_of_vector(vector))  # Expected: 3.0
print("Size of Matrix:", m1.size())  # Expected: (3, 3)
print("Rank of Matrix:", m1.rank())  # Expected: 3
print("Nullity of Matrix:", m1.nullity())  # Expected: 0

# (b)
print("RREF of Matrix m1:")
rref_matrix = m1.rref(show_operations=True)
print(rref_matrix)

# (c) 
vectors = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
print("Linearly Independent:", m1.linearly_independent(vectors))  # Expected: False

# (d) 
print("Dimension of Span:", m1.dimension_of_span(vectors))  # Expected: 2
print("Basis for Span:", m1.basis_for_span(vectors))

# (e) 
U, V = m1.rank_factorization()
print("Rank Factorization U:", U)
print("Rank Factorization V:", V)

# (f) LU decomposition
print("LU Decomposition:")
L, U = m1.lu_decomposition()
print("L:", L)
print("U:", U)

# (g) PLU decomposition
print("PLU Decomposition:")
P, L, U = m1.plu_decomposition()
print("P:", P)
print("L:", L)
print("U:", U)