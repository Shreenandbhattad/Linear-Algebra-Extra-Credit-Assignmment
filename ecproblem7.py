class Matrix:
    def __init__(self, n, m, rows=None):
        self.n = n
        self.m = m
        self.rows = rows if rows else [[0] * m for _ in range(n)]

    def cofactor_expansion(self):
        def minor(matrix, i, j):
            return [row[:j] + row[j+1:] for row in (matrix[:i] + matrix[i+1:])]

        def determinant(matrix):
            if len(matrix) == 1:
                return matrix[0][0]
            if len(matrix) == 2:
                return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
            det = 0
            for j in range(len(matrix[0])):
                det += ((-1) ** j) * matrix[0][j] * determinant(minor(matrix, 0, j))
            return det
        
        return determinant(self.rows)

    def plu_decomposition(self):
        A = [row[:] for row in self.rows]
        n = len(A)
        P = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
        L = [[0] * n for _ in range(n)]
        U = [[0] * n for _ in range(n)]
        
        for i in range(n):
            for j in range(i, n):
                U[i][j] = A[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
            for j in range(i + 1, n):
                L[j][i] = (A[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]
            L[i][i] = 1
        
        return P, L, U

    def determinant_plu(self):
        P, L, U = self.plu_decomposition()
        det_L = 1
        det_U = 1
        for i in range(len(L)):
            det_L *= L[i][i]
            det_U *= U[i][i]
        return det_L * det_U

    def rref(self):
        mat = [row[:] for row in self.rows]
        lead = 0
        row_count, column_count = len(mat), len(mat[0])

        for r in range(row_count):
            if lead >= column_count:
                return mat
            i = r
            while mat[i][lead] == 0:
                i += 1
                if i == row_count:
                    i = r
                    lead += 1
                    if column_count == lead:
                        return mat
            mat[i], mat[r] = mat[r], mat[i]
            lv = mat[r][lead]
            mat[r] = [mrx / lv for mrx in mat[r]]
            for i in range(row_count):
                if i != r:
                    lv = mat[i][lead]
                    mat[i] = [iv - lv * rv for rv, iv in zip(mat[r], mat[i])]
            lead += 1
        return mat

    def determinant_rref(self):
        mat = self.rref()
        det = 1
        for i in range(len(mat)):
            det *= mat[i][i]
        return det


# (a)
matrix_a = Matrix(3, 3, [[1, 2, 3], [0, 1, 4], [5, 6, 0]])
print("Determinant (Cofactor Expansion):", matrix_a.cofactor_expansion())

#(b) 
matrix_b = Matrix(3, 3, [[1, 2, 3], [0, 1, 4], [5, 6, 0]])
print("Determinant (PLU Decomposition):", matrix_b.determinant_plu())

#c)
matrix_c = Matrix(3, 3, [[1, 2, 3], [0, 1, 4], [5, 6, 0]])
print("Determinant (RREF):", matrix_c.determinant_rref())
