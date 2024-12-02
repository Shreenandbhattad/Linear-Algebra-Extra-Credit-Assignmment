class Matrix:
    def __init__(self, n, m, rows=None):
        self.n = n
        self.m = m
        self.rows = rows if rows else [[0] * m for _ in range(n)]

    def is_square(self):
        return self.n == self.m

    def is_invertible(self):
        if not self.is_square():
            return False
        return self.rank() == self.n 

    def rank(self):
        mat = [row[:] for row in self.rows]
        lead = 0
        row_count, column_count = len(mat), len(mat[0])
        for r in range(row_count):
            if lead >= column_count:
                return sum(any(cell != 0 for cell in row) for row in mat)
            i = r
            while mat[i][lead] == 0:
                i += 1
                if i == row_count:
                    i = r
                    lead += 1
                    if column_count == lead:
                        return sum(any(cell != 0 for cell in row) for row in mat)
            mat[i], mat[r] = mat[r], mat[i]
            lv = mat[r][lead]
            mat[r] = [mrx / lv for mrx in mat[r]]
            for i in range(row_count):
                if i != r:
                    lv = mat[i][lead]
                    mat[i] = [iv - lv * rv for rv, iv in zip(mat[r], mat[i])]
            lead += 1
        return sum(any(cell != 0 for cell in row) for row in mat)

    def inverse_row_reduction(self):
        if not self.is_square() or not self.is_invertible():
            return "Matrix is not invertible"
        mat = [row[:] + [1 if i == j else 0 for j in range(self.n)] for i, row in enumerate(self.rows)]
        lead = 0
        row_count, column_count = len(mat), len(mat[0])
        for r in range(row_count):
            if lead >= column_count:
                return "Matrix is not invertible"
            i = r
            while mat[i][lead] == 0:
                i += 1
                if i == row_count:
                    i = r
                    lead += 1
                    if column_count == lead:
                        return "Matrix is not invertible"
            mat[i], mat[r] = mat[r], mat[i]
            lv = mat[r][lead]
            mat[r] = [mrx / lv for mrx in mat[r]]
            for i in range(row_count):
                if i != r:
                    lv = mat[i][lead]
                    mat[i] = [iv - lv * rv for rv, iv in zip(mat[r], mat[i])]
            lead += 1
        return Matrix(self.n, self.n, [row[self.n:] for row in mat])

    def minor(self, row, col):
        return [r[:col] + r[col+1:] for r in (self.rows[:row] + self.rows[row+1:])]

    def determinant(self):
        if not self.is_square():
            return "Matrix is not square"
        if self.n == 1:
            return self.rows[0][0]
        if self.n == 2:
            return self.rows[0][0] * self.rows[1][1] - self.rows[0][1] * self.rows[1][0]
        det = 0
        for c in range(self.n):
            det += ((-1) ** c) * self.rows[0][c] * Matrix(self.n - 1, self.n - 1, self.minor(0, c)).determinant()
        return det

    def adjoint(self):
        if not self.is_square() or self.determinant() == 0:
            return "Matrix is not invertible"
        adj = []
        for i in range(self.n):
            adj_row = []
            for j in range(self.n):
                minor_matrix = Matrix(self.n - 1, self.n - 1, self.minor(i, j))
                adj_row.append(((-1) ** (i + j)) * minor_matrix.determinant())
            adj.append(adj_row)
        return Matrix(self.n, self.n, adj).transpose()

    def transpose(self):
        return Matrix(self.n, self.m, [[self.rows[j][i] for j in range(self.n)] for i in range(self.m)])

    def inverse_adjoint(self):
        if not self.is_square() or self.determinant() == 0:
            return "Matrix is not invertible"
        adj = self.adjoint()
        det = self.determinant()
        return Matrix(self.n, self.n, [[adj.rows[i][j] / det for j in range(self.n)] for i in range(self.n)])

# Test cases
A = [[2, 1], [1, 3]]
mat = Matrix(2, 2, A)

print("Is matrix square:", mat.is_square())  
print("Is matrix invertible:", mat.is_invertible())  
print("Inverse by row reduction:", mat.inverse_row_reduction())  
print("Inverse by adjoint:", mat.inverse_adjoint())  
