class Matrix:
    def __init__(self, n, m, rows=None):
        self.n = n
        self.m = m
        self.rows = rows if rows else [[0] * m for _ in range(n)]
    
    def is_in_span(self, S, v):
        augmented_matrix = [row + [v[i]] for i, row in enumerate(S)]
        rref_matrix = self.rref(augmented_matrix)
        return rref_matrix[-1][-1] == 0

    def linear_combination(self, S, v):
        augmented_matrix = [row + [v[i]] for i, row in enumerate(S)]
        rref_matrix = self.rref(augmented_matrix)
        return [r[-1] for r in rref_matrix]

    def same_span(self, S1, S2):
        return self.rank(S1) == self.rank(S2)

    def rank(self, matrix):
        rref_matrix = self.rref(matrix)
        return sum(1 for row in rref_matrix if any(cell != 0 for cell in row))

    def rref(self, mat):
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
    
    def compute_coordinates(self, B, v):
        return self.linear_combination(B, v)

    def vector_from_coordinates(self, B, coords):
        return [sum(b * c for b, c in zip(row, coords)) for row in B]

    def change_of_basis(self, B1, B2):
        return [[sum(B2[i][j] * B1[i][k] for i in range(len(B1))) for j in range(len(B2))] for k in range(len(B1))]

    def coordinates_in_new_basis(self, B1, B2, coords_B1):
        change_of_basis_matrix = self.change_of_basis(B1, B2)
        return [sum(change_of_basis_matrix[i][j] * coords_B1[j] for j in range(len(coords_B1))) for i in range(len(B2))]
#a)
S = [[1, 2], [3, 4]]
v = [5, 6]
mat = Matrix(2, 2, S)
print("v in span of S:", mat.is_in_span(S, v))
#b0
S = [[1, 2], [3, 4]]
v = [5, 6]
mat = Matrix(2, 2, S)
print("Linear combination of v:", mat.linear_combination(S, v))
#c)
S1 = [[1, 2], [3, 4]]
S2 = [[5, 6], [7, 8]]
mat = Matrix(2, 2, S1)
print("S1 and S2 span the same subspace:", mat.same_span(S1, S2))
#d)
B = [[1, 0], [0, 1]]
v = [3, 4]
mat = Matrix(2, 2, B)
print("Coordinates of v:", mat.compute_coordinates(B, v))

coords = [3, 4]
print("Vector from coordinates:", mat.vector_from_coordinates(B, coords))
#e)
B1 = [[1, 0], [0, 1]]
B2 = [[2, 1], [1, 3]]
mat = Matrix(2, 2, B1)
print("Change of Basis Matrix:", mat.change_of_basis(B1, B2))
#f)
coords_B1 = [3, 4]
mat = Matrix(2, 2, B1)
print("Coordinates of v in B2:", mat.coordinates_in_new_basis(B1, B2, coords_B1))
