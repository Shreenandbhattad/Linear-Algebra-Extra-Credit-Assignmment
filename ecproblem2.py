class ComplexNumber:
    def __init__(self, real, imag):
        self.real = real
        self.imag = imag

    def __add__(self, other):
        return ComplexNumber(self.real + other.real, self.imag + other.imag)

    def __mul__(self, other):
        real = self.real * other.real - self.imag * other.imag
        imag = self.real * other.imag + self.imag * other.real
        return ComplexNumber(real, imag)

    def __truediv__(self, other):
        denom = other.real ** 2 + other.imag ** 2
        if denom == 0:
            raise ZeroDivisionError("Division by zero is not allowed.")
        real = (self.real * other.real + self.imag * other.imag) / denom
        imag = (self.imag * other.real - self.real * other.imag) / denom
        return ComplexNumber(real, imag)

    def conjugate(self):
        return ComplexNumber(self.real, -self.imag)

    def __eq__(self, other):
        return self.real == other.real and self.imag == other.imag

    def __repr__(self):
        return f"{self.real} + {self.imag}i"


class Matrix:
    def __init__(self, n, m, field="real", rows=None):
        self.n = n
        self.m = m
        self.field = field.lower()
        self.rows = rows if rows else [[0] * m for _ in range(n)]

    def is_zero(self):
        return all(all(x == 0 for x in row) for row in self.rows)

    def is_square(self):
        return self.n == self.m

    def is_symmetric(self):
        if not self.is_square():
            return False
        return all(self.rows[i][j] == self.rows[j][i] for i in range(self.n) for j in range(self.n))

    def is_hermitian(self):
        if not self.is_square():
            return False
        return all(self.rows[i][j] == self.rows[j][i].conjugate() if isinstance(self.rows[j][i], ComplexNumber) else False for i in range(self.n) for j in range(self.n))

    def is_identity(self):
        if not self.is_square():
            return False
        return all(self.rows[i][j] == (1 if i == j else 0) for i in range(self.n) for j in range(self.n))

    def is_singular(self):
        return not self.is_invertible()

    def is_invertible(self):
        return self.determinant() != 0 if self.is_square() else False

    def is_scalar(self):
        if not self.is_square():
            return False
        d = self.rows[0][0]
        return all(self.rows[i][i] == d for i in range(self.n)) and all(self.rows[i][j] == 0 for i in range(self.n) for j in range(self.n) if i != j)

    def is_nilpotent(self):
        if not self.is_square():
            return False
        temp = self
        for _ in range(1, self.n + 1):
            temp = temp * self
            if temp.is_zero():
                return True
        return False

    def is_orthogonal(self):
        if not self.is_square():
            return False
        return self * self.transpose() == Matrix.identity(self.n, self.field)

    def is_unitary(self):
        if not self.is_square():
            return False
        return self * self.transpose_conjugate() == Matrix.identity(self.n, self.field)

    def is_diagonalizable(self):
        if not self.is_square():
            return False
        return self.is_symmetric() or self.is_hermitian()

    def has_lu_decomposition(self):
        if not self.is_square():
            return False
        return not self.is_singular()

    def determinant(self):
        if not self.is_square():
            raise ValueError("Determinant only defined for square matrices.")
        if self.n == 1:
            return self.rows[0][0]
        if self.n == 2:
            return self.rows[0][0] * self.rows[1][1] - self.rows[0][1] * self.rows[1][0]
        det = 0
        for c in range(self.n):
            sub_matrix = Matrix(self.n - 1, self.n - 1, self.field, [
                row[:c] + row[c + 1:] for row in self.rows[1:]
            ])
            det += ((-1) ** c) * self.rows[0][c] * sub_matrix.determinant()
        return det

    def transpose(self):
        return Matrix(self.m, self.n, self.field, [[self.rows[j][i] for j in range(self.n)] for i in range(self.m)])

    def transpose_conjugate(self):
        return self.transpose().conjugate()

    def conjugate(self):
        return Matrix(self.n, self.m, self.field, [
            [(x.conjugate() if isinstance(x, ComplexNumber) else x) for x in row] for row in self.rows
        ])

    @staticmethod
    def identity(size, field="real"):
        return Matrix(size, size, field, [[1 if i == j else 0 for j in range(size)] for i in range(size)])


# Test cases
m1 = Matrix(3, 3, "real", [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
print("Identity:", m1.is_identity())

m2 = Matrix(3, 3, "real", [[0, 0, 0], [0, 0, 0], [0, 0, 0]])
print("Zero:", m2.is_zero())

m3 = Matrix(2, 2, "real", [[1, 2], [2, 1]])
print("Symmetric:", m3.is_symmetric())

m4 = Matrix(2, 2, "real", [[1, 2], [3, 4]])
print("Singular:", m4.is_singular())

m5 = Matrix(2, 2, "complex", [[ComplexNumber(1, 0), ComplexNumber(0, -1)], [ComplexNumber(0, 1), ComplexNumber(1, 0)]])
print("Hermitian:", m5.is_hermitian())
