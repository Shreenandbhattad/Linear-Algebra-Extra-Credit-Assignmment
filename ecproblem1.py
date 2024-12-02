class ComplexNumber:
    def __init__(self, real, imag):
        self.real = real
        self.imag = imag

    def __add__(self, other):
        return ComplexNumber(self.real + other.real, self.imag + other.imag)

    def __mul__(self, other):
        return ComplexNumber(
            self.real * other.real - self.imag * other.imag,
            self.real * other.imag + self.imag * other.real
        )

    def __truediv__(self, other):
        denom = other.real**2 + other.imag**2
        if denom == 0:
            raise ZeroDivisionError("Division by zero.")
        return ComplexNumber(
            (self.real * other.real + self.imag * other.imag) / denom,
            (self.imag * other.real - self.real * other.imag) / denom
        )

    def abs(self):
        return (self.real**2 + self.imag**2)**0.5

    def conjugate(self):
        return ComplexNumber(self.real, -self.imag)

    def __repr__(self):
        return f"{self.real} + {self.imag}i"


class Vector:
    def __init__(self, n, field="real", elements=None):
        self.n, self.field = n, field.lower()
        if elements:
            if len(elements) != n:
                raise ValueError("Invalid element length.")
            self.elements = elements
        else:
            self.elements = [self._input_value() for _ in range(n)]

    def _input_value(self):
        return (
            ComplexNumber(float(input("Real: ")), float(input("Imag: ")))
            if self.field == "complex"
            else float(input("Value: "))
        )

    def __add__(self, other):
        return Vector(self.n, self.field, [a + b for a, b in zip(self.elements, other.elements)])

    def __mul__(self, other):
        return sum(a * b for a, b in zip(self.elements, other.elements))

    def __repr__(self):
        return f"Vector({self.elements})"


class Matrix:
    def __init__(self, n, m, field="real", vectors=None, rows=None):
        self.n, self.m, self.field = n, m, field.lower()
        self.rows = (
            rows if rows else [[vec.elements[i] for vec in vectors] for i in range(n)]
            if vectors
            else [[self._input_value() for _ in range(m)] for _ in range(n)]
        )

    def _input_value(self):
        return (
            ComplexNumber(float(input("Real: ")), float(input("Imag: ")))
            if self.field == "complex"
            else float(input("Value: "))
        )

    def __add__(self, other):
        if self.n != other.n or self.m != other.m:
            raise ValueError("Matrix dimensions must match for addition.")
        return Matrix(self.n, self.m, self.field, rows=[
            [self.rows[i][j] + other.rows[i][j] for j in range(self.m)]
            for i in range(self.n)
        ])

    def __mul__(self, other):
        if not isinstance(other, Matrix):
            raise TypeError("Matrix multiplication is only defined for matrices.")
        if self.m != other.n:
            raise ValueError("Incompatible dimensions for multiplication.")
        return Matrix(self.n, other.m, self.field, rows=[
            [sum(self.rows[i][k] * other.rows[k][j] for k in range(self.m)) for j in range(other.m)]
            for i in range(self.n)
        ])

    def transpose(self):
        return Matrix(self.m, self.n, self.field, rows=[[self.rows[j][i] for j in range(self.n)] for i in range(self.m)])

    def conjugate(self):
        return Matrix(self.n, self.m, self.field, rows=[
            [elem.conjugate() if isinstance(elem, ComplexNumber) else elem for elem in row]
            for row in self.rows
        ])

    def transpose_conjugate(self):
        return self.transpose().conjugate()

    def get_row(self, index):
        if index < 0 or index >= self.n:
            raise IndexError("Row index out of bounds.")
        return self.rows[index]

    def get_column(self, index):
        if index < 0 or index >= self.m:
            raise IndexError("Column index out of bounds.")
        return [self.rows[i][index] for i in range(self.n)]

    def __repr__(self):
        return "\n".join(" ".join(map(str, row)) for row in self.rows)


# Test Cases
#a)
c1 = ComplexNumber(3, 4)  # 3 + 4i
c2 = ComplexNumber(1, 2)  # 1 + 2i
print("Addition (a):", c1 + c2)  # 4 + 6i
print("Multiplication (a):", c1 * c2)  # -5 + 10i
print("Division (a):", c1 / c2)  # 2.2 - 0.4i
print("Absolute value (a):", c1.abs())  # 5.0
print("Conjugate (a):", c1.conjugate())  # 3 - 4i
#b)
v1 = Vector(3, "real", [1, 2, 3])  # Real vector
v2 = Vector(3, "complex", [ComplexNumber(1, 2), ComplexNumber(3, 4), ComplexNumber(5, 6)])  # Complex vector
v3 = Vector(3, "real", [4, 5, 6])  # Real vector

print("Vector Addition (b):", v1 + v3)  # [5, 7, 9]
print("Dot Product (b):", v1 * v3)  # 1*4 + 2*5 + 3*6 = 32
#c)
m1 = Matrix(2, 3, "real", rows=[[1, 2, 3], [4, 5, 6]])  # Real matrix
v1 = Vector(2, "complex", [ComplexNumber(1, 1), ComplexNumber(2, -1)])  # Complex column vector
v2 = Vector(2, "complex", [ComplexNumber(3, 0), ComplexNumber(-1, 4)])  # Complex column vector
m2 = Matrix(2, 2, "complex", vectors=[v1, v2])  # Complex matrix from vectors

print("Matrix m1 (c):\n", m1)
print("Matrix m2 (c):\n", m2)
#d)
v1 = Vector(3, "real", [1, 2, 3])  # Real column vector
v2 = Vector(3, "real", [4, 5, 6])  # Real column vector
m3 = Matrix(3, 2, "real", vectors=[v1, v2])  # Matrix from column vectors

print("Matrix from Vectors (d):\n", m3)
#e)
m4 = Matrix(2, 2, "real", rows=[[1, 2], [3, 4]])
m5 = Matrix(2, 2, "real", rows=[[5, 6], [7, 8]])

print("Matrix Addition (e):\n", m4 + m5)  # [[6, 8], [10, 12]]

m6 = Matrix(2, 3, "real", rows=[[1, 2, 3], [4, 5, 6]])
m7 = Matrix(3, 2, "real", rows=[[7, 8], [9, 10], [11, 12]])

print("Matrix Multiplication (e):\n", m6 * m7)  # [[58, 64], [139, 154]]
#f)
m8 = Matrix(3, 3, "real", rows=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
print("Row 1 (f):", m8.get_row(0))  # [1, 2, 3]
print("Column 1 (f):", m8.get_column(0))  # [1, 4, 7]
#g)
m9 = Matrix(2, 2, "complex", rows=[[ComplexNumber(1, 1), ComplexNumber(2, -1)], [ComplexNumber(3, 0), ComplexNumber(4, 4)]])
print("Transpose (g):\n", m9.transpose())  # Transpose of the matrix
print("Conjugate (g):\n", m9.conjugate())  # Conjugate of each entry
print("Transpose-Conjugate (g):\n", m9.transpose_conjugate())  # Combined transpose and conjugate

