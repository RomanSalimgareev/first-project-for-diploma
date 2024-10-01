#pragma once

// Матрица L по Халецкому
Matrix<double> matrixL(const Matrix<double>& matrix)
{
	size_t n = matrix.sizeRows();
	Matrix<double> L(n, n);

	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j <= i; ++j) {
			double sum = 0;
			if (j == i) {
				for (size_t k = 0; k < j; ++k) {
					sum += pow(L[j][k], 2);
				}
				L[j][j] = sqrt(matrix[j][j] - sum);
			}
			else {
				for (size_t k = 0; k < j; ++k) {
					sum += (L[i][k] * L[j][k]);
				}
				L[i][j] = (matrix[i][j] - sum) / L[j][j];
			}
		}
	}

	return L;
}

// Транспонирование
template <typename T>
Matrix<T> transpose (const Matrix<T>& non_transpose)
{
	size_t r = non_transpose.sizeRows();
	size_t c = non_transpose.sizeColumns();

	Matrix<T> new_matrix(c, r);
	for (size_t i = 0; i < c; ++i)
	{
		for (size_t j = 0; j < r; ++j)
		{
			new_matrix[i][j] = non_transpose[j][i];
		}
	}
	return new_matrix;
}

// Метод Гаусса с выбором ведущего элемента
std::vector<double> solverGauss (const Matrix<double>& matrixCoefficients,
								const std::vector<double>& ColumnsOfFreeMembers)
{
	const size_t rows = matrixCoefficients.sizeRows();
	const size_t columns = matrixCoefficients.sizeColumns();
	Matrix<double> A = matrixCoefficients;
	std::vector<double> B = ColumnsOfFreeMembers;
	std::vector<double> result(rows);

	for (size_t i = 0; i < columns; ++i)
	{
		size_t maxIndex = i;
		for (size_t j = i; j < rows; ++j)
		{
			if (abs(A[j][i]) > abs(A[maxIndex][i]))
			{
				maxIndex = j;
			}
		}

		std::swap(A[i], A[maxIndex]);
		std::swap(B[i], B[maxIndex]);
		
		for (size_t k = i; k < rows; ++k)
		{
			double firstCoefficient = A[k][i];
			if (k == i)
				B[k] = B[k] / firstCoefficient;
			else
				B[k] -= B[i] * firstCoefficient;

			for (size_t m = i; m < columns; ++m)
			{
				if (k == i)
					A[k][m] = A[k][m] / firstCoefficient;
				else
					A[k][m] -= A[i][m] * firstCoefficient;
			}
		}
	}

	for (size_t i = rows; i-- > 0; )
	{
		double sum = 0;
		for (size_t j = i + 1; j < rows; ++j)
		{
			sum += result[j] * A[i][j];
		}
		result[i] = B[i] - sum;
	}

	return result;
}