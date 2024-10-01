#pragma once

template <typename T>
std::ostream& operator << (std::ostream& out, const Matrix<T>& matix)
{
	const size_t rows = matix.sizeRows();
	const size_t columns = matix.sizeColumns();

	for (size_t i = 0; i < rows; ++i)
	{
		for (size_t j = 0; j < columns; ++j)
		{
			std::cout << matix[i][j] << " ";
		}
		std::cout << "\n\n";
	}
	return out;
}

template <typename T>
std::ostream& operator << (std::ostream& out, const std::vector<T>& vectorOutput)
{
	const size_t rows = vectorOutput.size();

	for (size_t i = 0; i < rows; ++i)
	{
		std::cout << vectorOutput[i] << " ";
	}
	return out;
}

std::vector<double> operator* (const double& first, const std::vector<double>& second)
{
	std::vector<double> result = second;
	const size_t sizeSecond = second.size();
	for (size_t i = 0; i < sizeSecond; ++i)
	{
		result[i] = first * second[i];
	}
	return result;
}

std::vector<double> operator* (const std::vector<double>& first, const double& second)
{
	size_t sizeFirst = first.size();
	std::vector<double> result(sizeFirst);
	for (size_t i = 0; i < sizeFirst; ++i)
	{
		result[i] = second * first[i];
	}
	return result;
}

double operator* (const std::vector<double>& first,
				  const std::vector<double>& second)
{
	double result = 0;
	const size_t sizeFirst = first.size();
	for (size_t i = 0; i < sizeFirst; ++i)
	{
		result += first[i] * second[i];
	}
	return result;
}

std::vector<double>& operator+= (std::vector<double>& first,
								 const std::vector<double>& second)
{
	const size_t sizeFirst = first.size();
	for (size_t i = 0; i < sizeFirst; ++i)
	{
		first[i] += second[i];
	}
	return first;
}

std::vector<double> operator+ (const std::vector<double>& first,
							   const std::vector<double>& second)
{
	std::vector<double> result = first;
	const size_t sizeFirst = first.size();
	for (size_t i = 0; i < sizeFirst; ++i)
	{
		result[i] += second[i];
	}
	return result;
}

std::vector<double> operator- (const std::vector<double>& first, const std::vector<double>& second)
{
	std::vector<double> result = first;
	const size_t sizeFirst = first.size();
	for (size_t i = 0; i < sizeFirst; ++i)
	{
		result[i] -= second[i];
	}
	return result;
}

auto operator* (const double& first, const Matrix<double>& second)
{
	const size_t rowsSecond = second.sizeRows();
	const size_t columnsSecond = second.sizeColumns();
	Matrix<double> product(rowsSecond, columnsSecond);

	for (size_t i = 0; i < rowsSecond; ++i)
	{
		for (size_t j = 0; j < columnsSecond; ++j)
		{
			product[i][j] = second[i][j] * first;
		}
	}
	return product;
}

std::vector<double> operator* (const Matrix<double>& first, const std::vector<double>& second)
{
	const size_t rowsFirst = first.sizeRows();
	const size_t columnsFirst = first.sizeColumns();
	std::vector <double> result(rowsFirst);

	for (size_t i = 0; i < rowsFirst; ++i)
	{
		for (size_t j = 0; j < columnsFirst; ++j)
		{
			result[i] += first[i][j] * second[j];
		}
	}
	return result;
}

Matrix<double> operator+ (const Matrix<double>& first, const Matrix<double>& second)
{
	Matrix<double> result = first;
	const size_t rowsFirst = first.sizeRows();
	const size_t columnsFirst = first.sizeColumns();

	for (size_t i = 0; i < rowsFirst; ++i)
	{
		for (size_t j = 0; j < columnsFirst; ++j)
		{
			result[i][j] += second[i][j];
		}
	}

	return result;
}

Matrix<double> operator- (const Matrix<double>& first, const Matrix<double>& second)
{
	Matrix<double> result = first;
	const size_t rowsFirst = first.sizeRows();
	const size_t columnsFirst = first.sizeColumns();

	for (size_t i = 0; i < rowsFirst; ++i)
	{
		for (size_t j = 0; j < columnsFirst; ++j)
		{
			result[i][j] -= second[i][j];
		}
	}

	return result;
}