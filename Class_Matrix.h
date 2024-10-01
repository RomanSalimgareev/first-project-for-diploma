// Класс Матрица
#pragma once
#include <vector>

template <typename T>
class Matrix
{
public:

	Matrix () = default;

	Matrix (const size_t& r, const size_t& c)
	{
		m_matrix.resize(r);
		for (size_t i = 0; i < r; i++)
		{
			m_matrix[i].resize(c);
		}
	}

	Matrix(std::initializer_list<std::vector<T>> values) : m_matrix(values) {}

	const std::vector<T>& operator[] (size_t i) const
	{
		return m_matrix.at(i);
	}

	std::vector<T>& operator[] (size_t i)
	{
		return m_matrix.at(i);
	}

	size_t sizeRows() const
	{
		return m_matrix.size();
	}

	size_t sizeColumns() const
	{
		return m_matrix[0].size();
	}

	void eraseMatrix (const size_t& i)
	{
		m_matrix.erase(m_matrix.begin() + i);
	}

private:

	std::vector<std::vector<T>> m_matrix;

};
