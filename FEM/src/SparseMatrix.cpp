#include "SparseMatrix.h"
#include "functions.h"
#include <iomanip>
#include <iostream>

void SparseMatrix::ZeroOutRowAndCol(uint32_t row)
{
	auto ibeg = m_Ig[row];
	auto iend = m_Ig[row + 1];
	for (auto i = ibeg; i < iend; i++)
	{
		m_Ggl[i] = 0.0;
	}
	for (auto j = row + 1; j < m_N; j++)
	{
		auto jbeg = m_Ig[j];
		auto jend = m_Ig[j + 1];
		auto index = binarySearch(m_Jg, row, jbeg, jend - 1);
		if (index != -1)
		{
			m_Ggl[index] = 0.0;
		}
	}
}

std::vector<double> SparseMatrix::MultiplyByVector(const std::vector<double>& vector)
{
	std::vector<double> result(vector.size(), 0.0);
	for (auto i = 0u; i < vector.size(); i++)
	{
		result[i] = m_Di[i] * vector[i];
	}

	for (auto i = 0u; i < vector.size(); i++)
	{
		for (auto j = m_Ig[i]; j < m_Ig[i + 1]; j++)
		{
			result[i] += m_Ggl[j] * vector[m_Jg[j]];
			result[m_Jg[j]] += m_Ggl[j] * vector[i];
		}
	}
	return result;
}


void SparseMatrix::Print(std::ostream& out)
{
	for (uint32_t i = 0; i < m_N; i++)
	{
		auto beg = m_Ig[i];
		auto end = m_Ig[i + 1];
		auto length = end - beg;
		if (length == 0)
		{
			for (uint32_t j = 0; j < i; j++)
			{
				out << std::setw(9) << std::setprecision(3) << 0.0 << " ";
			}
			out << std::setw(9) << std::setprecision(3) << m_Di[i] << "\n";
		}
		else
		{
			for (uint32_t j = 0, k = 0; j < i; j++)
			{
				if (binarySearch(m_Jg, j, beg, end - 1) != -1)
				{
					out << std::setw(9) << std::setprecision(3) << m_Ggl[beg + k] << " ";
					k++;
				}
				else
				{
					out << std::setw(9) << std::setprecision(3) << 0.0 << " ";
				}
			}
			out << std::setw(9) << std::setprecision(3) << m_Di[i] << "\n";
		}
	}
}
