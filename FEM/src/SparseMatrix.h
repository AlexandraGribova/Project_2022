#pragma once
#include <vector>
#include <ostream>

class SparseMatrix
{
public:

	SparseMatrix() = default;
	SparseMatrix(std::vector<uint32_t>&ig, std::vector<uint32_t>&jg, std::vector<double>&di, std::vector<double>&ggl) :
		m_Ig(std::move(ig)), m_Jg(std::move(jg)), m_Di(std::move(di)), m_Ggl(std::move(ggl)), m_N(m_Di.size())
	{}

	auto& GetIg() { return m_Ig; }
	auto& GetJg() { return m_Jg; }
	auto& GetDi() { return m_Di; }
	auto& GetGgl() { return m_Ggl; }

	void ZeroOutRowAndCol(uint32_t row);
	std::vector<double> MultiplyByVector(const std::vector<double>& vector);
	void Print(std::ostream& out);
	
	std::vector<double> operator *(const std::vector<double>& vector) { return MultiplyByVector(vector); }

private:
	std::vector<uint32_t> m_Ig;
	std::vector<uint32_t> m_Jg;
	std::vector<double> m_Di;
	std::vector<double> m_Ggl;
	uint32_t m_N;
};

