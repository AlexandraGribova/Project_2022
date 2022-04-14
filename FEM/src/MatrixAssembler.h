#pragma once

#include <array>

#include "ApproximationInfo.h"
#include "DataStructures.h"
#include "GridData.h"
#include "functions.h"

class MatrixAssembler
{
	using Fn2d = std::function<double(double, double)>;
public:
	MatrixAssembler() = default;
	MatrixAssembler(const ApproximationInfo& info) :
		m_ApproximationInfo(info)
	{}

	std::unique_ptr<SparseMatrix> AssembleMatrix(const GridData& gridData, std::vector<uint32_t>& ig, std::vector<uint32_t>& jg, Fn2d additionalCoeffs = [](double, double) { return 1.0; });
	std::vector<double> AssembleLoadVector(const GridData& gridData);
	void ApplyDirichlet(const GridData& gridData, std::vector<double>& loadVector);
	void ApplyNeumann(const GridData& gridData, std::vector<double>& loadVector);

	auto& GetApproximationInfo() { return m_ApproximationInfo; }
	void SetApproximationInfo(const ApproximationInfo& info) { m_ApproximationInfo = info; }


protected:
	ApproximationInfo m_ApproximationInfo;
};