#pragma once

#include <functional>

#include "DataStructures.h"
#include "SparseMatrix.h"
#include "ApproximationInfo.h"
#include "MatrixAssembler.h"

enum class FEMMode 
{
	None = 0,
	Bilinear, Biquadratic,
	Linear, Quadratic
};

class FEMSolver
{
public:
	using Fn2d = std::function<double(double, double)>;

	FEMSolver(const std::string& path);

	auto& GetSolution() const { return m_Solution; }
	auto& GetLoadVector() const { return m_LoadVector; }
	auto& GetStiffnessMatrix() const { return m_StiffnessMatrix; }

	auto& GetMode() const { return m_Mode; }
	auto& GetPlastPressure() const { return m_PlastPressure; }

	auto& GetGridData() const { return m_GridData; }

	void SetGridData(GridData& gridData);
	void SetGlobalMatrixAssembler(const MatrixAssembler& assembler) { m_MatrixAssembler = assembler; }

	void LoadViscosity(const std::string& path);
	void LoadDomainData(const std::string& path);
	void LoadMode(const std::string& path);
	void LoadPlastPressure(const std::string& path);

	void CalculateStiffnesMatrix();
	void CalculateLoadVector();
	void ApplyDirichlet();
	void ApplyNeumann();

	void CalculateSolution();

private:
	FEMMode m_Mode;
	std::unique_ptr<SparseMatrix> m_StiffnessMatrix;
	std::vector<double> m_LoadVector;

	std::vector<double> m_Solution;
	
	GridData m_GridData;

	MatrixAssembler m_MatrixAssembler;

	std::vector<double> m_Viscosity;
	std::vector<std::pair<double, double>> m_Permeabilities;
	double m_PlastPressure;

	double m_Phi;
	double m_Saturation;

private:
	void PerformGaussianReduction();
	void ReduceRow(uint32_t row);

};
