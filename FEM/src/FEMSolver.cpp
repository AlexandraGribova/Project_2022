#include <iostream>
#include <utility>
#include <ostream>
#include <fstream>

#include "FEMSolver.h"
#include "functions.h"
#include "GridData.h"

FEMSolver::FEMSolver(const std::string& path) : 
	m_Permeabilities(1, { 1.0, INFINITY }), m_Phi(1.0), m_Saturation(1.0)
{
	LoadMode(path);
	switch (m_Mode)
	{
	case FEMMode::None:
		break;
	case FEMMode::Bilinear:
		m_MatrixAssembler.SetApproximationInfo(ApprInfo::BilinearInfo);
		break;
	case FEMMode::Biquadratic:
		m_MatrixAssembler.SetApproximationInfo(ApprInfo::BiquadraticInfo);
		break;
	case FEMMode::Linear:
		m_MatrixAssembler.SetApproximationInfo(ApprInfo::LinearInfo);
		break;
	case FEMMode::Quadratic:
		m_MatrixAssembler.SetApproximationInfo(ApprInfo::QuadraticInfo);
		break;
	default:
		break;
	}
}

void FEMSolver::SetGridData(GridData& gridData)
{
	
	if (m_Mode == FEMMode::Biquadratic)
		gridData.ModifyGrid();
	if (m_Mode == FEMMode::Quadratic)
		gridData.ModifyTriangleGrid();

	if (m_Mode == FEMMode::Linear || m_Mode == FEMMode::Quadratic)
		gridData.Triangulate();
	m_GridData = std::move(gridData);
}

void FEMSolver::LoadViscosity(const std::string& path)
{
	std::ifstream in(path);
	uint32_t N;
	in >> N;
	m_Viscosity.resize(N);
	for (uint32_t i = 0; i < N; i++)
	{
		double v; 
		in >> v;
		m_Viscosity[i] = v;
	}
}

void FEMSolver::LoadDomainData(const std::string& path)
{
	std::ifstream in(path);
	in >> m_Saturation;
	in >> m_Phi;
	uint32_t numberOfPermeabilities;
	in >> numberOfPermeabilities;
	m_Permeabilities.clear();
	for (uint32_t i = 0; i < numberOfPermeabilities; i++)
	{
		double permeability = 0.0, boundary = 0.0;
		in >> permeability >> boundary;
		m_Permeabilities.push_back({ permeability, boundary });
	}
}

void FEMSolver::LoadMode(const std::string& path)
{
	std::ifstream in(path);
	std::string mode;
	in >> mode;
	if (mode == "BL")
		m_Mode = FEMMode::Bilinear;
	else if (mode == "BQ")
		m_Mode = FEMMode::Biquadratic;
	else if (mode == "L")
		m_Mode = FEMMode::Linear;
	else if (mode == "Q")
		m_Mode = FEMMode::Quadratic;
}

void FEMSolver::LoadPlastPressure(const std::string& path)
{
	std::ifstream in(path);
	in >> m_PlastPressure;
}

void FEMSolver::CalculateStiffnesMatrix()
{
	if (m_GridData.Elements.size() == 0 || m_GridData.Nodes.size() == 0)
	{
		std::cout << "You have to set nodes and elements.\n";
		return;
	}

	auto&& [ig, jg] = generatePortrait(m_GridData.Nodes.size(), m_GridData.Elements);

	// множитель перед матрицей
	auto multiply = m_Viscosity.size() > 0 ? 0.0 : 1.0;
	for (const auto& vis : m_Viscosity)
		multiply += 1.0 / vis;
	m_MatrixAssembler.GetApproximationInfo().SetStiffnessCoefficient(multiply * m_MatrixAssembler.GetApproximationInfo().GetStiffnessCoefficient());

	auto permeability = [&](double x, double y) {
		// m_Permeabilities должно быть отсортировано по возрастанию по второму параметру.
		for (const auto& p : m_Permeabilities)
		{
			auto&& [value, boundary] = p;
			if (x < boundary) return value;
		}
		return 1.0;
	};

	m_StiffnessMatrix = std::move(m_MatrixAssembler.AssembleMatrix(m_GridData, ig, jg, permeability));
}

void FEMSolver::CalculateLoadVector()
{
	if (m_GridData.Elements.size() == 0 || m_GridData.Nodes.size() == 0)
	{
		std::cout << "You have to set nodes and elements.\n";
		return;
	}
	// В нашем случае это равно нулю.
	m_LoadVector = std::move(m_MatrixAssembler.AssembleLoadVector(m_GridData));
}

void FEMSolver::ApplyDirichlet()
{
	m_MatrixAssembler.ApplyDirichlet(m_GridData, m_LoadVector);
	PerformGaussianReduction();
}


void FEMSolver::ApplyNeumann()
{
	m_MatrixAssembler.ApplyNeumann(m_GridData, m_LoadVector);
}


void FEMSolver::CalculateSolution()
{
	m_Solution.resize(m_LoadVector.size(), 0.0);
	auto eps = 1e-21;
	auto maxIter = static_cast<uint32_t>(1e+7);
	auto loadVectorNorm = sqrt(dotProduct(m_LoadVector, m_LoadVector));
	auto z = m_LoadVector;
	auto r = m_LoadVector;

	auto p = *m_StiffnessMatrix.get() * z;
	auto rDotR = dotProduct(r, r);
	auto residual = sqrt(rDotR) / loadVectorNorm;

	for (uint32_t k = 1; k < maxIter && residual > eps; k++)
	{
		auto pDotP = dotProduct(p, p);
		auto a = dotProduct(p, r) / pDotP;

		for (uint32_t i = 0; i < m_LoadVector.size(); i++)
		{
			m_Solution[i] += a * z[i];
			r[i] -= a * p[i];
		}
		auto Ar = *m_StiffnessMatrix.get() * r;
		auto b = -dotProduct(p, Ar) / pDotP;
		for (uint32_t i = 0; i < m_LoadVector.size(); i++)
		{
			z[i] = r[i] + b * z[i];
			p[i] = Ar[i] + b * p[i];
		}
		rDotR = dotProduct(r, r);
		residual = sqrt(rDotR) / loadVectorNorm;
	}
}


void FEMSolver::PerformGaussianReduction()
{
	for (const auto& dirichletData : m_GridData.DirichletConditions)
	{
		auto row = dirichletData.Node;
		ReduceRow(row);
	}
	for (const auto& dirichletData : m_GridData.DirichletConditions)
	{
		auto row = dirichletData.Node;
		m_StiffnessMatrix->ZeroOutRowAndCol(row);
	}
}

void FEMSolver::ReduceRow(uint32_t row)
{
	auto& di = m_StiffnessMatrix->GetDi();
	auto& ggl = m_StiffnessMatrix->GetGgl();
	auto& ig = m_StiffnessMatrix->GetIg();
	auto& jg = m_StiffnessMatrix->GetJg();

	di[row] = 1.0;
	auto ibeg = ig[row];
	auto iend = ig[row + 1];
	for (uint32_t j = ibeg; j < iend; j++)
	{
		if (binarySearch(m_GridData.DirichletConditions, jg[j], 0, m_GridData.DirichletConditions.size() - 1) != -1) continue;
		auto multiplier = ggl[j];
		m_LoadVector[jg[j]] -= m_LoadVector[row] * multiplier;
	}
	for (uint32_t j = row + 1; j < m_LoadVector.size(); j++)
	{
		if (binarySearch(m_GridData.DirichletConditions, j, 0, m_GridData.DirichletConditions.size() - 1) != -1) continue;
		auto jbeg = ig[j];
		auto jend = ig[j + 1];
		auto index = binarySearch(jg, row, jbeg, jend - 1);
		if (index != -1)
		{
			auto multiplier = ggl[index];
			m_LoadVector[j] -= m_LoadVector[row] * multiplier;
		}
	}
}
