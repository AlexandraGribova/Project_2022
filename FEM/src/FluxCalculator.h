#pragma once

#include "GridData.h"
#include "BasisInfo.h"


class FluxCalculator
{
public:
	FluxCalculator();
	
	void Init();

	void CalculateFlux(const BasisInfo& basisInfo, const GridData& grid, const std::vector<double>& solution);

	void SetMode(BasisType mode) { m_Mode = mode; Init(); }

	auto& GetFlux() const { return m_Fluxes; }

private:

	std::vector<std::vector<Vec2D>> GetPrecalculatedGradXiIntegral(const BasisInfo& basisInfo, const auto& origins, const auto& directions);
	// Συ
	std::vector<std::vector<std::vector<Vec2D>>> GetPrecalculatedGradXiValues(const BasisInfo& basisInfo, const auto& origins, const auto& directions);

	std::vector<std::vector<Vec2D>> GetNormals();

	std::vector<std::vector<Point2D>> GetEdgesOrigins();
	std::vector<std::vector<Vec2D>> GetEdgesDirections();

	std::vector<double> JacobianModifier();

	double GetEdgeLength(const Element& elem, uint32_t elemNum, uint32_t edgeNum, double width, double height);

private:
	std::vector<double> m_Fluxes;

	std::function<void(Vec2D&, const std::vector<double>&)> m_NormalModifier;

	BasisType m_Mode;
};

