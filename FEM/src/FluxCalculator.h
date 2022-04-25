#pragma once

#include "GridData.h"
#include "BasisInfo.h"


class FluxCalculator
{
public:
	FluxCalculator() = default;
	
	void CalculateFlux(const BasisInfo& basisInfo, const GridData& grid, const std::vector<double>& solution);

	void SetMode(BasisType mode) { m_Mode = mode; }

	auto& GetFlux() const { return m_Fluxes; }

private:
	std::vector<std::vector<Vec2D>> GetNormals();

	std::vector<std::vector<Point2D>> GetEdgesOrigins();
	std::vector<std::vector<Vec2D>> GetEdgesDirections();

	std::vector<double> JacobianModifier();

private:
	std::vector<double> m_Fluxes;

	BasisType m_Mode;
};

