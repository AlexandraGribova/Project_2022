#include <algorithm>

#include "FluxCalculator.h"
#include "functions.h"
#include "Math.h"

void FluxCalculator::CalculateFlux(const BasisInfo& basisInfo, const GridData& grid, const std::vector<double>& solution)
{
	m_Fluxes = std::vector<double>(grid.Elements.back().Edges.back() + 1, 0.0);

	// Здесь и далее нормали и т.д. хранятся для четных и нечетных элементов.
	// Для прямоугольников они совпадают, но мне кажется что это лучше, чем 
	// каждый раз проверять тип элемента в if. Потом лучше переделать как в indexingFunction (или нет). <-------------------------
	static auto normals = GetNormals();
	static auto origins = GetEdgesOrigins();
	static auto directions = GetEdgesDirections();

	std::vector<std::vector<std::vector<Vec2D>>> precalculatedIntegrals;
	precalculatedIntegrals.emplace_back(GetPrecalculatedGradXiIntegral(basisInfo, origins[0], directions[0]));
	precalculatedIntegrals.emplace_back(GetPrecalculatedGradXiIntegral(basisInfo, origins[1], directions[1]));

	static auto jacobiMod = JacobianModifier();

	static auto integrationPoints = std::vector<double>{ 0.0, 0.5, 1.0 };
	static auto integrationCoefficients = std::vector<double>{ 1.0, 4.0, 1.0 };

	static auto indF = grid.IndexingFunction;

	// Для треугольников будет сложнее, потому что они разные, пока не буду специально запрашивать грани.
	// Еще мы считаем интеграл по "локальной грани" и у "перевернутого" треугольника это приводит к 
	// отрицательному якобиану.
	for (auto elemNum = 0u; elemNum < grid.Elements.size(); elemNum++)
	{
		auto& elem = grid.Elements[elemNum];
		auto elemDisparity = elemNum & 1; // Нечетность номера элемента для якобиана перевернутого треугольника (disparity - это не нечетность, но мне показалось смешно)
		auto&& [width, height] = getDimensions(grid.Nodes, elem.Nodes);
		auto J = width * height;
		auto JInv = std::vector<double>{ 1.0 / width, 1.0 / height }; // Храним только диагональ.
		for (auto edgeNum = 0u; edgeNum < elem.Edges.size(); edgeNum++)
		{
			auto& normal = normals[elemDisparity][edgeNum];
			auto edgeFlux = 0.0;
			// xi - это типо кси
			for (auto xiNum = 0; xiNum < basisInfo.GetBasisFunctions().size(); xiNum++)
			{
				auto gradIntegral = Vec2D{ 0.0, 0.0 };
				gradIntegral.X = precalculatedIntegrals[elemDisparity][indF(elemNum, xiNum)][edgeNum].X * JInv[0] * jacobiMod[elemDisparity];
				gradIntegral.Y = precalculatedIntegrals[elemDisparity][indF(elemNum, xiNum)][edgeNum].Y * JInv[1] * jacobiMod[elemDisparity];
				edgeFlux += solution[elem.Nodes[xiNum]] * (gradIntegral.X * normal.X + gradIntegral.Y * normal.Y);
			}
			edgeFlux *= elem.EdgeDirections[edgeNum] * J / 6.0;
			// Наверное это не лучший способ.
			m_Fluxes[elem.Edges[edgeNum]] = m_Fluxes[elem.Edges[edgeNum]] == 0.0 ? edgeFlux : (edgeFlux + m_Fluxes[elem.Edges[edgeNum]]) / 2.0;
		}
	}
}

std::vector<std::vector<Vec2D>> FluxCalculator::GetPrecalculatedGradXiIntegral(const BasisInfo& basisInfo, const auto& origins, const auto& directions)
{
	static auto integrationCoefficients = std::vector<double>{ 1.0, 4.0, 1.0 };

	auto precalculatedGradXiIntegral = std::vector<std::vector<Vec2D>>(basisInfo.GetBasisFunctions().size());
	std::for_each(precalculatedGradXiIntegral.begin(), precalculatedGradXiIntegral.end(), [&](auto& val) {val.resize(origins.size()); });

	auto precalculatedGradXiValues = GetPrecalculatedGradXiValues(basisInfo, origins, directions);
	for (auto xiNum = 0u; xiNum < precalculatedGradXiValues.size(); xiNum++)
	{
		for (auto edgeNum = 0u; edgeNum <  precalculatedGradXiValues[xiNum].size(); edgeNum++)
		{
			auto edge = precalculatedGradXiValues[xiNum][edgeNum];
			auto gradIntegral = Vec2D{ 0.0, 0.0 };
			for (auto pointNum = 0u; pointNum < edge.size(); pointNum++)
			{
				gradIntegral += integrationCoefficients[pointNum] * edge[pointNum];
			}
			precalculatedGradXiIntegral[xiNum][edgeNum] = gradIntegral;
		}
	}
	return precalculatedGradXiIntegral;
}

std::vector<std::vector<std::vector<Vec2D>>> FluxCalculator::GetPrecalculatedGradXiValues(const BasisInfo& basisInfo, const auto& origins, const auto& directions)
{
	static auto integrationPoints = std::vector<double>{ 0.0, 0.5, 1.0 };
	
	auto precalculatedGradXiValues = std::vector<std::vector<std::vector<Vec2D>>>(basisInfo.GetBasisFunctions().size());
	std::for_each(precalculatedGradXiValues.begin(), precalculatedGradXiValues.end(), [&](auto& val) {val.resize(origins.size()); });

	for (auto xiNum = 0u; xiNum < basisInfo.GetBasisFunctions().size(); xiNum++)
	{
		auto xi = basisInfo.GetBasisFunctions()[xiNum];
		for (auto edgeNum = 0u; edgeNum < origins.size(); edgeNum++)
		{
			auto& origin = origins[edgeNum];
			auto& direction = directions[edgeNum];

			for (auto intPointNum = 0u; intPointNum < integrationPoints.size(); intPointNum++)
			{
				auto& integralPoint = integrationPoints[intPointNum];
				auto point = origin + integralPoint * direction;
				auto gradXiAtPoint = gradient(xi, point);
				precalculatedGradXiValues[xiNum][edgeNum].push_back(gradXiAtPoint);
			}
		}
	}

	return precalculatedGradXiValues;
}

std::vector<std::vector<Vec2D>> FluxCalculator::GetNormals()
{
	std::vector<std::vector<Vec2D>> normals;
	switch (m_Mode)
	{
	case BasisType::None:
		break;
	case BasisType::Bilinear:
	case BasisType::Biquadratic:
		normals = { 
			{
				{0.0, 1.0},
				{1.0, 0.0},
				{1.0, 0.0},
				{0.0, 1.0} 
			},
			{
				{0.0, 1.0},
				{1.0, 0.0},
				{1.0, 0.0},
				{0.0, 1.0}
			}
		};
		break;

	case BasisType::Linear:
	case BasisType::Quadratic:
		normals = {
			{
				{0.0, 1.0},
				{1.0, 0.0},
				{1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0)},
			},
			{
				{1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0)},
				{1.0, 0.0},
				{0.0, 1.0}
			}
		};
		break;
	default:
		break;
	}
	return normals;
}

std::vector<std::vector<Point2D>> FluxCalculator::GetEdgesOrigins()
{
	switch (m_Mode)
	{
	case BasisType::None:
		break;
	case BasisType::Bilinear:
	case BasisType::Biquadratic:
		return {
			{ {0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}, {1.0, 0.0} },
			{ {0.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}, {1.0, 0.0} }
		};
		break;

	case BasisType::Linear:
	case BasisType::Quadratic:
		return {
			{ {0.0, 0.0}, {0.0, 0.0}, {0.0, 1.0} },
			{ {0.0, 1.0}, {0.0, 0.0}, {0.0, 0.0} } // Интергируем по локальному треугольнику
		};
		break;

	default:
		break;
	}

}

std::vector<std::vector<Vec2D>> FluxCalculator::GetEdgesDirections()
{
	switch (m_Mode)
	{
	case BasisType::None:
		break;
	case BasisType::Bilinear:
	case BasisType::Biquadratic:
		return {
			{ {1.0, 0.0}, {0.0, 1.0}, {0.0, 1.0}, {1.0, 0.0} },
			{ {1.0, 0.0}, {0.0, 1.0}, {0.0, 1.0}, {1.0, 0.0} }
		};
		break;

	case BasisType::Linear:
	case BasisType::Quadratic:
		return {
			{ {0.0,  1.0}, {0.0, 1.0}, {1.0, -1.0} },
			{ {1.0, -1.0}, {0.0, 1.0}, {1.0,  0.0} }
		};
		break;

	default:
		break;
	}
}

std::vector<double> FluxCalculator::JacobianModifier()
{
	switch (m_Mode)
	{
	case BasisType::None:
		break;
	case BasisType::Bilinear:
	case BasisType::Biquadratic:
		return { 1.0, 1.0 };
		break;

	case BasisType::Linear:
	case BasisType::Quadratic:
		return { 1.0, -1.0 };
		break;

	default:
		break;
	}
}
