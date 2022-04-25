#include "FluxCalculator.h"
#include "functions.h"
#include "Math.h"

void FluxCalculator::CalculateFlux(const BasisInfo& basisInfo, const GridData& grid, const std::vector<double>& solution)
{
	m_Fluxes = std::vector<double>(grid.Elements.back().Edges.back() + 1, 0.0);

	static auto normals = GetNormals();
	static auto origins = GetEdgesOrigins();
	static auto directions = GetEdgesDirections();

	static auto jacobiMod = JacobianModifier();

	static auto integrationPoints = std::vector<double>{ 0.0, 0.5, 1.0 };
	static auto integrationCoefficients = std::vector<double>{ 1.0, 4.0, 1.0 };

	static auto indF = grid.IndexingFunction;

	// ƒл€ треугольников будет сложнее, потому что они разные, пока не буду специально запрашивать грани.
	// ≈ще мы считаем интеграл по "локальной грани" и у "перевернутого" треугольника это приводит к 
	// отрицательному €кобиану.
	for (auto elemNum = 0u; elemNum < grid.Elements.size(); elemNum++)
	{
		auto& elem = grid.Elements[elemNum];
		auto elemDisparity = elemNum & 1; // Ќечетность номера элемента дл€ €кобиана перевернутого треугольника (disparity - это не нечетность, но мне показалось смешно)
		auto&& [width, height] = getDimensions(grid.Nodes, elem.Nodes);
		auto J = width * height;
		auto JInv = std::vector<double>{ 1.0 / width, 1.0 / height }; // ’раним только диагональ.
		for (auto edgeNum = 0u; edgeNum < elem.Edges.size(); edgeNum++)
		{
			auto origin = origins[elemDisparity][edgeNum];
			auto direction = directions[elemDisparity][edgeNum];
			auto normal = normals[elemDisparity][edgeNum];
			auto edgeFlux = 0.0;
			// xi - это типо кси
			for (auto xiNum = 0; xiNum < basisInfo.GetBasisFunctions().size(); xiNum++)
			{
				auto gradIntegral = Vec2D{ 0.0, 0.0 };
				auto xi = basisInfo.GetBasisFunctions()[indF(elemNum, xiNum)];
				for (auto coeffNum = 0u; coeffNum < integrationCoefficients.size(); coeffNum++)
				{
					auto point = origin + integrationPoints[coeffNum] * direction;
					auto gradXi = gradient(xi, point);
					// Ёто можно очень сильно оптимизировать, т.к. градиент считаетс€ все врем€ в одних точках.
					// » на JInv можно тоже в самом конце домножать, оставлю пока так, чтобы пон€тнее было.
					gradIntegral.X += integrationCoefficients[coeffNum] * gradXi.X * JInv[0] * jacobiMod[elemDisparity];
					gradIntegral.Y += integrationCoefficients[coeffNum] * gradXi.Y * JInv[1] * jacobiMod[elemDisparity];
				}
				edgeFlux += solution[elem.Nodes[xiNum]] * (gradIntegral.X * normal.X + gradIntegral.Y * normal.Y);
			}
			edgeFlux *= elem.EdgeDirections[edgeNum] * J / 6.0;
			// Ќаверное это не лучший способ.
			m_Fluxes[elem.Edges[edgeNum]] = m_Fluxes[elem.Edges[edgeNum]] == 0.0 ? edgeFlux : (edgeFlux + m_Fluxes[elem.Edges[edgeNum]]) / 2.0;
		}
	}
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
			{ {0.0, 1.0}, {0.0, 0.0}, {0.0, 0.0} } // »нтергируем по локальному треугольнику
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
