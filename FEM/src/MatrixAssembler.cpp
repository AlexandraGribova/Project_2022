#include "MatrixAssembler.h"
#include <iostream>

std::unique_ptr<SparseMatrix> MatrixAssembler::AssembleMatrix(const GridData& gridData, std::vector<uint32_t>& ig, std::vector<uint32_t>& jg, Fn2d additionalCoeffs)
{
	std::vector<double> di(gridData.Nodes.size(), 0.0);
	std::vector<double> ggl(jg.size(), 0.0);
	auto localMatrixSize = m_ApproximationInfo.GetLocalStiffnessMatrices()[0].size();
	auto indF = gridData.IndexingFunction;
	for (uint32_t elN = 0; elN < gridData.Elements.size(); elN++)
	{
		auto element = gridData.Elements[elN];
		auto&& [x, y] = getCoords(gridData.Nodes, element);
		//std::cout << elN << ": " << "{ " << gridData.Nodes[element[0]].X << " " << gridData.Nodes[element[0]].Y << ", "
		//	<< gridData.Nodes[element[1]].X << " " << gridData.Nodes[element[1]].Y << ", "
		//	<< gridData.Nodes[element[2]].X << " " << gridData.Nodes[element[2]].Y << ", "
		//	<< gridData.Nodes[element[3]].X << " " << gridData.Nodes[element[3]].Y << " }\n";
		auto&& [hx, hy] = getDimensions(gridData.Nodes, element);
		double coeffs[2] = { hy / hx, hx / hy };

		auto localA = m_ApproximationInfo.GetLocalStiffnessMatrices()[0];
		auto localB = m_ApproximationInfo.GetLocalStiffnessMatrices()[1];

		for (uint32_t i = 0; i < localMatrixSize; i++)
		{
			auto value = localA[indF(elN, i)][indF(elN, i)] * coeffs[0] + localB[indF(elN, i)][indF(elN, i)] * coeffs[1];
			value *= m_ApproximationInfo.GetStiffnessCoefficient() * additionalCoeffs(x, y);
			di[element[i]] += value;
		}

		for (uint32_t i = 0; i < localMatrixSize; i++)
		{
			auto ibeg = ig[element[i]];
			auto iend = ig[element[i] + 1] - 1;
			for (uint32_t j = 0; j < i; j++)
			{
				auto index = binarySearch(jg, element[j], ibeg, iend);
				auto value = localA[indF(elN, i)][indF(elN, j)] * coeffs[0] + localB[indF(elN, i)][indF(elN, j)] * coeffs[1];
				value *= m_ApproximationInfo.GetStiffnessCoefficient() * additionalCoeffs(x, y);
				ggl[index] += value;
				ibeg++;
			}
		}
	}
	return std::make_unique<SparseMatrix>(ig, jg, di, ggl);
}

std::vector<double> MatrixAssembler::AssembleLoadVector(const GridData& gridData)
{
	std::vector<double> loadVector(gridData.Nodes.size(), 0.0);
	auto indF = gridData.IndexingFunction;

	for (uint32_t elN = 0; elN < gridData.Elements.size(); elN++)
	{
		auto element = gridData.Elements[elN];
		auto area = getArea(gridData.Nodes, element); // на самом деле это не площадь, а определитель

		auto numberOfVertices = element.size();
		std::vector<double> rhsFunctionAtNodes(numberOfVertices, 0.0);

		for (uint32_t i = 0; i < numberOfVertices; i++)
		{
			auto index = indF(elN, i);
			rhsFunctionAtNodes[index] = gridData.RhsFunction(gridData.Nodes[element[i]].X, gridData.Nodes[element[i]].Y);
		}

		auto massMatrix = m_ApproximationInfo.GetLocalMassMatrix();

		for (uint32_t i = 0; i < numberOfVertices; i++)
		{
			auto value = -area * dotProduct(rhsFunctionAtNodes, massMatrix[indF(elN,i)]);
			value *= m_ApproximationInfo.GetMassCoefficient();
			loadVector[element[i]] += value;
		}
	}
	return loadVector;
}

void MatrixAssembler::ApplyDirichlet(const GridData& gridData, std::vector<double>& loadVector)
{
	for (const auto& dirichletData : gridData.DirichletConditions)
	{
		auto&& [index, value] = dirichletData;
		loadVector[index] = value(gridData.Nodes[index].X, gridData.Nodes[index].Y);
	}
}

void MatrixAssembler::ApplyNeumann(const GridData& gridData, std::vector<double>& loadVector)
{
	for (const auto& neumannData : gridData.NeumannConditions)
	{
		auto element = gridData.Elements[neumannData.Element];
		auto nodes = neumannData.Nodes;
		auto first = element[nodes[0]];
		auto last = element[nodes[nodes.size() - 1]];
		double length = 0.0;

		length = abs((gridData.Nodes[last].X - gridData.Nodes[first].X) + (gridData.Nodes[last].Y - gridData.Nodes[first].Y));

		auto index = 0;
		for (const auto& node : nodes)
		{
			auto value = length * neumannData.Theta * m_ApproximationInfo.GetNeumannCoefficients()[index++];
			loadVector[element[node]] += value;
		}		
	}
}
