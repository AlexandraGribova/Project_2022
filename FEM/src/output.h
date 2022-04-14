#pragma once
#include <fstream>
#include <vector>

#include "DataStructures.h"
#include "functions.h"
#include "GridData.h"

template <typename T>
inline void outputToFile(const std::string& path, const std::vector<T>& vector)
{
	std::ofstream out(path);
	for (const auto& el : vector)
	{
		out << el << "\n";
	}
	out.close();
}

template <typename T>
inline void outputToFile(const std::string& path, const std::vector<std::vector<T>>& matrix)
{
	std::ofstream out(path);
	for (const auto& row : matrix)
	{
		for (const auto& el : row)
		{
			out << std::setw(9) << std::setprecision(3) << el << " ";
		}
		out << "\n";
	}
	out.close();
}

inline void prettyPrint(std::ostream& outputStream,  const std::vector<double>& res, const GridData& gridData, const std::function<double(double, double)>& function = {})
{
	std::vector<Point2D> referenceGrid;
	if (gridData.HasReferenceGrid)
	{
		std::ifstream in("input/referenceGrid.txt");
		uint32_t N;
		in >> N;
		referenceGrid = std::vector<Point2D>(N, { 0, 0 });
		for (uint32_t i = 0; i < N; i++)
		{
			auto x = 0.0;
			auto y = 0.0;
			in >> x >> y;
			referenceGrid[i] = { x, y };
		}
	}


	auto nodes = gridData.Nodes;
	if (function)
	{
		outputStream << "+----------------------------------------------------------------------------------------+\n";
		outputStream << "|    X    |    Y    |           P          |           T          |         P - T        |\n";
		outputStream << "+---------+---------+----------------------+----------------------+----------------------+\n";
	}
	else
	{
		outputStream << "+------------------------------------------+\n";
		outputStream << "|    X    |    Y    |           P          |\n";
		outputStream << "+---------+---------+----------------------+\n";
	}
	auto nx = gridData.ModifiedNx;
	auto step = nx == -1 ? 1 : 2;

	for (uint32_t i = 0; i < nodes.size(); i += step)
	{
		outputStream << "|";
		outputStream << std::fixed << std::setw(9) << std::setprecision(3) << nodes[i].X << "|";
		outputStream << std::fixed << std::setw(9) << std::setprecision(3) << nodes[i].Y << "|";
		outputStream << std::scientific << std::setw(22) << std::setprecision(DBL_DIG) << res[i] << "|";
		if (function)
		{
			outputStream << std::scientific << std::setw(22) << std::setprecision(DBL_DIG) << function(nodes[i].X, nodes[i].Y) << "|";
			outputStream << std::scientific << std::setw(22) << std::setprecision(DBL_DIG) << function(nodes[i].X, nodes[i].Y) - res[i] << "|\n";
		}
		else
		{
			outputStream << "\n";
		}
		if (i > 0 && nx != -1 && (i + 1) % nx == 0)
		{
			i += nx - 1;
		}
	}
	uint32_t j = 0;
	if (function)
	{
		auto sum1 = 0.0;
		auto sum2 = 0.0;
		for (uint32_t k = 0; k < nodes.size(); k += step)
		{
			if (gridData.HasReferenceGrid)
			{
				if (nodes[k].X == referenceGrid[j].X && nodes[k].Y == referenceGrid[j].Y) j++;
				else continue;
			}
			
			sum1 += (function(nodes[k].X, nodes[k].Y) - res[k]) * (function(nodes[k].X, nodes[k].Y) - res[k]);
			sum2 += function(nodes[k].X, nodes[k].Y) * function(nodes[k].X, nodes[k].Y);
			if (k > 0 && nx != -1 && (k + 1) % nx == 0)
			{
				k += nx - 1;
			}
		}
		auto rel = sqrt(sum1 / sum2);
		rel = rel == rel ? rel : 0.0; // rel != rel --> rel == nan.
		outputStream << "+----------------------------------------------------------------------------------------+\n";
		outputStream << "|   ||P-T||/||T||   |" << std::scientific << std::setw(68) << std::setprecision(DBL_DIG) << rel << "|\n";
		outputStream << "+----------------------------------------------------------------------------------------+\n";
	}
	else 
	{
		outputStream << "+------------------------------------------+\n";
	}
}