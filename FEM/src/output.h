#pragma once
#include <fstream>
#include <vector>
#include <format>

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
		outputStream << std::format("+{:-^88}+\n", "");
		outputStream << std::format("|{:^9}|{:^9}|{:^22}|{:^22}|{:^22}|\n", "X", "Y", "P", "T", "P - T");
		outputStream << std::format("+{:-^88}+\n", "");
	}
	else
	{
		outputStream << std::format("+{:-^42}+\n", "");
		outputStream << std::format("|{:^9}|{:^9}|{:^22}|\n", "X", "Y", "P");
		outputStream << std::format("+{:-^42}+\n", "");
	}
	auto nx = gridData.ModifiedNx;
	auto step = nx == -1 ? 1 : 2;

	for (uint32_t i = 0; i < nodes.size(); i += step)
	{
		outputStream << "|";
		outputStream << std::format("{:>9.3f}|{:>9.3f}|{:>22.15e}|", nodes[i].X, nodes[i].Y, res[i]);
		if (function)
		{
			outputStream << std::format("{:>22.15e}|{:>22.15e}|\n", function(nodes[i].X, nodes[i].Y), function(nodes[i].X, nodes[i].Y) - res[i]);
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
		outputStream << std::format("+{:-^88}+\n", "");
		outputStream << std::format("|{:^19}|{:>68.15e}|\n", "||P-T||/||T||", rel);
		outputStream << std::format("+{:-^88}+\n", "");
	}
	else 
	{
		outputStream << std::format("+{:-^42}+\n", "");
	}
}