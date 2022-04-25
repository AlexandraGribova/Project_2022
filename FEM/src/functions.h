#pragma once

#include "DataStructures.h"
#include "SparseMatrix.h"

std::pair<std::vector<uint32_t>, std::vector<uint32_t>> generatePortrait(uint32_t N, const std::vector<Element>& elements);

inline double getWidth(const std::vector<Point2D>& nodes, const std::vector<uint32_t>& element) 
{
	if (element.size() == 4)
		return nodes[element[1]].X - nodes[element[0]].X;
	if (element.size() == 6)
		return abs(nodes[element[2]].X - nodes[element[0]].X + nodes[element[3]].X - nodes[element[5]].X);
	if (element.size() == 3)
		return abs(nodes[element[1]].X - nodes[element[0]].X);
	return nodes[element[2]].X - nodes[element[0]].X;
}
inline double getHeight(const std::vector<Point2D>& nodes, const std::vector<uint32_t>& element)
{
	if (element.size() == 4)
		return nodes[element[2]].Y - nodes[element[0]].Y;
	if (element.size() == 6)
		return abs(nodes[element[5]].Y - nodes[element[0]].Y);
	if (element.size() == 3)
		return abs(nodes[element[2]].Y - nodes[element[0]].Y);
	return nodes[element[6]].Y - nodes[element[0]].Y;
}
inline std::pair<double, double> getDimensions(const std::vector<Point2D>& nodes, const std::vector<uint32_t>& element)
{
	return { getWidth(nodes, element), getHeight(nodes, element) };
}
inline double getArea(const std::vector<Point2D>& nodes, const std::vector<uint32_t>& element)
{
	return getWidth(nodes, element) * getHeight(nodes, element);
}
inline std::pair<double, double> getCoords(const std::vector<Point2D>& nodes, const std::vector<uint32_t>& element)
{
	auto x = std::min(nodes[element[0]].X, nodes[element[1]].X);
	auto y = nodes[element[0]].Y;
	return { x, y };
}

template <typename T>
inline int32_t binarySearch(const std::vector<T>& values, const T& value, uint32_t left, uint32_t right)
{
	while (left != right)
	{
		auto mid = (left + right) / 2 + 1;
		if (values[mid] > value)
			right = mid - 1;
		else
			left = mid;
	}
	if (values[left] == value)
		return static_cast<int32_t>(left);
	return -1;
}

inline int32_t binarySearch(const std::vector<DirichletData>& values, uint32_t value, uint32_t left, uint32_t right)
{
	while (left != right)
	{
		auto mid = (left + right) / 2 + 1;
		if (values[mid].Node > value)
			right = mid - 1;
		else
			left = mid;
	}
	if (values[left].Node == value)
		return static_cast<int32_t>(left);
	return -1;
}

inline int32_t binarySearchAny(const std::vector<DirichletData>& values, const std::vector<uint32_t>& vals, uint32_t left, uint32_t right)
{
	for (const auto& val : vals)
	{
		if (binarySearch(values, val, left, right) != -1) return -1;
	}
	return 0;
}