#pragma once
#include <cinttypes>
#include <vector>
#include <functional>

typedef std::vector<uint32_t> Element;

struct Point2D
{
	Point2D() = default;
	Point2D(double x, double y) :
		X(x), Y(y)
	{}
	double X, Y;
};

struct DirichletData
{
	using Fn2d = std::function<double(double, double)>;
	DirichletData() = default;
	DirichletData(uint32_t node, Fn2d value) :
		Node(node), Value(value)
	{}
	uint32_t Node;
	Fn2d Value; 
};

struct NeumannData
{
	NeumannData() = default;
	NeumannData(uint32_t element, const std::vector<uint32_t>& nodes, double theta) :
		Element(element), Nodes(nodes), Theta(theta)
	{}
	uint32_t Element;
	std::vector<uint32_t> Nodes;
	double Theta;
};
