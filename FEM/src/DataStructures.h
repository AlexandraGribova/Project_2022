#pragma once
#include <cinttypes>
#include <vector>
#include <functional>

using NodeIDs = std::vector<uint32_t>;
using Fn2d = std::function<double(double, double)>;

struct Element
{
	std::vector<uint32_t> Nodes;
	std::vector<uint32_t> Edges;
	std::vector<int32_t> EdgeDirections;
};

struct Vec2D
{
	Vec2D() = default;
	Vec2D(double x, double y) :
		X(x), Y(y)
	{}

	Vec2D operator-() const { return Vec2D{ -X, -Y }; };
	Vec2D operator*(double val) const { return { X * val, Y * val }; };
	Vec2D operator/(double val) const { return *this * (1.0 / val); };
	Vec2D operator+(const Vec2D& other) const { return { X + other.X, Y + other.Y }; };
	Vec2D operator-(const Vec2D& other) const { return { X - other.X, Y - other.Y }; };

	Vec2D& operator+=(const Vec2D& other) { X += other.X, Y += other.Y; return *this; }
	Vec2D& operator-=(const Vec2D& other) { X -= other.X, Y -= other.Y; return *this; }

	Vec2D& operator*=(double val) { X *= val, Y *= val; return *this; }
	Vec2D& operator/=(double val) { return *this *= (1.0 / val); }

	double X, Y;
};

inline Vec2D operator*(double c, const Vec2D& vec2) { return vec2 * c; };
inline Vec2D operator/(double c, const Vec2D& vec2) { return vec2 / c; };

using Point2D = Vec2D;

struct DirichletData
{
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
