#pragma once

#include "DataStructures.h"

inline double dfdx(Fn2d function, const Point2D& point)
{
	auto step = sqrt(std::numeric_limits<double>::epsilon());
	// In speculative execution we trust.
	if (point.X != 0.0) step *= std::abs(point.X);
	return (function(point.X + step, point.Y) - function(point.X - step, point.Y)) / (2.0 * step);
}

inline double dfdy(Fn2d function, const Point2D& point)
{
	auto step = sqrt(std::numeric_limits<double>::epsilon());
	// In speculative execution we trust.
	if (point.Y != 0.0) step *= std::abs(point.Y);
	return (function(point.X, point.Y + step) - function(point.X, point.Y - step)) / (2.0 * step);
}

inline Vec2D gradient(Fn2d function, const Point2D& point)
{
	return Vec2D{ dfdx(function, point), dfdy(function, point) };
}

inline double dotProduct(const std::vector<double> a, const std::vector<double> b)
{
	double res = 0.0;
	for (auto i = 0u; i < a.size(); i++)
	{
		res += a[i] * b[i];
	}
	return res;
}