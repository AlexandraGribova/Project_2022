#pragma once

#include <vector>
#include <functional>

enum class BasisType
{
	None = 0,
	Bilinear, Biquadratic,
	Linear, Quadratic
};

class BasisInfo
{
public:
	using Fn2d = std::function<double(double, double)>;
	BasisInfo() = default;
	BasisInfo(const std::vector<Fn2d>& basisFunctions) :
		m_BasisFunctions(basisFunctions) {}
	
	auto GetBasisFunctions() const { return m_BasisFunctions; }
private:
	std::vector<Fn2d> m_BasisFunctions;
};

namespace BsInfo
{
	using Fn2d = std::function<double(double, double)>;
	const BasisInfo Billinear(
		[]() constexpr -> std::vector<Fn2d>
		{
			// y1(y) = x1(y), y2(y) = x2(y)
			auto x1 = [](double x) { return 1.0 - x; };
			auto x2 = [](double x) { return x; };

			auto xi1 = [&](double x, double y) { return x1(x) * x1(y); };
			auto xi2 = [&](double x, double y) { return x2(x) * x1(y); };
			auto xi3 = [&](double x, double y) { return x1(x) * x2(y); };
			auto xi4 = [&](double x, double y) { return x2(x) * x2(y); };

			return { xi1, xi2, xi3, xi4 };
		}()
	);
	const BasisInfo Biquadratic(
		[]() constexpr -> std::vector<Fn2d>
		{
			// y1(y) = x1(y), y2(y) = x2(y), y3(y) = x3(y)
			auto x1 = [](double x) { return 2.0 * x * x - 3.0 * x + 1; };
			auto x2 = [](double x) { return -4.0 * x * x + 4.0 * x; };
			auto x3 = [](double x) { return 2.0 * x * x - x; };

			auto xi1 = [&](double x, double y) { return x1(x) * x1(y); };
			auto xi2 = [&](double x, double y) { return x2(x) * x1(y); };
			auto xi3 = [&](double x, double y) { return x3(x) * x1(y); };
			auto xi4 = [&](double x, double y) { return x1(x) * x2(y); };
			auto xi5 = [&](double x, double y) { return x2(x) * x2(y); };
			auto xi6 = [&](double x, double y) { return x3(x) * x2(y); };
			auto xi7 = [&](double x, double y) { return x1(x) * x3(y); };
			auto xi8 = [&](double x, double y) { return x2(x) * x3(y); };
			auto xi9 = [&](double x, double y) { return x3(x) * x3(y); };

			return { xi1, xi2, xi3, xi4, xi5, xi6, xi7, xi8, xi9 };
		}()
	);
	const BasisInfo Quadratic(
		[]() constexpr -> std::vector<Fn2d>
		{
			auto l1 = [&](double x, double y) { return 1.0 - x - y; };
			auto l2 = [&](double x, double y) { return x; };
			auto l3 = [&](double x, double y) { return y; };

			auto xi1 = [&](double x, double y) { auto val = l1(x, y); return val * (2.0 * val - 1.0); };
			auto xi2 = [&](double x, double y) { auto val = l3(x, y); return val * (2.0 * val - 1.0);  };
			auto xi3 = [&](double x, double y) { auto val = l2(x, y); return val * (2.0 * val - 1.0);  };
			auto xi4 = [&](double x, double y) { return 4.0 * l1(x, y) * l3(x, y); };
			auto xi5 = [&](double x, double y) { return 4.0 * l2(x, y) * l3(x, y); };
			auto xi6 = [&](double x, double y) { return 4.0 * l1(x, y) * l2(x, y); };

			return { xi1, xi2, xi3, xi4, xi5, xi6 };
		}()
	);
	const BasisInfo Linear(
		[]() constexpr -> std::vector<Fn2d>
		{
			auto xi1 = [&](double x, double y) { return 1.0 - x - y; };
			auto xi2 = [&](double x, double y) { return x; };
			auto xi3 = [&](double x, double y) { return y; };
			return { xi1, xi2, xi3 };
		}()
	);
}

