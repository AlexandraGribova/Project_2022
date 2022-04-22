#pragma once

#include <vector>

class ApproximationInfo
{
	using Matrix = std::vector<std::vector<double>>;
public:

	ApproximationInfo() = default;
	ApproximationInfo(const std::vector<Matrix>& stfMatrices, const Matrix& msMatrices, double stfCoeff, double msCoeff, const std::vector<double>& nmCoeffs);
	~ApproximationInfo() {}

	void AddLocalStiffnessMatrix(const Matrix& stiffnessMatrix) { m_LocalStiffnessMatrices.push_back(stiffnessMatrix); }
	void SetLocalMassMatrix(const Matrix& massMatrix) { m_LocalMassMatrix = massMatrix; }
	auto& GetLocalStiffnessMatrices() { return m_LocalStiffnessMatrices; }
	auto& GetLocalMassMatrix() { return m_LocalMassMatrix; }

	void SetStiffnessCoefficient(double coeff) { m_StiffnessCoefficient = coeff; }
	auto GetStiffnessCoefficient() { return m_StiffnessCoefficient; }
	auto GetMassCoefficient() { return m_MassCoefficient; }
	auto& GetNeumannCoefficients() { return m_NeumannCoefficients; }

private:
	std::vector<Matrix> m_LocalStiffnessMatrices;
	Matrix m_LocalMassMatrix;

	double m_StiffnessCoefficient = 1.0;
	double m_MassCoefficient = 1.0;
	std::vector<double> m_NeumannCoefficients;

};

namespace ApprInfo
{
	using Matrix = std::vector<std::vector<double>>;
	const ApproximationInfo BilinearInfo(
		std::vector<Matrix>
		{
			{
				{  2.0, -2.0,  1.0, -1.0 },
				{ -2.0,  2.0, -1.0,  1.0 },
				{  1.0, -1.0,  2.0, -2.0 },
				{ -1.0,  1.0, -2.0,  2.0 }
			},
			{
				{  2.0,  1.0, -2.0, -1.0 },
				{  1.0,  2.0, -1.0, -2.0 },
				{ -2.0, -1.0,  2.0,  1.0 },
				{ -1.0, -2.0,  1.0,  2.0 }
			}							 
		},
		Matrix
		{
			{ 4.0, 2.0, 2.0, 1.0 },
			{ 2.0, 4.0, 1.0, 2.0 },
			{ 2.0, 1.0, 4.0, 2.0 },
			{ 1.0, 2.0, 2.0, 4.0 },
		},
		1.0 / 6.0, 1.0 / 36.0, {3.0 / 6.0, 3.0 / 6.0}
	);
	const ApproximationInfo BiquadraticInfo(
		std::vector<Matrix>
		{
			{
				{  28.0,  -32.0,   4.0,   14.0,  -16.0,    2.0,   -7.0,    8.0,   -1.0 },
				{ -32.0,   64.0, -32.0,  -16.0,   32.0,  -16.0,    8.0,  -16.0,    8.0 },
				{   4.0,  -32.0,  28.0,    2.0,  -16.0,   14.0,   -1.0,    8.0,   -7.0 },
				{  14.0,  -16.0,   2.0,  112.0, -128.0,   16.0,   14.0,  -16.0,    2.0 },
				{ -16.0,   32.0, -16.0, -128.0,  256.0, -128.0,  -16.0,   32.0,  -16.0 },
				{   2.0,  -16.0,  14.0,   16.0, -128.0,  112.0,    2.0,  -16.0,   14.0 },
				{  -7.0,    8.0,  -1.0,   14.0,  -16.0,    2.0,   28.0,  -32.0,    4.0 },
				{   8.0,  -16.0,   8.0,  -16.0,   32.0,  -16.0,  -32.0,   64.0,  -32.0 },
				{  -1.0,    8.0,  -7.0,    2.0,  -16.0,   14.0,    4.0,  -32.0,   28.0 }
			},
			{
				{  28.0,   14.0,  -7.0,  -32.0,  -16.0,    8.0,    4.0,    2.0,   -1.0 },
				{  14.0,  112.0,  14.0,  -16.0, -128.0,  -16.0,    2.0,   16.0,    2.0 },
				{  -7.0,   14.0,  28.0,    8.0,  -16.0,  -32.0,   -1.0,    2.0,    4.0 },
				{ -32.0,  -16.0,   8.0,   64.0,   32.0,  -16.0,  -32.0,  -16.0,    8.0 },
				{ -16.0, -128.0, -16.0,   32.0,  256.0,   32.0,  -16.0, -128.0,  -16.0 },
				{   8.0,  -16.0, -32.0,  -16.0,   32.0,   64.0,    8.0,  -16.0,  -32.0 },
				{   4.0,    2.0,  -1.0,  -32.0,  -16.0,    8.0,   28.0,   14.0,   -7.0 },
				{   2.0,   16.0,   2.0,  -16.0, -128.0,  -16.0,   14.0,  112.0,   14.0 },
				{  -1.0,    2.0,   4.0,    8.0,  -16.0,  -32.0,   -7.0,   14.0,   28.0 },
			}
		},
		Matrix
		{
			{ 16.0, 8.0, -4.0, 8.0, 4.0, -2.0, -4.0, -2.0, 1 },
			{ 8.0, 64.0, 8.0, 4.0, 32.0, 4.0, -2.0, -16.0, -2 },
			{ -4.0, 8.0, 16.0, -2.0, 4.0, 8.0, 1.0, -2.0, -4 },
			{ 8.0, 4.0, -2.0, 64.0, 32.0, -16.0, 8.0, 4.0, -2 },
			{ 4.0, 32.0, 4.0, 32.0, 256.0, 32.0, 4.0, 32.0, 4 },
			{ -2.0, 4.0, 8.0, -16.0, 32.0, 64.0, -2.0, 4.0, 8 },
			{ -4.0, -2.0, 1.0, 8.0, 4.0, -2.0, 16.0, 8.0, -4 },
			{ -2.0, -16.0, -2.0, 4.0, 32.0, 4.0, 8.0, 64.0, 8 },
			{ 1.0, -2.0, -4.0, -2.0, 4.0, 8.0, -4.0, 8.0, 16 }
		},
		1.0 / 90.0, 1.0 / 900.0, { 5.0 / 30.0, 20.0 / 30.0, 5.0 / 30.0 }
	);
	const ApproximationInfo QuadraticInfo(
		std::vector<Matrix>
		{
			{
				{  3.0,  0.0,  1.0,  0.0,  0.0, -4.0 },
				{  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
				{  1.0,  0.0,  3.0,  0.0,  0.0, -4.0 },
				{  0.0,  0.0,  0.0,  8.0, -8.0,  0.0 },
				{  0.0,  0.0,  0.0, -8.0,  8.0,  0.0 },
				{ -4.0,  0.0, -4.0,  0.0,  0.0,  8.0 },
			},
			{
				{  3.0,  1.0,  0.0, -4.0,  0.0,  0.0 },
				{  1.0,  3.0,  0.0, -4.0,  0.0,  0.0 },
				{  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
				{ -4.0, -4.0,  0.0,  8.0,  0.0,  0.0 },
				{  0.0,  0.0,  0.0,  0.0,  8.0, -8.0 },
				{  0.0,  0.0,  0.0,  0.0, -8.0,  8.0 },
			}
		},
		Matrix
		{
			{  6.0, -1.0, -1.0,  0.0,  -4.0,   0.0 },
			{ -1.0,  6.0, -1.0,  0.0,   0.0,  -4.0 },
			{ -1.0, -1.0,  6.0, -4.0,   0.0,   0.0 },
			{  0.0,  0.0, -4.0, 32.0,  16.0,  16.0 },
			{ -4.0,  0.0,  0.0, 16.0,  32.0,  16.0 },
			{  0.0, -4.0,  0.0, 16.0,  16.0,  32.0 }
		},
		1.0 / 6.0, 1.0 / 360.0, { 5.0 / 30.0, 20.0 / 30.0, 5.0 / 30.0}
	);
	const ApproximationInfo LinearInfo(
		std::vector<Matrix>
		{
			{
				{  1.0, -1.0,  0.0 },
				{ -1.0,  1.0,  0.0 },
				{  0.0,  0.0,  0.0 }
			},
			{
				{  1.0,  0.0, -1.0 },
				{  0.0,  0.0,  0.0 },
				{ -1.0,  0.0,  1.0 }
			},
		},
		Matrix
		{
			{ 2.0, 1.0, 1.0 },
			{ 1.0, 2.0, 1.0 },
			{ 1.0, 1.0, 2.0 }
		},
		1.0 / 2.0, 1.0 / 24.0, { 3.0 / 6.0, 3.0 / 6.0 }
	);
}


