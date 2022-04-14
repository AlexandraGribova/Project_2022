#include "ApproximationInfo.h"

ApproximationInfo::ApproximationInfo(const std::vector<Matrix>& stfMatrices, const Matrix& msMatrices, double stfCoeff, double msCoeff, const std::vector<double>& nmCoeffs) :
	m_LocalStiffnessMatrices(stfMatrices), m_LocalMassMatrix(msMatrices),
	m_StiffnessCoefficient(stfCoeff), m_MassCoefficient(msCoeff), m_NeumannCoefficients(nmCoeffs)
{
}
