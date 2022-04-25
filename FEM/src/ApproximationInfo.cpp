#include "ApproximationInfo.h"

ApproximationInfo::ApproximationInfo(const BasisInfo& basisInfo, const std::vector<Matrix>& stfMatrices, const Matrix& msMatrices, double stfCoeff, double msCoeff, const std::vector<double>& nmCoeffs) :
	m_BasisInfo(basisInfo),
	m_LocalStiffnessMatrices(stfMatrices), m_LocalMassMatrix(msMatrices),
	m_StiffnessCoefficient(stfCoeff), m_MassCoefficient(msCoeff), m_NeumannCoefficients(nmCoeffs)
{
}
