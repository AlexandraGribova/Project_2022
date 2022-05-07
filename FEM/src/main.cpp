#include <algorithm>
#include <iostream>
#include <functional>
#include <array>
#include <utility>
#include <iomanip>
#include <fstream>
#include <cassert>

#include "output.h"
#include "DataStructures.h"
#include "functions.h"
#include "SparseMatrix.h"
#include "GridData.h"
#include "FEMSolver.h"

double testFunction(double x, double y) { return x+y; }
// Ёто еЄ лапласиан P
double rhsFunction(double x, double y) { return 0; }


int main()
{
	FEMSolver FemSolver("input/FEMMode.txt"); // BL - билинейные, BQ - биквадратические, L - линейные, Q - квадратические.
	FemSolver.LoadPlastPressure("input/Plast.txt");


	GridData gridData;
	gridData.BoundaryFunction = [&](double x, double y) { return testFunction(x, y); };
	//gridData.BoundaryFunction = [&](double x, double y) { return FemSolver.GetPlastPressure(); };
	gridData.RhsFunction = [](double x, double y) { return rhsFunction(x, y); };
	gridData.LoadNodes("input/grid.txt");
	gridData.LoadDirichletConditions("input/BC1.txt");
	if (FemSolver.GetMode() == BasisType::Linear || FemSolver.GetMode() == BasisType::Quadratic)
	{
		gridData.LoadElements("input/outTriangle.txt");
		gridData.LoadNeumannConditions("input/BC2_Tri.txt");
	}
	else
	{
		gridData.LoadElements("input/out.txt");
		gridData.LoadNeumannConditions("input/BC2.txt");
	}
	gridData.SetupEdges();

	FemSolver.SetGridData(gridData);
	//FemSolver.LoadDomainData("input/mat.txt");
	//FemSolver.LoadViscosity("input/phaseprop.txt");
	FemSolver.CalculateStiffnesMatrix();
	FemSolver.CalculateLoadVector();
	FemSolver.ApplyNeumann();
	FemSolver.ApplyDirichlet();
	FemSolver.CalculateSolution();
	auto& solution = FemSolver.GetSolution();
	FemSolver.CalculateFlux();
	std::vector<double> flux = FemSolver.GetFlux();

	gridData.flux_balancer(flux);///

	std::ofstream resout("output/result.txt");

	// ≈сли не указать 4 агрумент (известную функцию), то просто выведетс€ ответ.
	prettyPrint(resout, solution, FemSolver.GetGridData());
	//prettyPrint(std::cout, solution, FemSolver.GetGridData());
}

