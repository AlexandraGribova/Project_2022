#pragma once

#include "DataStructures.h"

class GridData
{
	using Fn2d = std::function<double(double, double)>;
	using FnInd = std::function<uint32_t(uint32_t i, uint32_t j)>;
public:
	GridData();
	void LoadElements(const std::string& path);
	void LoadNodes(const std::string& path);
	void LoadDirichletConditions(const std::string& path);
	void LoadNeumannConditions(const std::string& path);

	void ModifyTriangleGrid();
	void ModifyGrid();
	void Triangulate();
	void SetupEdges();


	void ig_jg_generation();
	int ig_creation(int elem1, int elem2, int edge, std::vector<int> &jg);
	std::vector<int32_t> GetNumberEdge(uint32_t numedge);

	std::vector<Element> Elements;
	std::vector<NodeIDs> LL;
	std::vector<Point2D> Nodes;

	Fn2d RhsFunction;
	Fn2d BoundaryFunction;

	FnInd IndexingFunction;

	std::vector<DirichletData> DirichletConditions;
	std::vector<NeumannData> NeumannConditions;

	int32_t ModifiedNx = -1;
	bool HasReferenceGrid = false;

private:
	void UpdateDirichlet(const std::vector<Point2D>& newNodes);
	void UpdateNeumann();

};
