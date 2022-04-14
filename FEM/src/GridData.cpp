#include "GridData.h"

#include <fstream>
#include <array>
#include <iostream>

GridData::GridData()
{
	IndexingFunction = [](uint32_t i, uint32_t j) { return j; };
}

void GridData::LoadElements(const std::string& path)
{
	std::ifstream in(path);
	uint32_t numberOfEls;
	uint32_t numberOfNodes;
	in >> numberOfEls;
	in >> numberOfNodes;
	Elements = std::vector<Element>(numberOfEls, std::vector<uint32_t>(numberOfNodes, {}));
	for (uint32_t i = 0; i < numberOfEls; i++)
	{
		Element el;
		el.reserve(numberOfNodes);
		for (uint32_t j = 0; j < numberOfNodes; j++)
		{
			uint32_t node;
			in >> node;
			el.push_back(node);
		}
		Elements[i] = el;
	}
}

void GridData::LoadNodes(const std::string& path)
{
	std::ifstream in(path);
	uint32_t N;
	bool ref = false;
	in >> N;
	in >> ref;
	HasReferenceGrid = ref;
	Nodes.resize(N);
	for (uint32_t i = 0; i < N; i++)
	{
		auto x = 0.0, y = 0.0;
		in >> x >> y;
		Nodes[i] = { x, y };
	}
}

void GridData::LoadDirichletConditions(const std::string& path)
{
	std::ifstream in(path);
	uint32_t N;
	in >> N;
	DirichletConditions.resize(N);
	for (uint32_t i = 0; i < N; i++)
	{
		uint32_t node;
		in >> node;
		DirichletConditions[i] = { node, BoundaryFunction };
	}
}

void GridData::LoadNeumannConditions(const std::string& path)
{
	std::ifstream in(path);
	uint32_t N;
	in >> N;
	NeumannConditions.resize(N);
	for (uint32_t i = 0; i < N; i++)
	{
		uint32_t el;
		in >> el;
		uint32_t beg, end;
		in >> beg >> end;
		double theta;
		in >> theta;
		NeumannConditions[i] = { el, {beg, end}, theta };
	}
}

void GridData::ModifyTriangleGrid()
{
	/*
	3-------4-------5          10--11--12--13--14 			2-------3         6---7---8
	| \	    | \     |			| \	    | \     |	  		| \	    |		  | \	  |
	|	\   |	\	|	---->	5	6   7	8	9	;		|	\   |  ---->  3	  4   5
	|	  \	|	  \ |			|	  \	|	  \ | 			|	  \	|		  |	    \ |
	0-------1-------2           0---1---2---3---4 			0-------1		  0---1---2
	*/

	auto Nx = Elements[0][2];
	ModifiedNx = 2 * (Elements[0][0] / Nx) * (2 * Nx - 1) + 2 * (Elements[0][0] % Nx) + 2 * Nx - 1;
	auto numberOfNewNodes = 2 * (Elements[Elements.size() - 2][0] / Nx) * (2 * Nx - 1) + 2 * (Elements[Elements.size() - 2][0] % Nx) + 2 * (2 * Nx - 1) + 2 + 1;
	std::vector<Point2D> newNodes(numberOfNewNodes, { 0.0, 0.0 });
	std::vector<Element> newElements(Elements.size(), Element(6, 0));

	for (uint32_t i = 0; i < Elements.size(); i++)
	{
		auto element = Elements[i];
		if (!(i & 1))
		{
			auto Ki = 2 * (element[0] / Nx) * (2 * Nx - 1) + 2 * (element[0] % Nx);
			auto midBottom  = Point2D{ (Nodes[element[1]].X + Nodes[element[0]].X) / 2.0,  Nodes[element[0]].Y };
			auto midLeft	= Point2D{  Nodes[element[0]].X,                              (Nodes[element[2]].Y + Nodes[element[0]].Y) / 2.0 };
			auto center		= Point2D{ (Nodes[element[1]].X + Nodes[element[0]].X) / 2.0, (Nodes[element[2]].Y + Nodes[element[0]].Y) / 2.0 };
			newNodes[Ki] = Nodes[element[0]];
			newNodes[Ki + 1] = midBottom;
			newNodes[Ki + 2] = Nodes[element[1]];
			newNodes[Ki + 2 * Nx - 1] = midLeft;
			newNodes[Ki + 2 * Nx] = center;
			newNodes[Ki + 2 * (2 * Nx - 1)] = Nodes[element[2]];
			newElements[i][0] = Ki;
			newElements[i][1] = Ki + 1;
			newElements[i][2] = Ki + 2;
			newElements[i][3] = Ki + 2 * Nx - 1;
			newElements[i][4] = Ki + 2 * Nx;
			newElements[i][5] = Ki + 2 * (2 * Nx - 1);

		}
		else
		{
			auto Ki = 2 * ((element[0] - 1) / Nx) * (2 * Nx - 1) + 2 * ((element[0] - 1) % Nx);
			auto midRight   = Point2D{  Nodes[element[0]].X,                              (Nodes[element[2]].Y + Nodes[element[0]].Y) / 2.0 };
			auto midTop		= Point2D{ (Nodes[element[1]].X + Nodes[element[2]].X) / 2.0,  Nodes[element[2]].Y };
			newNodes[Ki + 2 * Nx + 1] = midRight;
			newNodes[Ki + 2 * (2 * Nx - 1) + 1] = midTop;
			newNodes[Ki + 2 * (2 * Nx - 1) + 2] = Nodes[element[2]];
			newElements[i][0] = Ki + 2;
			newElements[i][1] = Ki + 2 * Nx;
			newElements[i][2] = Ki + 2 * Nx + 1;
			newElements[i][3] = Ki + 2 * (2 * Nx - 1);
			newElements[i][4] = Ki + 2 * (2 * Nx - 1) + 1;
			newElements[i][5] = Ki + 2 * (2 * Nx - 1) + 2;
		}
	}
	
	UpdateDirichlet(newNodes);

	UpdateNeumann();

	Elements = std::move(newElements);
	Nodes = std::move(newNodes);
}

void GridData::ModifyGrid()
{
	/*
	3-------4-------5             10---11--12---13--14				  2-------3              6---7---8
	|		|		|              |		|		 |				  |		  |				 |		 |
	|		|		|	---->	   5	6	7	 8	 9		;		  |		  |	   ---->	 3	 4	 5
	|		|		|		       |		|		 |				  |		  |				 |		 |
	0-------1-------2              0----1---2----3---4				  0-------1				 0---1---2
	*/

	auto Nx = Elements[0][2];
	ModifiedNx = 2 * (Elements[0][0] / Nx) * (2 * Nx - 1) + 2 * (Elements[0][0] % Nx) + 2 * Nx - 1;
	auto numberOfNewNodes = 2 * (Elements[Elements.size() - 1][0] / Nx) * (2 * Nx - 1) + 2 * (Elements[Elements.size() - 1][0] % Nx) + 2 * (2 * Nx - 1) + 2 + 1;
	std::vector<Point2D> newNodes(numberOfNewNodes, { 0.0, 0.0 });
	std::vector<Element> newElements(Elements.size(), Element(9, 0));

	for (uint32_t i = 0; i < Elements.size(); i++)
	{
		auto element = Elements[i];
		auto midBottom  = Point2D{ (Nodes[element[1]].X + Nodes[element[0]].X) / 2.0,  Nodes[element[0]].Y };
		auto midLeft	= Point2D{  Nodes[element[0]].X,                              (Nodes[element[2]].Y + Nodes[element[0]].Y) / 2.0 };
		auto center		= Point2D{ (Nodes[element[1]].X + Nodes[element[0]].X) / 2.0, (Nodes[element[2]].Y + Nodes[element[0]].Y) / 2.0 };
		auto midRight	= Point2D{  Nodes[element[1]].X,                              (Nodes[element[2]].Y + Nodes[element[0]].Y) / 2.0 };
		auto midTop		= Point2D{ (Nodes[element[1]].X + Nodes[element[0]].X) / 2.0,  Nodes[element[2]].Y };
		auto Ki = 2 * (element[0] / Nx) * (2 * Nx - 1) + 2 * (element[0] % Nx);
		newNodes[Ki] = Nodes[element[0]];
		newNodes[Ki + 1] = midBottom;
		newNodes[Ki + 2] = Nodes[element[1]];
		newNodes[Ki + 2 * Nx - 1] = midLeft;
		newNodes[Ki + 2 * Nx] = center;
		newNodes[Ki + 2 * Nx + 1] = midRight;
		newNodes[Ki + 2 * (2 * Nx - 1)] = Nodes[element[2]];
		newNodes[Ki + 2 * (2 * Nx - 1) + 1] = midTop;
		newNodes[Ki + 2 * (2 * Nx - 1) + 2] = Nodes[element[3]];
		for (uint32_t j = 0; j < 9; j++)
		{
			std::array<uint32_t, 3> offsets = { Ki, Ki + 2 * Nx - 1, Ki + 2 * (2 * Nx - 1) };
			auto index = j / 3;
			newElements[i][j] = offsets[index] + j % 3;
		}
	}

	UpdateDirichlet(newNodes);

	UpdateNeumann();

	Elements = std::move(newElements);
	Nodes = std::move(newNodes);
}

void GridData::Triangulate()
{
	if (Elements[0].size() == 3) {
		LL = {
			{0, 1, 2},
			{2, 1, 0}
		};
	}
	else if (Elements[0].size() == 6)
	{
		LL = {
			{ 0, 5, 2, 3, 4, 1 },
			{ 1, 4, 3, 2, 5, 0 },
		};
	}

	IndexingFunction = [&](uint32_t i, uint32_t j) { return LL[i & 1][j]; };
}

void GridData::UpdateDirichlet(const std::vector<Point2D>& newNodes)
{
	Point2D bottomLeft = { Nodes[DirichletConditions[0].Node].X, Nodes[DirichletConditions[0].Node].Y };
	auto width = Nodes[DirichletConditions[DirichletConditions.size() - 1].Node].X - bottomLeft.X;
	auto height = Nodes[DirichletConditions[DirichletConditions.size() - 1].Node].Y - bottomLeft.Y;
	DirichletConditions.clear();

	for (uint32_t i = 0; i < newNodes.size(); i++)
	{
		auto y = newNodes[i].Y;
		auto x = newNodes[i].X;
		if (y == bottomLeft.Y)
			DirichletConditions.push_back({ i, BoundaryFunction });
		else if (y == bottomLeft.Y + height)
			DirichletConditions.push_back({ i, BoundaryFunction });
		else if (x == bottomLeft.X)
			DirichletConditions.push_back({ i, BoundaryFunction });
		else if (x == bottomLeft.X + width)
			DirichletConditions.push_back({ i, BoundaryFunction });
	}
}

void GridData::UpdateNeumann()
{
	static std::vector<Element> rectangleNodes = {
		{ 6, 3, 0 },
		{ 0, 1, 2 },
		{ 8, 7, 6 },
		{ 2, 5, 8 }
	};

	static std::vector<Element> triangleNodes = {
		{ 5, 3, 0 },
		{ 0, 1, 2 },
		{ 5, 4, 3 },
		{ 0, 2, 5 }
	};

	std::vector<Element>& nodes = rectangleNodes;
	if (Elements[0].size() == 3)
		nodes = triangleNodes;

	for (auto& nd : NeumannConditions)
	{
		auto lastNode = nd.Nodes[nd.Nodes.size() - 1];
		switch (lastNode)
		{
		case 0: // Vertical - left
			nd.Nodes = nodes[0];
			break;
		case 1: // Horizontal - bottom
			nd.Nodes = nodes[1];
			break;
		case 2: // Horizontal - top
			nd.Nodes = nodes[2];
			break;
		case 3: // Vertical - right
			nd.Nodes = nodes[3];
			break;
		default:
			std::cout << "Incorrect mesh data\n";
			return;
		}
	}
}