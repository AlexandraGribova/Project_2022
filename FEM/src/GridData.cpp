#include "GridData.h"

#include <fstream>
#include <array>
#include <iostream>
#include <format>

GridData::GridData()
{
	IndexingFunction = [](uint32_t i, uint32_t j) { return j; };
}

void GridData::LoadElements(const std::string& path)
{
	std::ifstream in(path);
	if (!in)
	{
		std::cout << std::format("Was unable to open {}\n", path);
		return;
	}
	uint32_t numberOfEls;
	uint32_t numberOfNodes;
	in >> numberOfEls;
	in >> numberOfNodes;
	Elements = std::vector<Element>(numberOfEls);

	for (uint32_t i = 0; i < numberOfEls; i++)
	{
		Element el;
		el.Nodes.reserve(numberOfNodes);
		el.Edges.reserve(numberOfNodes);
		for (uint32_t j = 0; j < numberOfNodes; j++)
		{
			uint32_t node;
			in >> node;
			el.Nodes.push_back(node);
		}
		Elements[i] = el;
	}
}

void GridData::LoadNodes(const std::string& path)
{
	std::ifstream in(path);
	if (!in)
	{
		std::cout << std::format("Was unable to open {}\n", path);
		return;
	}
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
	if (!in)
	{
		std::cout << std::format("Was unable to open {}\n", path);
		return;
	}
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
	if (!in)
	{
		std::cout << std::format("Was unable to open {}\n", path);
		return;
	}
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

	auto nX = Elements[0].Nodes[2];
	ModifiedNx = 2 * (Elements[0].Nodes[0] / nX) * (2 * nX - 1) + 2 * (Elements[0].Nodes[0] % nX) + 2 * nX - 1;
	auto numberOfNewNodes = 2 * (Elements[Elements.size() - 2].Nodes[0] / nX) * (2 * nX - 1) + 2 * (Elements[Elements.size() - 2].Nodes[0] % nX) + 2 * (2 * nX - 1) + 2 + 1;
	std::vector<Point2D> newNodes(numberOfNewNodes, { 0.0, 0.0 });
	std::vector<Element> newElements(Elements.size()); // TODO opt

	for (uint32_t i = 0; i < Elements.size(); i++)
	{
		newElements[i].Nodes.resize(6);
		newElements[i].Edges = Elements[i].Edges;
		newElements[i].EdgeDirections = Elements[i].EdgeDirections;
		auto elementNodes = Elements[i].Nodes;
		if (!(i & 1))
		{
			auto Ki = 2 * (elementNodes[0] / nX) * (2 * nX - 1) + 2 * (elementNodes[0] % nX);
			auto midBottom = Point2D{ (Nodes[elementNodes[1]].X + Nodes[elementNodes[0]].X) / 2.0,  Nodes[elementNodes[0]].Y };
			auto midLeft = Point2D{ Nodes[elementNodes[0]].X, (Nodes[elementNodes[2]].Y + Nodes[elementNodes[0]].Y) / 2.0 };
			auto center = Point2D{ (Nodes[elementNodes[1]].X + Nodes[elementNodes[0]].X) / 2.0, (Nodes[elementNodes[2]].Y + Nodes[elementNodes[0]].Y) / 2.0 };
			newNodes[Ki] = Nodes[elementNodes[0]];
			newNodes[Ki + 1] = midBottom;
			newNodes[Ki + 2] = Nodes[elementNodes[1]];
			newNodes[Ki + 2 * nX - 1] = midLeft;
			newNodes[Ki + 2 * nX] = center;
			newNodes[Ki + 2 * (2 * nX - 1)] = Nodes[elementNodes[2]];
			newElements[i].Nodes[0] = Ki;
			newElements[i].Nodes[1] = Ki + 1;
			newElements[i].Nodes[2] = Ki + 2;
			newElements[i].Nodes[3] = Ki + 2 * nX - 1;
			newElements[i].Nodes[4] = Ki + 2 * nX;
			newElements[i].Nodes[5] = Ki + 2 * (2 * nX - 1);

		}
		else
		{
			auto Ki = 2 * ((elementNodes[0] - 1) / nX) * (2 * nX - 1) + 2 * ((elementNodes[0] - 1) % nX);
			auto midRight = Point2D{ Nodes[elementNodes[0]].X, (Nodes[elementNodes[2]].Y + Nodes[elementNodes[0]].Y) / 2.0 };
			auto midTop = Point2D{ (Nodes[elementNodes[1]].X + Nodes[elementNodes[2]].X) / 2.0, Nodes[elementNodes[2]].Y };
			newNodes[Ki + 2 * nX + 1] = midRight;
			newNodes[Ki + 2 * (2 * nX - 1) + 1] = midTop;
			newNodes[Ki + 2 * (2 * nX - 1) + 2] = Nodes[elementNodes[2]];
			newElements[i].Nodes[0] = Ki + 2;
			newElements[i].Nodes[1] = Ki + 2 * nX;
			newElements[i].Nodes[2] = Ki + 2 * nX + 1;
			newElements[i].Nodes[3] = Ki + 2 * (2 * nX - 1);
			newElements[i].Nodes[4] = Ki + 2 * (2 * nX - 1) + 1;
			newElements[i].Nodes[5] = Ki + 2 * (2 * nX - 1) + 2;
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

	auto nX = Elements[0].Nodes[2];
	ModifiedNx = 2 * (Elements[0].Nodes[0] / nX) * (2 * nX - 1) + 2 * (Elements[0].Nodes[0] % nX) + 2 * nX - 1;
	auto numberOfNewNodes = 2 * (Elements[Elements.size() - 1].Nodes[0] / nX) * (2 * nX - 1) + 2 * (Elements[Elements.size() - 1].Nodes[0] % nX) + 2 * (2 * nX - 1) + 2 + 1;
	std::vector<Point2D> newNodes(numberOfNewNodes, { 0.0, 0.0 });
	std::vector<Element> newElements(Elements.size());

	for (uint32_t i = 0; i < Elements.size(); i++)
	{
		newElements[i].Nodes.resize(9);
		newElements[i].Edges = Elements[i].Edges;
		newElements[i].EdgeDirections = Elements[i].EdgeDirections;
		auto elementNodes = Elements[i].Nodes;
		auto midBottom = Point2D{ (Nodes[elementNodes[1]].X + Nodes[elementNodes[0]].X) / 2.0, Nodes[elementNodes[0]].Y };
		auto midLeft = Point2D{ Nodes[elementNodes[0]].X, (Nodes[elementNodes[2]].Y + Nodes[elementNodes[0]].Y) / 2.0 };
		auto center = Point2D{ (Nodes[elementNodes[1]].X + Nodes[elementNodes[0]].X) / 2.0, (Nodes[elementNodes[2]].Y + Nodes[elementNodes[0]].Y) / 2.0 };
		auto midRight = Point2D{ Nodes[elementNodes[1]].X, (Nodes[elementNodes[2]].Y + Nodes[elementNodes[0]].Y) / 2.0 };
		auto midTop = Point2D{ (Nodes[elementNodes[1]].X + Nodes[elementNodes[0]].X) / 2.0, Nodes[elementNodes[2]].Y };
		auto Ki = 2 * (elementNodes[0] / nX) * (2 * nX - 1) + 2 * (elementNodes[0] % nX);
		newNodes[Ki] = Nodes[elementNodes[0]];
		newNodes[Ki + 1] = midBottom;
		newNodes[Ki + 2] = Nodes[elementNodes[1]];
		newNodes[Ki + 2 * nX - 1] = midLeft;
		newNodes[Ki + 2 * nX] = center;
		newNodes[Ki + 2 * nX + 1] = midRight;
		newNodes[Ki + 2 * (2 * nX - 1)] = Nodes[elementNodes[2]];
		newNodes[Ki + 2 * (2 * nX - 1) + 1] = midTop;
		newNodes[Ki + 2 * (2 * nX - 1) + 2] = Nodes[elementNodes[3]];
		for (uint32_t j = 0; j < 9; j++)
		{
			std::array<uint32_t, 3> offsets = { Ki, Ki + 2 * nX - 1, Ki + 2 * (2 * nX - 1) };
			auto index = j / 3;
			newElements[i].Nodes[j] = offsets[index] + j % 3;
		}
	}

	UpdateDirichlet(newNodes);

	UpdateNeumann();

	Elements = std::move(newElements);
	Nodes = std::move(newNodes);
}

void GridData::Triangulate()
{
	if (Elements[0].Nodes.size() == 3) {
		LL = {
			{0, 1, 2},
			{2, 1, 0}
		};
	}
	else if (Elements[0].Nodes.size() == 6)
	{
		LL = {
			{ 0, 5, 2, 3, 4, 1 },
			{ 1, 4, 3, 2, 5, 0 },
		};
	}

	IndexingFunction = [&](uint32_t i, uint32_t j) { return LL[i & 1][j]; };
}

void GridData::SetupEdges()
{
	auto nX = Elements[0].Nodes[2];
	if (Elements[0].Nodes.size() == 4)
	{
		for (auto elemNum = 0u; elemNum < Elements.size(); elemNum++)
		{
			auto& elem = Elements[elemNum];
			auto index = elemNum / (nX - 1) * (2 * nX - 1) + elemNum % (nX - 1);
			elem.Edges.push_back(index);
			elem.Edges.push_back(index + nX - 1);
			elem.Edges.push_back(index + nX);
			elem.Edges.push_back(index + 2 * nX - 1);
			elem.EdgeDirections = { -1, -1, 1, 1 };
		}
	}
	else if (Elements[0].Nodes.size() == 3)
	{
		for (auto elemNum = 0u; elemNum < Elements.size(); elemNum++)
		{
			auto& elem = Elements[elemNum];
			if (!(elemNum & 1)) // |\ 
			{
				auto row = (elemNum / 2) / (nX - 1);
				auto col = (elemNum / 2) % (nX - 1);
				auto index = row * (3 * nX - 2) + col;

				elem.Edges.push_back(index);
				elem.Edges.push_back(index + nX - 1 + col);
				elem.Edges.push_back(index + nX + col);
				elem.EdgeDirections = { -1, -1, 1 };
			}
			else				// \|
			{
				auto row = ((elemNum - 1) / 2) / (nX - 1);
				auto col = ((elemNum - 1) / 2) % (nX - 1);
				auto index = row * (3 * nX - 2) + col;
				elem.Edges.push_back(index + nX + col);
				elem.Edges.push_back(index + nX + col + 1);
				elem.Edges.push_back((row + 1) * (3 * nX - 2) + col); // Все это лучше проверить
				elem.EdgeDirections = { -1, 1, 1 };
			}
		}
	}
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
	static std::vector<NodeIDs> rectangleNodes = {
		{ 6, 3, 0 },
		{ 0, 1, 2 },
		{ 8, 7, 6 },
		{ 2, 5, 8 }
	};

	static std::vector<NodeIDs> triangleNodes = {
		{ 5, 3, 0 },
		{ 0, 1, 2 },
		{ 5, 4, 3 },
		{ 0, 2, 5 }
	};

	std::vector<NodeIDs>& nodes = rectangleNodes;
	if (Elements[0].Nodes.size() == 3)
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