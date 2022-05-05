#include "FluxBalancer.h"
#include "FEMSolver.h"
#include "FluxCalculator.h"
#include "GridData.h"
#include "Solver.cpp"
#include <algorithm>
#include<vector>
#include "Math.h"
#include "SparseMatrix.h"
#include "functions.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <math.h>
using namespace std;

vector<int32_t> GridData::get_Sg(int elem, vector<double> flux)
{
	vector<int32_t> Sg;
	double flux_edge;
	for (int i = 0; i < 4; i++)
	{
		flux_edge = flux[Elements[elem].Edges[i]];//ïîòîê íà ãðàíè i ýëåìåíòà elem
		if (flux_edge < 0) Sg.push_back((-1) * Elements[elem].EdgeDirections[i]);
		if (flux_edge > 0) Sg.push_back(Elements[elem].EdgeDirections[i]);
		if (flux_edge == 0) Sg.push_back(0);
	}
	return Sg;
}
int GridData::get_number(int elem, int edge)
{
	int number;
	for (number = 0; number < 4; number++)
		if (Elements[elem].Edges[number] == edge) return number;

	return -1;
}
int GridData::GetMaxFluxElement(int elem, vector<double> flux)
{
	int number;
	int out = 0;
	for (number = 0; number < 4; number++)
		if (abs(flux[Elements[elem].Edges[out]]) < abs(flux[Elements[elem].Edges[number]]))
		{
			out = number;
		}

	return out;
}

/*double GridData::AveragedFluxElement(int elem)
{

}*/

vector<double> GridData::vectorD(vector<double> flux, vector<double> betta)
{

	uint32_t n = Elements.size();
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //edge

	vector<double> d(quantity, 0);
	vector<int32_t> finit_elem(2, 0);
	vector<int32_t> edgeSg(4, 0);

	double averagedflux = 0;

	for (uint32_t i = 0; i < quantity; i++)  // edge
	{
		finit_elem = GetNumberEdge(i);
		for (uint32_t j = 0; j < 2; j++) // elements
		{
			if (finit_elem[j] != -1)
			{
				edgeSg = get_Sg(j, flux);

				//d[i] -= betta[j] * (edgeSg[0] * flux[Elements[finit_elem[j]].Edges[0]] + edgeSg[1] * flux[Elements[finit_elem[j]].Edges[1]] + edgeSg[2] * flux[Elements[finit_elem[j]].Edges[2]] + edgeSg[3] * flux[Elements[finit_elem[j]].Edges[3]]) / 4;
			}
		}

	}

	return d;
}
vector<int32_t> GridData::GetNumberEdge(uint32_t numedge)
{
	uint32_t n = Elements.size();
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //êîëè÷åñòâî ãðàíåé
	vector<vector<int32_t> > edge(quantity, vector<int32_t>(2, -1));



	bool flag = true;
	uint32_t sum = 0;
	for (uint32_t i = 0; i < quantity; i++) // öèêë ïî âñåì ãðàíÿì
	{
		for (uint32_t j = 0; j < n && flag; j++) // öèêë ïî âñåì ýëåìåíòàì
		{
			for (uint32_t k = 0; k < 4; k++) // öèêë ïî âñåì ãðàíÿì ýëåìåíòà
			{
				if (i == Elements[j].Edges[k]) {

					edge[i][sum] = j;

					sum++;
					if (sum == 2)
					{
						sum = 0;
						flag = false;
					}
				}
			}
		}
		flag = true; sum = 0;
	}


	return { edge[numedge][0],edge[numedge][1] };
}
int GridData::ig_creation(int elem1, int elem2, int edge, vector<int>& jg)//íà âõîä ïîäàëè äâà âåêòîðà ñ íîìåðàìè ãðàíåé è ñàìó ãðàíü íà êîòîðîé ìû òóñèì
{
	vector<uint32_t> vec;
	int number = 0;
	vec = Elements[elem1].Edges;
	if (elem2 != -1) vec.insert(vec.end(), Elements[elem2].Edges.begin(), Elements[elem2].Edges.end());
	sort(vec.begin(), vec.end());
	for (int i = 0; i < vec.size(); i++)
		if (vec[i] < edge)
		{
			number++;
			jg.push_back(vec[i] + 1);
		}
	return number;
}
void GridData::ig_jg_generation(vector<int>& ig, vector<int>& jg)
{
	int ig_elem;
	int i;
	vector<int32_t> finit_elem;
	finit_elem.resize(2);
	auto nX = Elements[0].Nodes[2];//÷èñëî óçëîâ íà íèæíåé ñòîðîíå îáëàñòè
	for (i = 0; i < nX; i++)//÷èñëî 1 â ìàññèâå ig ñîâïàäàåò ñ ÷èñëîì óçëîâ íà íèæíåé ñòîðîíå îáëàñòè
		ig.push_back(1);
	uint32_t n = Elements.size();
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //êîëè÷åñòâî ãðàíåé
	for (i -= 1; i < quantity; i++)//ãåíåðàöèÿ âåêòîðà ig è jg
	{
		finit_elem = GetNumberEdge(i);//âåðíóò äâà ýëåìåíòà äëÿ êîòîðûõ ãðàíü ñìåæíàÿ
		ig_elem = ig[i] + ig_creation(finit_elem[0], finit_elem[1], i, jg);
		ig.push_back(ig_elem);
	}
}
int GridData::find_elem(int edge, int elem1, int elem2)
{
	for (int i = 0; i < 4; i++)
	{
		if (Elements[elem1].Edges[i] == edge) return elem1;
		if (elem2 != -1 && Elements[elem2].Edges[i] == edge) return elem2;
	}
	return -1;
}
void GridData::b_matrix_init(std::vector<double>& gg, std::vector<double> betta, std::vector<double> flux)
{
	gg.clear();
	uint32_t n = Elements.size();
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //êîëè÷åñòâî ãðàíåé
	vector<int> elements(2, 0);
	double b_elem;
	int this_elem;
	for (int j = 0; j < quantity; j++)
	{
		elements = GetNumberEdge(j);
		vector<uint32_t> vec;
		int number = 0;
		vector<int32_t> Sg(0, 4);
		vec = Elements[elements[0]].Edges;
		if (elements[1] != -1) vec.insert(vec.end(), Elements[elements[1]].Edges.begin(), Elements[elements[1]].Edges.end());
		sort(vec.begin(), vec.end());
		for (int i = 0; i < vec.size(); i++)
			if (vec[i] < j)
			{
				this_elem = find_elem(vec[i], elements[0], elements[1]);
				Sg = get_Sg(this_elem, flux);
				b_elem = betta[this_elem] * Sg[get_number(this_elem, j)] * Sg[get_number(this_elem, vec[i])];
				gg.push_back(b_elem);
			}
	}
}
void GridData::flux_balancer(vector<double> flux)
{
	uint32_t n = Elements.size();//êîëè÷åñòâî ýëåìåíòîâ
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //êîëè÷åñòâî ãðàíåé
	double eps_balance = 1e-3;
	vector<double> d(quantity, 0);
	vector<double> betta(n, 1.0);
	vector<double> q(quantity, 0);
	vector<int> ig;
	vector<int> jg;
	vector<double> gg, diag(quantity, 2);
	bool Flag=false;
	ig_jg_generation(ig, jg);
	double sum = 0, global_sum=5;
	vector<int32_t> edgeSg(4, 0);
	int	maxiter = 50;
	LOS_ los;
	int i = 0;
	while (global_sum > eps_balance && Flag==false && i < maxiter)
	{
		i++;
		Flag = true;
		global_sum = 0;
		d = vectorD(flux, betta);
		if (los.get_DotProduct(d, quantity)< 0.0001)
		{
			break;
		}
		b_matrix_init(gg, betta, flux);
		los.solve(ig, jg, gg, diag, d, quantity);
		q = los.get_q();
		for (int i = 0; i < n; i++)
		{
			sum = 0;
			edgeSg = get_Sg(i, flux);
			int numberMax = GetMaxFluxElement(i, flux);
			for (int k = 0; k < 4; k++) 
				sum += edgeSg[k] * (abs(flux[Elements[i].Edges[k]]) + q[Elements[i].Edges[k]]);//небаланс
			global_sum += betta[i]*abs(sum);
			sum /=  abs(flux[Elements[i].Edges[numberMax]]);
			if (sum > eps_balance)
			{
				betta[i] *= 2;
				Flag = false;
			}

		}
	}


	vector<double> fflux(n, 0);
	vector<double> edgefflux(quantity, 0);
	for (int j = 0; j < n; j++)
	{
		edgeSg = get_Sg(j, flux);
		for (int k = 0; k < 4; k++) 
			fflux[j] += edgeSg[k] * (abs(flux[Elements[j].Edges[k]]) + q[Elements[j].Edges[k]]);
		
	}
	// пересчет граней 
	for (int j = 0; j < n; j++)
	{
		for (int k = 0; k < 4; k++)
			edgefflux[Elements[j].Edges[k]] =  (abs(flux[Elements[j].Edges[k]]) + q[Elements[j].Edges[k]]);
	}

	

}