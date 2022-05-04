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
		flux_edge = flux[Elements[elem].Edges[i]];//поток на грани i элемента elem
		if (flux_edge < 0) Sg.push_back( (-1) * Elements[elem].EdgeDirections[i]);
		if (flux_edge > 0) Sg.push_back(Elements[elem].EdgeDirections[i]);
		if (flux_edge==0) Sg.push_back(0);
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
		if (flux[Elements[elem].Edges[out]] < flux[Elements[elem].Edges[number]])
		{
			out = number ;
		}

	return out;
}

vector<double> GridData::vectorD(vector<double> flux, vector<double> betta)
{
	
	uint32_t n = Elements.size();
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //количество граней

	vector<double> d(quantity, 0);
	vector<int32_t> finit_elem(2,0);
	vector<int32_t> edgeSg(4, 0);

	for (uint32_t i = 0; i < quantity; i++) // цикл по всем граням
	{
		finit_elem = GetNumberEdge(i);// два элемента для которых грань смежная
		for (uint32_t j = 0; j < 2; j++)// пройтись по этим двум элементам
		{
			if (finit_elem[j] != -1)
			{
				edgeSg = get_Sg(j, flux);
				d[i] += betta[j] * (edgeSg[0]*flux[Elements[finit_elem[j]].Edges[0]]+edgeSg[1] * flux[Elements[finit_elem[j]].Edges[1]] +edgeSg[2] * flux[Elements[finit_elem[j]].Edges[2]] +edgeSg[3] * flux[Elements[finit_elem[j]].Edges[3]]);
			}
		}

	}

	return d;
}
vector<int32_t> GridData::GetNumberEdge(uint32_t numedge)
{
	uint32_t n = Elements.size();
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //количество граней
	vector<vector<int32_t> > edge(quantity, vector<int32_t>(2, -1));



	bool flag = true;
	uint32_t sum = 0;
	for (uint32_t i = 0; i < quantity; i++) // цикл по всем граням
	{
		for (uint32_t j = 0; j < n && flag; j++) // цикл по всем элементам
		{
			for (uint32_t k = 0; k < 4; k++) // цикл по всем граням элемента
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
int GridData::ig_creation(int elem1, int elem2, int edge, vector<int> &jg)//на вход подали два вектора с номерами граней и саму грань на которой мы тусим
{
	vector<uint32_t> vec;
	int number=0;
	vec=Elements[elem1].Edges;
	if(elem2!=-1) vec.insert(vec.end(), Elements[elem2].Edges.begin(), Elements[elem2].Edges.end());
	sort(vec.begin(), vec.end());
	for (int i = 0; i < vec.size(); i++)
		if (vec[i] < edge)
		{
			number++;
			jg.push_back(vec[i] + 1);
		}
	return number;
}
void GridData::ig_jg_generation(vector<int>&ig, vector<int>& jg)
{	
	int ig_elem;
	int i;
	vector<int32_t> finit_elem;
	finit_elem.resize(2);
	auto nX = Elements[0].Nodes[2];//число узлов на нижней стороне области
	for (i = 0; i < nX; i++)//число 1 в массиве ig совпадает с числом узлов на нижней стороне области
		ig.push_back(1);
	uint32_t n = Elements.size();
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //количество граней
	for (i-=1; i < quantity; i++)//генерация вектора ig и jg
	{
		finit_elem = GetNumberEdge(i);//вернут два элемента для которых грань смежная
		ig_elem = ig[i] + ig_creation(finit_elem[0], finit_elem[1], i, jg);
		ig.push_back(ig_elem);
	}
}
int GridData::find_elem(int edge, int elem1, int elem2)
{
	for (int i = 0; i < 4; i++)
	{
		if (Elements[elem1].Edges[i] == edge) return elem1;
		if (elem2!=-1 && Elements[elem2].Edges[i] == edge) return elem2;
	}
	return -1;
}
void GridData::b_matrix_init(std::vector<double> &gg, std::vector<double> betta, std::vector<double> flux)
{
	uint32_t n = Elements.size();
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //количество граней
	vector<int> elements(2,0);
	double b_elem;
	int this_elem;
	for (int j = 0; j < quantity; j++)
	{
		elements = GetNumberEdge(j);
		vector<uint32_t> vec;
		int number = 0;
		vector<int32_t> Sg(0,4);
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
	uint32_t n = Elements.size();//количество элементов
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //количество граней
	double eps_balance = 1e-1;
	vector<double> d(quantity, 0);
	vector<double> betta(n, 1.0);
	vector<double> q;
	vector<int> ig;
	vector<int> jg;
	vector<double> gg, diag(quantity, 2);
	ig_jg_generation(ig, jg);
	double sum = 0;
	int k = 1;
	vector<int32_t> edgeSg(4, 0);
	
	while (k <100)
	{
	    //k = 0;
		d = vectorD(flux, betta);
		b_matrix_init(gg, betta, flux);
		LOS_ los(ig, jg, gg, diag, d, quantity);
		k++;
		q = los.get_q();
		for (int i = 0; i < n; i++)
		{
			
				edgeSg = get_Sg(i, flux);
				int numberMax = GetMaxFluxElement(i,flux);
				for (int k = 0; k < 4; k++) // 0 2 3 5 -- 0 6 10 0 -- 0 -1 1 0 -- 2 8 12 3
				{
					sum += edgeSg[k] * abs(flux[Elements[i].Edges[k]]) + q[Elements[i].Edges[k]];
				}
				sum /= flux[Elements[i].Edges[numberMax]];
				
				//sum = ((edgeSg[0] * abs(flux[Elements[i].Edges[0]]) + q[Elements[i].Edges[0]]) + (edgeSg[1] * abs(flux[Elements[i].Edges[1]]) + q[Elements[i].Edges[1]]) + (edgeSg[2] * abs(flux[Elements[i].Edges[2]])+q[Elements[i].Edges[2]]) + (edgeSg[3] * abs(flux[Elements[i].Edges[3]]))/(flux[Elements[i].Edges[numberMax]]+q[Elements[i].Edges[3]]));
				
				if (sum > eps_balance)
				{
					betta[i] /= 2;
					//k++;
				}
			
		}
	}


	vector<double> fflux(quantity, 0);
	for (int j = 0; j < n; j++)
	{
		edgeSg = get_Sg(j, flux);
		
		for (int k = 0; k < 4; k++) // 0 2 3 5 -- 0 6 10 0 -- 0 -1 1 0 -- 2 8 12 3
		{
			fflux[Elements[j].Edges[k]] += edgeSg[k] * abs(flux[Elements[j].Edges[k]]) + q[Elements[j].Edges[k]];
		}
	}


}