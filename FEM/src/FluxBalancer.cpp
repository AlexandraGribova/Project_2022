#include "FluxBalancer.h"
#include "FEMSolver.h"
#include "FluxCalculator.h"
#include "GridData.h"
#include <algorithm>
#include<vector>
#include "Math.h"
#include "SparseMatrix.h"
#include "functions.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
using namespace std;

vector<int32_t> GridData::get_Sg(int elem, vector<double> flux)
{
	vector<int32_t> Sg;
	double flux_edge;
	for (int i = 0; i < 4; i++)
	{		
		flux_edge = flux[Elements[elem].Edges[i]];//????? ?? ????? i ???????? elem
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

vector<int32_t> GridData::vectorD(vector<double> flux, vector<double> betta)
{
	
	uint32_t n = Elements.size();
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //?????????? ??????

	vector<int32_t> d(quantity, 0);
	vector<int32_t> finit_elem(2,0);
	vector<int32_t> edgeSg(4, 0);

	for (uint32_t i = 0; i < quantity; i++) // ???? ?? ???? ??????
	{
		finit_elem = GetNumberEdge(i);// ??? ???????? ??? ??????? ????? ???????
		for (uint32_t j = 0; j < 2; j++)// ???????? ?? ???? ???? ?????????
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
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //?????????? ??????
	vector<vector<int32_t> > edge(quantity, vector<int32_t>(2, -1));



	bool flag = true;
	uint32_t sum = 0;
	for (uint32_t i = 0; i < quantity; i++) // ???? ?? ???? ??????
	{
		for (uint32_t j = 0; j < n && flag; j++) // ???? ?? ???? ?????????
		{
			for (uint32_t k = 0; k < 4; k++) // ???? ?? ???? ?????? ????????
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

int GridData::ig_creation(int elem1, int elem2, int edge, vector<int> &jg)//?? ???? ?????? ??? ??????? ? ???????? ?????? ? ???? ????? ?? ??????? ?? ?????
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
	auto nX = Elements[0].Nodes[2];//????? ????? ?? ?????? ??????? ???????
	for (i = 0; i < nX; i++)//????? 1 ? ??????? ig ????????? ? ?????? ????? ?? ?????? ??????? ???????
		ig.push_back(1);
	uint32_t n = Elements.size();
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //?????????? ??????
	for (i-=1; i < quantity; i++)//????????? ??????? ig ? jg
	{
		finit_elem = GetNumberEdge(i);//?????? ??? ???????? ??? ??????? ????? ???????
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
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //?????????? ??????
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
	uint32_t n = Elements.size();//?????????? ?????????
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //?????????? ??????

	vector<int32_t> d(quantity, 0);
	vector<double> betta(n, 1.0);
	vector<int> ig;
	vector<int> jg;
	vector<double> gg;
	ig_jg_generation(ig, jg);
	d = vectorD(flux, betta);
	b_matrix_init(gg, betta, flux);
}