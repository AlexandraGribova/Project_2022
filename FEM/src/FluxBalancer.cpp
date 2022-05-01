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


vector<int32_t> GridData::GetNumberEdge(uint32_t numedge)
{
	uint32_t n = Elements.size();
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //���������� ������
	vector<vector<int32_t> > edge(quantity, vector<int32_t>(2, -1));

	bool flag = true;
	uint32_t sum = 0;
	for (uint32_t i = 0; i < quantity; i++) // ���� �� ���� ������
	{
		for (uint32_t j = 0; j < n && flag; j++) // ���� �� ���� ���������
		{
			for (uint32_t k = 0; k < 4; k++) // ���� �� ���� ������ ��������
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

int GridData::ig_creation(int elem1, int elem2, int edge, vector<int> &jg)//�� ���� ������ ��� ������� � �������� ������ � ���� ����� �� ������� �� �����
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
	int finit_elem1, finit_elem2;
	auto nX = Elements[0].Nodes[2];//����� ����� �� ������ ������� �������
	for (i = 0; i < nX; i++)//����� 1 � ������� ig ��������� � ������ ����� �� ������ ������� �������
		ig.push_back(1);
	uint32_t n = Elements.size();
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //���������� ������
	for (i-=1; i < quantity; i++)//��������� ������� ig � jg
	{
		finit_elem1= GetNumberEdge(i)[0];//������ ��� �������� ��� ������� ���� �������
		finit_elem2 = GetNumberEdge(i)[1];
		ig_elem = ig[i] + ig_creation(finit_elem1, finit_elem2, i, jg);
		ig.push_back(ig_elem);
	}
}

void GridData::flux_balancer(vector<double> flux)
{
	vector<int> ig;
	vector<int> jg;
	vector<double> betta;
	uint32_t n = Elements.size();
	betta.resize(n);
	ig_jg_generation(ig, jg);
}