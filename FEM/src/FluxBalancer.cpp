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
		if (abs(flux_edge) < 1e-6) Sg.push_back(0);
		else
		{
			if (flux_edge < 0) Sg.push_back((-1) * Elements[elem].EdgeDirections[i]);
			if (flux_edge > 0) Sg.push_back(Elements[elem].EdgeDirections[i]);
		}
		
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
		if (abs(flux[Elements[elem].Edges[out]] ) < abs(flux[Elements[elem].Edges[number]]))
		{
			out = number;
		}

	return out;
}

vector<double> GridData::vectorD(vector<double> flux, vector<double> betta)
{

	uint32_t n = Elements.size();
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //êîëè÷åñòâî ãðàíåé
	int Sg_num;
	vector<double> d(quantity, 0);
	vector<int32_t> finit_elem(2, 0);
	vector<int32_t> edgeSg(4, 0);

	for (uint32_t i = 0; i < quantity; i++) // öèêë ïî âñåì ãðàíÿì
	{
		finit_elem = GetNumberEdge(i);// äâà ýëåìåíòà äëÿ êîòîðûõ ãðàíü ñìåæíàÿ
		for (uint32_t j = 0; j < 2; j++)// ïðîéòèñü ïî ýòèì äâóì ýëåìåíòàì
		{
			if (finit_elem[j] != -1)
			{
				Sg_num = get_number(finit_elem[j], i);
				edgeSg = get_Sg(j, flux);
				d[i]-= betta[j] * edgeSg[Sg_num] * (edgeSg[0] * flux[Elements[finit_elem[j]].Edges[0]] + edgeSg[1] * flux[Elements[finit_elem[j]].Edges[1]] + edgeSg[2] * flux[Elements[finit_elem[j]].Edges[2]] + edgeSg[3] * flux[Elements[finit_elem[j]].Edges[3]]);
			}
		}

	}

	return d;
}
void GridData::Change_q(vector<double>& q)
{
	vector<int> border;
	int elem = 0;
	auto nX = Elements[0].Nodes[2];
	uint32_t n = Elements.size();
	uint32_t quantity = Elements[n - 1].Edges[3] + 1;
	for (elem; elem < nX; elem++)
		border.push_back(elem);
	for (elem--; elem < quantity - nX + 1;)
	{
		
		if (elem + 2 * nX - 1 >= quantity) break;
		elem += 2 * nX - 1;
		border.push_back(elem);
		/*elem += nX - 1;
		border.push_back(elem);
		if (elem+nX == quantity) break; //Это для всех фиксированных грвниц
		elem += nX;
		border.push_back(elem);*/
	}
	for (elem+=nX; elem < quantity; elem++)
		border.push_back(elem);
	/*for (elem++; elem < quantity;elem++)
		border.push_back(elem);*/

	for (int i = 0; i < quantity; i++)
		for(int j=0; j<border.size(); j++)
			if (i == border[j])
			{
				q[i] = 0;
				break;
			}
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
			jg.push_back(vec[i]);
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
		ig.push_back(0);
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

void  GridData::get_diag(vector<double>& diag, vector<double> betta)
{
	uint32_t n = Elements.size();//êîëè÷åñòâî ýëåìåíòîâ
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //êîëè÷åñòâî ãðàíåé
	vector<int32_t> finit_elem;
	finit_elem.resize(2);
	for (int i = 0; i < quantity; i++)
	{
		finit_elem = GetNumberEdge(i);
		if (finit_elem[1] == -1) diag[i] = 2 * betta[finit_elem[0]];
		else diag[i] = betta[finit_elem[0]]+ betta[finit_elem[1]];
	}
}


void GridData::flux_balancer(vector<double> flux)
{
	uint32_t n = Elements.size();//êîëè÷åñòâî ýëåìåíòîâ
	uint32_t quantity = Elements[n - 1].Edges[3] + 1; //êîëè÷åñòâî ãðàíåé
	double eps_balance = 1e-6;
	vector<double> d(quantity, 0);
	vector<double> betta(n, 0.5);
	vector<double> q;
	vector<int> border;
	vector<int> ig;
	vector<int> jg;
	vector<double> gg, diag(quantity, 2);
	bool Flag=false;
	ig_jg_generation(ig, jg);
	double sum = 0, global_sum=1;
	vector<int32_t> edgeSg(4, 0);
	LOS_ los;



	vector<double> neb(n, 0);
	flux[15] += 50;
	/*flux[35] += 0.1;
	flux[48] += 0.1;
	flux[63] += 0.1;
	flux[104] += 0.1;*/
	for (int j = 0; j < n; j++)
	{
		edgeSg = get_Sg(j, flux);
		for (int k = 0; k < 4; k++)
			neb[j] += edgeSg[k] * (abs(flux[Elements[j].Edges[k]]));
		//cout << std::uppercase << std::scientific << neb[j] << endl;
	}
	int iter = 0;
	while (global_sum> eps_balance && Flag==false)
	{
		iter++;
		Flag = true;
		global_sum = 0;
		d = vectorD(flux, betta);
		if (los.get_DotProduct(d, quantity) < 0.0001)
		{
			q.clear();//надо ли?
			for (int h = 0; h < quantity; h++) q.push_back(0);
			break;
		}
		b_matrix_init(gg, betta, flux);
		
		for (int i = 0; i < NeumannConditions.size(); i+=4)
		{
			int global_edge;
			for (int k = 0; k < 4; k++)
			{
				global_edge=Elements[NeumannConditions[i].Element].Edges[k];
				diag[global_edge] = 1;
				d[global_edge] = 0;
			}
		
		}
		get_diag(diag, betta);
		los.solve(ig, jg, gg, diag, d, quantity);
		q = los.get_q();
		Change_q(q);//приравниваем надбавку к 0 на границе области
		for (int i = 0; i < NeumannConditions.size(); i+=4)
		{
			int global_edge;
			for (int k = 0; k < 4; k++)
			{
				global_edge = Elements[NeumannConditions[i].Element].Edges[k];
				q[global_edge] = 0;
			}

		}
		
		for (int i = 0; i < n; i++)
		{
			sum = 0;
			edgeSg = get_Sg(i, flux);
			int numberMax = GetMaxFluxElement(i, flux);
			for (int k = 0; k < 4; k++) 
				sum += edgeSg[k] * (abs(flux[Elements[i].Edges[k]]) + q[Elements[i].Edges[k]]);//небаланс
			global_sum += betta[i]*abs(sum);
			sum /= abs(flux[Elements[i].Edges[numberMax]]);

			if (abs(sum) > eps_balance)
			{
				betta[i] = sqrt(betta[i]);
				Flag = false;
			}

		}
		if (iter == 100) break;
	}

	vector<double> fflux(n, 0);
	for (int j = 0; j < n; j++)
	{
		edgeSg = get_Sg(j, flux);
		for (int k = 0; k < 4; k++)
			fflux[j] += edgeSg[k] * (abs(flux[Elements[j].Edges[k]]) + q[Elements[j].Edges[k]]);
		
		//cout << std::uppercase << std::scientific << fflux[j] << endl;
	}

	vector<double> ffluxedge(quantity, 0);
	for (int j = 0; j < quantity; j++)
	{
		
			ffluxedge[j] = (abs(flux[j]) + q[j]);
			cout << std::uppercase << std::scientific << q[j] << endl;
	}

}