#include "functions.h"
#include "SparseMatrix.h"

std::pair<std::vector<uint32_t>, std::vector<uint32_t>> generatePortrait(uint32_t N, const std::vector<Element>& elements)
{
	uint32_t bigNumberWhichIDontKnowButMoreIThinkAboutItMoreItSeemsToBeLessThanNTimes8 = N * 8;
	std::vector<std::vector<uint32_t>> list(2, std::vector<uint32_t>(bigNumberWhichIDontKnowButMoreIThinkAboutItMoreItSeemsToBeLessThanNTimes8, 0));
	std::vector<uint32_t> ig, jg;
	ig.reserve(N);
	jg.reserve(bigNumberWhichIDontKnowButMoreIThinkAboutItMoreItSeemsToBeLessThanNTimes8);
	std::vector<uint32_t> listbeg(bigNumberWhichIDontKnowButMoreIThinkAboutItMoreItSeemsToBeLessThanNTimes8 * 2, 0);
	uint32_t listSize = 0;

	for (const auto elem : elements)
	{
		auto& elemNodes = elem.Nodes;
		for (uint32_t i = 0; i < elemNodes.size(); i++)
		{
			auto k = elemNodes[i];
			for (uint32_t j = i + 1; j < elemNodes.size(); j++)
			{
				auto ind1 = k;
				auto ind2 = elemNodes[j];
				if (ind2 < ind1)
				{
					ind1 = ind2;
					ind2 = k;
				}
				auto iaddr = listbeg[ind2];
				if (!iaddr)
				{
					listSize++;
					listbeg[ind2] = listSize;
					list[0][listSize] = ind1;
					list[1][listSize] = 0;
				}
				else
				{
					while (list[0][iaddr] < ind1 && list[1][iaddr] > 0)
					{
						iaddr = list[1][iaddr];
					}
					if (list[0][iaddr] > ind1)
					{
						listSize++;
						list[0][listSize] = list[0][iaddr];
						list[1][listSize] = list[1][iaddr];
						list[0][iaddr] = ind1;
						list[1][iaddr] = listSize;
					}
					else
					{
						if (list[0][iaddr] < ind1)
						{
							listSize++;
							list[1][iaddr] = listSize;
							list[0][listSize] = ind1;
							list[1][listSize] = 0;
						}
					}
				}
			}
		}
	}
	ig.push_back(0);
	for (uint32_t i = 0; i < N; i++)
	{
		ig.push_back(ig[i]);
		auto iaddr = listbeg[i];
		while (iaddr != 0)
		{
			jg.push_back(list[0][iaddr]);
			ig[i + 1]++;
			iaddr = list[1][iaddr];
		}
	}
	return { std::move(ig), std::move(jg) };
}
