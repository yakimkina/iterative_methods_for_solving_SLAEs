#include "main.h"

vector<vector<VALUE_TYPE>>	transpose(vector<vector<VALUE_TYPE>> &slae)
{
	vector<vector<VALUE_TYPE>> tr;

	int	size = slae.size();

	for (int j = 0; j < size; j++)
	{
		vector<VALUE_TYPE> tmp;
		for (int i = 0; i < size; i++)
			tmp.push_back(slae[i][j]);

		tr.push_back(tmp);
	}

	return (tr);
}

