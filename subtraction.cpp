#include "main.h"

vector<VALUE_TYPE>	subtract_vectors(vector<VALUE_TYPE> a, vector<VALUE_TYPE> b, int n)
{
	vector<VALUE_TYPE> tmp;

	for (int i = 0; i < n; i++)
		tmp.push_back(a[i] - b[i]);

	return tmp;
}

vector<vector<VALUE_TYPE>>	subtract_matrices(vector<vector<VALUE_TYPE>> a, vector<vector<VALUE_TYPE>> b, int m, int n)
{
	vector<vector<VALUE_TYPE>>	slae;

	for (int i = 0; i < m; i++)
	{
		vector<VALUE_TYPE>	tmp;

		for (int j = 0; j < n; j++)
			tmp.push_back(a[i][j] - b[i][j]);

		slae.push_back(tmp);
	}

	return (slae);
}


