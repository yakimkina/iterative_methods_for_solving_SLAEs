#include "main.h"

VALUE_TYPE	sum(vector<VALUE_TYPE> a_i, vector<vector<VALUE_TYPE>> &b, int j, int size)
{
	VALUE_TYPE s = 0;

	for (int i = 0; i < size; i++)
		s += a_i[i] * b[i][j];

	return s;
}

VALUE_TYPE	sum(vector<VALUE_TYPE> a_i, vector<VALUE_TYPE> &b, int size)
{
	VALUE_TYPE s = 0;

	for (int i = 0; i < size; i++)
		s += a_i[i] * b[i];

	return s;
}

vector<vector<VALUE_TYPE>>	multiply(vector<vector<VALUE_TYPE>> &a, vector<vector<VALUE_TYPE>> &b)
{
	vector<vector<VALUE_TYPE>> slae;

	/* a is [m x n] size, b is [n x k] size */
	int		m = a.size() - A__ID_COLUMN;
	int		n = a[0].size() - A__B_COLUMN;
	int		k = b[0].size() - B__B_COLUMN;

//	cout << m << "|" << n << "|" << k << "|" << b.size() - B__ID_COLUMN;
	if (n != b.size() - B__ID_COLUMN)
	{
		cerr << "[ERROR] Умножение матриц прервано: неправильный размер матриц." << endl;
		exit(0);
	}

	for (int i = 0; i < m; i++)
	{
		vector<VALUE_TYPE> tmp;
		for (int j = 0; j < k; j++)
			tmp.push_back(sum(a[i], b, j, n));

		slae.push_back(tmp);
	}

	return slae;
}

vector<VALUE_TYPE>	multiply_with_add(vector<vector<VALUE_TYPE>> slae, vector<VALUE_TYPE> x, int m, int n)
{
	vector<VALUE_TYPE>	mult;

	for (int i = 0; i < m; i++)
		mult.push_back(sum(slae[i], x, n) + slae[i][n - 1]);

	return mult;
}