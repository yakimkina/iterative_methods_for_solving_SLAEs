#include "QR_method.h"

vector<vector<VALUE_TYPE>>	create_T(int size)
{
	vector<vector<VALUE_TYPE>> T;

	for (int i = 0; i < size; i++)
	{
		vector<VALUE_TYPE> tmp;
		for (int j = 0; j < size; j++)
			(i == j) ? tmp.push_back(1) : tmp.push_back(0);

		T.push_back(tmp);
	}

	return (T);
}

void	build_R(vector<vector<VALUE_TYPE>> &matR, vector<vector<VALUE_TYPE>> &matT, int size)
{
	for (int i = 0; i < size; i++)
	{
//		if (matR[i][i] == 0) continue;  // ???????

		for (int j = i + 1; j < size; j++)
		{
			if (matR[j][i] == 0) continue; // ???
			
			VALUE_TYPE x_main = matR[i][i];
			VALUE_TYPE x_j = matR[j][i];
			VALUE_TYPE sq = sqrt(x_main * x_main + x_j * x_j);
			VALUE_TYPE c = x_main / sq;
			VALUE_TYPE s = x_j / sq;

			for (int k = 0; k <= size; k++)
			{
				VALUE_TYPE x1 = matR[i][k];
				VALUE_TYPE x2 = matR[j][k];
				matR[i][k] = c * x1 + s * x2;
				matR[j][k] = (-1) * s * x1 + c * x2;

				if (k == size) break;
				
				VALUE_TYPE t1 = matT[i][k];
				VALUE_TYPE t2 = matT[j][k];
				matT[i][k] = c * t1 + s * t2;
				matT[j][k] = (-1) * s * t1 + c * t2;
			}
		}
	}
}

