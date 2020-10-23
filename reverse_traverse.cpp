//#include "Gauss_method.h"
#include "main.h"
static	void	subtraction(vector<vector<VALUE_TYPE>> &slae, int k, int size)
{
	for (int i = k - 1; i >= 0; i--)
	{
		VALUE_TYPE a = slae[i][k] / slae[k][k];
		slae[i][k] = 0; // same as: slae[i][k] -= a * slae[k][k];
		slae[i][size] -= a * slae[k][size];
	}
}

bool	all_zero_up(vector<vector<VALUE_TYPE>> &slae, int k)
{
	for (int i = k; i >= 0; i--)
		if (slae[i][k] != 0) return (false);
	return (true);
}

void	reverse_traverse(vector<vector<VALUE_TYPE>> &slae, int size)
{
	for (int i = size - 1; i >= 0; i--)
	{
		if (all_zero_up(slae, i))
			continue;
		
		subtraction(slae, i, size);
	}
}

