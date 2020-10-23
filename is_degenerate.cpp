#include "main.h"

bool	is_degenerate(vector<vector<VALUE_TYPE>> &slae, int size)
{
	for (int i = 0; i < size; i++)
	{
		if (abs(slae[i][i]) < PREC_ZERO)
			return true;
	}

	return false;
}

