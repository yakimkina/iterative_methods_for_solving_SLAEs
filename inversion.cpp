#include "main.h"

void	insert_x(vector<vector<VALUE_TYPE>> &slae, int size, vector<VALUE_TYPE> x, int pos)
{
	for (int i = 0; i < size; i++)
		slae[i][pos] = x[i];
}

void	replace_b(vector<vector<VALUE_TYPE>> &slae, int size, int pos)
{
	for (int i = 0; i < size; i++)
		slae[i][size] = (i == pos) ? 1 : 0;
}

vector<vector<VALUE_TYPE>>	inverse(vector<vector<VALUE_TYPE>> slae, int size, vector<VALUE_TYPE> (*method)(vector<vector<VALUE_TYPE>>))
{
	vector<vector<VALUE_TYPE>>	inverse_slae(size, vector<VALUE_TYPE>(size));

	/* Ax_i = e_i,
	 * where x_i - column_i A^(-1) */

	for (int i = 0; i < size; i++)
	{
		/* b replace on e_i*/
		replace_b(slae, size, i);

		vector<VALUE_TYPE>	x_i = method(slae);

		/* make inverse slae: insert into column_i  solution x_i */
		insert_x(inverse_slae, size, x_i, i);
	}

	return inverse_slae;
}
