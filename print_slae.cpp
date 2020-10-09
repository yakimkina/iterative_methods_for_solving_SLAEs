#include "main.h"

void	print_vector(vector<VALUE_TYPE> x)
{
	cout.precision(PREC);

	for (VALUE_TYPE i : x)
		cout << i << "  ";
	cout << endl;
}

/*void 	print_slae(vector<vector<VALUE_TYPE>> slae)
{
	cout.precision(PREC);

	int	size = slae.size();
	for (int i = 0; i < size - WITHOUT_IDS; i++)
	{
		if (i == size - 1) cout << endl;

		for (int j = 0; j < size; j++)
		{
			if (j == size - 1) cout << "| ";
			cout.width(WIDTH);
			cout << fixed << slae[i][j] << "  ";
		}
		cout << endl;
	}
}*/

void 	print_slae(vector<vector<VALUE_TYPE>> slae, int m, int n)
{
	cout.precision(PREC);

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (j == n - 1) cout << "| ";
			cout.width(WIDTH);
			cout << fixed << slae[i][j] << "  ";
		}
		cout << endl;
	}
}

void	print_slae(vector<vector<VALUE_TYPE>> slae, int size)
{
	cout.precision(PREC);
	
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			cout.width(WIDTH);
			cout << fixed << slae[i][j] << "  ";
		}
		cout << endl;
	}
}

