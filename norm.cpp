#include "main.h"

VALUE_TYPE	vector_norm(vector<VALUE_TYPE> x, int size)
{
	VALUE_TYPE	sum = 0;
	for (int i = 0; i < size; i++)
		sum += x[i] * x[i];

	return sqrt(sum);
}

VALUE_TYPE	norm_1(vector<vector<VALUE_TYPE>> slae, int size)
{
	VALUE_TYPE	max = 0;

	for (int j = 0; j < size; j++)
	{
		VALUE_TYPE	sum = 0;
		for (int i = 0; i < size; i++)
			sum += abs(slae[i][j]);

		if (sum > max) max = sum;
	}

	return max;
}

VALUE_TYPE	vector_norm_1(vector<VALUE_TYPE> x, int size)
{
	VALUE_TYPE	sum = abs(x[0]);

	for (int i = 1; i < size; i++)
		sum += abs(x[i]);

	cout << RESET;
	return sum;
}

VALUE_TYPE	norm_inf(vector<vector<VALUE_TYPE>> slae, int size)
{
	VALUE_TYPE	max = 0;

	for (int i = 0; i < size; i++)
	{
		VALUE_TYPE	sum = 0;
		for (int j = 0; j < size; j++)
			sum += abs(slae[i][j]);

		if (sum > max) max = sum;
	}

	return max;
}

VALUE_TYPE	vector_norm_inf(vector<VALUE_TYPE> x, int size)
{
	VALUE_TYPE	max = abs(x[0]);

	for (int i = 1; i < size; i++)
		if (abs(x[i]) > max) max = abs(x[i]);

	return max;
}

