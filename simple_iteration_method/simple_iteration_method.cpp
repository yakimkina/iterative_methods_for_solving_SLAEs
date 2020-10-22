#include "simple_iteration_method.h"

vector<vector<VALUE_TYPE>>	create_C(vector<vector<VALUE_TYPE>> slae, VALUE_TYPE tau, int m, int n)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n - 1; j++)
			slae[i][j] = (-1) * (tau * slae[i][j] - ((i == j) ? 1 : 0));

		slae[i][n - 1] *= tau;
	}

//	cout << endl << "create C: " << endl;
//	print_slae(slae, m, n);

	return slae;
}

vector<VALUE_TYPE>	simple_iteration_method(vector<vector<VALUE_TYPE>> slae)
{
	int	m = slae.size();
	int	n = slae[0].size();

	VALUE_TYPE	tau = 0.05;

	vector<vector<VALUE_TYPE>>	matC = create_C(slae, tau, m, n);
	VALUE_TYPE	norm_C = norm_inf(matC, m);
	VALUE_TYPE	accuracy = (1 - norm_C) / norm_C * EPSILON;

	vector<VALUE_TYPE>	xk = {53, -90, 107, 68};
	vector<VALUE_TYPE>	xk_1 = multiply_with_add(matC, xk, m, n);
	vector<VALUE_TYPE>	delta = subtract_vectors(xk_1, xk, m);

	int k = 0; // количесвто итераций
	while (vector_norm_inf(delta, m) > accuracy)
	{
		xk = xk_1;
		xk_1 = multiply_with_add(matC, xk, m, n);
		delta = subtract_vectors(xk_1, xk, m);

		k += 1;
//		cout << endl;
//		cout << vector_norm_inf(delta, m);
//		cout << endl;

	}

	cout << endl << "[МПИ] количество итераций = " << k << endl;
	return xk_1;
}

