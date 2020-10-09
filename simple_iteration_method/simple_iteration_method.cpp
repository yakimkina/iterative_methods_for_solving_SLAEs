#include "simple_iteration_method.h"

vector<vector<VALUE_TYPE>>	create_C(vector<vector<VALUE_TYPE>>	slae, VALUE_TYPE tau, int m, int n)
{
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n - 1; j++)
			slae[i][j] = (-1) * (tau * slae[i][j] - ((i == j) ? 1 : 0));

		slae[i][n - 1] *= tau;
	}

	cout << endl << "create C: " << endl;
	print_slae(slae, m, n);

	return slae;
}

vector<VALUE_TYPE>	simple_iteration_method(vector<vector<VALUE_TYPE>> slae)
{
	int	m = slae.size();
	int	n = slae[0].size();

	VALUE_TYPE	tau = 0.05; //VALUE_TYPE(1) / 38;
//	cout << "tau = " << tau << " " << endl;

	vector<vector<VALUE_TYPE>>	matC = create_C(slae, tau, m, n);
	VALUE_TYPE	norm_C = norm_inf(matC, m);
	VALUE_TYPE	accuracy = (1 - norm_C) / norm_C * EPSILON;

	cout << "norm C = " << norm_C << endl;
	vector<VALUE_TYPE>	xk = {53, -90, 107, 68};
	vector<VALUE_TYPE>	xk_1 = multiply_with_add(matC, xk, m, n);

	cout << "xk1 = ";
	print_vector(xk_1);
	cout << endl;
	vector<VALUE_TYPE>	delta = subtract_vectors(xk_1, xk, m);
	cout.precision(10);
	cout << "acc = " << accuracy << endl;
//	cout << vector_norm_inf(subtract_vectors(xk_1, xk, m), m) << endl;
	cout << vector_norm_inf(delta, m) << endl;
	int i = 100;
	while (vector_norm_inf(delta, m) > accuracy && i--)
	{
		cout.precision(10);
		xk = xk_1;
//		cout << "hear" << endl;
//		print_vector(xk);
		vector<VALUE_TYPE>	xk_1 = multiply_with_add(matC, xk, m, n);
//		cout << endl;
		vector<VALUE_TYPE> delta = subtract_vectors(xk_1, xk, m);
//		print_vector(sub);
		cout << endl;
		cout << vector_norm_inf(delta, m);
		cout << endl;
//		cout << endl << "MULT AFTER " << vector_norm_inf(subtract_vectors(xk_1, xk, m), m) << endl;
	}

//	cout.precision(10);
//	cout << "p = " << p << endl;
	return xk_1;
}

