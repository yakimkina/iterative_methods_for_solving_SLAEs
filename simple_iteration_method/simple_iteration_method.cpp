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
	cout << RED << "##########################" << endl << "[МЕТОД ПРОСТОЙ ИТЕРАЦИИ]" << RESET << endl;

	int	m = slae.size();
	int	n = slae[0].size();

	VALUE_TYPE	tau = TAU;

	vector<vector<VALUE_TYPE>>	matC = create_C(slae, tau, m, n);
	VALUE_TYPE	norm_C = norm_inf(matC, m);
	VALUE_TYPE	accuracy = (1 - norm_C) / norm_C * EPSILON;

	vector<VALUE_TYPE>	xk = get_x0(slae, m, n);
	vector<VALUE_TYPE>	x0 = xk;
//	cout << "x0 = ";
//	print_vector(xk);
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
	cout << "количество итераций = " << k << endl;
	cout.precision(6);
	cout << "априорная = " << pow(norm_C, m) * vector_norm_inf(subtract_vectors(x0, xk_1, m), m) << endl;
	cout.precision(10);
	cout << "апостериорная = " << norm_C / (1 - norm_C) * vector_norm_inf(delta, m) << endl;
	cout << "точность = " << accuracy << endl;
	cout.precision(4);
	cout << "норма C = " << norm_C << endl;



	return xk_1;
}

