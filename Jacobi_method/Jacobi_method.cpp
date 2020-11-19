#include "Jacobi_method.h"

vector<VALUE_TYPE>	get_x0(vector<vector<VALUE_TYPE>> slae, int m, int n)
{
	vector<VALUE_TYPE>	x0;

	for (int i = 0; i < m; i++)
		x0.push_back(slae[i][n - 1]);

	return x0;
}

vector<vector<VALUE_TYPE>>	create_C(vector<vector<VALUE_TYPE>> slae, int m, int n)
{
	for (int i = 0; i < m; i++)
	{
		VALUE_TYPE	a_ii = slae[i][i];

		for (int j = 0; j < n - 1; j++)
			slae[i][j] = (i == j) ? 0 : ((-1) * slae[i][j] / a_ii);

		slae[i][n - 1] /= a_ii;
	}

//	cout << endl << "create C: " << endl;
//	print_slae(slae, m, n);

	return slae;
}

vector<VALUE_TYPE>	Jacobi_method(vector<vector<VALUE_TYPE>> slae)
{
	cout << RED << "##########################" << endl << "[МЕТОД ЯКОБИ]" << RESET << endl;
	int	m = slae.size();
	int	n = slae[0].size();
	
	vector<vector<VALUE_TYPE>>	matC = create_C(slae, m, n);
	VALUE_TYPE	norm_C = norm_inf(matC, m);
	VALUE_TYPE	accuracy = (1 - norm_C) / norm_C * EPSILON;

	vector<VALUE_TYPE>	xk = get_x0(slae, m, n);
	vector<VALUE_TYPE>	x0 = xk;
//	cout << "x0 = ";
//	print_vector(xk);
//	vector<VALUE_TYPE>	xk = {3.1, -1.1};
	vector<VALUE_TYPE>	xk_1 = multiply_with_add(matC, xk, m, n);
	vector<VALUE_TYPE>	delta = subtract_vectors(xk_1, xk, m);

	int k = 0; // количесвто итераций
	int ap = 0; // апостериорная оценка
	while (vector_norm_inf(delta, m) > accuracy)
	{
//		cout << "xk = " << endl;
//		print_vector(xk);
		xk = xk_1;
		xk_1 = multiply_with_add(matC, xk, m, n);
		delta = subtract_vectors(xk_1, xk, m);

		if (vector_norm_inf(delta, m) > EPSILON) ap += 1;


		k += 1;
//		cout << endl;
//		cout << vector_norm_inf(delta, m);
//		cout << endl;
	}

	cout << "количество итераций = " << k << endl;
	cout.precision(6);
	cout<< BLUE << "априорная = " << log(EPSILON / vector_norm_inf(subtract_vectors(x0, xk_1, m), m)) / log(norm_C) << endl;
	cout.precision(10);
	cout << "апостериорная = " << ap  << RESET<< endl;
	cout << "точность = " << accuracy << endl;
	cout.precision(4);

	cout << "норма C = " << norm_C << endl;

	return xk_1;
}
