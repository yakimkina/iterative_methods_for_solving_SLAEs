#include "Jacobi_method.h"

vector<vector<VALUE_TYPE>>	create_B(vector<vector<VALUE_TYPE>> slae, int m, int n)
{
	for (int i = 0; i < m; i++)
	{
		VALUE_TYPE	a_ii = slae[i][i];

		for (int j = 0; j < n - 1; j++)
			slae[i][j] = (i == j) ? 0 : ((-1) * slae[i][j] / a_ii);

		slae[i][n - 1] /= a_ii;
	}

//	cout << endl << "create B: " << endl;
//	print_slae(slae, m, n);

	return slae;
}

vector<VALUE_TYPE>	Jacobi_method(vector<vector<VALUE_TYPE>> slae)
{
	int	m = slae.size();
	int	n = slae[0].size();


	vector<vector<VALUE_TYPE>>	matB = create_B(slae, m, n);
	VALUE_TYPE	norm_B = norm_inf(matB, m);
	VALUE_TYPE	accuracy = (1 - norm_B) / norm_B * EPSILON;

	vector<VALUE_TYPE>	xk = {53, -90, 107, 68};
	vector<VALUE_TYPE>	xk_1 = multiply_with_add(matB, xk, m, n);
	vector<VALUE_TYPE>	delta = subtract_vectors(xk_1, xk, m);

	int k = 0; // количесвто итераций
	while (vector_norm_inf(delta, m) > accuracy)
	{
		xk = xk_1;
		xk_1 = multiply_with_add(matB, xk, m, n);
		delta = subtract_vectors(xk_1, xk, m);

		k += 1;
//		cout << endl;
//		cout << vector_norm_inf(delta, m);
//		cout << endl;
	}

	cout << endl << "[МЕТОД ЯКОБИ] количество итераций = " << k << endl;
	return xk_1;
}
