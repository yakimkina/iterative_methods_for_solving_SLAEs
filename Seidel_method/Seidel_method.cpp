#include "Seidel_method.h"

vector<VALUE_TYPE>	Seidel_method(vector<vector<VALUE_TYPE>> slae)
{
	int	m = slae.size();
	int	n = slae[0].size();

	vector<vector<VALUE_TYPE>>	matC = create_B(slae, m, n);
	VALUE_TYPE	norm_C = norm_inf(matC, m);
	VALUE_TYPE	norm_C_u = norm_inf_u(matC, m);
	VALUE_TYPE	accuracy = (1 - norm_C) / norm_C_u * EPSILON;

	vector<VALUE_TYPE>	xk = {53, -90, 107, 68};
	vector<VALUE_TYPE>	xk_1 = multiply_with_add_iter(matC, xk, m, n);
	vector<VALUE_TYPE>	delta = subtract_vectors(xk_1, xk, m);

	int k = 0; // количесвто итераций
	while (vector_norm_inf(delta, m) > accuracy)
	{
		xk = xk_1;
		xk_1 = multiply_with_add_iter(matC, xk, m, n);
		delta = subtract_vectors(xk_1, xk, m);

		k += 1;
//		cout << endl;
//		cout << vector_norm_inf(delta, m);
//		cout << endl;
	}
	
	cout << endl << "[МЕТОД ЗЕЙДЕЛЯ] количество итераций = " << k << endl;
	return xk_1;
}
