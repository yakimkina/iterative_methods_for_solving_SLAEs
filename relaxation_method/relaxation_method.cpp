#include "relaxation_method.h"

// проверить правильность подсчета матрицы С - теоретическая


vector<vector<VALUE_TYPE>>	create_E_with_w(VALUE_TYPE w, int size)
{
	vector<vector<VALUE_TYPE>>	E;

	for (int i = 0; i < size; i++)
	{
		vector<VALUE_TYPE>	tmp;
		for (int j = 0; j < size; j++)
			tmp.push_back((i == j) ? (1 - w) : 0);

		E.push_back(tmp);
	}

	return E;
}

vector<vector<VALUE_TYPE>>	create_C_theoretic(vector<vector<VALUE_TYPE>> slae, VALUE_TYPE w, int m, int n)
{
	vector<vector<VALUE_TYPE>>	matD = create_D(slae, m);
	vector<vector<VALUE_TYPE>>	matL = create_L(slae, m);
	vector<vector<VALUE_TYPE>>	matU = create_U(slae, m);

	vector<vector<VALUE_TYPE>>	invD = inverse(matD, m, QR_method); // D^-1
	vector<vector<VALUE_TYPE>>	DL = multiply(invD, matL, m, w); // w * D^-1 * L
	add_E(DL, m); // E + D^-1 * L
	vector<vector<VALUE_TYPE>>	inv = inverse(DL, m, QR_method); // (E + w * D^-1 * L)^-1
	vector<vector<VALUE_TYPE>>	DU = multiply(invD, matU, m, w); // w * D^-1 * U
	vector<vector<VALUE_TYPE>>	matE_w = create_E_with_w(w, m); // (1 - w) * E
	vector<vector<VALUE_TYPE>>	tmp = subtract_matrices(matE_w, DU, m, m); // (1 - w) * E - w * D^-1 * U

	return multiply(inv, tmp, m);
}


vector<VALUE_TYPE>	relaxation_method(vector<vector<VALUE_TYPE>> slae)
{
	cout << RED << "##########################" << endl << "[МЕТОД РЕЛАКСАЦИИ]" << RESET << endl;
	
	int	m = slae.size();
	int	n = slae[0].size();

	VALUE_TYPE	w = RELAX;

	vector<vector<VALUE_TYPE>>	matC_theoretic = create_C_theoretic(slae, w, m, n);
	cout << "Matrix C:" << endl;
	print_slae(matC_theoretic, m);

	vector<vector<VALUE_TYPE>>	matC = create_C(slae, m, n);
	VALUE_TYPE	norm_C = norm_inf(matC_theoretic, m);
//	VALUE_TYPE	norm_C_u = norm_inf_u(matC, m);
	VALUE_TYPE	accuracy = (1 - norm_C) / norm_C * EPSILON;

	vector<VALUE_TYPE>	xk = get_x0(slae, m, n);
	vector<VALUE_TYPE>	x0 = xk;
//	cout << "x0 = ";
//	print_vector(xk);
	vector<VALUE_TYPE>	xk_1 = multiply_with_add_iter_relax(matC, xk, w, m, n);
	vector<VALUE_TYPE>	delta = subtract_vectors(xk_1, xk, m);

	int k = 0; // количесвто итераций
	int ap = 0; // апостериорная оценка
	while (vector_norm_inf(delta, m) > accuracy)
	{
		xk = xk_1;
		xk_1 = multiply_with_add_iter_relax(matC, xk, w, m, n);
		delta = subtract_vectors(xk_1, xk, m);

		if (vector_norm_inf(delta, m) > EPSILON) ap += 1;

		k += 1;
//		cout << endl;
//		cout << vector_norm_inf(delta, m);
//		cout << endl;
	}

	cout << "количество итераций = " << k << endl;
	cout.precision(6);
	cout << BLUE<< "априорная = " << log(EPSILON / vector_norm_inf(subtract_vectors(x0, xk_1, m), m)) / log(norm_C) << endl;
	cout.precision(10);
	cout << "апостериорная = " << ap  << RESET<< endl;
	cout << "точность = " << accuracy << endl;
	cout.precision(4);

	cout << "норма C = " << norm_C << endl;

	return xk_1;
}


