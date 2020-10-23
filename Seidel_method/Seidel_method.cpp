#include "Seidel_method.h"

tridiagonal_C	create_C(tridiagonal_slae slae, int n)
{
	tridiagonal_C	tmp;

	vector<VALUE_TYPE>	a;
	vector<VALUE_TYPE>	c;
	vector<VALUE_TYPE>	d;

/*	for (int i = 0; i < n; i++)
	{
		VALUE_TYPE	a_ii = slae[i][i];

		for (int j = 0; j < n - 1; j++)
			slae[i][j] = (i == j) ? 0 : ((-1) * slae[i][j] / a_ii);

		slae[i][n - 1] /= a_ii;
	}*/

//	cout << endl << "create C: " << endl;
//	print_slae(slae, m, n);

	return tmp;
}

vector<VALUE_TYPE>	Seidel_method(tridiagonal_slae slae, int n)
{
	
	vector<VALUE_TYPE> xk_1;
//	vector<vector<VALUE_TYPE>>	matC = create_B(slae, m, n);
//	VALUE_TYPE	norm_C = norm_inf(matC, m);
//	VALUE_TYPE	norm_C_u = norm_inf_u(matC, m);
//	VALUE_TYPE	accuracy = (1 - norm_C) / norm_C_u * EPSILON;
//
//	vector<VALUE_TYPE>	xk = get_x0(slae, m, n);
//	cout << "x0 = ";
//	print_vector(xk);
//	vector<VALUE_TYPE>	xk_1 = multiply_with_add_iter(matC, xk, m, n);
//	vector<VALUE_TYPE>	delta = subtract_vectors(xk_1, xk, m);
//
//	int k = 0; // количесвто итераций
//	while (vector_norm_inf(delta, m) > accuracy)
//	{
//		xk = xk_1;
//		xk_1 = multiply_with_add_iter(matC, xk, m, n);
//		delta = subtract_vectors(xk_1, xk, m);
//
//		k += 1;
//		cout << endl;
//		cout << vector_norm_inf(delta, m);
//		cout << endl;
//	}
//
//	cout << endl << "[МЕТОД ЗЕЙДЕЛЯ] количество итераций = " << k << endl;
	return xk_1;
}

static	void	add_x_id(vector<vector<VALUE_TYPE>> &slae)
{
	vector<VALUE_TYPE> x;

	for (int i = 1; i <= slae.size(); i++)
		x.push_back(i);

	// id = 0 relates to b column
	x.push_back(0);

	slae.push_back(x);
}

vector<vector<VALUE_TYPE>> create_D(vector<vector<VALUE_TYPE>> slae, int size)
{
	vector<vector<VALUE_TYPE>>	D;

	for (int i = 0; i < size; i++)
	{
		vector<VALUE_TYPE>	tmp;
		for (int j = 0; j < size + 1; j++)
			tmp.push_back((j == i) ? slae[i][j] : 0);

		D.push_back(tmp);
	}

	add_x_id(D);

//	cout << "Matrix D:" << endl;
//	print_slae(D, size, size + 1);

	return D;
}


vector<vector<VALUE_TYPE>> create_L(vector<vector<VALUE_TYPE>> slae, int size)
{
	vector<vector<VALUE_TYPE>>	L;

	for (int i = 0; i < size; i++)
	{
		vector<VALUE_TYPE>	tmp;
		for (int j = 0; j < size; j++)
			tmp.push_back((j < i) ? slae[i][j] : 0);

		L.push_back(tmp);
	}

	return L;
}

vector<vector<VALUE_TYPE>> create_U(vector<vector<VALUE_TYPE>> slae, int size)
{
	vector<vector<VALUE_TYPE>>	U;

	for (int i = 0; i < size; i++)
	{
		vector<VALUE_TYPE>	tmp;
		for (int j = 0; j < size; j++)
			tmp.push_back((j > i) ? slae[i][j] : 0);

		U.push_back(tmp);
	}

	return U;
}

void	add_E(vector<vector<VALUE_TYPE>> &slae, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			slae[i][j] += (i == j) ? 1 : 0;
		
		slae[i].push_back(0); // столбец b
	}

	add_x_id(slae);
}

// объединить с методом релаксации где w = 1
vector<vector<VALUE_TYPE>>	create_C_theoretic(vector<vector<VALUE_TYPE>> slae, int m, int n)
{
	vector<vector<VALUE_TYPE>>	matD = create_D(slae, m);
	vector<vector<VALUE_TYPE>>	matL = create_L(slae, m);
	vector<vector<VALUE_TYPE>>	matU = create_U(slae, m);

	vector<vector<VALUE_TYPE>>	invD = inverse(matD, m, QR_method); // D^-1
	vector<vector<VALUE_TYPE>>	DL = multiply(invD, matL, m); // D^-1 * L
	add_E(DL, m); // E + D^-1 * L
	vector<vector<VALUE_TYPE>>	inv = inverse(DL, m, QR_method); // (E + D^-1 * L)^-1
	vector<vector<VALUE_TYPE>>	DU = multiply(invD, matU, m, -1); // (-1) * D^-1 * U

	return multiply(inv, DU, m);
}

vector<VALUE_TYPE>	Seidel_method(vector<vector<VALUE_TYPE>> slae)
{
	cout << RED << "##########################" << endl << "[МЕТОД ЗЕЙДЕЛЯ]" << RESET << endl;
	int	m = slae.size();
	int	n = slae[0].size();

	vector<vector<VALUE_TYPE>>	matC_theoretic = create_C_theoretic(slae, m, n);
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
