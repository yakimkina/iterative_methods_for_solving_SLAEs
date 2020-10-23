#include "main.h"

// какое t для МПИ

// как проверить на вырожденность

// какую точность 10^-3 или рассчитать, как в учебнике

// какое х0 приближение использовать?

// добавить в парсер табы

void 	solve_tridiagonal_slae()
{
	int n = 200 + VAR; // номер варианта
	tridiagonal_slae	slae;
	vector<VALUE_TYPE>	a(n);
	vector<VALUE_TYPE>	b(n);
	vector<VALUE_TYPE>	c(n);
	vector<VALUE_TYPE>	d(n);


	fill(a.begin(), a.end(), 1);
	fill(b.begin(), b.end(), 4);
	fill(c.begin(), c.end(), 1);
	d[0] = 6;
	for (int i = 1; i < n - 1; ++i)
		fill(d.begin(), d.end() - 1, 10 - 2 * (i % 2));
	d[n - 1] = 9 - 3 * (n % 2);

	slae.a = a;
	slae.b = b;
	slae.c = c;
	slae.d = d;

	vector<VALUE_TYPE>	Seidel_solution = Seidel_method(slae, n);
	if (!Seidel_solution.empty())
	{
		cout << endl << GREEN << "[МЕТОД ЗЕЙДЕЛЯ] " << MAGENTA << "РЕШЕНИЕ:" << RESET << endl;
		print_vector(Seidel_solution);
	}

/*	vector<VALUE_TYPE>	relaxation_solution = relaxation_method(slae, n);
	if (!relaxation_solution.empty())
	{
		cout << endl << GREEN << "[МЕТОД РЕЛАКСАЦИИ] " << MAGENTA << "РЕШЕНИЕ:" << RESET << endl;
		print_vector(relaxation_solution);
	}*/
}

int main()
{
	if (TRIDIAGONAL_SLAE)
		solve_tridiagonal_slae();
	else
	{
		vector<vector<VALUE_TYPE>>	slae = parsing_file();

		cout << endl << MAGENTA << "СЛАУ:" << RESET << endl;
		print_slae(slae, slae.size(), slae[0].size());

		vector<VALUE_TYPE>	simple_iteration_solution = simple_iteration_method(slae);
		if (!simple_iteration_solution.empty())
		{
			cout << endl << GREEN << "[МЕТОД ПРОСТОЙ ИТЕРАЦИИ] " << MAGENTA << "РЕШЕНИЕ:" << RESET << endl;
			print_vector(simple_iteration_solution);
		}

		vector<VALUE_TYPE>	Jacobi_solution = Jacobi_method(slae);
		if (!Jacobi_solution.empty())
		{
			cout << endl << GREEN << "[МЕТОД ЯКОБИ] " << MAGENTA << "РЕШЕНИЕ:" << RESET << endl;
			print_vector(Jacobi_solution);
		}

		vector<VALUE_TYPE>	Seidel_solution = Seidel_method(slae);
		if (!Seidel_solution.empty())
		{
			cout << endl << GREEN << "[МЕТОД ЗЕЙДЕЛЯ] " << MAGENTA << "РЕШЕНИЕ:" << RESET << endl;
			print_vector(Seidel_solution);
		}

		vector<VALUE_TYPE>	relaxation_solution = relaxation_method(slae);
		if (!relaxation_solution.empty())
		{
			cout << endl << GREEN << "[МЕТОД РЕЛАКСАЦИИ] " << MAGENTA << "РЕШЕНИЕ:" << RESET << endl;
			print_vector(relaxation_solution);
		}
	}
	return 0;
}
