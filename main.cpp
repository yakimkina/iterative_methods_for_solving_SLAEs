#include "main.h"

// какое t для МПИ

// как проверить на вырожденность

// какую точность 10^-3 или рассчитать, как в учебнике

// какое х0 приближение использовать?

int main()
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

	return 0;
}
