#include "QR_method.h"

static	vector<VALUE_TYPE>	get_solution(vector<vector<VALUE_TYPE>> &slae, int size)
{
	vector<VALUE_TYPE> x;

	for (int i = 0; i < size; i++)
		x.push_back(slae[i][size] / slae[i][i]);

	return x;
}

vector<VALUE_TYPE>	QR_method(vector<vector<VALUE_TYPE>> slae)
{
	int size = slae.size() - 1;
	
	if (PRINT_STEPS)
		cout << endl << GREEN << "[МЕТОД QR]" << RESET << endl;

	vector<vector<VALUE_TYPE>> matR = slae;
	vector<vector<VALUE_TYPE>> matT = create_T(size);
	build_R(matR, matT, size);

	/* check if matrix is degenerate */
	if (is_degenerate(matR, size))
	{
		cout << RED << "[ERROR] Матрица вырождена!" << RESET << endl;
		return vector<VALUE_TYPE>(0);
	}

	if (PRINT_STEPS)
	{
		cout << endl << MAGENTA << "[МЕТОД QR] Матрица R:" << RESET << endl;
		print_slae(matR, size, size + 1);
	}

	vector<vector<VALUE_TYPE>>	matQ = transpose(matT);

	if (PRINT_STEPS)
	{
		cout << endl << MAGENTA << "[МЕТОД QR] Матрица Q:" << RESET << endl;
		print_slae(matQ, matQ.size());
	}

	if (PRINT_CHECK)
	{
		/* check if QR == A */

//		make in class: if (A == multiplication(matQ, matR))

		cout << endl << YELLOW << "[проверка] A - QR:" << endl << "A:" << RESET << endl;
		vector<vector<VALUE_TYPE>>	matA = multiply(matQ, matR);
		vector<vector<VALUE_TYPE>>	delta = subtract_matrices(slae, matA, size, size);
		print_slae(delta, delta.size());
	}

	reverse_traverse(matR, size);

	return get_solution(matR, size);
}

