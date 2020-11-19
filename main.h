#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#define FILE_NAME	"../TEST/D0.TXT"
//#define FILE_NAME	"../_Lab1,2/SYS2/P_DAT30.TXT"
//#define FILE_NAME	"../_Lab1,2/SYS1/DATA30.TXT"
#define WIDTH	8
#define PREC	4
#define PREC_ZERO	0.000001 // проверка деления на 0 - добавить (проверка на вырожденность?)
#define WITHOUT_IDS	1

/* метод простой итерации */
//#define TAU
#define	TRIDIAGONAL_SLAE	false
#define VAR		19
//#define EPSILON	0.0001
#define EPSILON	0.0000001
#define	TAU		0.02
#define RELAX	0.5

/* colorful output */
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */

/* multiplication */
#define A__ID_COLUMN	0
#define A__B_COLUMN		0
#define B__ID_COLUMN	1
#define B__B_COLUMN		1

using namespace std;

typedef double VALUE_TYPE;

struct	tridiagonal_slae
{
	vector<VALUE_TYPE>	a;
	vector<VALUE_TYPE>	b;
	vector<VALUE_TYPE>	c;
	vector<VALUE_TYPE>	d;
};

struct	tridiagonal_C
{
	vector<VALUE_TYPE>	a;
	vector<VALUE_TYPE>	c;
	vector<VALUE_TYPE>	d;
};

vector<vector<VALUE_TYPE>> 	parsing_file();

vector<VALUE_TYPE>	simple_iteration_method(vector<vector<VALUE_TYPE>> slae);
vector<VALUE_TYPE>	Jacobi_method(vector<vector<VALUE_TYPE>> slae);
vector<VALUE_TYPE>	Seidel_method(vector<vector<VALUE_TYPE>> slae);
vector<VALUE_TYPE>	relaxation_method(vector<vector<VALUE_TYPE>> slae);

vector<VALUE_TYPE>	Seidel_method(tridiagonal_slae slae, int n);


void	print_vector(vector<VALUE_TYPE> x);
void 	print_slae(vector<vector<VALUE_TYPE>> slae, int m, int n); /* print SLAE with b vector */
void	print_slae(vector<vector<VALUE_TYPE>> slae, int size); /* print SLAE without b vector */

VALUE_TYPE	norm_1(vector<vector<VALUE_TYPE>> slae, int size);
VALUE_TYPE	vector_norm_1(vector<VALUE_TYPE> x, int size);
VALUE_TYPE	norm_1_u(vector<vector<VALUE_TYPE>> slae, int size);
VALUE_TYPE	norm_inf(vector<vector<VALUE_TYPE>> slae, int size);
VALUE_TYPE	vector_norm_inf(vector<VALUE_TYPE> x, int size);
VALUE_TYPE	norm_inf_u(vector<vector<VALUE_TYPE>> slae, int size);

vector<VALUE_TYPE>	subtract_vectors(vector<VALUE_TYPE> a, vector<VALUE_TYPE> b, int n);
vector<vector<VALUE_TYPE>>	multiply(vector<vector<VALUE_TYPE>> &a, vector<vector<VALUE_TYPE>> &b);
vector<vector<VALUE_TYPE>>	multiply(vector<vector<VALUE_TYPE>> &a, vector<vector<VALUE_TYPE>> &b, int size, VALUE_TYPE value = 1);
vector<VALUE_TYPE>	multiply_with_add(vector<vector<VALUE_TYPE>> slae, vector<VALUE_TYPE> x, int m, int n);
vector<VALUE_TYPE>	multiply_with_add_iter(vector<vector<VALUE_TYPE>> slae, vector<VALUE_TYPE> x, int m, int n);
vector<VALUE_TYPE>	multiply_with_add_iter_relax(vector<vector<VALUE_TYPE>> slae, vector<VALUE_TYPE> x, VALUE_TYPE w, int m, int n);


vector<vector<VALUE_TYPE>>	create_C(vector<vector<VALUE_TYPE>> slae, int m, int n);
vector<VALUE_TYPE>	get_x0(vector<vector<VALUE_TYPE>> slae, int m, int n);

vector<vector<VALUE_TYPE>>	inverse(vector<vector<VALUE_TYPE>> slae, int size, vector<VALUE_TYPE> (*method)(vector<vector<VALUE_TYPE>>));
vector<VALUE_TYPE>	QR_method(vector<vector<VALUE_TYPE>> slae);
vector<vector<VALUE_TYPE>>	subtract_matrices(vector<vector<VALUE_TYPE>> a, vector<vector<VALUE_TYPE>> b, int m, int n);
vector<vector<VALUE_TYPE>>	transpose(vector<vector<VALUE_TYPE>> &slae);
bool	is_degenerate(vector<vector<VALUE_TYPE>> &slae, int size);
void	reverse_traverse(vector<vector<VALUE_TYPE>> &slae, int size);


vector<vector<VALUE_TYPE>> create_D(vector<vector<VALUE_TYPE>> slae, int size);
vector<vector<VALUE_TYPE>> create_L(vector<vector<VALUE_TYPE>> slae, int size);
vector<vector<VALUE_TYPE>> create_U(vector<vector<VALUE_TYPE>> slae, int size);
void	add_E(vector<vector<VALUE_TYPE>> &slae, int size);

#endif
