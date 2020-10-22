cmake_minimum_required(VERSION 3.17)
project(iterative_methods_for_solving_SLAEs)

set(CMAKE_CXX_STANDARD 14)

add_executable(iterative_methods_for_solving_SLAEs
		main.cpp print_slae.cpp norm.cpp subtraction.cpp multiplication.cpp main.h
		parser/parser.cpp parser/parser.h
		simple_iteration_method/simple_iteration_method.cpp simple_iteration_method/simple_iteration_method.h Jacobi_method/Jacobi_method.cpp Jacobi_method/Jacobi_method.h Seidel_method/Seidel_method.cpp Seidel_method/Seidel_method.h relaxation_method/relaxation_method.cpp relaxation_method/relaxation_method.h)