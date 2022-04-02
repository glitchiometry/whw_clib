#ifndef MATRIX_H
#define MATRIX_H
#include "stdlib.h"
#include "stdio.h"
#include "basics.h"
// points or matrices, collections thereof
typedef struct
{
	int m; int n;
	double **e;
} matrix_real;

typedef struct
{
	int m; int n;
	array_voidstar *data;
} sparse_matrix_real;

// RESUME: include basic sparse matrix operations.

// A variable length collection of matrices of a given dimension
typedef struct
{
	int len;
	int mem;
	int m; int n;
	double ***e;
} matrices_real;

// A variable length collection of matrices of arbitrary dimension
typedef struct
{
	int len;
	int mem;
	matrix_real *e;
} array_matrix_real;

typedef struct
{
	int len;
	int mem;
	matrices_real *e;
} array_matrices_real;

// A variable length collection of variable length collections of matrices
//	of a given dimension
typedef struct
{
	int len;
	int mem;
	array_matrix_real *e;
} aarray_matrix_real;

matrix_real matrix_real_init(int m, int n);
matrices_real matrices_real_init(int mem, int m, int n);
void prep_matrices_real(matrices_real *vrs);
void add_mem_matrices_real(matrices_real *vrs);
void add_mem_matrices_real_until(matrices_real *vrs, int lim);
void add2matrices_real(matrices_real *vrs, double **nmatrix);
void remove_matrices_real(matrices_real *vrs, int m_index);
void free_matrix_real(matrix_real *vr);
void free_matrices_real(matrices_real *vrs);
array_matrix_real array_matrix_real_init(int mem);
void add_mem_array_matrix_real(array_matrix_real *avr);
void add_mem_array_matrix_real_until(array_matrix_real *avr, int lim);
void extend_array_matrix_real(array_matrix_real *avr, int m, int n);
void free_array_matrix_real(array_matrix_real *avr);
matrix_real plus_matrix_real(matrix_real *v1, matrix_real *v2);
matrix_real plus_matrices_real(matrices_real *vs, matrix_real *v2);
matrix_real minus_matrix_real(matrix_real *v1, matrix_real *v2);
matrix_real multiply_matrix_real(matrix_real *v, double s);
double matrix_real_dot(matrix_real *v1, matrix_real *v2);
double matrices_real_dot(matrices_real *vs, matrix_real *v2, int i);
double matrices_reals_dot(matrices_real *v1s, matrices_real *v2s, int i, int j);

array_matrices_real array_matrices_real_init(int mem);
void add_mem_array_matrices_real(array_matrices_real *avr);
void add_mem_array_matrices_real_until(array_matrices_real *avr, int lim);
void extend_array_matrices_real(array_matrices_real *avr, int mem, int m, int n);
void free_array_matrices_real(array_matrices_real *avr);

// matrices and arrays thereof of integers 
typedef struct
{
	int m; int n;
	int **e;
} matrix_int;

// A variable length collection of matrices of a given dimension
typedef struct
{
	int len;
	int mem;
	int m; int n;
	int ***e;
} matrices_int;

// A variable length collection of matrices of arbitrary dimension
typedef struct
{
	int len;
	int mem;
	matrix_int *e;
} array_matrix_int;

typedef struct
{
	int len;
	int mem;
	matrices_int *e;
} array_matrices_int;


matrix_int matrix_int_init(int m, int n);
matrices_int matrices_int_init(int mem, int m, int n);
void prep_matrices_int(matrices_int *vrs);
void add_mem_matrices_int(matrices_int *vrs);
void add_mem_matrices_int_until(matrices_int *vrs, int lim);
void add2matrices_int(matrices_int *vrs, int **pt);
void remove_matrices_int(matrices_int *vrs, int vertex);
array_matrix_int array_matrix_int_init(int mem);
void free_matrix_int(matrix_int *vr);
void free_matrices_int(matrices_int *vrs);
void add_mem_array_matrix_int(array_matrix_int *avr);
void add_mem_array_matrix_int_until(array_matrix_int *avr, int lim);
void extend_array_matrix_int(array_matrix_int *avr, int m, int n);
matrix_int plus_matrix_int(matrix_int *v1, matrix_int *v2);
matrix_int plus_matrices_int(matrices_int *vs, matrix_int *v2);
matrix_int minus_matrix_int(matrix_int *v1, matrix_int *v2);
matrix_int multiply_matrix_int(matrix_int *v, int s);
int matrix_int_dot(matrix_int *v1, matrix_int *v2);
int matrices_int_dot(matrices_int *vs, matrix_int *v2, int i);
int matrices_ints_dot(matrices_int *v1s, matrices_int *v2s, int i, int j);

array_matrices_int array_matrices_int_init(int mem);
void add_mem_array_matrices_int(array_matrices_int *avr);
void add_mem_array_matrices_int_until(array_matrices_int *avr, int lim);
void extend_array_matrices_int(array_matrices_int *avr, int mem, int m, int n);
void free_array_matrices_int(array_matrices_int *avr);

#endif
