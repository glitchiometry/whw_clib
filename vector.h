#ifndef VECTOR_H
#define VECTOR_H
#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#define VSRGET vectors_real_get
#define VSIGET vectors_int_get
#define VSCGET vectors_char_get
#define VSRSET vectors_real_set
#define VSISET vectors_int_set
#define VSCSET vectors_char_set

// points or vectors, collections thereof
typedef struct
{
	int dim;
	double *e;
} vector_real;

// A variable length collection of vectors of a given dimension
typedef struct
{
	int len;
	int mem;
	int dim;
	double ***e;
} vectors_real;

// A variable length collection of vectors of arbitrary dimension
typedef struct
{
	int len;
	int mem;
	vector_real *e;
} array_vector_real;

// A variable length array of collections of vectors
typedef struct
{
	int len;
	int mem;
	vectors_real *e;
} array_vectors_real;

// A variable length collection of variable length collections of vectors
//	of a given dimension
typedef struct
{
	int len;
	int mem;
	array_vector_real *e;
} aarray_vector_real;

void aarray_vector_real_init(aarray_vector_real *aavr, int mem); 
void free_aarray_vector_real(aarray_vector_real *aavr);
// RESUME: define these
void add2aarray_vector_real(aarray_vector_real *aavr, array_vector_real elem);
void reset_mem_aarray_vector_real(aarray_vector_real *aavr);
void add_mem_aarray_vector_real(aarray_vector_real *aavr);
void add_mem_aarray_vector_real_until(aarray_vector_real *aavr, int len);
void prep_aarray_vector_real(aarray_vector_real *aavr);
void extend_aarray_vector_real(aarray_vector_real *aavr);

vector_real vector_real_init(int dim);
void set_equal_double(double *v1, double *v2, int len);

void vectors_real_init(vectors_real *vrs, int mem, int dim);
double vectors_real_get(vectors_real *vrs, int elem, int coord);
void vectors_real_set(vectors_real *vrs, int elem, int coord, double val);
void prep_vectors_real(vectors_real *vrs);
void add_mem_vectors_real(vectors_real *vrs);
void add_mem_vectors_real_until(vectors_real *vrs, int lim);
void add2vectors_real(vectors_real *vrs, double *pt);
void add2vectors_real_zeros(vectors_real *vrs);
void remove_vectors_real(vectors_real *vrs, int vertex);
void free_vector_real(vector_real *vr);
void free_vectors_real(vectors_real *vrs);
void array_vector_real_init(array_vector_real *avr, int mem);
void add_mem_array_vector_real(array_vector_real *avr);
void add_mem_array_vector_real_until(array_vector_real *avr, int lim);
void extend_array_vector_real(array_vector_real *avr, int dim);
void free_array_vector_real(array_vector_real *avr);
void add2array_vector_real(array_vector_real *avr, vector_real elem);


void xform_p_vector_real(vector_real *v1, vector_real *v3);
vector_real plus_vector_real(vector_real *v1, vector_real *v2);
vector_real plus_vectors_real(vectors_real *vs, vector_real *v2);
vector_real minus_vector_real(vector_real *v1, vector_real *v2);
vector_real multiply_vector_real(vector_real *v, double s);

double vector_real_dot(vector_real *v1, vector_real *v2);
double vectors_real_dot(vectors_real *vs, vector_real *v2, int i);
double vectors_reals_dot(vectors_real *v1s, vectors_real *v2s, int i, int j);

array_vectors_real array_vectors_real_init(int mem);
void add_mem_array_vectors_real(array_vectors_real *avr);
void add_mem_array_vectors_real_until(array_vectors_real *avr, int lim);
void extend_array_vectors_real(array_vectors_real *avr, int mem, int dim);
void free_array_vectors_real(array_vectors_real *avr);

// vectors and arrays thereof of integers 
typedef struct
{
	int dim;
	int *e;
} vector_int;

// A variable length collection of vectors of a given dimension
typedef struct
{
	int len;
	int mem;
	int dim;
	int ***e;
} vectors_int;

// A variable length collection of vectors of arbitrary dimension
typedef struct
{
	int len;
	int mem;
	vector_int *e;
} array_vector_int;

typedef struct
{
	int len;
	int mem;
	vectors_int *e;
} array_vectors_int;

typedef struct
{
	int mem;
	int len;
	array_vector_int *e;
} aarray_vector_int;

void aarray_vector_int_init(aarray_vector_int *aavr, int mem); 
void free_aarray_vector_int(aarray_vector_int *aavr);

// RESUME: define these
void add2aarray_vector_int(aarray_vector_int *aavr, array_vector_int elem);
void reset_mem_aarray_vector_int(aarray_vector_int *aavr);
void add_mem_aarray_vector_int(aarray_vector_int *aavr);
void add_mem_aarray_vector_int_until(aarray_vector_int *aavr, int len);
void prep_aarray_vector_int(aarray_vector_int *aavr);
void extend_aarray_vector_int(aarray_vector_int *aavr);


void set_equal_int(int *a1, int *a2, int len);
vector_int vector_int_init(int dim);
void vectors_int_init(vectors_int *vrs, int mem, int dim);
int vectors_int_get(vectors_int *vrs, int elem, int coord);
void vectors_int_set(vectors_int *vrs, int elem, int coord, int val);

void prep_vectors_int(vectors_int *vrs);
void add_mem_vectors_int(vectors_int *vrs);
void add_mem_vectors_int_until(vectors_int *vrs, int lim);
void add2vectors_int(vectors_int *vrs, int *pt);
void add2vectors_int_zeros(vectors_int *vrs);
void remove_vectors_int(vectors_int *vrs, int vertex);
void array_vector_int_init(array_vector_int *avr, int mem);
void free_array_vector_int(array_vector_int *avr);  
void free_vector_int(vector_int *vr);
void free_vectors_int(vectors_int *vrs);
void add_mem_array_vector_int(array_vector_int *avr);
void add_mem_array_vector_int_until(array_vector_int *avr, int lim);
void extend_array_vector_int(array_vector_int *avr, int dim);
void add2array_vector_int(array_vector_int *avr, vector_int elem); // RESUME: define this

void xform_p_vector_int(vector_int *v1, vector_int *v3);

vector_int plus_vector_int(vector_int *v1, vector_int *v2);
vector_int plus_vectors_int(vectors_int *vs, vector_int *v2);
vector_int minus_vector_int(vector_int *v1, vector_int *v2);
vector_int multiply_vector_int(vector_int *v, int s);
int vector_int_dot(vector_int *v1, vector_int *v2);
int vectors_int_dot(vectors_int *vs, vector_int *v2, int i);
int vectors_ints_dot(vectors_int *v1s, vectors_int *v2s, int i, int j);

array_vectors_int array_vectors_int_init(int mem);
void add_mem_array_vectors_int(array_vectors_int *avr);
void add_mem_array_vectors_int_until(array_vectors_int *avr, int lim);
void extend_array_vectors_int(array_vectors_int *avr, int mem, int dim);
void free_array_vectors_int(array_vectors_int *avr);


// points or vectors, collections thereof
typedef struct
{
	int dim;
	char *e;
} vector_char;

// A variable length collection of vectors of a given dimension
typedef struct
{
	int len;
	int mem;
	int dim;
	char ***e;
} vectors_char;

// A variable length collection of vectors of arbitrary dimension
typedef struct
{
	int len;
	int mem;
	vector_char *e;
} array_vector_char;

typedef struct
{
	int len;
	int mem;
	vectors_char *e;
} array_vectors_char;

// A variable length collection of variable length collections of vectors
//	of a given dimension
typedef struct
{
	int len;
	int mem;
	array_vector_char *e;
} aarray_vector_char;

void aarray_vector_char_init(aarray_vector_char *aavr, int mem);
void free_aarray_vector_char(aarray_vector_char *aavr);
void add2aarray_vector_char(aarray_vector_char *aavr, array_vector_char elem);
// RESUME: define these
void reset_mem_aarray_vector_char(aarray_vector_char *aavr);
void add_mem_aarray_vector_char(aarray_vector_char *aavr);
void add_mem_aarray_vector_char_until(aarray_vector_char *aavr, int len);
void prep_aarray_vector_char(aarray_vector_char *aavr);
void extend_aarray_vector_char(aarray_vector_char *aavr);

vector_char vector_char_init(int dim);
void set_equal_char(char *v1, char *v2, int len);

void vectors_char_init(vectors_char *vrs, int mem, int dim);
char vectors_char_get(vectors_char *vrs, int elem, int coord);
void vectors_char_set(vectors_char *vrs, int elem, int coord, char val);
void prep_vectors_char(vectors_char *vrs);
void add_mem_vectors_char(vectors_char *vrs);
void add_mem_vectors_char_until(vectors_char *vrs, int lim);
void add2vectors_char(vectors_char *vrs, char *pt);
void add2vectors_char_zeros(vectors_char *vrs);
void remove_vectors_char(vectors_char *vrs, int vertex);
void free_vector_char(vector_char *vr);
void free_vectors_char(vectors_char *vrs);
void array_vector_char_init(array_vector_char *avr, int mem); 
void add_mem_array_vector_char(array_vector_char *avr);
void add_mem_array_vector_char_until(array_vector_char *avr, int lim);
void extend_array_vector_char(array_vector_char *avr, int dim);
void free_array_vector_char(array_vector_char *avr);
void add2array_vector_char(array_vector_char *avr, vector_char elem); // RESUME: define this

void xform_p_vector_char(vector_char *v1, vector_char *v3);
vector_char plus_vector_char(vector_char *v1, vector_char *v2);
vector_char plus_vectors_char(vectors_char *vs, vector_char *v2);
vector_char minus_vector_char(vector_char *v1, vector_char *v2);
vector_char multiply_vector_char(vector_char *v, char s);

char vector_char_dot(vector_char *v1, vector_char *v2);
char vectors_char_dot(vectors_char *vs, vector_char *v2, int i);
char vectors_chars_dot(vectors_char *v1s, vectors_char *v2s, int i, int j);

array_vectors_char array_vectors_char_init(int mem);
void add_mem_array_vectors_char(array_vectors_char *avr);
void add_mem_array_vectors_char_until(array_vectors_char *avr, int lim);
void extend_array_vectors_char(array_vectors_char *avr, int mem, int dim);
void free_array_vectors_char(array_vectors_char *avr);

double vector_real_norm(vector_real *vr);
double vector_int_norm(vector_int *vi);

#endif
