#ifndef NUMERICAL_H
#define NUMERICAL_H
#include "stdlib.h"
#include "string.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_permutation.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_poly.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_interp.h"
#include "gsl/gsl_rng.h"
#include "basics.h"
#include "vector.h"
#include "limits.h"
#include "float.h"
#define MAX_STEP_COUNT (int) 1e6
#define FIND_INTERIOR_POINT_MAX_ATTEMPTS 1000000000
#define INT_POWER_CUTOFF 2
#define MAX_AV_ITER 10
#define MAX_BALL_OVERLAP_ITER 5
#define INIT_SAMPLE_NUMBER 4
#define MAX_APPROX_VOLUME_SEARCH 10
#define N_EVAL_LIM 1e6
#define DK_ITER_MAX 10000
#define P1DS_E polynomial1D_static_evaluator
#define P3DS_E polynomial3D_static_evaluator
#define POLY1DS polynomial1D_static
#define STANDARD_TOL 1e-9
#define GSL_RNG_SEED 23

void display_c(_Complex double a);

gsl_vector *zero3d; 
extern double third_pi;

int int_power(double a, int n);
double cuberoot(double x);
// Note: the following routines can be incorporated in a separate program that produces parameter values in this file.
typedef struct box_parti3D
{
	struct box_parti3D *desc;
	double llim[3];
	double ulim[3];
	char state;
} box_partition3D;

void box_partition3D_init(double *llim, double *ulim, box_partition3D *bpart);

void box_partition3D_extend(box_partition3D *bpart);

double box_partition3D_type_measure(box_partition3D *bpart, char state);

void box_partition3D_split(box_partition3D *bpart, char state);
void refine_boundary_box_partition3D_ntimes(box_partition3D *bpart, char (*interior)(double *, void *), void *pars, int depth);
void refine_boundary_box_partition3D(box_partition3D *bpart, char (*interior)(double *, void *), void *pars);

void free_box_partition3D(box_partition3D *bpart);

typedef struct
{
	double val;
	double err;
} approx_double;

void find_interior_point(char (*interior)(double *x, void *), void *interior_pars, double *llim, double *ulim, double *guess);

approx_double approx_volume(char (*interior)(double *x, void *pars), void *pars, double epsilon, double *llim, double *ulim, double *int_pt); // TEST!

typedef struct pnode
{
	struct pnode *next;
	int data[3];
	double c;
} mnmial3D_node;

typedef struct
{
	// Hashing function
	int (*hf)(int *, double *pars);
	// Parameters
	double *hf_pars;
	// Hash table
	mnmial3D_node **table;
	int mem;
} mnmial3D_hashtable;

mnmial3D_hashtable mnmial3D_hashtable_init(int mem, int (*hf)(int *, double *), double *pars);
void free_mnmial3D_hashtable(mnmial3D_hashtable *ht);
void add2mnmial3D_hashtable(mnmial3D_hashtable *ht, int data[3], double val);

typedef struct
{
	double ***c;
	int degree[3];
	int **degree_z_xy;
	int *degree_y_x;
} polynomial3D_static;

typedef struct
{
	polynomial3D_static *p3d;
	double c;
	array_int subset;
} sample_constraint_polynomial3D_static;

typedef struct
{
	double *c;
	int degree;
} polynomial1D_static;

typedef struct
{
	vectors_int sequence;
	vectors_int monomial_generators;
	vector_real values;
} polynomial3D_static_evaluator; // P3DS_E

void polynomial3D_static_evaluator_init(polynomial3D_static_evaluator *p3dse, polynomial3D_static *p3d);

typedef struct
{
	vectors_int sequence;
	vector_int monomial_generators;
	vector_real values;
} polynomial1D_static_evaluator;

void polynomial1D_static_evaluator_init(polynomial1D_static_evaluator *p1dse, polynomial1D_static *p1d);

void polynomial3D_static_init(int *degree, polynomial3D_static *p3d);
void polynomial3D_static_set(polynomial3D_static *p3d, int i, int j, int k, double val);
double polynomial3D_static_get(polynomial3D_static *p3d, int i, int j, int k);
void polynomial3D_static_set_incr(polynomial3D_static *p3d, int i, int j, int k, double incr);
void free_polynomial3D_static(polynomial3D_static *p3d);
void read_polynomial3D_static(char *ifname, polynomial3D_static *p3d); // TEST
double polynomial3D_static_ev(polynomial3D_static *p3d, double *x); // TEST
void polynomial3D_static_grad(polynomial3D_static *p3d, polynomial3D_static **grad_p3d); // TEST
void polynomial3D_static_div(polynomial3D_static *v_p3d, polynomial3D_static *div_p3d); // TEST
char interior_polynomial3D_static_domain(double *x, void *pars); // TEST
char interior_multiple_polynomial3D_static_domain(double *x, void *pars); // TEST!
double gsl_vector_normsq(gsl_vector *v); // TEST

void Newton_multivar(void (*func)(double *, double *, void *), void *func_pars, void (*Dfunc)(double *, double **, void *), void *Dfunc_pars, double *guess, double tol, int dim);

void Newton_polynomial3D_static(polynomial3D_static *v_p3d, double *guess, double tol); // TEST
void polynomial3D_static_multiply(polynomial3D_static *A, polynomial3D_static *B, polynomial3D_static *C); // TEST
void power_linear_factor_polynomial3D_static(double *coeffs, int i, polynomial3D_static *p3d); // TEST
void scale_polynomial3D_static(polynomial3D_static *p3d, double factor); // TEST!
void add2polynomial3D_static(polynomial3D_static *p3d, polynomial3D_static *p3d_); // TEST!
void polynomial3D_static_affine_xform(polynomial3D_static *p3d, gsl_matrix *G, polynomial3D_static *xfp3d); // TEST

void display_polynomial3D_static(polynomial3D_static *p3d); // TEST
void transcribe_polynomial3D_static(polynomial3D_static *A, polynomial3D_static *B);
void recenter_polynomial3D_static(polynomial3D_static *p3d, double *x, polynomial3D_static *p3d_);
void truncate_polynomial3D_static(polynomial3D_static *p3d, int *multi_degree);

void polynomial3D_static_at_curve(polynomial3D_static *A, polynomial1D_static **curve, polynomial1D_static *B);
void polynomial3D_static_at_line(polynomial3D_static *A, double *coeffs, polynomial1D_static *B);
void polynomial3D_static_at2coords(polynomial3D_static *A, int odim, double v1, double v2, polynomial1D_static *B);

void polynomial3D_static_fit(polynomial3D_static *A, double **sample_pts, double *sample_vals, int N_samples, int max_degree, double *degree_weights);
void polynomial3D_static_lvl_set_fit(polynomial3D_static *A, double **sample_pts, int N_samples, int max_degree, double *degree_weights);


void polynomial1D_static_init(int degree, polynomial1D_static *p1d);
void free_polynomial1D_static(polynomial1D_static *p1d);
void polynomial1D_static_free(polynomial1D_static *p1d);
void transcribe_polynomial1D_static(polynomial1D_static *A, polynomial1D_static *B);
void polynomial1D_static_multiply(polynomial1D_static *A, polynomial1D_static *B, polynomial1D_static *C);
void polynomial1D_static_3D_static_multiply(polynomial1D_static *A, polynomial3D_static *B, polynomial3D_static *C, int dim);
void transcribe_polynomial1D_static_3D_static(polynomial1D_static *A, polynomial3D_static *B, int dim);
void transcribe_polynomial3D_static(polynomial3D_static *A, polynomial3D_static *B);
double polynomial1D_static_ev(polynomial1D_static *p1d, double x);
_Complex double polynomial1D_static_evc(polynomial1D_static *p1d, _Complex double x);
void polynomial1D_static_deriv(polynomial1D_static *p1d, polynomial1D_static *dp1d);

void polynomial1D_static_fit(polynomial1D_static *A, double *sample_pts, double *sample_vals, int N_samples, int max_degree);
void display_polynomial1D_static(polynomial1D_static *p1d); // TEST

void polynomial1D_static_DK_solve(polynomial1D_static *p1d, _Complex double *roots, double tol);

void Newton_polynomial1D_static(polynomial1D_static *f, double *guess, double tol);

// Methods for splines representing solids of revolution
void epsilon_nbhd_spline_rev(double (*inner)(double, void *), void *pars, double *x_samples, int N_samples, double epsilon, gsl_spline **outer, double tolerance);
void ball_overlap_spline_rev(double (*shape)(double, void *), void *pars, double ball_radius, double (*radial_dist)(double, void *), void *rd_pars, double *x_samples, int N_samples, gsl_spline **ovl, double tolerance);
char interior_multiple_domain(double *x, void *pars);

double d1polynomial1D_static_ev(double x, void *pars);
double d1spline_ev(double x, void *pars);
double d1func_offset_ev(double x, void *pars);

char interior_solid_revolution(double *x, void *pars);
char interior_solid_revolution_spline(double *x, void *pars);
char interior_solid_revolution_polynomial1D_static(double *x, void *pars);

typedef struct
{
	int dim;
	gsl_spline **spline;
} vector_spline;

void vector_spline_alloc(int dim, vector_spline *vspl);
void vector_spline_init(vector_spline *vspl, double *x_samples, double **v_samples, int N_samples, gsl_interp_type *gslt);
void vector_spline_free(vector_spline *vspl);

typedef struct
{
	double x[3];
	double r;
	double rsq;
} sphere;

char interior_sphere(double *x, void *pars);
void linspace(double *samples, double llim, double ulim, int N_samples);
void alloc_double(double **array, int len);

typedef struct multi_int
{
	struct multi_int *parts;
	double estimate;
	double confidence;
	double tol;
	int dim;
	array_bit 
} multivar_integrator;

// Operations between splines
void product_spline(gsl_spline *s1, gsl_spline *s2, gsl_spline **p12);
void quotient_spline(gsl_spline *s1, gsl_spline *s2, gsl_spline **q12);
void inverse_spline(gsl_spline *spl, gsl_spline **inv_spl);
void rescale_spline(gsl_spline **spl, double coeff);
void difference_spline(gsl_spline *s1, gsl_spline *s2, gsl_spline **s1ms2);
void divide_spline(gsl_spline **numerator, gsl_spline *denominator);
void derivative_spline(gsl_spline *f, gsl_spline **dfdx);
void norm_vector_spline(gsl_spline *curve, gsl_spline **norm_x, gsl_spline **norm_y);
void add_spline(gsl_spline *s1, gsl_spline *s2, gsl_spline **s1ps2);
void shared_domain_tally(gsl_spline *s1, gsl_spline *s2, double *xmin, double *xmax, int *t1, int *t2);
void apply_op_spline_coarse_fine(gsl_spline *s1, gsl_spline *s2, gsl_spline **s1ms2, double (*op)(double, double), double xmin, double xmax, int Nsamples);
void apply_op_spline_fine_coarse(gsl_spline *s2, gsl_spline *s1, gsl_spline **s2ms1, double (*op)(double, double), double xmin, double xmax, int Nsamples); 

char seek_domain_interior(char (*interior)(double *, void *), void *pars, double *llim, double *ulim, double *interior_pt, int depth); // TEST
approx_double integrate_2D(double (*f)(double, double, void *), void *pars, double (*domain)(double, double, void *), void *domain_pars, double *xbnds, double *ybnds, double tol); // TEST
double integrate_2D_fubini_fast(double (*f)(double, double, void *), void *fpars, double *xbnds, double (*y_min)(double, void *), void *y_min_pars, double (*y_max)(double, void *), void *y_max_pars, double tol);
double integrate_2D_fubini(double (*f)(double, double, void *), void *fpars, double *xbnds, double (*y_min)(double, void *), void *y_min_pars, double (*y_max)(double, void *), void *y_max_pars, double tol);

// Elementary functions in 1D and 2D
double test_func2D(double, double, void *);
double test_func2D_area(double y, double z, void *pars);
double quadratic_form_2D(double, double, void *);
double sqrt1D(double, void *);
double composite_2D(double, double, void *);
double composite_1D(double, void *);
double d2func_difference(double, double, void *);
double d2func_add(double, double, void *);
double d2func_linear(double, double, void *);
double d1func_linear(double, void *);
double d2func_product(double, double, void *);

// Random number functions
void parse_gslr_type(const gsl_rng_type **gslr_type, char *gslr_mode); // RESUME: test this function!

double monomial_eval_basic(double *pt, int *deg, int dim);
double monomial_eval1D(double x, int degree);

// Inversion

double invert_bisection_method(double a, void *pars);

// Interfaces to apply various gsl routines


#endif
