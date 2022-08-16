#include "numerical.h"
#define NEWTON_ITER_MAX 100000

void display_c(_Complex double a)
{
	if (__real__ a != 0)
	{
		printf("%g", __real__ a);
		if (__imag__ a > 0)
		{
			printf(" + %g i", __imag__ a);
		} 
		else if (__imag__ a < 0)
		{
			printf(" - %g i", - __imag__ a);
		}
		printf("\n");
	}
	else if (__imag__ a != 0)
	{
		printf("%g i\n", __imag__ a);
	}
}

char probably_interior = 1;
char probably_boundary = 2;
char probably_exterior = 3;
double third_pi = M_PI / 3.;
char axes[3] = {'x', 'y', 'z'};

int corner_case_count[3];

void alloc_double(double **array, int len)
{
	(*array) = (double *) calloc(len, sizeof(double));
}

void linspace(double *samples, double llim, double ulim, int N_samples)
{
	double delz = (ulim - llim) / (N_samples - 1);
	samples[0] = llim;
	int im1 = 0;
	for (int i = 1; i < N_samples; i++)
	{
		samples[i] = samples[im1] + delz;
		im1 = i;
	}
}

void box_partition3D_init(double *llim, double *ulim, box_partition3D *bpart)
{
	for (int i = 0; i < 3; i++)
	{
		(*bpart).llim[i] = llim[i];
		(*bpart).ulim[i] = ulim[i];
	}
	(*bpart).state = 0;
	(*bpart).desc = NULL;
	return;
}

void box_partition3D_extend(box_partition3D *bpart)
{
	(*bpart).desc = (box_partition3D *) calloc(8, sizeof(box_partition3D));
	double center[3];
	for (int i = 0; i < 3; i++) center[i] = ((*bpart).llim[i] + (*bpart).ulim[i]) * 0.5;
	//printf("Center: %g %g %g\n", center[0], center[1], center[2]);
	double llims[8][3];
	double ulims[8][3];
	// double vols[8];
	for (int i = 0; i < 8; i++)
	{
		llims[i][0] = !((i >> 2) & 1)? (*bpart).llim[0]: center[0];
		llims[i][1] = !((i >> 1) & 1)? (*bpart).llim[1]: center[1];
		llims[i][2] = !(i & 1)? (*bpart).llim[2]: center[2];
		ulims[i][0] = !((i >> 2) & 1) ? center[0]: (*bpart).ulim[0];
		ulims[i][1] = !((i >> 1) & 1)? center[1]: (*bpart).ulim[1];
		ulims[i][2] = !(i & 1)? center[2]: (*bpart).ulim[2];
		//vols[i] = 1;
		//for (int ii = 0; ii < 3; ii++) vols[i] *= (ulims[i][ii] - llims[i][ii]);
		// printf("%g ", vols[i]);
	}
	// printf("\n");
	for (int i = 0; i < 8; i++)
	{
		box_partition3D_init(llims[i], ulims[i], &(*bpart).desc[i]);
		(*bpart).desc[i].state = (*bpart).state;
	}
	(*bpart).state = -1;
}

int box_partition3D_tally(box_partition3D *bpart)
{
	if ((*bpart).state == -1)
	{
		int subtally = 0;
		for (int i = 0; i < 8; i++)
		{
			subtally += box_partition3D_tally(&(*bpart).desc[i]);
		}
		return subtally;
	}
	else
	{
		return 1;
	}
}

double box_partition3D_type_measure(box_partition3D *bpart, char state)
{
	if ((*bpart).state == state)
	{
		double vol = (*bpart).ulim[0] - (*bpart).llim[0];
		for (int i = 1; i < 3; i++) vol *= (*bpart).ulim[i] - (*bpart).llim[i];
		return vol;
	}
	else 
	{
		double vol = 0;
		if ((*bpart).state == -1)
		{
			for (int i = 0; i < 8; i++)
			{
				vol += box_partition3D_type_measure(&(*bpart).desc[i], state);
			}
		}
		return vol;
	}
}

void box_partition3D_split(box_partition3D *bpart, char state) // TEST
{
	if ((*bpart).state == -1)
	{
		for (int i = 0; i < 8; i++)
		{
			box_partition3D_split(&(*bpart).desc[i], state);
		}
	}
	else if ((*bpart).state == state)
	{
		box_partition3D_extend(bpart);
	}
}

void refine_boundary_box_partition3D_ntimes(box_partition3D *bpart, char (*interior)(double *, void *), void *pars_, int depth)
{
	if (depth > 0)
	{
		refine_boundary_box_partition3D(bpart, interior, pars_);
		int depthm1 = depth - 1;
		for (int i = 0; i < 8; i++)
		{
			if ((*bpart).desc[i].state == probably_boundary)
			{
				refine_boundary_box_partition3D_ntimes(&(*bpart).desc[i], interior, pars_, depthm1);
			}
		}
	}
}

// This method is somewhat inefficient, but gets the job done
void refine_boundary_box_partition3D(box_partition3D *bpart, char (*interior)(double *, void *), void *pars_)
{
	if ((*bpart).state == -1)
	{
		for (int i = 0; i < 8; i++)
		{
			refine_boundary_box_partition3D(&((*bpart).desc[i]), interior, pars_);
		}
	}
	else if ((*bpart).state == probably_boundary)
	{
		box_partition3D_extend(bpart);
		for (int i = 0; i < 8; i++)
		{
			unsigned char corner_status = 0;
			for (int ii = 0; ii < 8; ii++)
			{
				double corner[3];
				for (int iii = 0; iii < 3; iii++) corner[iii] = (ii >> iii) & 1? (*bpart).desc[i].ulim[iii]: (*bpart).desc[i].llim[iii];
				if (interior(corner, pars_))
				{
					corner_status |= (1 << ii);
				}	
			}
			if ((corner_status != 0) && (corner_status != 255))
			{
				// printf("corner status: %d\n", (int) corner_status);
				(*bpart).desc[i].state = probably_boundary;
				corner_case_count[1] += 1;
			}
			else 
			{
				if (corner_status == 255) 
				{
					corner_case_count[2] += 1;
					(*bpart).desc[i].state = probably_interior;
				} 
				else 
				{
					corner_case_count[0] += 1;
					(*bpart).desc[i].state = probably_exterior;
				}
			}
		}
	}
}

void write_corner_case_count(char *diagfname, int iter_num)
{
	FILE *ofile = fopen(diagfname, "a");
	if (ofile != NULL)
	{
		fprintf(ofile, "%d ", iter_num);
		for (int i = 0; i < 3; i++) fprintf(ofile, "%d ", corner_case_count[i]);
		fprintf(ofile, "\n");
		fclose(ofile);
	}
}

char seek_domain_interior(char (*interior)(double *, void *), void *pars, double *llim, double *ulim, double *interior_pt, int depth) // TEST
{
	double center[3];
	for (int i = 0; i < 3; i++) center[i] = 0.5 * (llim[i] + ulim[i]);
	//printf("%g %g %g\n", center[0], center[1], center[2]);
	char status = interior(center, pars);
	if (!status)
	{
		if (depth > 0)
		{
			double rllim[3];
			double rulim[3];
			for (int i = 0; i < 8; i++)
			{
				for (int dim = 0; dim < 3; dim++)
				{
					rllim[dim] = (i >> dim) & 1? llim[i]: center[i];
					rulim[dim] = (i >> dim) & 1? center[i]: ulim[i];
				}
				status = seek_domain_interior(interior, pars, rllim, rulim, interior_pt, depth - 1);
				if (status) break;
			}
		}
	}
	else
	{
		interior_pt[0] = center[0];
		interior_pt[1] = center[1];
		interior_pt[2] = center[2];
	}
	return status;
}

// RESUME: check this!
void find_interior_point_random(char (*interior)(double *, void *), void *interior_pars, double *llim, double *ulim, double *guess)
{
	int N_attempts = 0;
	gsl_rng *gslr = gsl_rng_alloc(gsl_rng_taus);
	gsl_rng_set(gslr, GSL_RNG_SEED); 
	//void find_interior_point(char (*interior)(double *x, void *), void *interior_pars, double *llim, double *ulim, double *guess);
	double lens[3] = {ulim[0] - llim[0], ulim[1] - llim[1], ulim[2] - llim[2]};
	while (N_attempts < FIND_INTERIOR_POINT_MAX_ATTEMPTS)
	{
		double x[3];
		x[0] = llim[0] + lens[0] * gsl_rng_uniform(gslr);
		x[1] = llim[1] + lens[1] * gsl_rng_uniform(gslr);
		x[2] = llim[2] + lens[2] * gsl_rng_uniform(gslr);
		if (!interior(x, interior_pars)) {}
		else
		{
			guess[0] = x[0];
			guess[1] = x[1];
			guess[2] = x[2];
			gsl_rng_free(gslr);
			return;
		}
		N_attempts += 1;
	}
	guess[0] = llim[0] - 1;
	guess[1] = llim[1];
	guess[2] = llim[2];
}

// INCOMPLETE: test and improve performance
// Optional diagnostic info: corner tally to improve algorithm (e.g. to anticipate structure)
//	Consider changing this to incorporate boundary structure, and maybe using simplices instead of or in addition to blocks.
//	Also, check the interior point search function.
approx_double approx_volume(char (*interior)(double *, void *), void *pars_, double epsilon, double *llim, double *ulim, double *int_pt) 
{
	char interior_corner = 0;
	double corner[3];
	for (int i = 0; i < 8; i++)
	{
		printf("Checking if corner %d is interior\n", i);
		corner[0] = (i >> 2) & 1? llim[0]: ulim[0];
		corner[1] = (i >> 1) & 1? llim[1]: ulim[1];
		corner[2] = i & 1? llim[2]: ulim[2];
		if (interior(corner, pars_))
		{
			interior_corner = 1;
			printf("interior corner found\n");
			break;
		}
	}
	printf("(done)\n");
	if (interior_corner) 
	{
		printf("Found interior corner! Commencing with volume calculation\n");
	}
	else 
	{
		double interior_pt[3];
		if (int_pt == NULL)
		{
			// Find a single interior point
			double rect_volume = 1;
			for (int i = 0; i < 3; i++) rect_volume *= (ulim[i] - llim[i]);
			printf("Volume of enclosing rectangle: %g\n", rect_volume);
			int search_scale = (int) (0.8 * log(cuberoot(rect_volume / (epsilon * epsilon * epsilon))))+ 1;
			int search_count = 0;
			char interior_found = 0;
			printf("Seeking interior point at scale 2^-%d\n", search_scale);
			printf("%g %g %g %g %g %g\n", llim[0], llim[1], llim[2], ulim[0], ulim[1], ulim[2]);
			interior_found = seek_domain_interior(interior, pars_, llim, ulim, interior_pt, search_scale); 
			printf("(done)\n");
			if (interior_found) 
			{}
			else
			{
				printf("Unable to find interior of domain in specified rectangle\n");
				approx_double est;
				est.val = 0;
				est.err = 1;
				return est;
			}
		}
		else
		{
			interior_pt[0] = int_pt[0];
			interior_pt[1] = int_pt[1];
			interior_pt[2] = int_pt[2];
		}
		// Partition the original rectangle into subdomains 
		double rllims[3];
		double rulims[3];
		approx_double approx_vol;
		approx_vol.val = 0;
		approx_vol.err = 0;
		for (int i = 0; i < 8; i++)
		{
			rllims[0] = (i >> 2) & 1? llim[0]: interior_pt[0];
			rulims[0] = (i >> 2) & 1? interior_pt[0]: ulim[0];
			rllims[1] = (i >> 1) & 1? llim[1]: interior_pt[1];
			rulims[1] = (i >> 1) & 1? interior_pt[1]: ulim[1];
			rllims[2] = i & 1? llim[2]: interior_pt[2];
			rulims[2] = i & 1? interior_pt[2]: ulim[2];
			// Add the approximate volumes from each subdomain
			printf("Computing domain volume of octant %d\n", i);
			approx_double incr = approx_volume(interior, pars_, epsilon, rllims, rulims, interior_pt);
			printf("(done: vol = %g)\n", incr.val);
			approx_vol.val += incr.val;
			approx_vol.err += incr.err * incr.val;
		}
		approx_vol.err = (approx_vol.val > 0? approx_vol.err / approx_vol.val: 0);
		return approx_vol;
	}
	void **pars = (void **) pars_;
	box_partition3D bpart;
	box_partition3D_init(llim, ulim, &bpart);
	bpart.state = probably_boundary;
	//printf("Upper and lower bounds: %g %g %g %g %g %g\n", bpart.llim[0], bpart.llim[1], bpart.llim[2], bpart.ulim[0], bpart.ulim[1], bpart.ulim[2]);
	double encl_vol = 1;
	for (int i = 0; i < 3; i++) encl_vol *= (bpart.ulim[i] - bpart.llim[i]);
	double lower_est = 0;
	double bndry_est = box_partition3D_type_measure(&bpart, probably_boundary);
	// FILE *ccc_ofile = fopen(diagfname, "w");
	// if (ccc_ofile != NULL) fclose(ccc_ofile);
	int iter_num = 0;
	while (bndry_est >= lower_est * epsilon && iter_num < MAX_AV_ITER)
	{
		for (int i = 0; i < 256; i++) corner_case_count[i] = 0;
		int depth = (int) (2 * log(bndry_est / (epsilon * encl_vol)));
		depth = (depth < 1? 1: depth);
		printf("Refining boundary to scale 2^(-%d)\n", depth);
		refine_boundary_box_partition3D_ntimes(&bpart, interior, pars_, depth);
		// write_corner_case_count(diagfname, iter_num);
		printf("Computing measures\n");
		lower_est = box_partition3D_type_measure(&bpart, probably_interior);
		bndry_est = box_partition3D_type_measure(&bpart, probably_boundary);
		int tally = box_partition3D_tally(&bpart);
		//printf("%g %g %g %d\n", lower_est, bndry_est, lower_est * epsilon, tally);
		iter_num += 1;
	}
	//if (iter_num == MAX_AV_ITER) printf("Maximum volume estimation iteration exceeded\n");
	//printf("final: %g %g\n", lower_est, bndry_est);
	approx_double result;
	result.val = lower_est > 0? lower_est + 0.5 * bndry_est: 0;
	result.err = lower_est > 0? bndry_est / result.val: bndry_est;
	free_box_partition3D(&bpart);
	return result;
}

void free_box_partition3D(box_partition3D *bpart)
{
	if ((*bpart).state == -1)
	{
		for (int i = 0; i < 8; i++)
		{
			free_box_partition3D(&(*bpart).desc[i]);
		}
		free((*bpart).desc);
		(*bpart).desc = NULL;
		(*bpart).state = 0;
	}
}

// Check this!
void polynomial3D_static_init(int *degree, polynomial3D_static *p3d)
{
	(*p3d).degree[0] = degree[0];
	(*p3d).degree[1] = degree[1];
	(*p3d).degree[2] = degree[2];
	int dxp1 = degree[0] + 1;
	int dyp1 = degree[1] + 1;
	int dzp1 = degree[2] + 1;
	(*p3d).c = (double ***) calloc(dxp1, sizeof(double **));
	(*p3d).degree_z_xy = (int **) calloc(dxp1, sizeof(int *));
	(*p3d).degree_y_x = (int *) calloc(dxp1, sizeof(int));
	for (int i = 0; i <= degree[0]; i++)
	{
		(*p3d).c[i] = (double **) calloc(dyp1, sizeof(double *));
		(*p3d).degree_y_x[i] = 0;
		(*p3d).degree_z_xy[i] = (int *) calloc(dyp1, sizeof(int));
		for (int ii = 0; ii <= degree[1]; ii++)
		{
			(*p3d).c[i][ii] = (double *) calloc(dzp1, sizeof(double));
			(*p3d).degree_z_xy[i][ii] = 0;
			for (int iii = 0; iii <= degree[2]; iii++) (*p3d).c[i][ii][iii] = 0;
		}
	}
	return;
}

void free_polynomial3D_static(polynomial3D_static *p3d)
{
	for (int i = 0; i <= (*p3d).degree[0]; i++)
	{
		for (int ii = 0; ii <= (*p3d).degree[1]; ii++)
		{
			free((*p3d).c[i][ii]);
		}
		free((*p3d).degree_z_xy[i]);
		free((*p3d).c[i]);
	}
	free((*p3d).degree_z_xy);
	free((*p3d).degree_y_x);
	free((*p3d).c);
}

void polynomial3D_static_set_incr(polynomial3D_static *p3d, int i, int j, int k, double incr)
{
	if (incr != 0)
	{
		if (k <= (*p3d).degree_z_xy[i][j]) {}
		else
		{
			(*p3d).degree_z_xy[i][j] = k;
		}
		if (j <= (*p3d).degree_y_x[i]) {}
		else
		{
			(*p3d).degree_y_x[i] = j;
		}
	}
	(*p3d).c[i][j][k] += incr;
}

void polynomial3D_static_set(polynomial3D_static *p3d, int i, int j, int k, double val)
{
	if (val != 0)
	{
		if (k <= (*p3d).degree_z_xy[i][j]) {}
		else
		{
			(*p3d).degree_z_xy[i][j] = k;
		}
		if (j <= (*p3d).degree_y_x[i]) {}
		else
		{
			(*p3d).degree_y_x[i] = j;
		}
	}
	else
	{
		if (k != (*p3d).degree_z_xy[i][j]) {}
		else if ((*p3d).degree_z_xy[i][j] > 0)
		{
			do
			{
				(*p3d).degree_z_xy[i][j] -= 1;
			} while ((*p3d).c[i][j][(*p3d).degree_z_xy[i][j]] == 0 && (*p3d).degree_z_xy[i][j] > 0);
		}
		if (j != (*p3d).degree_y_x[i]) {}
		else if ((*p3d).degree_z_xy[i][j] == 0)
		{
			while ((*p3d).degree_z_xy[i][(*p3d).degree_y_x[i]] == 0 && (*p3d).degree_y_x[i] > 0)
			{
				(*p3d).degree_y_x[i] -= 1;
			}
		}
	}
	(*p3d).c[i][j][k] = val;
}

double polynomial3D_static_get(polynomial3D_static *p3d, int i, int j, int k)
{
	return (*p3d).c[i][j][k];
}



void polynomial3D_static_at_line(polynomial3D_static *A, double *coeffs, polynomial1D_static *B)
{
	polynomial1D_static **line = (polynomial1D_static **) calloc(3, sizeof(polynomial1D_static *));
	for (int i = 0; i < 3; i++)
	{
		line[i] = (polynomial1D_static *) calloc(1, sizeof(polynomial1D_static));
		polynomial1D_static_init(1, line[i]);
		(*line[i]).c[0] = coeffs[i];
		(*line[i]).c[1] = coeffs[i + 3];
	}
	polynomial3D_static_at_curve(A, line, B);
	for (int i = 0; i < 3; i++) 
	{
		free_polynomial1D_static(line[i]);
		free(line[i]);
	}
	free(line);
}

// RESUME: implement this after settling on a good algorithm to optimize the evaluation sequence.
//		At present, a sort of incremental 'greedy search' seems like it would be the simplest
//		approach (i.e. iteratively building the evaluation sequence by adding 
//		the monomial that requires the fewest additional steps to evaluate.)
//		There are probably many cases where this is approach is suboptimal, but it should 
//			be reasonably effective in low-degree settings (e.g. with polynomials whose
//			degrees are no greater than 10 or so.)
void polynomial3D_static_evaluator_init(polynomial3D_static_evaluator *p3dse, polynomial3D_static *p3d)
{
	// First determine common differences in multi degree among nonzero coefficients
	//	multi-degree differences that exceed half the degree along each dimension are ignored,
	//	as these can only generate at most one new term (and can probably be added at the end if necessary)
	// Find a minimal set of differences and starting monomials that produce a fully connected 'graph' of nonzero coefficients
	//	Loop over sets of differences that produce a connected graph
	//	Reduce the complexity of each set by finding relationships between differences.
	//		Express 'high norm' differences in terms of lower norm ones.
	//		Determine which differences are represented in the set of nonzero monomials
	//		Look for monomials that are multiples of generators
	//	Each pairwise difference, at an instance level, is 'essentially equivalent' to one additional term.
	//	Together with 'root' monomials, a polynomial should be able to be reconstructed with a set of linearly
	//		independent differences whose tallies add up to the total number of nonzero monomials. Even better
	//		if these differences include monomials themselves.
	//	A difference value can be assigned an 'importance' based on the number of monomial terms that it can produce,
	//		as well as the computational cost of evaluating it (a sort of 'opportunity cost').
	//		Each 'importance framework' may produce a different proposed evaluation scheme.
	//	Monomial differences are necessary; but the order in which they are evaluated can vary (as can auxiliary monomials.)
	//		e.g.: a polynomial with degree support 1, 2, 4, 8, 16.  1st differences: 1, 2, 4, 8, 3, 7, 15*, 6, 14*, 12*,
	//									2nd differences: 1, 2, 4, 5, 3, 9*, 7, 6, 13*, 11*, 8, 10*, 14*, 12*
	//	idea: starting from the trivial evaluation sequence, where each monomial is evaluated independently of the others, apply elementary 
	//		transformations until an efficient evaluation sequence has been found, reducing the total number of independent generators
	//		(while preserving the set of monomials that can be generated.)  e.g. if m1 = m2 + m3, then m1 can be replaced with an evaluation 
	//		step involving m2 and m3.  More generally, if m1 = m2 + a12, then m1 can be replaced with a12 and an evaluation
	//		step involving m2 and a12.  There are therefore a large number of possible substitions; are all reduction methods guaranteed
	//		to eventually reach a common set of generators and evaluations?  Only if evaluation sequences can be simplified concurrently with
	//		generators.  Consider for example the degree sequence 1, 2, 4, 16 -> 1, 2, 4, 12, (4 + 12) -> 1, 2, 4, 8, (4 + 8), (4 + (4 + 8))
	//		where in the last sequence the evaluation step (4 + 8) can be eliminated, as (4 + (4 + 8)) can be expressed as (8 + 8).  Let us 
	//		consider reducing the degree of terms further in a ridiculous way, just to make the necessity of simplification clear: 
	//		'' -> 1, 2, 3, (1 + 3), 8, ((1 + 3) + 8), ((1 + 3) + ((1 + 3) + 8)) -> 1, 2, 3, (1 + 3), 5, (5 + 3), ((1 + 3) + (5 + 3)), 
	//											((1 + 3) + ((1 + 3) + (5 + 3))),
	//		'' -> 1, 2, 3, (1 + 3), (2 + 3), ((2 + 3) + 3), ((1 + 3) + ((2 + 3) + 3)), ((1 + 3) + ((1 + 3) + ((2 + 3) + 3))),
	//		'' -> etc.  
	//		We would like to be able to recognize that some of the terms are redundant or unnecessary.  For example, that (4 + (4 + 8)) is
	//		(8 + 8), and so the evaluation sequence (4 + 8) can be dropped.  This is mainly to avoid having to rely on exhaustive enumeration,
	//		which could get out of hand in the higher dimensional case.  The resolution could be that, periodically at least, monomial terms 
	//		that are already expressed as part of an evaluation sequence can be 'revisited' in terms of new generators, or even evaluation series.
	//		-> 1, 2, 3, (2 + 2), (2 + 3), ((2 + 3) + 3), ((2 + 2) + ((2 + 3) + 3)), (((2 + 3) + 3) + ((2 + 3) + 3)),
	//	Update: (1/8/2022) It seems likely that a 'perfectly optimal' approach to evaluating polynomials is impractical (or more perhaps more accurately,
	//		it simply eludes me.)  At least in the multivariable case (the single variable version of the problem is a bit simpler, but still
	//		has its own subtleties.)  The evaluation sequence must traverse the M non-vanishing monomial terms in some order, and auxiliary eval. 
	//		terms can be rearranged so that they appear immediately before the first monomial that depends on them; in this way, one can interpret
	//		an evaluation sequence as a process of iteratively adding new monomials sequentially to a growing list.  However, it is not immediately
	//		obvious that the problem 'localizes', or that the global optimum involves the shortest sequence with each new monomial term that is 
	//		added (except at the very end.)  The problem is 'frustrated', in the sense that at the beginning there are benefits to both minimizing
	//		the number of auxiliary evaluations as well as using the initial constructive process to generate auxiliar evaluations that are broadly
	//		useful later, in evaluating later monomials in the fewest number of steps (one may think of these as auxiliary evaluations that are
	//		in any case necessary in the later stages of the optimal sequence, and which might as well be used at the beginning as well.)
	//		There are even simple one-dimensional cases with several global minima.  The full complexity of the instance space has yet to be fully
	//		mapped out, but I suspect that it is a somewhat 'glassy' problem, especially in higher dimensions.  That being said, there is probably
	//		no shortage of 'scaffold-based' algorithms with optimal scaling, but which 'miss' the global minimum.  Still, it would be nice to have
	//		more insight into this question... Is it as frustrated as it seems, or is there more structure to it that I'm missing?
	//			
	int dp1[3];
	int half_degree[3];
	for (int i = 0; i < 3; i++) 
	{
		dp1[i] = (*p3d).degree[i] + 1;
		half_degree[i] = (*p3d).degree[i] / 2;
	}
	int counts[dp1[0]][dp1[1]][dp1[2]];
	aarray_int nonzero;
	aarray_int_init(&nonzero, 1);
	for (int i = 0; i <= (*p3d).degree[0]; i++) 
	{
		for (int ii = 0; ii <= (*p3d).degree[1]; ii++)
		{
			for (int iii = 0; iii <= (*p3d).degree[2]; iii++)
			{
				counts[i][ii][iii] = 0;
				if ((*p3d).c[i][ii][iii] != 0)
				{
					int len = nonzero.len;
					extend_aarray_int(&nonzero);
					add2array_int(&(nonzero.e[len]), i);
					add2array_int(&(nonzero.e[len]), ii);
					add2array_int(&(nonzero.e[len]), iii);
				}
			}
		}
	}
	for (int i = 0; i < nonzero.len; i++)
	{
		for (int ii = 0; ii < i; ii++)
		{
			int diff[3];
			char sat = 1;
			for (int j = 0; j < 3; j++) 
			{
				diff[j] = nonzero.e[i].e[j] - nonzero.e[ii].e[j];
				if (diff[j] > -half_degree[j] && diff[j] < half_degree[j]) {}
				else
				{
					sat = 0;
					break;
				}
			}
			if (sat)
			{
				for (int j = 0; j < 3; j++) diff[j] += half_degree[j];
				counts[diff[0]][diff[1]][diff[2]] += 1;
			}
		}
	}
	free_aarray_int(&nonzero);
}

// RESUME: Check this!
void polynomial3D_static_at_curve(polynomial3D_static *A, polynomial1D_static **curve, polynomial1D_static *B)
{
	// B[deg_x_i + deg_y_ii + deg_z_iii] += (*A).c[i][ii][iii] * (deg_x, i) * (deg_y, ii) * (deg_z, iii);
	// B += (*A).c[i][ii][iii] * (curve_x)^i * (curve_y)^ii * (curve_z)^iii
	int deg_B = (*A).degree[0] * (*curve[0]).degree + (*A).degree[1] * (*curve[1]).degree + (*A).degree[2] * (*curve[2]).degree;
	polynomial1D_static_init(deg_B, B);
	polynomial1D_static *factor_x, *factor_xy, *factor_xyz, *factor_x_, *factor_xy_, *factor_xyz_;
	factor_x = (polynomial1D_static *) calloc(1, sizeof(polynomial1D_static));
	factor_xy = (polynomial1D_static *) calloc(1, sizeof(polynomial1D_static));
	factor_xyz = (polynomial1D_static *) calloc(1, sizeof(polynomial1D_static));
	factor_x_ = (polynomial1D_static *) calloc(1, sizeof(polynomial1D_static));
	factor_xy_ = (polynomial1D_static *) calloc(1, sizeof(polynomial1D_static));
	factor_xyz_ = (polynomial1D_static *) calloc(1, sizeof(polynomial1D_static));
	polynomial1D_static_init(0, factor_x);
	(*factor_x).c[0] = 1;
	for (int i = 0; i <= (*A).degree[0]; i++)
	{
		transcribe_polynomial1D_static(factor_x, factor_xy);
		for (int ii = 0; ii <= (*A).degree_y_x[i]; ii++)
		{
			transcribe_polynomial1D_static(factor_xy, factor_xyz);
			for (int iii = 0; iii <= (*A).degree_z_xy[i][ii]; iii++)
			{
				for (int ci = 0; ci <= (*factor_xyz).degree; ci++)
				{
					(*B).c[ci] += (*factor_xyz).c[ci] * (*A).c[i][ii][iii];
				}
				polynomial1D_static *aux_xyz;
				polynomial1D_static_multiply(factor_xyz, curve[2], factor_xyz_);
				polynomial1D_static_free(factor_xyz);
				aux_xyz = factor_xyz;
				factor_xyz = factor_xyz_;
				factor_xyz_ = aux_xyz;
			}
			polynomial1D_static_free(factor_xyz);
			polynomial1D_static *aux_xy;
			polynomial1D_static_multiply(factor_xy, curve[1], factor_xy_);
			polynomial1D_static_free(factor_xy);
			aux_xy = factor_xy;
			factor_xy = factor_xy_;
			factor_xy_ = aux_xy;
		}
		polynomial1D_static_free(factor_xy);
		polynomial1D_static *aux_x;
		polynomial1D_static_multiply(factor_x, curve[0], factor_x_);
		polynomial1D_static_free(factor_x);
		aux_x = factor_x;
		factor_x = factor_x_;
		factor_x_ = aux_x;
	}
	//printf("Attempting to free auxiliary factors\n");
	polynomial1D_static_free(factor_x);
	//printf("freed factor_x\n");
	//polynomial1D_static_free(factor_xy);
	//printf("freed factor_xy\n");
	//polynomial1D_static_free(factor_xyz);
	//printf("freed factor_xyz.  Attempting to free auxiliary factor pointers\n");
	free(factor_x);
	free(factor_xy);
	free(factor_xyz);
	free(factor_x_);
	free(factor_xy_);
	free(factor_xyz_);
	// Determine the degree of 'B'
	while ((*B).c[(*B).degree] == 0)
	{
		(*B).degree -= 1;
	}
}

void polynomial3D_static_at2coords(polynomial3D_static *A, int odim, double v1, double v2, polynomial1D_static *B)
{
	int dim1, dim2;
	switch (odim)
	{
		case 0:
			dim1 = 1;
			dim2 = 2;
			break;
		case 1:
			dim1 = 0;
			dim2 = 2;
			break;
		case 2:
			dim1 = 0;
			dim2 = 1;
			break;
	}
	polynomial1D_static_init((*A).degree[odim], B);
	for (int i = 0; i <= (*A).degree[odim]; i++) 
	{
		double v1_ii = 1;
		for (int ii = 0; ii <= (*A).degree[dim1]; ii++) 
		{
			double v1_ii_v2_iii = v1_ii;
			for (int iii = 0; iii <= (*A).degree[dim2]; iii++)
			{
				switch (odim)
				{
					case 0:
						(*B).c[i] += (*A).c[i][ii][iii] * v1_ii_v2_iii;
						break;
					case 1:
						(*B).c[i] += (*A).c[ii][i][iii] * v1_ii_v2_iii;
						break;

					case 2:
						(*B).c[i] += (*A).c[ii][iii][i] * v1_ii_v2_iii;
						break;
				}
				v1_ii_v2_iii *= v2;
			}
			v1_ii *= v1;
		}
	}
}

void transcribe_polynomial3D_static(polynomial3D_static *A, polynomial3D_static *B)
{
	polynomial3D_static_init((*A).degree, B);
	for (int i = 0; i <= (*A).degree[0]; i++) for (int ii = 0; ii <= (*A).degree[1]; ii++) for (int iii = 0; iii <= (*A).degree[2]; iii++)
	{
		polynomial3D_static_set(B, i, ii, iii, polynomial3D_static_get(A, i, ii, iii));
	}
}

void read_polynomial3D_static(char *ifname, polynomial3D_static *p3d)
{
	FILE *ifile = fopen(ifname, "r");
	if (ifile != NULL)
	{
		// Seek the polynomial header (with degree info)
		char buf[256];
		int EOF_ = 0;
		while (strcmp(buf, "POLYNOMIAL3D") != 0 && EOF_ != EOF)
		{
			EOF_ = fscanf(ifile, "%s\n", buf);
		}
		if (EOF_ == EOF)
		{
			fclose(ifile);
			return;
		}
		int degree[3];
		fscanf(ifile, "%d %d %d", &degree[0], &degree[1], &degree[2]);
		polynomial3D_static_init(degree, p3d);
		// Read coefficients term by term
		while (strcmp(buf, "END") != 0 && EOF_ != EOF)
		{
			EOF_ = fscanf(ifile, "%s", buf);
			if (EOF_ == EOF || strcmp(buf, "END") == 0) break;
			if (strcmp(buf, "+") == 0) {}
			else
			{
				double c;
				int px, py, pz;
				sscanf(buf, "%lgx^%dy^%dz^%d", &c, &px, &py, &pz);
				polynomial3D_static_set(p3d, px, py, pz, c);
			}
		}
		fclose(ifile);
		return;
	}
	else
	{
		return;
	}
}

double polynomial3D_static_ev(polynomial3D_static *p3d, double *x)
{
	double result = 0;
	double xfactor = 1;
	for (int i = 0; i <= (*p3d).degree[0]; i++)
	{
		double yfactor = 1;
		for (int ii = 0; ii <= (*p3d).degree_y_x[i]; ii++)
		{
			double zfactor = 1, xyfactor = xfactor * yfactor;
			for (int iii = 0; iii <= (*p3d).degree_z_xy[i][ii]; iii++)
			{
				if ((*p3d).c[i][ii][iii] != 0) result += xyfactor * zfactor * (*p3d).c[i][ii][iii]; 
				zfactor *= x[2];
			}
			yfactor *= x[1];
		}
		xfactor *= x[0];
	}
	return result;
}

double monomial_eval_basic(double *pt, int *deg, int dim)
{
	double r = 1.0;
	for (int i = 0; i < dim; i++)
	{
		int rem = deg[i];
		double factor = 1.0;
		double elem = pt[i];
		while (rem > 0)
		{
			if (rem & 1)
			{
				factor *= elem;
			}
			elem *= elem;
			rem >>= 1;
		}
		r *= factor;
	}
	return r;
}

void polynomial3D_static_fit(polynomial3D_static *A, double **sample_pts, double *sample_vals, int N_samples, int max_degree, double *degree_weights)
{
	// Requires linear algebra: sum_i (f(x_i; C) - y_i)^2 + -> sum_i (f(x_i; C) - y_i) x_i^{\alpha} = 0, x_i^{\alpha} \equiv \prod_j x_{ij}^{\alpha_j}
	//	i.e. sum_\beta C_\beta sum_i x_i^{\beta + \alpha} = sum_i y_i x_i^{\alpha}, \forall \alpha
	//		M C = b; M_{\alpha \beta} \equiv sum_i x_i^{\beta + \alpha}, b \equiv sum_i y_i x_i^{\alpha}
	//	LU C = b; U C = Linv b, 
	int max_degrees[3];
	for (int i = 0; i < 3; i++)
	{
		max_degrees[i] = (int) ((double) max_degree / degree_weights[i]);
	}
	polynomial3D_static_init(max_degrees, A);
	for (int i = 0; i < 3; i++) max_degrees[i] += 1;
	int max_degree_x = max_degrees[0];
	int term_map[max_degrees[0]][max_degrees[1]][max_degrees[2]];
	for (int i = 0; i < max_degrees[0]; i++) for (int ii = 0; ii < max_degrees[1]; ii++) for (int iii = 0; iii < max_degrees[2]; iii++)
	{
		term_map[i][ii][iii] = -1;
	}
	aarray_int terms;
	aarray_int_init_precise(&terms, 16, 3);
	int N_terms = 0;
	double wtx = 0;
	for (int i = 0; i < max_degree_x; i++)
	{
		double wty = wtx;
		int max_degree_y = (int) (((double) max_degree - wtx) / degree_weights[1]) + 1;
		for (int ii = 0; ii < max_degree_y; ii++)
		{
			double wtz = wty;
			int max_degree_z = (int) (((double) max_degree - wty) / degree_weights[2]) + 1;
			for (int iii = 0; iii < max_degree_z; iii++)
			{
				if (wtz < max_degree)
				{
					term_map[i][ii][iii] = N_terms;
					extend_aarray_int(&terms);
					add_mem_array_int_until(&(terms.e[N_terms]), 3);
					terms.e[N_terms].e[0] = i;
					terms.e[N_terms].e[1] = ii;
					terms.e[N_terms].e[2] = iii;
					terms.e[N_terms].len = 3;
					N_terms = terms.len;
				}
				wtz += degree_weights[2];
			}
			wty += degree_weights[1];
		}
		wtx += degree_weights[0];
	}
	gsl_matrix *gslm = gsl_matrix_calloc(N_terms, N_terms);
	gsl_vector *b = gsl_vector_calloc(N_terms);
	for (int i = 0; i < N_samples; i++)
	{
		double pow_x_beta = 1.;
		for (int ti = 0; ti < N_terms; ti++)
		{
			double pow_x_alpha = monomial_eval_basic(sample_pts[i], terms.e[ti].e, 3);
			for (int tii = 0; tii <= ti; tii++)
			{
				double pow_x_beta = monomial_eval_basic(sample_pts[i], terms.e[tii].e, 3);
				gsl_matrix_set(gslm, ti, tii, gsl_matrix_get(gslm, ti, tii) + pow_x_beta * pow_x_alpha);
			}
			gsl_vector_set(b, ti, gsl_vector_get(b, ti) + pow_x_alpha * sample_vals[i]);
		}
	}
	for (int ti = 0; ti < N_terms; ti++)
	{
		for (int tii = 0; tii < ti; tii++)
		{
			gsl_matrix_set(gslm, tii, ti, gsl_matrix_get(gslm, ti, tii));
		}
	}
	gsl_vector *x = gsl_vector_calloc(N_terms);
	gsl_permutation *pi = gsl_permutation_calloc(N_terms);
	int sign;
	gsl_linalg_LU_decomp(gslm, pi, &sign);
	gsl_linalg_LU_solve(gslm, pi, b, x);
	for (int i = 0; i < N_terms; i++)
	{
		(*A).c[terms.e[i].e[0]][terms.e[i].e[1]][terms.e[i].e[2]] = gsl_vector_get(x, i);
	}
	gsl_permutation_free(pi);
	gsl_vector_free(x);
	gsl_vector_free(b);
	free_aarray_int(&terms);
	gsl_matrix_free(gslm);
}

//	RESUME: (polynomial3D_static_lvl_set_fit) This function needs some kind of regularizer, which may require making several assumptions 
//		about the structure of the points in question.
//			- Connectivity of sample points (i.e. a 'mesh' structure rather than scattered samples)
//			- Connected components of the level set
//		Unfortunately, the structure that seems most relevant would also make this a highly non-linear problem.
//		Instead, it might be possible to restrict the class of polynomials that are used to construct level sets,
//			as well as the sample surfaces to which the level set algorithm is applicable.  
//			For example, one might consider approximating shapes or surfaces with reflection symmetry, and 
//			limiting polynomials to those with even degree terms, and require a "positive definite" 
//			result (with negative constant) within a specific degree range.  It isn't obvious to me that this
//			approach would converge, however (restricting the solution space to positive coefficients of
//			positive definite monomials could be overly constrictive for general convex shapes (almost surely 
//			for non-convex shapes). It would be interesting to learn what shapes or surfaces can be approximated
//			to arbitrary precision this way, however.)
void polynomial3D_static_lvl_set_fit(polynomial3D_static *A, double **sample_pts, int N_samples, int max_degree, double *degree_weights)
{
	int max_degrees[3];
	for (int i = 0; i < 3; i++)
	{
		max_degrees[i] = (int) ((double) max_degree / degree_weights[i]);
	}
	polynomial3D_static_init(max_degrees, A);
	for (int i = 0; i < 3; i++) max_degrees[i] += 1;
	int max_degree_x = max_degrees[0];
	int term_map[max_degrees[0]][max_degrees[1]][max_degrees[2]];
	for (int i = 0; i < max_degrees[0]; i++) for (int ii = 0; ii < max_degrees[1]; ii++) for (int iii = 0; iii < max_degrees[2]; iii++)
	{
		term_map[i][ii][iii] = -1;
	}
	aarray_int terms;
	aarray_int_init_precise(&terms, 16, 3);
	int N_terms = 0;
	double wtx = 0;
	for (int i = 0; i < max_degree_x; i++)
	{
		double wty = wtx;
		int max_degree_y = (int) (((double) max_degree - wtx) / degree_weights[1]) + 1;
		for (int ii = 0; ii < max_degree_y; ii++)
		{
			double wtz = wty;
			int max_degree_z = (int) (((double) max_degree - wty) / degree_weights[2]) + 1;
			for (int iii = 0; iii < max_degree_z; iii++)
			{
				if (wtz < max_degree)
				{
					term_map[i][ii][iii] = N_terms;
					extend_aarray_int(&terms);
					add_mem_array_int_until(&(terms.e[N_terms]), 3);
					terms.e[N_terms].e[0] = i;
					terms.e[N_terms].e[1] = ii;
					terms.e[N_terms].e[2] = iii;
					terms.e[N_terms].len = 3;
					N_terms = terms.len;
				}
				wtz += degree_weights[2];
			}
			wty += degree_weights[1];
		}
		wtx += degree_weights[0];
	}
	gsl_matrix *gslm = gsl_matrix_calloc(N_terms, N_terms);
	gsl_vector *b = gsl_vector_calloc(N_terms);
	int ti0 = term_map[0][0][0];
	gsl_matrix_set(gslm, ti0, ti0, 1.0);
	gsl_vector_set(b, ti0, -1.);
	for (int i = 0; i < N_samples; i++)
	{
		double pow_x_beta = 1.;
		for (int ti = 0; ti < N_terms; ti++)
		{
			double pow_x_alpha = monomial_eval_basic(sample_pts[i], terms.e[ti].e, 3);
			for (int tii = 0; tii <= ti; tii++)
			{
				double pow_x_beta = monomial_eval_basic(sample_pts[i], terms.e[tii].e, 3);
				gsl_matrix_set(gslm, ti, tii, gsl_matrix_get(gslm, ti, tii) + pow_x_beta * pow_x_alpha);
			}
		}
	}
	for (int ti = 0; ti < N_terms; ti++)
	{
		for (int tii = 0; tii < ti; tii++)
		{
			gsl_matrix_set(gslm, tii, ti, gsl_matrix_get(gslm, ti, tii));
		}
	}
	gsl_vector *x = gsl_vector_calloc(N_terms);
	gsl_permutation *pi = gsl_permutation_calloc(N_terms);
	int sign;
	gsl_linalg_LU_decomp(gslm, pi, &sign);
	gsl_linalg_LU_solve(gslm, pi, b, x);
	for (int i = 0; i < N_terms; i++)
	{
		polynomial3D_static_set(A, terms.e[i].e[0], terms.e[i].e[1], terms.e[i].e[2], gsl_vector_get(x, i));
	}
	gsl_permutation_free(pi);
	gsl_vector_free(x);
	gsl_vector_free(b);
	free_aarray_int(&terms);
	gsl_matrix_free(gslm);

}

// NEWISH
void polynomial3D_static_grad(polynomial3D_static *p3d, polynomial3D_static **grad_p3d)
{
	(*grad_p3d) = (polynomial3D_static *) calloc(3, sizeof(polynomial3D_static));
	int gdegree[3] = {(*p3d).degree[0], (*p3d).degree[1], (*p3d).degree[2]};
	gdegree[0] -= 1;
	polynomial3D_static_init(gdegree, &((*grad_p3d)[0]));
	gdegree[0] += 1;
	gdegree[1] -= 1;
	polynomial3D_static_init(gdegree, &((*grad_p3d)[1]));
	gdegree[1] += 1;
	gdegree[2] -= 1;
	polynomial3D_static_init(gdegree, &((*grad_p3d)[2]));
	for (int i = 0; i <= (*p3d).degree[0]; i++)
	{
		(*grad_p3d)[1].degree_y_x[i] = (*p3d).degree_y_x[i] - 1;
		for (int ii = 0; ii <= (*p3d).degree_y_x[i]; ii++)
		{
			for (int iii = 0; iii <= (*p3d).degree_z_xy[i][ii]; iii++)
			{
				if (i < (*p3d).degree[0])
				{
					polynomial3D_static_set(&(*grad_p3d)[0], i, ii, iii, (*p3d).c[i + 1][ii][iii] * (i + 1));
				}
				if (ii < (*p3d).degree_y_x[i])
				{
					polynomial3D_static_set(&(*grad_p3d)[1], i, ii, iii, (*p3d).c[i][ii + 1][iii] * (ii + 1));
				}
				if (iii < (*p3d).degree_z_xy[i][ii])
				{
					polynomial3D_static_set(&(*grad_p3d)[2], i, ii, iii, (*p3d).c[i][ii][iii + 1] * (iii + 1));
				}
			}
		}
	}
}

// NEWISH
void polynomial3D_static_div(polynomial3D_static *v_p3d, polynomial3D_static *div_p3d)
{
	int div_degree[3];
	for (int i = 0; i < 3; i++)
	{
		div_degree[i] = (v_p3d[0]).degree[i];
		for (int ii = 1; ii < 3; ii++)
		{
			div_degree[i] = (div_degree[i] >= (v_p3d[ii]).degree[i]? div_degree[i]: (v_p3d[ii]).degree[i]);
		}
	}
	polynomial3D_static_init(div_degree, div_p3d);
	for (int dim = 0; dim < 3; dim++)
	{
		int dimp1m3 = (dim + 1) % 3;
		int dimp2m3 = (dim + 2) % 3;
		for (int i = 0; i < (v_p3d[dim]).degree[dim]; i++)
		{
			int ip1 = i + 1;
			for (int ii = 0; ii <= (v_p3d[dim]).degree[dimp1m3]; ii++)
			{
				for (int iii = 0; iii <= (v_p3d[dim]).degree[dimp2m3]; iii++)
				{
					switch (dim)
					{
						case 0:
							polynomial3D_static_set_incr(div_p3d, i, ii, iii, (v_p3d[0]).c[ip1][ii][iii] * ip1);
							break;
						case 1:
							polynomial3D_static_set_incr(div_p3d, iii, i, ii, (v_p3d[1]).c[iii][ip1][ii] * ip1);
							break;
						case 2:
							polynomial3D_static_set_incr(div_p3d, ii, iii, i, (v_p3d[2]).c[ii][iii][ip1] * ip1);
							break;
					}
				}
			}
		}
	}
}

// NEWISH
char interior_polynomial3D_static_domain(double *x, void *pars_)
{
	void **pars = (void **) pars_;
	polynomial3D_static *p3d;
	p3d = (polynomial3D_static *) pars[0];
	double C = pars[1] != NULL? *((double *) pars[1]): 0;
	double p3d_value = polynomial3D_static_ev(p3d, x);
	if (p3d_value <= C) return 1;
	else return 0;
}

char interior_multiple_polynomial3D_static_domain(double *x, void *pars_)
{
	void **pars = (void **) pars_;
	polynomial3D_static **p3d;
	p3d = (polynomial3D_static **) pars[0];
	double *C = (double *) pars[1];
	int N_p3d = *((int *) pars[2]);
	for (int i = 0; i < N_p3d; i++)
	{
		if (polynomial3D_static_ev(p3d[i], x) > C[i]) return 0;
	}
	return 1;
}

// NEW
// NOTE: depending on how vectors are accessed, this method may be unnecessarily slow
double gsl_vector_normsq(gsl_vector *v)
{
	int len = (*v).size;
	double normvsq = 0;
	for (int i = 0; i < (*v).size; i++)
	{
		double vi = gsl_vector_get(v, i);
		normvsq += vi * vi;
	}
	return normvsq;
}

// Newton's method:
//	f(x) \approx Df.(x - x0) -> x0 approx x - inv(Df).f(x)
void Newton_multivar(void (*func)(double *, double *, void *), void *func_pars, void (*Dfunc)(double *, double **, void *), void *Dfunc_pars, double *guess, double tol, int dim)
{
	gsl_matrix *Dfunc_eval = gsl_matrix_alloc(dim, dim);
	double **Dfunc_eval_array = bb_matrix_alloc(dim, dim);
	gsl_vector *func_eval = gsl_vector_alloc(dim);
	double func_eval_array[dim];
	double tolsq = tol * tol;
	double err = 2 * tolsq;
	gsl_vector *delx = gsl_vector_alloc(dim);
	gsl_permutation *pi = gsl_permutation_alloc(dim);
	gsl_permutation_init(pi);
	int Nsteps = 0;
	int sample_rate = 1;
	while (err > tolsq && Nsteps < MAX_STEP_COUNT) 
	{
		if (Nsteps % sample_rate == 0)
		{
			printf("Newton's method, step %d, err = %g\n", Nsteps, err);
		}
		Nsteps += 1;
		// evaluate the derivative and function at the current best guess
		func(guess, func_eval_array, func_pars);
		Dfunc(guess, Dfunc_eval_array, Dfunc_pars);
		for (int i = 0; i < dim; i++)
		{
			gsl_vector_set(func_eval, i, func_eval_array[i]);
			for (int ii = 0; ii < dim; ii++)
			{
				gsl_matrix_set(Dfunc_eval, i, ii, Dfunc_eval_array[i][ii]);
			}
		}
		// compute the displacement to the next estimated root
		int sign_pi;
		int status = gsl_linalg_LU_decomp(Dfunc_eval, pi, &sign_pi);
		gsl_linalg_LU_solve(Dfunc_eval, pi, func_eval, delx);
		for (int i = 0; i < dim; i++) guess[i] -= gsl_vector_get(delx, i);
		err = gsl_vector_normsq(delx) + gsl_vector_normsq(func_eval);
	}
	bb_matrix_free(Dfunc_eval_array, dim);
	gsl_matrix_free(Dfunc_eval);
	gsl_vector_free(func_eval);
	gsl_vector_free(delx);
	gsl_permutation_free(pi);
	printf("(done)\n");
}

// RESUME: Test this!
void Newton_polynomial3D_static(polynomial3D_static *v_p3d, double *guess, double tol)
{
	gsl_matrix *Dv_p3d_eval = gsl_matrix_alloc(3, 3); 
	gsl_vector *v_p3d_eval = gsl_vector_alloc(3);
	polynomial3D_static *Dv_p3d[3];
	for (int i = 0; i < 3; i++)
	{
		polynomial3D_static_grad(&(v_p3d[i]), &Dv_p3d[i]);
	}
	double tolsq = tol * tol;
	double err = 2 * tolsq;
	gsl_vector *delx = gsl_vector_alloc(3);
	gsl_permutation *pi = gsl_permutation_alloc(3);
	gsl_permutation_init(pi);
	int Nsteps = 0;
	int sample_rate = 1;
	while (err > tolsq && Nsteps < MAX_STEP_COUNT) 
	{
		Nsteps += 1;
		// evaluate the derivative and function at the current best guess
		for (int i = 0; i < 3; i++)
		{
			gsl_vector_set(v_p3d_eval, i, polynomial3D_static_ev(&v_p3d[i], guess));
			for (int ii = 0; ii < 3; ii++)
			{
				gsl_matrix_set(Dv_p3d_eval, i, ii, polynomial3D_static_ev(&(Dv_p3d[i][ii]), guess));
			}
		}
		if (Nsteps % sample_rate == 0)
		{
			printf("Newton's method, step %d, err = %g\n", Nsteps, err);
		}
		// compute the displacement to the next estimated root
		int sign_pi;
		int status = gsl_linalg_LU_decomp(Dv_p3d_eval, pi, &sign_pi);
		gsl_linalg_LU_solve(Dv_p3d_eval, pi, v_p3d_eval, delx);
		for (int i = 0; i < 3; i++) guess[i] -= gsl_vector_get(delx, i);
		err = gsl_vector_normsq(delx) + gsl_vector_normsq(v_p3d_eval);
	}
	for (int i = 0; i < 3; i++) 
	{
		for (int ii = 0; ii < 3; ii++) free_polynomial3D_static(&(Dv_p3d[i][ii]));
		free((Dv_p3d[i]));
	}
	gsl_matrix_free(Dv_p3d_eval);
	gsl_vector_free(v_p3d_eval);
	gsl_vector_free(delx);
	gsl_permutation_free(pi);
	printf("(done)\n");
}

// RESUME: test this!
void polynomial3D_static_multiply(polynomial3D_static *A, polynomial3D_static *B, polynomial3D_static *C)
{
	int cdegree[3] = {(*A).degree[0] + (*B).degree[0], (*A).degree[1] + (*B).degree[1], (*A).degree[2] + (*B).degree[2]};
	polynomial3D_static_init(cdegree, C);
	for (int i = 0; i <= (*A).degree[0]; i++)
	{
		for (int ii = 0; ii <= (*A).degree_y_x[i]; ii++)
		{
			for (int iii = 0; iii <= (*A).degree_z_xy[i][ii]; iii++)
			{
				if ((*A).c[i][ii][iii] != 0) {}
				else
				{
					continue;
				}
				for (int j = 0; j <= (*B).degree[0]; j++)
				{
					int ipj = i + j;
					for (int jj = 0; jj <= (*B).degree_y_x[i]; jj++)
					{
						int iipjj = ii + jj;
						for (int jjj = 0; jjj <= (*B).degree_z_xy[i][ii]; jjj++)
						{
							if ((*B).c[j][jj][jjj] != 0) {}
							else continue;
							int iiipjjj = iii + jjj;
							(*C).c[ipj][iipjj][iiipjjj] += (*A).c[i][ii][iii] * (*B).c[j][jj][jjj];
							(*C).degree_y_x[ipj] = iipjj <= (*C).degree_y_x[ipj]? (*C).degree_y_x[ipj]: iipjj;
							(*C).degree_z_xy[ipj][iipjj] = iiipjjj <= (*C).degree_z_xy[ipj][iipjj]? (*C).degree_z_xy[ipj][iipjj]: iiipjjj;
						}
					}
				}
			}
		}
	}
	return;
}

// RESUME: test this!
int int_power(double a, int n)
{
	double result;
	if (n < INT_POWER_CUTOFF)
	{
		result = (n > 0? a: 1);
		for (int i = 1; i < n; i++)
		{
			result *= a;
		}
		return result;
	}
	else
	{
		int remainder = n % 2;
		int hn = n / 2;
		result = int_power(a, hn);
		double resultsq = result * result;
		if (remainder != 0)
		{
			double rresult = int_power(a, remainder);
			return resultsq * rresult;
		} 
		else return resultsq;
	}
}

// RESUME: test this!
void power_linear_factor_polynomial3D_static(double *coeffs, int i, polynomial3D_static *xfx_i)
{
	if (coeffs[2] != 0) {}
	else
	{
		exit(EXIT_FAILURE);
		return;
	}	
	int degree[3] = {i, i, i};
	polynomial3D_static_init(degree, xfx_i);
	double cx = 1;
	double cz_max = int_power(coeffs[2], i); 
	int init_trinomial_coeff = 1; // i! / (ip!ipp!ippp!) -> i! / (ip! (i - ip)!)
	for (int ip = 0; ip <= i; ip++)
	{
		int imip = i - ip;
		double cy = 1;
		double cz = cz_max;
		int trinomial_coeff = init_trinomial_coeff;
		double cxyz = cz * cx;
		for (int ipp = 0; ipp <= imip; ipp++)
		{
			int ippp = imip - ipp;
			(*xfx_i).c[ip][ipp][ippp] = cxyz * trinomial_coeff;
			(*xfx_i).degree_z_xy[ip][ipp] = ippp;
			cxyz *= coeffs[1];
			trinomial_coeff *= ippp; // C / (ip! ipp! (i - ip - ipp)!) -> C / (ip! (ipp + 1)! (i - ip - ipp - 1)!);
			trinomial_coeff /= (ipp + 1);
		}
		(*xfx_i).degree_y_x[ip] = imip;
		cx *= coeffs[0];
		cz_max /= coeffs[2];
		init_trinomial_coeff *= imip; // C/((ip + 1)!(i - ip - 1)!);
		init_trinomial_coeff /= (ip + 1);
	}
	return;
}

// NEW
void scale_polynomial3D_static(polynomial3D_static *p3d, double factor)
{
	for (int i = 0; i <= (*p3d).degree[0]; i++) for (int ii = 0; ii <= (*p3d).degree[1]; ii++) for (int iii = 0; iii <= (*p3d).degree[2]; iii++)
	{
		(*p3d).c[i][ii][iii] *= factor;
	}
}

// NEW
void add2polynomial3D_static(polynomial3D_static *p3d, polynomial3D_static *p3d_)
{
	if ((*p3d).degree[0] >= (*p3d_).degree[0] && (*p3d).degree[1] >= (*p3d_).degree[1] && (*p3d).degree[2] >= (*p3d_).degree[2]) {}
	else
	{
		exit(EXIT_FAILURE);
		return;
	}
	for (int i = 0; i <= (*p3d_).degree[0]; i++) for (int ii = 0; ii <= (*p3d_).degree[1]; ii++) for (int iii = 0; iii <= (*p3d_).degree[2]; iii++)
	{
		(*p3d).degree_y_x[i] = (*p3d).degree_y_x[i] >= ii? (*p3d).degree_y_x[i]: ii;
		(*p3d).degree_z_xy[i][ii] = (*p3d).degree_z_xy[i][ii] >= iii? (*p3d).degree_z_xy[i][ii]: iii;
		(*p3d).c[i][ii][iii] += (*p3d_).c[i][ii][iii];
	}
}

// NEW (Note: the calculations of single variable powers seem a bit sketchy to me.)
void polynomial3D_static_affine_xform(polynomial3D_static *p3d, gsl_matrix *G, polynomial3D_static *xfp3d)
{
	int max_degree = -1;
	for (int i = 0; i < 3; i++) 
	{
		max_degree = ((*p3d).degree[i] <= max_degree? max_degree: (*p3d).degree[i]);
	}
	printf("Max degree = %d\n", max_degree);
	display_polynomial3D_static(p3d);
	printf("%d %d %d\n", (*p3d).degree[0], (*p3d).degree[1], (*p3d).degree[2]);
	int rdegree[3] = {max_degree, max_degree, max_degree};
	polynomial3D_static_init(rdegree, xfp3d);
	
	int init_degree[3] = {0, 0, 0};
	int init_gen_degree[3] = {1, 1, 1};
	
	polynomial3D_static xyzgens[3]; 
	for (int i = 0; i < 3; i++)
	{
		polynomial3D_static_init(init_gen_degree, &xyzgens[i]);
		for (int ii = 0; ii < 3; ii++)
		{
			(xyzgens[i]).c[0][0][0] = gsl_matrix_get(G, 3, i);
			(xyzgens[i]).c[1][0][0] = gsl_matrix_get(G, 0, i);
			(xyzgens[i]).c[0][1][0] = gsl_matrix_get(G, 1, i);
			(xyzgens[i]).c[0][0][1] = gsl_matrix_get(G, 2, i);
		}
		printf("%c generator:\n", axes[i]);
		display_polynomial3D_static(&xyzgens[i]);
	}
	printf("Performing rotation\n");
	polynomial3D_static *xfx_i = (polynomial3D_static *) calloc(1, sizeof(polynomial3D_static));
	polynomial3D_static *xfx_i_p = (polynomial3D_static *) calloc(1, sizeof(polynomial3D_static));
	polynomial3D_static_init(init_degree, xfx_i);
	(*xfx_i).c[0][0][0] = 1;
	for (int i = 0; i <= (*p3d).degree[0]; i++)
	{
		
		polynomial3D_static *xfy_ii = (polynomial3D_static *) calloc(1, sizeof(polynomial3D_static));
		polynomial3D_static *xfy_ii_p = (polynomial3D_static *) calloc(1, sizeof(polynomial3D_static));
		polynomial3D_static_init(init_degree, xfy_ii);
		(*xfy_ii).c[0][0][0] = 1;
		for (int ii = 0; ii <= (*p3d).degree_y_x[i]; ii++)
		{
			polynomial3D_static xfx_ixfy_ii;
			polynomial3D_static_multiply(xfx_i, xfy_ii, &xfx_ixfy_ii);
			polynomial3D_static *xfz_iii = (polynomial3D_static *) calloc(1, sizeof(polynomial3D_static));
			polynomial3D_static *xfz_iii_p = (polynomial3D_static *) calloc(1, sizeof(polynomial3D_static));
			polynomial3D_static_init(init_degree, xfz_iii);
			(*xfz_iii).c[0][0][0] = 1;
			for (int iii = 0; iii <= (*p3d).degree_z_xy[i][ii]; iii++)
			{
				if ((*p3d).c[i][ii][iii] != 0) {}
				else
				{
					continue;
				}
				polynomial3D_static incr;
				polynomial3D_static_multiply(&xfx_ixfy_ii, xfz_iii, &incr);
				scale_polynomial3D_static(&incr, (*p3d).c[i][ii][iii]); 
				add2polynomial3D_static(xfp3d, &incr); 
				if (iii < (*p3d).degree_z_xy[i][ii])
				{
					printf("Test: term %d %d %d\n", i, ii, iii);
					polynomial3D_static *aux_z;
					// polynomial3D_static *xfz_iii_p = (polynomial3D_static *) calloc(1, sizeof(polynomial3D_static));
					polynomial3D_static_multiply(xfz_iii, &xyzgens[2], xfz_iii_p);
					free_polynomial3D_static(xfz_iii);
					// free(xfz_iii);
					aux_z = xfz_iii;
					xfz_iii = xfz_iii_p;
					xfz_iii_p = aux_z;
				}
			}
			free_polynomial3D_static(xfz_iii);
			free(xfz_iii);
			free(xfz_iii_p);
			if (ii < (*p3d).degree[1])
			{
				//polynomial3D_static *xfy_ii_p = (polynomial3D_static *) calloc(1, sizeof(polynomial3D_static)); // a
				polynomial3D_static_multiply(xfy_ii, &xyzgens[1], xfy_ii_p);
				free_polynomial3D_static(xfy_ii);
				//free(xfy_ii); // a
				polynomial3D_static *aux_y;
				aux_y = xfy_ii;
				xfy_ii = xfy_ii_p;
				xfy_ii_p = aux_y;
			}
			free_polynomial3D_static(&xfx_ixfy_ii);
		}
		free_polynomial3D_static(xfy_ii);
		free(xfy_ii);
		free(xfy_ii_p); // r
		if (i < (*p3d).degree[0])
		{
			// polynomial3D_static *xfx_i_p = (polynomial3D_static *) calloc(1, sizeof(polynomial3D_static)); 
			polynomial3D_static_multiply(xfx_i, &xyzgens[0], xfx_i_p);
			free_polynomial3D_static(xfx_i);
			// free(xfx_i);
			polynomial3D_static *aux_x = xfx_i;
			xfx_i = xfx_i_p;
			xfx_i_p = aux_x;
		}
	}
	free_polynomial3D_static(xfx_i);
	free(xfx_i);
	free(xfx_i_p);
	for (int i = 0; i < 3; i++)
	{
		free_polynomial3D_static(&xyzgens[i]);
	}
	return;
}

void display_polynomial3D_static(polynomial3D_static *p3d)
{
	char lt_set = 0;
	char allzero = 1;
	for (int i = 0; i <= (*p3d).degree[0]; i++) for (int ii = 0; ii <= (*p3d).degree[1]; ii++) for (int iii = 0; iii <= (*p3d).degree[2]; iii++)
	{
		if ((*p3d).c[i][ii][iii] != 0)
		{
			allzero = 0;
			if (lt_set)
			{
				printf(" + ");
			}
			else
			{
				lt_set = 1;
			}
			printf("%g", (*p3d).c[i][ii][iii]);
			if (i > 1) printf("x^%d", i);
			else if (i == 1) printf("x");
			if (ii > 1) printf("y^%d", ii);
			else if (ii == 1) printf("y");
			if (iii > 1) printf("z^%d", iii);
			else if (iii == 1) printf("z");
		}
	}
	if (allzero) printf("0");
	printf("\n");
}

void recenter_polynomial3D_static(polynomial3D_static *p3d, double *x, polynomial3D_static *p3d_)
{
	polynomial3D_static_init((*p3d).degree, p3d_);
	polynomial1D_static xgen, ygen, zgen;
	polynomial1D_static_init(1, &xgen);
	polynomial1D_static_init(1, &ygen);
	polynomial1D_static_init(1, &zgen);
	xgen.c[0] = x[0];
	xgen.c[1] = 1;
	ygen.c[0] = x[1];
	ygen.c[1] = 1;
	zgen.c[0] = x[2];
	zgen.c[1] = 1;
	polynomial1D_static *xfactor = (polynomial1D_static *) calloc(1, sizeof(polynomial1D_static)); 
	polynomial1D_static_init(0, xfactor);
	(*xfactor).c[0] = 0;
	polynomial1D_static *xfactor_ = (polynomial1D_static *) calloc(1, sizeof(polynomial1D_static));
	polynomial3D_static *xyfactor = (polynomial3D_static *) calloc(1, sizeof(polynomial3D_static)); 
	polynomial3D_static *xyfactor_ = (polynomial3D_static *) calloc(1, sizeof(polynomial3D_static));
	polynomial3D_static *xyzfactor = (polynomial3D_static *) calloc(1, sizeof(polynomial3D_static));
	polynomial3D_static *xyzfactor_ = (polynomial3D_static *) calloc(1, sizeof(polynomial3D_static));
	for (int i = 0; i <= (*p3d).degree[0]; i++)
	{
		transcribe_polynomial1D_static_3D_static(xfactor, xyfactor, 0); 
		for (int ii = 0; ii <= (*p3d).degree_y_x[i]; ii++)
		{
			transcribe_polynomial3D_static(xyfactor, xyzfactor); 
			for (int iii = 0; iii <= (*p3d).degree_z_xy[i][ii]; iii++)
			{
				polynomial3D_static incr;
				transcribe_polynomial3D_static(xyzfactor, &incr);
				scale_polynomial3D_static(&incr, (*p3d).c[i][ii][iii]);
				add2polynomial3D_static(p3d_, &incr);
				if (iii < (*p3d).degree[2])
				{
					polynomial1D_static_3D_static_multiply(&zgen, xyzfactor, xyzfactor_, 2); 
					free_polynomial3D_static(xyzfactor);
					polynomial3D_static *aux = xyzfactor;
					xyzfactor = xyzfactor_;
					xyzfactor_ = aux;
				}
			}
			if (ii < (*p3d).degree_y_x[i])
			{
				polynomial1D_static_3D_static_multiply(&ygen, xyfactor, xyfactor_, 1); // 
				free_polynomial3D_static(xyfactor);
				polynomial3D_static *aux = xyfactor;
				xyfactor = xyfactor_;
				xyfactor_ = aux;
			}
		}
		if (i < (*p3d).degree[0])
		{
			polynomial1D_static_multiply(&xgen, xfactor, xfactor_);
			free_polynomial1D_static(xfactor);
			polynomial1D_static *aux = xfactor;
			xfactor = xfactor_;
			xfactor_ = aux;
		}
	}
	free_polynomial1D_static(xfactor);
	free_polynomial3D_static(xyfactor);
	free_polynomial3D_static(xyzfactor);
	free_polynomial1D_static(&xgen);
	free_polynomial1D_static(&ygen);
	free_polynomial1D_static(&zgen);
}

void truncate_polynomial3D_static(polynomial3D_static *p3d, int *multi_degree)
{
	double wts[3] = {1. / multi_degree[0], 1. / multi_degree[1], 1. / multi_degree[2]};
	double xyzwts[3] = {0, 0, 0};
	for (int i = 0; i <= (*p3d).degree[0]; i++) 
	{
		xyzwts[1] = xyzwts[0];
		for (int ii = 0; ii <= (*p3d).degree[1]; ii++) 
		{
			xyzwts[2] = xyzwts[1];
			for (int iii = 0; iii <= (*p3d).degree[2]; iii++)
			{
				if (xyzwts[2] < 1) {}
				else
				{
					(*p3d).c[i][ii][iii] = 0;
				}
				xyzwts[2] += wts[2];
			}
			xyzwts[1] += wts[1];
		}
		xyzwts[0] += wts[0];
	}
}

double sample_pt_energy_function(double **coords, int Nsamples, void **pars)
{
	sample_constraint_polynomial3D_static *constr = (sample_constraint_polynomial3D_static *) pars[0];
	int *Nconstraints = (int *) pars[1];
	double *lambdas = (double *) pars[2];
	nbrlist *mesh = (nbrlist *) pars[3];
	double E = 0;
	for (int i = 0; i < Nsamples; i++)
	{
		for (int nii = 0; nii < (*mesh).v.e[i].len; nii++)
		{
			int ii = (*mesh).v.e[i].e[nii];
			if (ii < i) {}
			else continue;
			double delx, dely, delz, delsq;
			delx = coords[i][0] - coords[ii][0];
			dely = coords[i][1] - coords[ii][1];
			delz = coords[i][2] - coords[ii][2];
			delsq = delx * delx + dely * dely + delz * delz;
			E += delsq;
		}
	}
	for (int ci = 0; ci < (*Nconstraints); ci++)
	{
		for (int ii = 0; ii < constr[ci].subset.len; ii++)
		{
			int i = constr[ci].subset.e[ii];
			E += lambdas[ci] * (polynomial3D_static_ev(constr[ci].p3d, coords[i]) - constr[ci].c);
		}
	}
	return E;
}

// TEST
double cuberoot(double x) 
{
	double r3x = x > 1? x: 1; // x - x^3 / 3y + 1/3 = x, 1 - x^2 / y, x^3 = y, (eta)^3 = 1 / y^2
	double err = 1;
	double invx = 1. / x;
	while (err > 1e-14 * x)
	{
		double r3xsq = r3x * r3x;
		double eval = r3xsq * r3x;
		double slope = 3 * r3xsq;
		err = eval - x; // f(x) approx (x - r) * f'(x)
		double delx = err / slope;
		r3x -= delx;
		err = err >= 0? err: -err;
		err += delx >= 0? delx: -delx;
	}
	return r3x;
}

// NOTE: before pursuing this further, devise a general principle for partitioning the surface into graphable components
//	RESUME:
double surface_area_polynomial3D_static_levelset(polynomial3D_static* p3d, double c, double *llim, double *ulim, double tol)
{
	double tolsq = tol * tol;
	// Define sample points on the surface within bounding box
	sample_constraint_polynomial3D_static *constr = (sample_constraint_polynomial3D_static *) calloc(7, sizeof(sample_constraint_polynomial3D_static));
	constr[6].p3d = p3d;
	constr[6].c = c;
	int N_samples_init  = 1. / tolsq;
	vectors_real samples;
	vectors_real_init(&samples, N_samples_init, 3);
	//	Solve for intersection of lines with the polynomial surface
	//		(here, complex arithmetic would help avoid unnecessary delays from the iteration counter.  
	//			Consider using the Durand-Kerner or Aberth methods.)
	//	Determine nearest neighbors along the surface
	//		
	//	Iteratively minimize energy function over sample point coordinates and refine the resulting mesh
	// Triangulate or interpolate in a nbhd of each point
	// Estimate the enclosed area
	// while ave_sqerror > tolsq
	//	refine and relax mesh
	//	compute area, 
	return 0;
}

void polynomial1D_static_init(int degree, polynomial1D_static *p1d)
{
	(*p1d).c = (double *) calloc(degree + 1, sizeof(double));
	for (int i = 0; i <= degree; i++)
	{
		(*p1d).c[i] = 0;
	}
	(*p1d).degree = degree;
}

_Complex double polynomial1D_static_evc(polynomial1D_static *p1d, _Complex double x)
{
	if ((*p1d).degree >= 0)
	{
		int term = (*p1d).degree;
		_Complex double y = (*p1d).c[term];
		while (term > 0)
		{
			y *= x;
			term -= 1;
			y += (*p1d).c[term];
		}
		return y;
	}
	return (_Complex double) 1e99;
}

// RESUME: test this!
void polynomial1D_static_DK_solve(polynomial1D_static *p1d, _Complex double *roots, double tol)
{
	double disc = (*p1d).c[(*p1d).degree] - 1;
	polynomial1D_static *aux = (polynomial1D_static *) calloc(1, sizeof(polynomial1D_static));
	(*aux).degree = -1;
	char reset = 0;
	if (-tol < disc && disc < tol) {}
	else
	{
		reset = 1;
		transcribe_polynomial1D_static(p1d, aux);
		double inv_cn = 1. / (*aux).c[(*aux).degree];
		for (int i = 0; i < (*aux).degree; i++)
		{
			(*aux).c[i] *= inv_cn;
		}
		(*aux).c[(*aux).degree] = 1.;
		polynomial1D_static *temp = p1d;
		p1d = aux;
		aux = temp;
		display_polynomial1D_static(p1d);
	}

	_Complex double incr[(*p1d).degree];
	double err = 1.0;
	double tolsq = tol * tol;
	int N_iter = 0;
	while (err > tolsq && N_iter < DK_ITER_MAX)
	{
		N_iter += 1;
		// Evaluate the increment
		for (int ri = 0; ri < (*p1d).degree; ri++)
		{
			_Complex double factor = 1.0;
			for (int rii = 0; rii < (*p1d).degree; rii++)
			{
				if (rii != ri)
				{
					factor *= roots[ri] - roots[rii];
				}
			}
			incr[ri] = - polynomial1D_static_evc(p1d, roots[ri]) / factor;
		}
		// Update the root estimate
		for (int ri = 0; ri < (*p1d).degree; ri++)
		{
			roots[ri] += incr[ri];
		}
		err = 0;
		for (int i = 0; i < (*p1d).degree; i++)
		{
			double re_incr = __real__ incr[i];
			double im_incr = __imag__ incr[i];
			err += re_incr * re_incr + im_incr * im_incr;
		}
	}
	if (reset)
	{
		polynomial1D_static *temp = aux;
		aux = p1d;
		p1d = temp;
		free_polynomial1D_static(aux);
	}
}


double polynomial1D_static_ev(polynomial1D_static *p1d, double x)
{
	if ((*p1d).degree >= 0)
	{
		int term = (*p1d).degree;
		double y = (*p1d).c[term];
		while (term > 0)
		{
			y *= x;
			term -= 1;
			y += (*p1d).c[term];
		}
		return y;
	}
	return 1e99;
}

void transcribe_polynomial1D_static(polynomial1D_static *A, polynomial1D_static *B)
{
	if ((*B).degree == -1) {}
	else if ((*B).degree > 0)
	{
		free((*B).c);
	}
	(*B).c = (double *) calloc((*A).degree + 1, sizeof(double));
	for (int i = 0; i <= (*A).degree; i++) (*B).c[i] = (*A).c[i];
	(*B).degree = (*A).degree;
}

void free_polynomial1D_static(polynomial1D_static *p1d)
{
	free((*p1d).c);
	(*p1d).degree = 0;
}

void polynomial1D_static_free(polynomial1D_static *A)
{
	free((*A).c);
	(*A).degree = 0;
}

void polynomial1D_static_deriv(polynomial1D_static *p1d, polynomial1D_static *dp1d)
{
	polynomial1D_static_init((*p1d).degree - 1, dp1d);
	int i = 0;
	for (int ip1 = 1; ip1 <= (*p1d).degree; ip1++)
	{
		(*dp1d).c[i] = ip1 * (*p1d).c[ip1];
		i = ip1;
	}
}

void polynomial1D_static_multiply(polynomial1D_static *A, polynomial1D_static *B, polynomial1D_static *C)
{
	polynomial1D_static_init((*A).degree + (*B).degree, C);
	for (int i = 0; i <= (*A).degree; i++) for (int ii = 0; ii <= (*B).degree; ii++)
	{
		(*C).c[i + ii] += (*A).c[i] * (*B).c[ii];
	}
}

void polynomial1D_static_3D_static_multiply(polynomial1D_static *A, polynomial3D_static *B, polynomial3D_static *C, int dim)
{
	int degreeC[3] = {(*B).degree[0], (*B).degree[1], (*B).degree[2]};
	degreeC[dim] += (*A).degree;
	polynomial3D_static_init(degreeC, C);
	for (int i = 0; i <= (*A).degree; i++) for (int j = 0; j <= (*B).degree[0]; j++) for (int jj = 0; jj <= (*B).degree[1]; jj++) for (int jjj = 0; jjj <= (*B).degree[2]; jjj++)
	{
		double incr_val = (*B).c[j][jj][jjj] * (*A).c[i];
		switch (dim)
		{
			case 0:
				(*C).c[i + j][jj][jjj] += incr_val;
				break;
			case 1:
				(*C).c[j][jj + i][jjj] += incr_val;
				break;
			case 2:
				(*C).c[j][jj][jjj + i] += incr_val; 
				break;
			break;
		}
	}
}

void Newton_polynomial1D_static(polynomial1D_static *f, double *guess, double tol)
{
	polynomial1D_static df;
	polynomial1D_static_deriv(f, &df);
	double err = 1;
	int iter_count = 0;
	while (err > tol && iter_count < NEWTON_ITER_MAX)
	{
		// f approx df * (x - r), r approx x - f/df
		double f_eval = polynomial1D_static_ev(f, *guess);
		double df_eval = polynomial1D_static_ev(&df, *guess);
		double incr = f_eval / df_eval;
		err = f_eval > 0? f_eval: -f_eval;
		err += incr > 0? incr: -incr;
		*guess -= incr;
		iter_count += 1;
	}
}

void transcribe_polynomial1D_static_3D_static(polynomial1D_static *A, polynomial3D_static *B, int dim)
{
	int degree[3] = {0, 0, 0};
	degree[dim] = (*A).degree;
	polynomial3D_static_init(degree, B);
	for (int i = 0; i <= (*A).degree; i++)
	{
		switch (dim)
		{
			case 0:
				(*B).c[i][0][0] = (*A).c[i];
				break;
			case 1:
				(*B).c[0][i][0] = (*A).c[i];
				break;
			case 2:
				(*B).c[0][0][i] = (*A).c[i];
				break;
			break;
		}
	}
}

double monomial_eval1D(double x, int degree)
{
	double result = 1.0;
	double incr = x;
	int pow = degree;
	while (pow > 0)
	{
		if (pow & 1) result *= incr;
		pow >>= 1;
		incr *= incr;
	}
	return result;
}

void display_polynomial1D_static(polynomial1D_static *p1d)
{
	int first_nonzero = 0;
	while ((*p1d).c[first_nonzero] == 0)
	{
		first_nonzero += 1;
	}
	if (first_nonzero == 0) printf("%g ", (*p1d).c[0]);
	else
	{
		printf("%gt^%d ", (*p1d).c[first_nonzero], first_nonzero);
	}
	if ((*p1d).degree > 0) 
	{
		for (int i = first_nonzero + 1; i <= (*p1d).degree; i++)
		{
			if ((*p1d).c[i] > 0) printf(" + %gt^%d", (*p1d).c[i], i);
			else if ((*p1d).c[i] < 0) printf(" - %gt^%d", -(*p1d).c[i], i);
		}
	}
	printf("\n");
}


void polynomial1D_static_fit(polynomial1D_static *A, double *sample_pts, double *sample_vals, int N_samples, int max_degree)
{
	// Requires linear algebra: sum_i (f(x_i; C) - y_i)^2 + -> sum_i (f(x_i; C) - y_i) x_i^{\alpha} = 0, x_i^{\alpha} \equiv \prod_j x_{ij}^{\alpha_j}
	//	i.e. sum_\beta C_\beta sum_i x_i^{\beta + \alpha} = sum_i y_i x_i^{\alpha}, \forall \alpha
	//		M C = b; M_{\alpha \beta} \equiv sum_i x_i^{\beta + \alpha}, b \equiv sum_i y_i x_i^{\alpha}
	//	LU C = b; U C = Linv b, 
	polynomial1D_static_init(max_degree, A);
	int N_terms = max_degree + 1;
	gsl_matrix *gslm = gsl_matrix_calloc(N_terms, N_terms);
	gsl_vector *b = gsl_vector_calloc(N_terms);
	for (int i = 0; i < N_samples; i++)
	{
		for (int ti = 0; ti < N_terms; ti++)
		{
			double pow_x_alpha = monomial_eval1D(sample_pts[i], ti);
			for (int tii = 0; tii <= ti; tii++)
			{
				double pow_x_beta = monomial_eval1D(sample_pts[i], tii);
				gsl_matrix_set(gslm, ti, tii, gsl_matrix_get(gslm, ti, tii) + pow_x_beta * pow_x_alpha);
			}
			gsl_vector_set(b, ti, gsl_vector_get(b, ti) + pow_x_alpha * sample_vals[i]);
		}
	}
	for (int ti = 0; ti < N_terms; ti++)
	{
		for (int tii = 0; tii < ti; tii++)
		{
			gsl_matrix_set(gslm, tii, ti, gsl_matrix_get(gslm, ti, tii));
		}
	}
	gsl_vector *x = gsl_vector_calloc(N_terms);
	gsl_permutation *pi = gsl_permutation_calloc(N_terms);
	int sign;
	gsl_linalg_LU_decomp(gslm, pi, &sign);
	gsl_linalg_LU_solve(gslm, pi, b, x);
	for (int i = 0; i < N_terms; i++)
	{
		(*A).c[i] = gsl_vector_get(x, i);
	}
	gsl_permutation_free(pi);
	gsl_vector_free(x);
	gsl_vector_free(b);
	gsl_matrix_free(gslm);
}

char measure_epsilon_nbhd_compliance(double (*inner)(double, void *), void *pars, gsl_spline *outer, double *x_samples, int N_samples, double epsilon, double tolerance)
{
	// Distance from points on outer curve to inner curve:
	//	solve the system of equations 
	//		(x - x0) + (y(x) - y0) * y'(x) = 0, 
	//		or (x - x0) + (y - y0(x0)) * y0'(x0) = 0 
	//		(the latter might be easier with gsl's built-in routines for splines)
	//	Taking the second: 
	//		f'(x0) = -1 - y0'(x0)^2 + (y - y0(x0)) * y0''(x0),
	// After finding the closest points, compute the distance, and check that it is within the tolerance of 
	//	epsilon
	double nmtol = 0.01 * tolerance; // An exact estimate of the nearest point isn't necessary
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	char sat = 1;
	double llimsq = epsilon - tolerance;
	double ulimsq = epsilon + tolerance;
	llimsq *= llimsq;
	ulimsq *= ulimsq;
	for (int i = 0; i < N_samples; i++)
	{
		double err = 1;
		double guess = x_samples[i];
		double x = x_samples[i];
		double y = inner(x, pars);
		double y0 = gsl_spline_eval(outer, guess, acc);
		while (err > nmtol)
		{
			double y0p = gsl_spline_eval_deriv(outer, guess, acc);
			double y0pp = gsl_spline_eval_deriv2(outer, guess, acc);
			double ymy0 = y - y0;
			double f_eval = (x - guess) + ymy0 * y0p;
			double fprime = -1 - y0p * y0p + ymy0 * y0pp;
			// f_eval = fprime * (x - r);
			err = f_eval / fprime;
			guess -= err;
			err = err > 0? err: -err;
			err += f_eval > 0? f_eval: -f_eval;
			y0 = gsl_spline_eval(outer, guess, acc);
		}
		double delx = x - guess;
		double dely = y - y0;
		double delsq = delx * delx + dely * dely;
		if (delsq > llimsq && delsq < ulimsq) {}
		else
		{
			sat = 0;
			break;
		}
	}
	gsl_interp_accel_free(acc);
	return sat;
}

// RESUME: update instances of approx_volume to be consistent with the new version.
// RESUME: check the original version of approx_volume to determine why it wasn't working as expected
// INCOMPLETE: apparent failure to detect overlaps? Slow convergence of integrals? Consider simplex methods to improve convergence.
void ball_overlap_spline_rev(double (*shape)(double, void *), void *shape_pars_, double ball_radius, double (*radial_dist)(double, void *), void *rd_pars_, double *x_samples, int N_samples, gsl_spline **ovl, double tolerance) // TEST
{
	void **shape_pars = (void **) shape_pars_;
	void **rd_pars = (void **) rd_pars_;

	double x_bnds[2];
	x_bnds[0] = x_samples[0];
	x_bnds[1] = x_samples[N_samples - 1];
	// try to minimize the number of samples: avoid recomputing.
	//  Initialize the interior function
	void **int1pars = (void **) calloc(2, sizeof(void *));
	int1pars[0] = (void *) shape;
	int1pars[1] = (void *) shape_pars;
	sphere sph;
	sph.r = ball_radius;
	sph.x[0] = 0;
	sph.x[1] = 0;
	sph.x[2] = 0;
	void *int2pars = (void *) &sph;
	void **intmpars = (void **) calloc(3, sizeof(void *));
	char (**domains)(double *, void *) = calloc(2, sizeof(char (**)(double *, void *)));
	domains[0] = interior_solid_revolution;
	domains[1] = interior_sphere;
	intmpars[0] = (void *) domains;
	void **domain_pars = (void **) calloc(2, sizeof(void *));
	domain_pars[0] = (void *) int1pars;
	domain_pars[1] = int2pars;
	intmpars[1] = (void *) domain_pars;
	int Ndomains = 2;
	intmpars[2] = (void *) &Ndomains;
	// Define the initial sample set
	int N_rsamples = N_samples;
	double *z_samples;
	double *ovl_samples;
	alloc_double(&z_samples, N_rsamples);
	alloc_double(&ovl_samples, N_rsamples);
	int i_min = 0, i_max = N_rsamples - 1;
	double max_ovl = 0;
	printf("Computing initial estimate with %d samples\n", N_rsamples);
	double delz = z_samples[1] - z_samples[0];
	for (int i = 0; i < N_rsamples; i++)
	{
		z_samples[i] = x_samples[i];
		double lbnds[3], ubnds[3];
		printf("Attempting to evaluate radial distance function\n");
		sph.x[0] = radial_dist(x_samples[i], rd_pars_);
		printf("\t(done)\n");
		sph.x[2] = x_samples[i];
		for (int ii = 0; ii < 3; ii++)
		{
			lbnds[ii] = sph.x[ii] - sph.r;
			ubnds[ii] = sph.x[ii] + sph.r;
		}
		lbnds[2] = lbnds[2] >= z_samples[0] + tolerance * delz? lbnds[2]: z_samples[0] + tolerance * delz;
		ubnds[2] = ubnds[2] <= x_bnds[1] - tolerance * delz? ubnds[2]: x_bnds[1] - tolerance * delz;
		printf("Computing overlap at position %g %g to solid of revolution of radius %g within box of dimensions %g %g %g %g %g %g\n", sph.x[0], sph.x[2], shape(sph.x[2], shape_pars_), lbnds[0], lbnds[1], lbnds[2], ubnds[0], ubnds[1], ubnds[2]);
		double int_pt[3] = {sph.x[0], 0, sph.x[2]};
		int_pt[0] = 0.5 * (shape(sph.x[2], shape_pars_) + sph.x[0] - sph.r);
		approx_double appr_ovl = approx_volume(interior_multiple_domain, (void *) intmpars, 0.5 * tolerance, lbnds, ubnds, int_pt);
		ovl_samples[i] = appr_ovl.val;
		printf("\tVolume = %g\n", ovl_samples[i]);
		max_ovl = max_ovl > ovl_samples[i]? max_ovl: ovl_samples[i];
	}
	// Define the spline
	(*ovl) = gsl_spline_alloc(gsl_interp_steffen, N_rsamples);
	gsl_spline_init(*ovl, z_samples, ovl_samples, N_rsamples);
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	// Refine samples:
	double err = 1;
	int iter_count = 0;
	while (err > tolerance && iter_count < MAX_BALL_OVERLAP_ITER)
	{
		iter_count += 1;
		printf("Performing iteration %d\n", iter_count);
		// Refine the sample set
		int nN_rsamples = N_rsamples << 1;
		double *nz_samples;
		alloc_double(&nz_samples, nN_rsamples); // freed
		int N_rsamplesm1 = N_rsamples - 1;
		linspace(nz_samples, z_samples[0], z_samples[N_rsamplesm1], nN_rsamples);
		double rdelz = nz_samples[1] - nz_samples[0];
		double *novl_samples;
		alloc_double(&novl_samples, nN_rsamples);
		for (int i = 0; i < N_rsamples; i++)
		{
			int i_ = 2 * i;
			int i_p = i_ + 1;
			double lbnds[3], ubnds[3];
			novl_samples[i_] = ovl_samples[i];
			sph.x[0] = radial_dist(nz_samples[i_p], rd_pars_);
			sph.x[2] = nz_samples[i_p];
			// RESUME: Change this!
			for (int ii = 0; ii < 3; ii++)
			{
				lbnds[ii] = sph.x[ii] - sph.r;
				ubnds[ii] = sph.x[ii] + sph.r;
			}
			printf("Computing overlap vol at position %g %g\n", sph.x[0], sph.x[2]);
			lbnds[2] = lbnds[2] >= z_samples[0] + tolerance * rdelz? lbnds[2]: z_samples[0] + tolerance * rdelz;
			ubnds[2] = ubnds[2] <= z_samples[N_rsamplesm1] - tolerance * rdelz? ubnds[2]: z_samples[N_rsamplesm1] - tolerance * rdelz;
			double int_pt[3] = {0.5 * (shape(sph.x[2], shape_pars_) + sph.x[0] - sph.r), 0, sph.x[2]};
			approx_double approx_ovl = approx_volume(interior_multiple_domain, (void *) intmpars, 0.5 * tolerance, lbnds, ubnds, int_pt);
			novl_samples[i_p] = approx_ovl.val;
			max_ovl = max_ovl > novl_samples[i_p]? max_ovl: novl_samples[i_p];
		}
		// Compute the average/maximum error 
		err = 0;
		int i_p = 1;
		for (int i = 0; i < N_samples; i++)
		{
			double diff = gsl_spline_eval(*ovl, nz_samples[i_p], acc);
			diff -= novl_samples[i_p];
			diff = diff >= 0? diff: -diff;
			err = err > diff? err: diff;
			i_p += 2;
		}
		if (max_ovl > 0) err /= max_ovl;
		// Replace the original samples with the refined set
		printf("Refining sample set\n");
		free(z_samples);
		free(ovl_samples);
		z_samples = nz_samples;
		ovl_samples = novl_samples;
		N_rsamples = nN_rsamples;
		// Update the spline
		printf("Updating the spline\n");
		gsl_spline_free(*ovl);
		*ovl = gsl_spline_alloc(gsl_interp_steffen, N_rsamples);
		gsl_spline_init(*ovl, z_samples, ovl_samples, N_rsamples);
		printf("(iteration complete)\n");
	}
	free(int1pars);
	free(intmpars);
	free(domains);
	free(domain_pars);
	free(ovl_samples);
	free(z_samples);
}


void epsilon_nbhd_spline_rev(double (*inner)(double, void *), void *pars, double *x_samples_, int N_samples, double epsilon, gsl_spline **outer, double tolerance) // TEST
{
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	double *x_samples = (double *) calloc(N_samples, sizeof(double)); // freed in first loop if initial spline is noncompliant, or post loop otherwise.
	for (int i = 0; i < N_samples; i++) 
	{
		x_samples[i] = x_samples_[i];
	}
	double abs_epsilon = fabs(epsilon);
	char noncompliant = 1;
	double x_bnds[2];
	x_bnds[0] = x_samples[0];
	x_bnds[1] = x_samples[N_samples - 1];
	while (noncompliant)
	{
		double *y_outer = (double *) calloc(N_samples, sizeof(double)); // freed, guaranteed at the end of each loop before end case is tested
		double *x_outer = (double *) calloc(N_samples, sizeof(double)); // freed (same as above)
		// Evaluate y_samples for the outer epsilon nbhd for each x sample
		double xm = x_samples[0];
		double xp = 0.5 * (x_samples[0] + x_samples[1]);
		double y_upper = inner(xp, pars);
		double y_lower = inner(xm, pars);
		for (int i = 0; i < N_samples; i++)
		{
			// Estimate the derivative of the inner curve near the i-th sample
			double h = xp - xm;
			double dy_inner = (y_upper - y_lower) / h;
			double y_inner = inner(x_samples[i], pars);
			double incr = epsilon / sqrt(1 + dy_inner * dy_inner);
			y_outer[i] = y_inner + incr;
			x_outer[i] = x_samples[i] - dy_inner * incr;
			y_lower = y_upper;
			xm = xp;
			if (i < N_samples - 2) xp = 0.5 * (x_samples[i + 1] + x_samples[i + 2]);
			else xp = x_samples[N_samples - 1];
			y_upper = inner(xp, pars);
		}
		double delx_a = x_outer[1] - x_outer[0];
		double delx_b = x_outer[N_samples - 1] - x_outer[N_samples - 2];
		x_bnds[0] = x_outer[0] > x_samples[0]? x_outer[0]: x_samples[0];
		x_bnds[1] = x_outer[N_samples - 1] < x_samples[N_samples - 1]? x_outer[N_samples - 1]: x_samples[N_samples - 1];
		x_bnds[0] += tolerance * delx_a;
		x_bnds[1] -= tolerance * delx_b;
		(*outer) = gsl_spline_alloc(gsl_interp_steffen, N_samples);
		gsl_spline_init(*outer, x_outer, y_outer, N_samples);
		int Nn_samples = N_samples << 1;
		double *x_nsamples = (double *) calloc(Nn_samples, sizeof(double)); // freed: both in the end case, and otherwise in the next iteration if continued
		linspace(x_nsamples, x_bnds[0], x_bnds[1], Nn_samples);
		printf("Testing consistency of approximation\n");
		char sat = measure_epsilon_nbhd_compliance(inner, pars, *outer, x_nsamples, Nn_samples, abs_epsilon, tolerance);
		printf("(done)\n");
		free(y_outer);
		free(x_outer);
		if (sat)
		{
			// end case
			noncompliant = 0;
			free(x_nsamples); 
			break;
		}
		else
		{
			// continue case
			gsl_spline_free(*outer);
			free(x_samples); // frees original allocation, or x_nsamples allocated in previous iteration
			x_samples = x_nsamples;
			N_samples = Nn_samples;
			linspace(x_samples, x_bnds[0], x_bnds[1], N_samples);
		}
	}
	gsl_interp_accel_free(acc);
	free(x_samples);
}

void interior_multiple_domain_unpack_pars(void *pars_, int *N_domains, char (***domains)(double *, void *), void ***dpars_)
{
	void **pars = (void **) pars_;
	(*N_domains) = *((int *) pars[2]);
	*domains = (char (**)(double *, void *)) pars[0];
	*dpars_ = (void **) pars[1];
}

char interior_multiple_domain(double *x, void *pars_)
{
	int N_domains;
	char (**domains)(double *, void *);
	void **dpars_;
	interior_multiple_domain_unpack_pars(pars_, &N_domains, &domains, &dpars_);
	/*void **pars = (void **) pars_;
	int N_domains = *((int *) pars[2]);
	char (**domains)(double *, void *);
	domains = (char (**)(double *, void *)) pars[0];*/
	for (int i = 0; i < N_domains; i++)
	{
		if (!domains[i](x, dpars_[i])) return 0;
	}
	return 1;
}

char interior_affine_xform_domain(double *x, void *pars_)
{
	void **pars = (void **) pars_;
	char (*interior)(double *, void *) = (char (*)(double *, void *)) pars[0];
	void *ipars = (void *) pars[1];
	gsl_matrix *G = (gsl_matrix *) pars[2];
	int dim = *((int *) pars[3]);
	double Gx[dim];
	for (int i = 0; i < dim; i++) 
	{
		Gx[i] = gsl_matrix_get(G, 3, i);
		for (int ii = 0; ii < dim; ii++)
		{
			Gx[i] += gsl_matrix_get(G, i, ii) * x[ii];
		}
	}
	return interior(Gx, ipars);
}

double d1polynomial1D_static_ev(double x, void *pars_)
{
	void **pars = (void **) pars_;
	polynomial1D_static *p1d = (polynomial1D_static *) pars[0];
	return polynomial1D_static_ev(p1d, x);
}

double d1spline_ev(double x, void *pars_)
{
	void **pars = (void **) pars_;
	gsl_spline *func = (gsl_spline *) pars[0];
	double xmin = (*(*func).interp).xmin;
	double xmax = (*(*func).interp).xmax;
	if (x < xmax && x > xmin) {}
	else
	{
		printf("Attempting to evaluate spline at %g, beyond its domain (%g %g)\n", x, xmin, xmax);
		exit(EXIT_FAILURE);
	}
	gsl_interp_accel *acc = (gsl_interp_accel *) pars[1];
	return gsl_spline_eval(func, x, acc);
}

double d1func_offset_ev(double x, void *pars_)
{
	void **pars = (void **) pars_;
	double (*func)(double, void *) = (double (*)(double, void*)) pars[0];
	void *func_pars = (void *) pars[1];
	double h = *((double *) pars[2]);
	return func(x, func_pars) + h;
}

char interior_solid_revolution(double *x, void *pars_)
{
	void **pars = (void **) pars_;
	double rsq = x[0] * x[0] + x[1] * x[1];
	double (*r_of_z)(double, void *) = (double (*)(double, void *)) pars[0];
	void *roz_pars = (void *) pars[1];
	double rthr = r_of_z(x[2], roz_pars);
	return rsq < rthr * rthr;
}

char interior_solid_revolution_spline(double *x, void *pars_)
{
	void **pars = (void **) pars_;
	double rsq = x[0] * x[0] + x[1] * x[1];
	gsl_spline *func = (gsl_spline *) pars[0];
	gsl_interp_accel *acc = (gsl_interp_accel *) pars[1];
	double rofz = gsl_spline_eval(func, x[2], acc);
	return rsq < (rofz * rofz);
}

char interior_solid_revolution_polynomial1D_static(double *x, void *pars_)
{
	void **pars = (void **) pars_;
	double rsq = x[0] * x[0] + x[1] * x[1];
	polynomial1D_static *p1d = (polynomial1D_static *) pars[0];
	double rofz = polynomial1D_static_ev(p1d, x[2]);
	return rsq < (rofz * rofz);
}

void vector_spline_alloc(int dim, vector_spline *vspl)
{
	(*vspl).spline = (gsl_spline **) calloc(dim, sizeof(gsl_spline *));
	(*vspl).dim = dim;
}

void vector_spline_init(vector_spline *vspl, double *x_samples, double **v_samples, int N_samples, gsl_interp_type *gslt)
{
	for (int i = 0; i < (*vspl).dim; i++)
	{
		(*vspl).spline[i] = gsl_spline_alloc(gslt, N_samples);
		gsl_spline_init((*vspl).spline[i], x_samples, v_samples[i], N_samples);
	}
}

void vector_spline_free(vector_spline *vspl)
{
	for (int i = 0; i < (*vspl).dim; i++) gsl_spline_free((*vspl).spline[i]);
	free((*vspl).spline);
}

char interior_sphere(double *x, void *pars_)
{
	sphere *sph = (sphere *) pars_;
	double delx, dely, delz;
	delx = x[0] - (*sph).x[0];
	dely = x[1] - (*sph).x[1];
	delz = x[2] - (*sph).x[2];
	return delx * delx + dely * dely + delz * delz < (*sph).r * (*sph).r;
}

void compare_precision(gsl_spline *s1, gsl_spline *s2, gsl_spline **finer, gsl_spline **coarser)
{
	if ((*s1).size >= (*s2).size)
	{
		*finer = s1;
		*coarser = s2;
	}
	else
	{
		*finer = s2;
		*coarser = s1;
	}
}

// Basic arithmetic operations
double op_add(double a, double b)
{
	return a + b;
}

double op_subtract(double a, double b)
{
	return a - b;
}

double op_divide(double a, double b)
{
	return a / b;
}

double op_multiply(double a, double b)
{
	return a * b;
}

// SPLINE OPERATIONS RESUME: check domains beforehand
void shared_domain_tally(gsl_spline *s1, gsl_spline *s2, double *xmin, double *xmax, int *t1, int *t2)
{
	*xmin = (*(*s1).interp).xmin >= (*(*s2).interp).xmin? (*(*s1).interp).xmin: (*(*s2).interp).xmin;
	*xmax = (*(*s1).interp).xmax <= (*(*s2).interp).xmax? (*(*s2).interp).xmax: (*(*s2).interp).xmax;
	(*t1) = 0;
	(*t2) = 0;
	for (int i = 0; i < (*s1).size; i++) (*t1) += (*s1).x[i] <= *xmax && (*s1).x[i] >= *xmin;
	for (int i = 0; i < (*s2).size; i++) (*t2) += (*s2).x[i] <= *xmax && (*s2).x[i] >= *xmin;
}

void apply_op_spline_coarse_fine(gsl_spline *s1, gsl_spline *s2, gsl_spline **s1ms2, double (*op)(double, double), double xmin, double xmax, int Nsamples)
{
	double ysamples[Nsamples];
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	(*s1ms2) = gsl_spline_alloc((*(*s2).interp).type, Nsamples);
	if (Nsamples == (*s2).size)
	{
		for (int i = 0; i < (*s2).size; i++)
		{
			ysamples[i] = op(gsl_spline_eval(s1, (*s2).x[i], acc), (*s2).y[i]);
		}
		gsl_spline_init(*s1ms2, (*s2).x, ysamples, (*s2).size);
	}
	else
	{
		double xsamples[Nsamples];
		linspace(xsamples, xmin, xmax, Nsamples);
		for (int i = 0; i < Nsamples; i++)
		{
			ysamples[i] = op(gsl_spline_eval(s1, xsamples[i], acc), gsl_spline_eval(s2, xsamples[i], acc));
		}
		gsl_spline_init(*s1ms2, xsamples, ysamples, Nsamples);
	}
	gsl_interp_accel_free(acc);
}

void apply_op_spline_fine_coarse(gsl_spline *s2, gsl_spline *s1, gsl_spline **s2ms1, double (*op)(double, double), double xmin, double xmax, int Nsamples) 
{
	double ysamples[Nsamples];
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	(*s2ms1) = gsl_spline_alloc((*(*s2).interp).type, Nsamples);
	if (Nsamples == (*s2).size)
	{
		for (int i = 0; i < (*s2).size; i++)
		{
			ysamples[i] = op((*s2).y[i],  gsl_spline_eval(s1, (*s2).x[i], acc)); 
		}
		gsl_spline_init(*s2ms1, (*s2).x, ysamples, (*s2).size);
	}
	else
	{
		double xsamples[Nsamples];
		linspace(xsamples, xmin, xmax, Nsamples);
		for (int i = 0; i < Nsamples; i++)
		{
			ysamples[i] = op(gsl_spline_eval(s2, xsamples[i], acc), gsl_spline_eval(s1, xsamples[i], acc));
		}
		gsl_spline_init(*s2ms1, xsamples, ysamples, Nsamples);
	}
	gsl_interp_accel_free(acc);
}

void product_spline(gsl_spline *s1, gsl_spline *s2, gsl_spline **p12)
{
	double xmin, xmax;
	int t1, t2, Nsamples;
	shared_domain_tally(s1, s2, &xmin, &xmax, &t1, &t2);
	Nsamples = t1 >= t2? t1: t2;
	if ((*s1).size >= (*s2).size) apply_op_spline_fine_coarse(s1, s2, p12, op_multiply, xmin, xmax, Nsamples);
	else apply_op_spline_coarse_fine(s1, s2, p12, op_multiply, xmin, xmax, Nsamples);
}

void quotient_spline(gsl_spline *s1, gsl_spline *s2, gsl_spline **q12)
{
	double xmin, xmax;
	int t1, t2, Nsamples;
	shared_domain_tally(s1, s2, &xmin, &xmax, &t1, &t2);
	Nsamples = t1 >= t2? t1: t2;
	if ((*s1).size >= (*s2).size) apply_op_spline_fine_coarse(s1, s2, q12, op_divide, xmin, xmax, Nsamples);
	else apply_op_spline_coarse_fine(s1, s2, q12, op_divide, xmin, xmax, Nsamples);

}

void inverse_spline(gsl_spline *spl, gsl_spline **inv_spl)
{
	double y_samples[(*spl).size];
	for (int i = 0; i < (*spl).size; i++) 
	{
		y_samples[i] = 1. / (*spl).y[i];
	}
	*inv_spl = gsl_spline_alloc((*(*spl).interp).type, (*spl).size);
	gsl_spline_init(*inv_spl, (*spl).x, y_samples, (*spl).size);
}

void rescale_spline(gsl_spline **spl, double coeff) // TEST THIS!
{
	for (int i = 0; i < (**spl).size; i++)
	{
		(**spl).y[i] *= coeff;
	}
	gsl_spline_init(*spl, (**spl).x, (**spl).y, (**spl).size);
}

void difference_spline(gsl_spline *s1, gsl_spline *s2, gsl_spline **s1ms2)
{
	double xmin, xmax;
	int t1, t2, Nsamples;
	shared_domain_tally(s1, s2, &xmin, &xmax, &t1, &t2);
	Nsamples = t1 >= t2? t1: t2;

	if ((*s1).size >= (*s2).size) apply_op_spline_fine_coarse(s1, s2, s1ms2, op_subtract, xmin, xmax, Nsamples);
	else apply_op_spline_coarse_fine(s1, s2, s1ms2, op_subtract, xmin, xmax, Nsamples);
}

void divide_spline(gsl_spline **numerator, gsl_spline *denominator) // TEST! (should follow from rescale)
{

	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	if ((*(**numerator).interp).xmin >= (*(*denominator).interp).xmin && (*(**numerator).interp).xmax <= (*(*denominator).interp).xmax)
	{
		for (int i = 0; i < (**numerator).size; i++)
		{
			(**numerator).y[i] /= gsl_spline_eval(denominator, (**numerator).x[i], acc);
		}
		gsl_spline_init(*numerator, (**numerator).x, (**numerator).y, (**numerator).size);
	}
	else
	{
		// Refine samples of the numerator
		double xmin = (*(**numerator).interp).xmin >= (*(*denominator).interp).xmin? (*(**numerator).interp).xmin: (*(*denominator).interp).xmin;
		double xmax = (*(**numerator).interp).xmax <= (*(*denominator).interp).xmax? (*(**numerator).interp).xmax: (*(*denominator).interp).xmax;
		double ysamples[(**numerator).size];
		double xsamples[(**numerator).size];
		linspace(xsamples, xmin, xmax, (**numerator).size);
		for (int i = 0; i < (**numerator).size; i++)
		{
			ysamples[i] = gsl_spline_eval(*numerator, xsamples[i], acc) / gsl_spline_eval(denominator, xsamples[i], acc);
		}
		gsl_spline_init(*numerator, xsamples, ysamples, (**numerator).size);
	}
	gsl_interp_accel_free(acc);
}

void derivative_spline(gsl_spline *f, gsl_spline **dfdx)
{
	*dfdx = gsl_spline_alloc((*(*f).interp).type, (*f).size);
	double ysamples[(*f).size];
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	for (int i = 0; i < (*f).size; i++)
	{
		ysamples[i] = gsl_spline_eval_deriv(f, (*f).x[i], acc);
	}
	gsl_spline_init(*dfdx, (*f).x, ysamples, (*f).size);
	gsl_interp_accel_free(acc);
}

void norm_vector_spline(gsl_spline *curve, gsl_spline **norm_x, gsl_spline **norm_y)
{
	(*norm_x) = gsl_spline_alloc((*(*curve).interp).type, (*curve).size);
	(*norm_y) = gsl_spline_alloc((*(*curve).interp).type, (*curve).size);
	double nx_samples[(*curve).size];
	double ny_samples[(*curve).size];
	gsl_interp_accel *acc = gsl_interp_accel_alloc();
	for (int i = 0; i < (*curve).size; i++)
	{
		double deriv = gsl_spline_eval_deriv(curve, (*curve).x[i], acc);
		ny_samples[i] = 1. / sqrt(1 + deriv * deriv);
		nx_samples[i] = -deriv * ny_samples[i];
	}
	gsl_spline_init(*norm_x, (*curve).x, nx_samples, (*curve).size);
	gsl_spline_init(*norm_y, (*curve).x, ny_samples, (*curve).size);
	gsl_interp_accel_free(acc);
}

void add_spline(gsl_spline *s1, gsl_spline *s2, gsl_spline **s1ps2)
{
	double xmin, xmax;
	int t1, t2, Nsamples;
	shared_domain_tally(s1, s2, &xmin, &xmax, &t1, &t2);
	Nsamples = t1 >= t2? t1: t2;
	if ((*s1).size >= (*s2).size) apply_op_spline_fine_coarse(s1, s2, s1ps2, op_add, xmin, xmax, Nsamples);
	else apply_op_spline_coarse_fine(s1, s2, s1ps2, op_add, xmin, xmax, Nsamples);
}

// RESUME
void closest_point(double (*surface)(double *, void *), void *surf_pars, void (*grad_surface)(double *, void *, double *), void *gsurf_pars, 
			double lvl_const, double *x, double *cp)
{
	// min (x - cp)^2, surface(cp, pars) = lvl_const
	// (x - cp) = lambda * grad_surface(cp, pars)
	// 	f(cp) = lambda * grad(cp, pars) + cp - x, surface(cp, pars) - lvl_const
	// Note: the NM minimization method probably fails to converge at inflection points
	if (grad_surface != NULL)
	{
		// Apply root-finding method to Lagrange-multiplier problem
		//	f(x) = 0, grad f(x) = lambda (x - cp), 
		
	}
	else
	{
		// Use the Nelder-Mead and/or interpolation-based methods
	}
}

double triangle_prism_vol_simple_intermediate(double h0, double h1, double h2)
{
	return (h0 + h1 + h2);
}

double triangle_prism_vol_intermediate(double h0, double h1, double h2, double *x0, double *x1, double *x2)
{
	double x01[2] = {x1[0] - x0[0], x1[1] - x0[1]};
	double x12[2] = {x2[0] - x1[0], x2[1] - x1[1]};
	double x20[2] = {x0[0] - x2[0], x0[1] - x2[1]};
	double c012 = x01[0] * x12[1] - x01[1] * x12[0];
	double c120 = x12[0] * x20[1] - x12[1] * x20[0];
	double c201 = x20[0] * x01[1] - x20[1] * x01[0];
	c012 = c012 > 0? c012: -c012;
	c120 = c120 > 0? c120: -c120;
	c201 = c201 > 0? c201: -c201;
	return (h0 * c201 + h1 * c012 + h2 * c120);
}

double triangle_prism_vol(double h0, double h1, double h2, double *x0, double *x1, double *x2)
{
	return triangle_prism_vol_intermediate(h0, h1, h2, x0, x1, x2) / 6;
}

double quad_prism_vol(double h0, double h1, double h2, double h3, double *x0, double *x1, double *x2, double *x3)
{
	return (triangle_prism_vol_intermediate(h0, h1, h2, x0, x1, x2) + triangle_prism_vol_intermediate(h3, h1, h2, x3, x1, x2)) / 6;
}

double pentagon_prism_vol(double h0, double h1, double h2, double bh0, double bh1, double *x0, double *x1, double *x2, double *bx0, double *bx1)
{
	return (triangle_prism_vol_intermediate(h0, h1, bh0, x0, x1, bx0) + triangle_prism_vol_intermediate(bh0, h1, h2, bx0, x1, x2) + triangle_prism_vol_intermediate(bh0, bh1, h2, bx0, bx1, x2)) / 6;
}

double hexagon_prism_vol(double h0, double bh0, double bh1, double h1, double bh2, double bh3, double *x0, double *bx0, double *bx1, double *x1, double *bx2, double *bx3)
{
	return (triangle_prism_vol_intermediate(bh0, bh1, bh2, bx0, bx1, bx2) + triangle_prism_vol_intermediate(bh0, bh2, bh3, bx0, bx2, bx3) + triangle_prism_vol_intermediate(bh0, h0, bh3, bx0, x0, bx3) + triangle_prism_vol_intermediate(bh1, h1, bh2, bx1, x1, bx2)) / 6;
}

#include "aux/numerical_fubini_case_by_thickness.h"

// Conventions for traversing plaquette vertices, as presented in plaq_xc and plaq_yc
// double plaq_xc[9] = {x0, x1, x2, x2, x2, x1, x0, x0, x1};
// double plaq_yc[9] = {y0, y0, y0, y1, y2, y2, y2, y1, y1};
int plaq_trav[4][4] = {{8, 3, 4, 5}, {8, 5, 6, 7}, {8, 7, 0, 1}, {8, 2, 3, 4}};

// NOTE: there is a more efficient way to approach this.  Instead of checking each lattice site for inclusion in the interior of the domain,
//		associate boundary-adjacent interior sites to upper and lower curves (y_max/y_min) by sweeping left to right, starting with all
//		lattice sites near the x bounds.  More specifically; after adding lattice sites between y_min(x0 + dx) and y_max(x0 + dx) (the 
//		first interior layer), select the interior site at x0 + dx = x1 (or xN) with minimal y coordinate y{a,N}, evaluate y_min(x{N+1}), determine
//		the site with minimal y coordinate at x{N+1}, and add all lattice sites with y coordinate between y{a,N} and y{a,N+1} to the 'interior 
//		boundary'.  Other lattice sites (between max(y{a,N}, y{a,N+1}) and min(y{b,N}, y{b,N+1}, where y{b,N} is the maximum lattice site) exclusive) 
//		qualify as interior sites, and 'boundary tiles' can be classified automatically by shape according to the number of increments in 'y' between
//		adjacent y minima (something like 1 pentagon, quadrangles/trapezoids, then a triangle at the very end when there is a change in altitude, 
//		or quadrangles when altitude remains unchanged, and a quadrangle or very singular pentagon at the very beginning and end (of the sweep over x.))
// RESUME: implement the more efficient version independently, and compare results to ensure correctness.
// RESUME: Check this!
double integrate_2D_fubini(double (*f)(double, double, void *), void *fpars, double *xbnds, double (*y_min)(double, void *), void *y_min_pars, double (*y_max)(double, void *), void *y_max_pars, double tol)
{
	// Scaling: expect quadratic error with linear function approximation: f(x) = f0(x) + 0.5 dx.H.dx + O(dx^3)
	//						int_dA f(x) = I0dA + int_{-hdx}^{hdx} int_{-hdy}^{hdy}  0.5 * dx.H.dx + O(dx^5)
	//	int_{-a}^a int_{-b}^b (A x^2 + B y^2 + C xy) = int_{-a}^a (2b A x^2 + B 2 b^3 / 3 ) dx
	//							= 4b A a^3 / 3 + 4 B b^3 a / 3 -> (dx^2 A + dy^2 B) dxdy / 12
	//	=> err goes (roughly) as dx^4 * A / (dx^2) = A * dx^2 ~ tol => dx ~ sqrt(tol / A);
	//	
	//	sum_<i, j> (f_i - f0_j - C_j.(x_i - x_j))^2 -> sum_<i, j> (f_i - f0_j - C_j.(x_i - x_j)) = 0, sum_<i, j> (f_i - f0_j - C_j.(x_i - x_j))(x_i - x_j) = 0
	//	
	//		x_i - x_j = dx * sqrt(2) / 2 * exp(I * pi (2i +1) / 4), rect lattice, 2 / sqrt(3) * dx / 2 * exp(2 * pi * I i / 3) or exp(2 * pi * I * (i + 0.5) / 3), triangular lattice, 
	//	local approx based on triangular tiling of surface: if (x - x_j)_x + (x - x_j)_y < 0.5 ell, f0 + ((f1 - f0) * (x - x_j)_x + (f2 - f0) * (x - x_j)_y) / ell, else if ... > 0.5, f3 + ((f3 - f1) * (x - x'_j)_x + (f3 - f2) * (x - x'_j)_y) / ell, etc.
	//		At boundaries: 
	//			- Interior points whose neighbors on square lattice are not contained in the region of integration
	//			- Determine neighboring boundary-crossing sites, solve for points of intersection with the boundary
	//			(using e.g. bisection method.)
	//			- Contract boundary crossing squares to boundary adjacent quadrangles or pentagons or triangles.
	//			- Evaluate the approximate integral of each domain:
	//				- Triangles: 
	//				|ell_{01} x ell_{02}| * f_0 / 3! + |ell_{10} x ell_{12}| * f_1 / 3! + |ell_{20} x ell_{21}| * f_2 / 3!
	//				- Quadrangles: (orientation dependent)
	//					Decompose into two triangles (possibly at random)
	//						- Choose a corner, then split along diagonal. 
	//				- Pentagons: 
	//					Decompose into three triangles (possibly at random)
	//						- Connect both corners adjacent to the boundary to the opposite interior corner
	//				- Hexagons (rare):	
	//					Decompose into four triangles (possibly at random)
	//						- Form a quadrangle from the four boundary adjacent points, and two triangles 
	//						associated with the two interior points.
	double est = 0;
	double dx = sqrt(tol);
	double len_x = xbnds[1] - xbnds[0];
	int Nxsamples = (int) (len_x / dx) + 1;
	dx = len_x / (Nxsamples - 1);
	double invdx = 1. / dx;
	double x0 = xbnds[0];
	double x1 = x0 + dx;
	double x2 = x1 + dx;
	double ymin0 = y_min(x0, y_min_pars);
	double ymax0 = y_max(x0, y_max_pars);
	double ymin1 = y_min(x1, y_min_pars);
	double ymax1 = y_max(x1, y_max_pars);
	double ymin2 = y_min(x2, y_min_pars);
	double ymax2 = y_max(x2, y_max_pars);
	aarray_int bdry_plaquettes;
	aarray_int_init(&bdry_plaquettes, INIT_A_MEM);
	array_int bdry_plaq_stat;
	array_int_init(&bdry_plaq_stat, INIT_A_MEM);
	int im1 = 1;
	int im2 = 0;
	double bdry_est = 0;
	for (int i = 2; i < Nxsamples; i++)
	{
		double est_y = 0;
		double ybnds[2] = {ymin1, ymax1};
		//...
		int iy = (int) (ybnds[0] * invdx);
		int Nysamples = (int) ((ymax1 - ymin1) * invdx) + 1;
		double y0 = iy * dx;
		double y1 = y0 + dx;
		double y2 = y1 + dx;
		vector_real bdry_pt = vector_real_init(3);
		bdry_pt.e[0] = x0;
		bdry_pt.e[1] = y0;
		int iim1 = 1;
		int iim2 = 0;
		for (int ii = 2; ii < Nysamples; ii++)
		{
			// Check that each adjacent plaquette is contained in the domain
			char plaq_status[8];
			double plaq_xc[9] = {x0, x1, x2, x2, x2, x1, x0, x0, x1};
			double plaq_yc[9] = {y0, y0, y0, y1, y2, y2, y2, y1, y1};
			double plaq_ycll[8] = {ymin0, ymin1, ymin2, ymin2, ymin2, ymin1, ymin0, ymin0};
			double plaq_ycul[8] = {ymax0, ymax1, ymax2, ymax2, ymax2, ymax1, ymax0, ymax0};
			char interior = 1;
			for (int iii = 0; iii < 8; iii++)
			{
				plaq_status[iii] = plaq_xc[iii] > xbnds[0] && plaq_xc[iii] < xbnds[1] && plaq_yc[iii] > plaq_ycll[iii] && plaq_yc[iii] < plaq_ycul[iii];
				interior = interior && plaq_status[iii];
			}
			if (interior)
			{
				est_y += f(x1, y1, fpars);
			}
			else
			{
				// Evaluate boundary terms according to the plaquette status
				// Check if each sub-plaquette has been included already
				int splaqs[4][2] = {{im1, iim1}, {im2, iim1}, {im2, iim2}, {im1, iim2}};
				for (int iii = 0; iii < 4; iii++)
				{
					if (aarray_int_contains(&bdry_plaquettes, splaqs[iii], 2)) {}
					else
					{
						array_int nplaquette;
						array_int_init(&nplaquette, 2);
						nplaquette.e[0] = splaqs[iii][0];
						nplaquette.e[1] = splaqs[iii][1];
						add2aarray_int(&bdry_plaquettes, nplaquette);
						int nplaq_stat = 0;
						// Check which corners of bdry plaquette are contained in the domain (besides the center point)
						nplaq_stat |= array_bit_int_masks[0];
						for (int iv = 1; iv < 4; iv++)
						{
							double xiv = plaq_xc[plaq_trav[iii][iv]];
							double yiv = plaq_yc[plaq_trav[iii][iv]];
							double yll = plaq_ycll[plaq_trav[iii][iv]];
							double yul = plaq_ycul[plaq_trav[iii][iv]];
							if (xiv > xbnds[0] && xiv < xbnds[1] && yiv > yll && yiv < yul)
							{
								nplaq_stat |= array_bit_int_masks[iv];
							}
						}
						// Try to evaluate the plaquette directly
						add2array_int(&bdry_plaq_stat, nplaq_stat);
						// Cases: 1000: triangle, 1100, 1001: quadrangle, 1010: hexagon, 1011, 1101, 1110: pentagon
						double xa[3] = {x1, y1, f(x1, y1, fpars)};
						if (nplaq_stat == 8)
						{
							// Triangle: 1000
							double xb[3] = {plaq_xc[plaq_trav[iii][1]], plaq_yc[plaq_trav[iii][1]], 0};
							double xc[3] = {plaq_xc[plaq_trav[iii][3]], plaq_yc[plaq_trav[iii][3]], 0};
							switch (iii)
							{
								void *inv_func_pars[4];
								case 0: // upper right
									xc[1] = plaq_ycul[plaq_trav[iii][3]];
									inv_func_pars[0] = (void *) y_max;
									inv_func_pars[1] = (void *) y_max_pars;
									inv_func_pars[2] = (void *) &xa[0];
									inv_func_pars[3] = (void *) &xb[0];
									xb[0] = invert_bisection_method(xb[1], (void *) inv_func_pars); 
									break;
								case 1: // upper left
									xb[1] = plaq_ycul[plaq_trav[iii][1]];
									inv_func_pars[0] = (void *) y_max;
									inv_func_pars[1] = (void *) y_max_pars;
									inv_func_pars[2] = (void *) &xc[0];
									inv_func_pars[3] = (void *) &xa[0];
									xc[0] = invert_bisection_method(xc[1], (void *) inv_func_pars);
									break;
								case 2: // lower left
									xc[1] = plaq_ycll[plaq_trav[iii][3]];
									inv_func_pars[0] = (void *) y_min;
									inv_func_pars[1] = (void *) y_min_pars;
									inv_func_pars[2] = (void *) &xb[0];
									inv_func_pars[3] = (void *) &xa[0];
									xb[0] = invert_bisection_method(xb[1], (void *) inv_func_pars);
									break;
								case 3: // lower right
									xb[1] = plaq_ycll[plaq_trav[iii][1]];
									inv_func_pars[0] = (void *) y_min;
									inv_func_pars[1] = (void *) y_min_pars;
									inv_func_pars[2] = (void *) &xa[0];
									inv_func_pars[3] = (void *) &xc[0];
									xc[0] = invert_bisection_method(xc[1], (void *) inv_func_pars);
									break;
							}
							xb[2] = f(xb[0], xb[1], fpars);
							xc[2] = f(xc[0], xc[1], fpars);
							bdry_est += triangle_prism_vol(xa[2], xb[2], xc[2], xa, xb, xc);
						}
						else if (nplaq_stat == 12 || nplaq_stat == 9)
						{
							// Quadrangle: 1100, 1001
							double xb[3] = {plaq_xc[plaq_trav[iii][1]], plaq_yc[plaq_trav[iii][1]], 0};
							double xc[3] = {plaq_xc[plaq_trav[iii][2]], plaq_yc[plaq_trav[iii][2]], 0};
							double xd[3] = {plaq_xc[plaq_trav[iii][3]], plaq_yc[plaq_trav[iii][2]], 0};
							switch (iii)
							{
								void *inv_pars[4];
								case 0:
									if (nplaq_stat == 12) // 1100: truncated top when iii = 0
									{
										xc[1] = plaq_ycul[plaq_trav[iii][2]];
										xd[1] = plaq_ycul[plaq_trav[iii][3]];
									}
									else // 1001: truncated right side when iii = 0
									{
										// Determine whether coordinates are above y_max or below y_min
										if (xb[1] > plaq_ycul[plaq_trav[iii][1]])
										{
											inv_pars[0] = (void *) y_max;
											inv_pars[1] = (void *) y_max_pars;
										}
										else if (xb[1] < plaq_ycll[plaq_trav[iii][1]])
										{
											inv_pars[0] = (void *) y_min;
											inv_pars[1] = (void *) y_min_pars;
										}
										else
										{
											printf("Error: inconsistency in evaluating boundary contribution to integral. Case %d%d%d%d doesn't match coordinate (%g %g) with associated bounds y_min, y_max = %g %g\n", (nplaq_stat >> 3) & 1, (nplaq_stat >> 2) & 1, (nplaq_stat >> 1) & 1, nplaq_stat & 1, xb[0], xb[1], plaq_ycll[plaq_trav[iii][1]], plaq_ycul[plaq_trav[iii][1]]);
											exit(EXIT_FAILURE);
										}
										inv_pars[2] = (void *) &xa[0];
										inv_pars[3] = (void *) &xb[0];
										xb[0] = invert_bisection_method(xb[1], (void *) inv_pars);
										if (xc[1] > plaq_ycul[plaq_trav[iii][2]])
										{
											inv_pars[0] = (void *) y_max;
											inv_pars[1] = (void *) y_max_pars;
										}
										else if (xc[1] < plaq_ycll[plaq_trav[iii][2]])
										{
											inv_pars[0] = (void *) y_min;
											inv_pars[1] = (void *) y_min_pars;
										}
										else
										{
											printf("Error: inconsistency in evaluating boundary contribution to integral. Case %d%d%d%d doesn't match coordinate (%g %g) with associated bounds y_min, y_max = %g %g\n", (nplaq_stat >> 3) & 1, (nplaq_stat >> 2) & 1, (nplaq_stat >> 1) & 1, nplaq_stat & 1, xc[0], xc[1], plaq_ycll[plaq_trav[iii][2]], plaq_ycul[plaq_trav[iii][2]]);
											exit(EXIT_FAILURE);
										}
										inv_pars[2] = (void *) &xa[0];
										inv_pars[3] = (void *) &xc[0];
										xc[0] = invert_bisection_method(xc[1], (void *) inv_pars);
									}
									break;
								case 1:
									if (nplaq_stat == 12) // 1100: truncated side left when iii = 1
									{
										// Determine whether coordinates are above y_max or below y_min
										if (xc[1] > plaq_ycul[plaq_trav[iii][2]])
										{
											inv_pars[0] = (void *) y_max;
											inv_pars[1] = (void *) y_max_pars;
										}
										else if (xc[1] < plaq_ycll[plaq_trav[iii][2]])
										{
											inv_pars[0] = (void *) y_min;
											inv_pars[1] = (void *) y_min_pars;
										}
										else
										{
											printf("Error: inconsistency in evaluating boundary contribution to integral. Case %d%d%d%d doesn't match coordinate (%g %g) with associated bounds y_min, y_max = %g %g\n", (nplaq_stat >> 3) & 1, (nplaq_stat >> 2) & 1, (nplaq_stat >> 1) & 1, nplaq_stat & 1, xc[0], xc[1], plaq_ycll[plaq_trav[iii][2]], plaq_ycul[plaq_trav[iii][2]]);
											exit(EXIT_FAILURE);
										}
										inv_pars[2] = (void *) &xc[0];
										inv_pars[3] = (void *) &xa[0];
										xc[0] = invert_bisection_method(xc[1], (void *) inv_pars);
										if (xd[1] > plaq_ycul[plaq_trav[iii][3]])
										{
											inv_pars[0] = (void *) y_max;
											inv_pars[1] = (void *) y_max_pars;
										}
										else if (xd[1] < plaq_ycll[plaq_trav[iii][3]])
										{
											inv_pars[0] = (void *) y_min;
											inv_pars[1] = (void *) y_min_pars;
										}
										else
										{
											printf("Error: inconsistency in evaluating boundary contribution to integral. Case %d%d%d%d doesn't match coordinate (%g %g) with associated bounds y_min, y_max = %g %g\n", (nplaq_stat >> 3) & 1, (nplaq_stat >> 2) & 1, (nplaq_stat >> 1) & 1, nplaq_stat & 1, xd[0], xd[1], plaq_ycll[plaq_trav[iii][3]], plaq_ycul[plaq_trav[iii][3]]);
											exit(EXIT_FAILURE);
										}
										inv_pars[2] = (void *) &xd[0];
										inv_pars[3] = (void *) &xa[0];
										xd[0] = invert_bisection_method(xd[1], (void *) inv_pars);
									}
									else // 1001: truncated top when iii = 1 (i.e. upper left quadrant)
									{
										xb[1] = plaq_ycul[plaq_trav[iii][1]];
										xc[1] = plaq_ycul[plaq_trav[iii][2]];
									}
									break;
								case 2:
									// Lower left quadrant
									if (nplaq_stat == 12) // 1100: truncated bottom when iii = 2
									{
										xc[1] = plaq_ycll[plaq_trav[iii][2]];
										xd[1] = plaq_ycll[plaq_trav[iii][3]];
									}
									else // 1001: truncated left when iii = 2 (b, c)
									{
										// Determine whether coordinates are above y_max or below y_min
										if (xb[1] > plaq_ycul[plaq_trav[iii][1]])
										{
											inv_pars[0] = (void *) y_max;
											inv_pars[1] = (void *) y_max_pars;
										}
										else if (xb[1] < plaq_ycll[plaq_trav[iii][1]])
										{
											inv_pars[0] = (void *) y_min;
											inv_pars[1] = (void *) y_min_pars;
										}
										else
										{
											printf("Error: inconsistency in evaluating boundary contribution to integral. Case %d%d%d%d doesn't match coordinate (%g %g) with associated bounds y_min, y_max = %g %g\n", (nplaq_stat >> 3) & 1, (nplaq_stat >> 2) & 1, (nplaq_stat >> 1) & 1, nplaq_stat & 1, xb[0], xb[1], plaq_ycll[plaq_trav[iii][1]], plaq_ycul[plaq_trav[iii][1]]);
											exit(EXIT_FAILURE);
										}
										inv_pars[2] = (void *) &xb[0];
										inv_pars[3] = (void *) &xa[0];
										xb[0] = invert_bisection_method(xb[1], (void *) inv_pars);
										if (xc[1] > plaq_ycul[plaq_trav[iii][2]])
										{
											inv_pars[0] = (void *) y_max;
											inv_pars[1] = (void *) y_max_pars;
										}
										else if (xc[1] < plaq_ycll[plaq_trav[iii][2]])
										{
											inv_pars[0] = (void *) y_min;
											inv_pars[1] = (void *) y_min_pars;
										}
										else
										{
											printf("Error: inconsistency in evaluating boundary contribution to integral. Case %d%d%d%d doesn't match coordinate (%g %g) with associated bounds y_min, y_max = %g %g\n", (nplaq_stat >> 3) & 1, (nplaq_stat >> 2) & 1, (nplaq_stat >> 1) & 1, nplaq_stat & 1, xc[0], xc[1], plaq_ycll[plaq_trav[iii][2]], plaq_ycul[plaq_trav[iii][2]]);
											exit(EXIT_FAILURE);
										}
										inv_pars[2] = (void *) &xc[0];
										inv_pars[3] = (void *) &xa[0];
										xc[0] = invert_bisection_method(xc[1], (void *) inv_pars);
									}
									break;
								case 3: // Lower right (RESUME HP: see if it might be possible to combine this with case 1, and likewise 
									//		for cases 0 and 2 UPDATE: very doable: requires bringing iii dependence within nplaq_stat branch.)
									if (nplaq_stat == 12) // 1100: truncated side right when iii = 3
									{
										// Determine whether coordinates are above y_max or below y_min
										if (xc[1] > plaq_ycul[plaq_trav[iii][2]])
										{
											inv_pars[0] = (void *) y_max;
											inv_pars[1] = (void *) y_max_pars;
										}
										else if (xc[1] < plaq_ycll[plaq_trav[iii][2]])
										{
											inv_pars[0] = (void *) y_min;
											inv_pars[1] = (void *) y_min_pars;
										}
										else
										{
											printf("Error: inconsistency in evaluating boundary contribution to integral. Case %d%d%d%d doesn't match coordinate (%g %g) with associated bounds y_min, y_max = %g %g\n", (nplaq_stat >> 3) & 1, (nplaq_stat >> 2) & 1, (nplaq_stat >> 1) & 1, nplaq_stat & 1, xc[0], xc[1], plaq_ycll[plaq_trav[iii][2]], plaq_ycul[plaq_trav[iii][2]]);
											exit(EXIT_FAILURE);
										}
										inv_pars[2] = (void *) &xa[0];
										inv_pars[3] = (void *) &xc[0];
										xc[0] = invert_bisection_method(xc[1], (void *) inv_pars);
										if (xd[1] > plaq_ycul[plaq_trav[iii][3]])
										{
											inv_pars[0] = (void *) y_max;
											inv_pars[1] = (void *) y_max_pars;
										}
										else if (xd[1] < plaq_ycll[plaq_trav[iii][3]])
										{
											inv_pars[0] = (void *) y_min;
											inv_pars[1] = (void *) y_min_pars;
										}
										else
										{
											printf("Error: inconsistency in evaluating boundary contribution to integral. Case %d%d%d%d doesn't match coordinate (%g %g) with associated bounds y_min, y_max = %g %g\n", (nplaq_stat >> 3) & 1, (nplaq_stat >> 2) & 1, (nplaq_stat >> 1) & 1, nplaq_stat & 1, xd[0], xd[1], plaq_ycll[plaq_trav[iii][3]], plaq_ycul[plaq_trav[iii][3]]);
											exit(EXIT_FAILURE);
										}
										inv_pars[2] = (void *) &xa[0];
										inv_pars[3] = (void *) &xd[0];
										xd[0] = invert_bisection_method(xd[1], (void *) inv_pars);
									}
									else // 1001: truncated bottom when iii = 3 (i.e. lower right quadrant)
									{
										xb[1] = plaq_ycll[plaq_trav[iii][1]];
										xc[1] = plaq_ycll[plaq_trav[iii][2]];
									}
									break;
							}
							xb[2] = f(xb[0], xb[1], fpars);
							xc[2] = f(xc[0], xc[1], fpars);
							xd[2] = f(xd[0], xd[1], fpars);
							bdry_est += quad_prism_vol(xa[2], xb[2], xc[2], xd[2], xa, xb, xc, xd);
						}
						else if (nplaq_stat == 11 || nplaq_stat == 13 || nplaq_stat == 14)
						{
							// Pentagon case: 1011, 1101, 1110,
							double xb[3], xc[3], xd[3], xe[3];
							switch (iii)
							{
								case 0:
									xb[1] = y1;
									xc[0] = x2;
									xd[0] = x2;
									xd[1] = y2;
									xe[0] = x1;
									xe[1] = y2;
									if (nplaq_stat == 11) // 1011 => missing lower right corner => y_min relevant
									{
										// Solve for xb[0]
										void *inv_pars[4];
										inv_pars[0] = (void *) y_min;
										inv_pars[1] = (void *) y_min_pars;
										inv_pars[2] = (void *) &xa[0];
										inv_pars[3] = (void *) &plaq_xc[plaq_trav[iii][1]];
										xb[0] = invert_bisection_method(xb[1], (void *) inv_pars);
										xc[1] = plaq_ycll[plaq_trav[iii][1]];
									}
									else if (nplaq_stat == 13) // 1101 => missing upper right corner => y_max relevant
									{
										xb[0] = plaq_xc[plaq_trav[iii][1]];
										xc[1] = plaq_ycul[plaq_trav[iii][1]];
										void *inv_pars[4];
										inv_pars[0] = (void *) y_max;
										inv_pars[1] = (void *) y_max_pars;
										inv_pars[2] = (void *) &xa[0];
										inv_pars[3] = (void *) &xd[0];
										xd[0] = invert_bisection_method(xd[1], (void *) inv_pars);
									}
									else if (nplaq_stat == 14) // 1110 => missing upper left corner => y_max relevant
									{
										xb[0] = x2;
										xc[1] = y2;
										void *inv_pars[4];
										xd[0] = x1;
										inv_pars[0] = (void *) y_max;
										inv_pars[1] = (void *) y_max_pars;
										inv_pars[2] = (void *) &xd[0];
										inv_pars[3] = (void *) &xb[0];
										xd[0] = invert_bisection_method(xd[1], (void *) inv_pars);
										xe[1] = plaq_ycul[plaq_trav[iii][3]];
									}
									break;
								case 1: // Upper left quadrant: (x1, x1, x0, x0), (y1, y2, y2, y1)
									if (nplaq_stat == 11) // 1011 => missing upper right corner => y_max relevant
									{
										xb[0] = x1;
										xb[1] = y2;
										// Solve for xc[0]
										void *inv_pars[4];
										inv_pars[0] = (void *) y_max;
										inv_pars[1] = (void *) y_max_pars;
										inv_pars[2] = (void *) &x0;
										inv_pars[3] = (void *) &x1;
										xc[1] = y2;
										xc[0] = invert_bisection_method(xc[1], (void *) inv_pars);
										xd[0] = x0;
										xd[1] = y2;
										xe[0] = x0;
										xe[1] = y1;
									}
									else if (nplaq_stat == 13) // 1101 => missing upper left corner => y_max relevant
									{
										xb[0] = x1;
										xb[1] = y2;
										void *inv_pars[4];
										inv_pars[0] = (void *) y_max;
										inv_pars[1] = (void *) y_max_pars;
										inv_pars[2] = (void *) &x0;
										inv_pars[3] = (void *) &x1;
										xc[1] = y2;
										xc[0] = invert_bisection_method(xc[1], (void *) inv_pars);
										xd[0] = x0;
										xd[1] = plaq_ycul[plaq_trav[iii][3]];
										xe[0] = x0;
										xe[1] = y1;
									}
									else if (nplaq_stat == 14) // 1110 => missing lower left corner => y_min relevant
									{
										xb[0] = x1;
										xb[1] = y2;
										xc[0] = x0;
										xc[1] = y2;
										void *inv_pars[4];
										xd[0] = x1;
										xd[1] = plaq_ycll[plaq_trav[iii][3]];
										inv_pars[0] = (void *) y_min;
										inv_pars[1] = (void *) y_min_pars;
										inv_pars[2] = (void *) &x0;
										inv_pars[3] = (void *) &x1;
										xe[1] = y1;
										xe[0] = invert_bisection_method(xe[1], (void *) inv_pars);
									}
									break;
								case 2: // Lower left quadrant: (x1, x0, x0, x1), (y1, y1, y0, y0)
									if (nplaq_stat == 11) // 1011 => missing upper left corner => y_max relevant
									{
										// Solve for xb[0]
										void *inv_pars[4];
										inv_pars[0] = (void *) y_max;
										inv_pars[1] = (void *) y_max_pars;
										inv_pars[2] = (void *) &x0;
										inv_pars[3] = (void *) &x1;
										xb[1] = y1;
										xb[0] = invert_bisection_method(xb[1], (void *) inv_pars);
										xc[0] = x0;
										xc[1] = plaq_ycul[plaq_trav[iii][2]];
										xd[0] = x0;
										xd[1] = y0;
										xe[0] = x1;
										xe[1] = y0;
									}
									else if (nplaq_stat == 13) // 1101 => missing lower left corner (xb normal)
									{
										xb[0] = x0;
										xb[1] = y1;
										xc[1] = plaq_ycll[plaq_trav[iii][2]];
										xc[0] = x0;
										void *inv_pars[4];
										inv_pars[0] = (void *) y_min;
										inv_pars[1] = (void *) y_min_pars;
										inv_pars[2] = (void *) &x0;
										inv_pars[3] = (void *) &x1;
										xd[1] = y0;
										xd[0] = invert_bisection_method(xd[1], (void *) inv_pars);
										xe[0] = x1;
										xe[1] = y1;
									}
									else if (nplaq_stat == 14) // 1110 => missing lower right corner => y_min relevant
									{
										xb[0] = x0;
										xb[1] = y1;
										xc[0] = x0;
										xc[1] = y0;
										void *inv_pars[4];
										inv_pars[0] = (void *) y_min;
										inv_pars[1] = (void *) y_min_pars;
										inv_pars[2] = (void *) &x0;
										inv_pars[3] = (void *) &x1;
										xd[1] = y0;
										xd[0] = invert_bisection_method(xd[1], (void *) inv_pars);
										xe[1] = plaq_ycul[plaq_trav[iii][3]];
									}
									break;
								case 3: // Lower right quadrant: (x1, x1, x2, x2), (y1, y0, y0, y1)
									if (nplaq_stat == 11) // 1011 => missing lower left corner => y_min relevant
									{
										xb[0] = x1;
										xb[1] = plaq_ycll[plaq_trav[iii][1]];
										// Solve for xc[0]
										void *inv_pars[4];
										inv_pars[0] = (void *) y_min;
										inv_pars[1] = (void *) y_min_pars;
										inv_pars[2] = (void *) &x1;
										inv_pars[3] = (void *) &x2;
										xc[1] = y0;
										xc[0] = invert_bisection_method(xc[1], (void *) inv_pars);
										xd[0] = x2;
										xd[1] = y0;
										xe[0] = x2;
										xe[1] = y1;
									}
									else if (nplaq_stat == 13) // 1101 => missing lower right corner => y_min relevant
									{
										xb[0] = x1;
										xb[1] = y0;
										void *inv_pars[4];
										inv_pars[0] = (void *) y_min;
										inv_pars[1] = (void *) y_min_pars;
										inv_pars[2] = (void *) &x1;
										inv_pars[3] = (void *) &x2;
										xc[1] = y0;
										xc[0] = invert_bisection_method(xc[1], (void *) inv_pars);
										xd[0] = x2;
										xd[1] = plaq_ycll[plaq_trav[iii][3]];
										xe[0] = x2;
										xe[1] = y1;
									}
									else if (nplaq_stat == 14) // 1110 => missing upper right corner => y_max relevant
									{
										xb[0] = x1;
										xb[1] = y0;
										xc[0] = x2;
										xc[1] = y0;
										xd[0] = x2;
										xd[1] = plaq_ycul[plaq_trav[iii][3]]; 
										void *inv_pars[4];
										inv_pars[0] = (void *) y_max;
										inv_pars[1] = (void *) y_max_pars;
										inv_pars[2] = (void *) &x1;
										inv_pars[3] = (void *) &x2;
										xe[1] = y1;
										xe[0] = invert_bisection_method(xe[1], (void *) inv_pars);
									}
									break;
							}
							xb[2] = f(xb[0], xb[1], fpars);
							xc[2] = f(xc[0], xc[1], fpars);
							xd[2] = f(xd[0], xd[1], fpars);
							xe[2] = f(xe[0], xe[1], fpars);
							bdry_est += pentagon_prism_vol(xa[2], xb[2], xc[2], xd[2], xe[2], xa, xb, xc, xd, xe);
						}
						else if (nplaq_stat == 10)
						{
							// Hexagon case: 1010
							double xb[3], xc[3], xd[3], xe[3], xf[3];
							void *inv_pars[4];
							switch (iii)
							{
								case 0: // Upper right: missing lower right and upper left: 
									// (x1, x2, x2, x1), (y1, y1, y2, y2)
									inv_pars[0] = (void *) y_min;
									inv_pars[1] = (void *) y_min_pars;
									inv_pars[2] = &x1;
									inv_pars[3] = &x2;
									xb[1] = y1;
									xb[0] = invert_bisection_method(xb[1], (void *) inv_pars);
									inv_pars[0] = (void *) y_max;
									inv_pars[1] = (void *) y_max_pars;
									xe[1] = y2;
									xe[0] = invert_bisection_method(xe[1], (void *) inv_pars);
									xc[0] = x2;
									xc[1] = plaq_ycll[plaq_trav[iii][1]];
									xd[0] = x2;
									xd[1] = y2;
									xf[0] = x1;
									xf[1] = plaq_ycul[plaq_trav[iii][3]];
									break;
								case 1: // Upper left: missing upper right, lower left 1010
									// (x1, x1, x0, x0), (y1, y2, y2, y1)
									xb[0] = x1;
									xb[1] = plaq_ycul[plaq_trav[iii][1]];
									inv_pars[0] = (void *) y_max;
									inv_pars[1] = (void *) y_max_pars;
									inv_pars[2] = (void *) &x0;
									inv_pars[3] = (void *) &x1;
									xc[1] = y2;
									xc[0] = invert_bisection_method(xc[1], (void *) inv_pars);
									xd[0] = x0;
									xd[1] = y2;
									xe[1] = plaq_ycll[plaq_trav[iii][3]];
									xe[0] = x0;
									inv_pars[0] = (void *) y_min;
									inv_pars[1] = (void *) y_min_pars;
									xf[1] = y1;
									xf[0] = invert_bisection_method(xf[1], (void *) inv_pars);
									break;
								case 2: // Lower left quadrant: missing upper left and lower right corners
									// (x1, x0, x0, x1), (y1, y1, y0, y0)
									inv_pars[0] = (void *) y_max;
									inv_pars[1] = (void *) y_max_pars;
									inv_pars[2] = (void *) &x0;
									inv_pars[3] = (void *) &x1;
									xb[1] = y1;
									xb[0] = invert_bisection_method(xb[1], (void *) inv_pars);
									xc[0] = x0;
									xc[1] = plaq_ycul[plaq_trav[iii][1]];
									xd[0] = x0;
									xd[1] = y0;
									inv_pars[0] = (void *) y_min;
									inv_pars[1] = (void *) y_min_pars;
									xe[1] = y0;
									xe[0] = invert_bisection_method(xe[1], (void *) inv_pars);
									xf[0] = x1;
									xf[1] = plaq_ycll[plaq_trav[iii][3]];
									break;
								case 3: // Lower right quadrant: missing lower left and upper right corners
									// (x1, x1, x2, x2), (y1, y0, y0, y1)
									xb[0] = x1;
									xb[1] = plaq_ycll[plaq_trav[iii][1]];
									inv_pars[0] = (void *) y_min;
									inv_pars[1] = (void *) y_min_pars;
									inv_pars[2] = (void *) &x1;
									inv_pars[3] = (void *) &x2;
									xc[1] = y0;
									xc[0] = invert_bisection_method(xc[1], (void *) inv_pars);
									xd[0] = x2;
									xd[1] = y0;
									xe[0] = x2;
									xe[1] = plaq_ycul[plaq_trav[iii][3]];
									inv_pars[0] = (void *) y_max;
									inv_pars[1] = (void *) y_max_pars;
									xf[1] = y1;
									xf[0] = invert_bisection_method(xf[1], (void *) inv_pars);
									break;
							}
							xb[2] = f(xb[0], xb[1], fpars);
							xc[2] = f(xc[0], xc[1], fpars);
							xd[2] = f(xd[0], xd[1], fpars);
							xe[2] = f(xe[0], xe[1], fpars);
							xf[2] = f(xf[0], xf[1], fpars);
							bdry_est += hexagon_prism_vol(xa[2], xb[2], xc[2], xd[2], xe[2], xf[2], xa, xb, xc, xd, xe, xf);
						}
					}
				}
			}
			y0 = y1;
			y1 = y2;
			y2 += dx;
			iim2 = iim1;
			iim1 = ii;
		}
		est += est_y * dx * dx;
		x0 = x1;
		x1 = x2;
		x2 += dx;
		ymin0 = ymin1;
		ymax0 = ymax1;
		ymin1 = ymin2;
		ymax1 = ymax2;
		ymin2 = y_min(x2, y_min_pars);
		ymax2 = y_max(x2, y_max_pars);
		im2 = im1;
		im1 = i;
	}
	free_aarray_int(&bdry_plaquettes);
	free_array_int(&bdry_plaq_stat);
	return est + bdry_est;
}

approx_double integrate_2D(double (*f)(double, double, void *), void *pars, double (*domain)(double, double, void *), void *domain_pars, double *xbnds, double *ybnds, double tol)
{
	double dx;
	double dy;
	// Define initial set of sample points
	int Nx = (int) ((xbnds[1] - xbnds[0]) / tol) + 1;
	int Ny = (int) ((ybnds[1] - ybnds[0]) / tol) + 1;
	dx = (xbnds[1] - xbnds[0]) / (Nx - 1);
	dy = (ybnds[1] - ybnds[0]) / (Ny - 1);
	double err = 1;
	double est = 0;
	int N_total = Nx * Ny;
	while (err > tol * est && N_total < N_EVAL_LIM) 
	{
		dx = (xbnds[1] - xbnds[0]) / (Nx - 1);
		dy = (ybnds[1] - ybnds[0]) / (Ny - 1);
		double dA = dx * dy;
		err = 0;
		est = 0;
		double x = xbnds[0];
		for (int i = 0; i < Nx; i++)
		{
			double y = ybnds[0];
			char interior = domain(x, y, domain_pars) > 0;
			for (int ii = 0; ii < Ny; ii++)
			{
				if (domain(x, y, domain_pars) > 0)
				{
					double incr = f(x, y, pars);
					est += incr;
					if (interior) {}
					else
					{
						interior = 1;
						err += incr;
					}
				}
				else
				{
					if (!interior) {}
					else
					{
						interior = 0;
						double incr = f(x, y - dy, pars);
						err += incr;
					}
				}
				y += dy;
			}
			x += dx;
		}
		est *= dA;
		err *= dA;
		Nx <<= 1;
		Ny <<= 1;
		N_total = Nx * Ny;
	}
	approx_double est_err;
	est_err.val = est + 0.5 * err;
	est_err.err = err;
	return est_err;
}

double test_func2D(double y, double z, void *pars)
{
	void **pars_ = (void **) pars;
	double Rsq = *((double *) pars_[0]);
	double rsq = *((double *) pars_[1]);
	double dist = *((double *) pars_[2]);
	double zsq = z * z;
	double ysq = y * y;
	double discr1 = Rsq - ysq;
	double discr2 = rsq - ysq - zsq;
	if (discr1 > 0 && discr2 > 0) return sqrt(discr1) + sqrt(discr2) - dist;
	else return -1;
}

double test_func2D_area(double y, double z, void *pars)
{
	void **pars_ = (void **) pars;
	double R_surfsq = *((double *) pars_[0]);
	double ysq = y * y;
	double arg = 1 + ysq / (R_surfsq - ysq);
	if (arg >= 0) return sqrt(arg);
	else return 0;
}

double quadratic_form_2D(double x, double y, void *pars)
{
	void **pars_ = (void **) pars;
	double Cxx, Cyy, Cxy;
	Cxx = *((double *) pars_[0]);
	Cyy = *((double *) pars_[1]);
	Cxy = *((double *) pars_[2]);
	return (Cxx * x + Cxy * y) * x + Cyy * y * y;
}

double sqrt1D(double x, void * pars)
{
	return sqrt(x);
}

double composite_2D(double x, double y, void *pars)
{
	void **pars_ = (void **) pars;
	double (*base_func)(double, double, void *) = (double (*)(double, double, void*)) pars_[0];
	void *base_pars = pars_[1];
	double (*comp_func)(double, void *) = (double (*)(double, void *)) pars_[2];
	void *comp_pars = pars_[3];
	return comp_func(base_func(x, y, base_pars), comp_pars);
}

double composite_1D(double x, void *pars)
{
	void **pars_ = (void **) pars;
	double (*base_func)(double, void *) = (double (*)(double, void*)) pars_[0];
	void *base_pars = pars_[1];
	double (*comp_func)(double, void *) = (double (*)(double, void*)) pars_[2];
	void *comp_pars = pars_[3];
	return comp_func(base_func(x, base_pars), comp_pars);
}

double d2func_difference(double x, double y, void *pars)
{
	void **pars_ = (void **) pars;
	double (*func1)(double, double, void *) = (double (*)(double, double, void*)) pars_[0];
	double (*func2)(double, double, void *) = (double (*)(double, double, void*)) pars_[2];
	return func1(x, y, pars_[1]) - func2(x, y, pars_[3]);
}

double d2func_add(double x, double y, void *pars)
{
	void **pars_ = (void **) pars;
	double (*func1)(double, double, void *) = (double (*)(double, double, void*)) pars_[0];
	double (*func2)(double, double, void *) = (double (*)(double, double, void*)) pars_[2];
	return func1(x, y, pars_[1]) + func2(x, y, pars_[3]);
}

double d2func_linear(double x, double y, void *pars)
{
	void **pars_ = (void **) pars;
	double cx = *((double *) pars_[1]);
	double cy = *((double *) pars_[2]);
	double c0 = *((double *) pars_[0]);
	return c0 + cx * x + cy * y;
}

double d1func_linear(double x, void *pars)
{
	void **pars_ = (void **) pars;
	double m = *((double *) pars_[1]);
	double b = *((double *) pars_[0]);
	return m * x + b;
}
double d2func_product(double x, double y, void *pars)
{
	void **pars_ = (void **) pars;
	double (*func1)(double, double, void *) = (double (*)(double, double, void*)) pars_[0];
	double (*func2)(double, double, void *) = (double (*)(double, double, void*)) pars_[2];
	return func1(x, y, pars_[1]) * func2(x, y, pars_[3]);
}

void parse_gslr_type(const gsl_rng_type **gslr_type, char *gslr_mode)
{
	toupper_string(gslr_mode);
	if (strcmp("GSL_RNG_MT19937", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_mt19937;
	}
	if (strcmp("GSL_RNG_RANLXS0", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_ranlxs0;
	}
	if (strcmp("GSL_RNG_RANLXS1", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_ranlxs1;
	}
	if (strcmp("GSL_RNG_RANLXS2", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_ranlxs2;
	}
	if (strcmp("GSL_RNG_RANLXD1", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_ranlxd1;
	}
	if (strcmp("GSL_RNG_RANLXD2", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_ranlxd2;
	}
	if (strcmp("GSL_RNG_RANLUX", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_ranlux;
	}
	if (strcmp("GSL_RNG_RANLUX389", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_ranlux389;
	}
	if (strcmp("GSL_RNG_CMRG", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_cmrg;
	}
	if (strcmp("GSL_RNG_MRG", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_mrg;
	}
	if (strcmp("GSL_RNG_TAUS", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_taus;
	}
	if (strcmp("GSL_RNG_TAUS2", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_taus2;
	}
	if (strcmp("GSL_RNG_GFSR4", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_gfsr4;
	}
	if (strcmp("GSL_RNG_RAND", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_rand;
	}
	if (strcmp("GSL_RNG_RANDOM_BSD", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_random_bsd;
	}
	if (strcmp("GSL_RNG_RANDOM_LIBC5", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_random_libc5;
	}
	if (strcmp("GSL_RNG_RANDOM_GLIBC2", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_random_glibc2;
	}
	if (strcmp("GSL_RNG_RAND48", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_rand48;
	}
	if (strcmp("GSL_RNG_RANF", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_ranf;
	}
	if (strcmp("GSL_RNG_RANMAR", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_ranmar;
	}
	if (strcmp("GSL_RNG_R250", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_r250;
	}
	if (strcmp("GSL_RNG_TT800", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_tt800;
	}
	if (strcmp("GSL_RNG_VAX", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_vax;
	}
	if (strcmp("GSL_RNG_TRANSPUTER", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_transputer;
	}
	if (strcmp("GSL_RNG_RANDU", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_randu;
	}
	if (strcmp("GSL_RNG_MINSTD", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_minstd;
	}
	if (strcmp("GSL_RNG_UNI", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_uni;
	}
	if (strcmp("GSL_RNG_UNI32", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_uni32;
	}
	if (strcmp("GSL_RNG_SLATEC", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_slatec;
	}
	if (strcmp("GSL_RNG_ZUF", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_zuf;
	}
	if (strcmp("GSL_RNG_KNUTHRAN2", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_knuthran2;
	}
	if (strcmp("GSL_RNG_KNUTHRAN2002", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_knuthran2002;
	}
	if (strcmp("GSL_RNG_KNUTHRAN", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_knuthran;
	}
	if (strcmp("GSL_RNG_BOROSH13", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_borosh13;
	}
	/*if (strcmp("GSL_RNG_BOROSH18", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_borosh18;
	}*/
	if (strcmp("GSL_RNG_FISHMAN20", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_fishman20;
	}
	if (strcmp("GSL_RNG_LECUYER21", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_lecuyer21;
	}
	if (strcmp("GSL_RNG_WATERMAN14", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_waterman14;
	}
	if (strcmp("GSL_RNG_FISHMAN2X", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_fishman2x;
	}
	if (strcmp("GSL_RNG_COVEYOU", gslr_mode) == 0)
	{
		*gslr_type = gsl_rng_coveyou;
	}
}

double invert_bisection_method(double a, void *pars)
{
	void **pars_ = (void **) pars;
	double (*func)(double, void *) = (double (*)(double, void *)) pars_[0];
	void *func_pars = (void *) pars_[1];
	double min_ = *((double *) pars_[2]);
	double max_ = *((double *) pars_[3]);
	double mid_;
	if (min_ < max_) {}
	else
	{
		mid_ = min_;
		min_ = max_;
		max_ = mid_;
	}
	mid_ = 0.5 * (min_ + max_);
	char lle = func(min_, func_pars) > 0;
	char ule = func(max_, func_pars) > 0;
	char mle = func(mid_, func_pars) > 0;
	double diff = max_ - min_;
	int N_attempts = 0;
	while (diff > STANDARD_TOL && N_attempts < N_EVAL_LIM) 
	{
		N_attempts += 1;
		if (mle == lle)
		{
			min_ = mid_;
		}
		else
		{
			max_ = mid_;
		}
		mid_ = 0.5 * (min_ + max_);
		mle = func(mid_, func_pars) > 0;
	}
	return mid_;
}

/*
	min sum_{ij} U(x_i - x_j), f(x_i) = 0, 1 <= i<= N
	sum_j nabla U(x_i - x_j) = lambda_i nabla f(x_i), 1 <= i <= N (3N constraints)
	f(x_i) = 0, 1 <= i <= N (N constraints)
	=> apply Newton's method (or Broyden q method) to F(x, lambda)
	Problem: not especially parallelizeable 
*/
