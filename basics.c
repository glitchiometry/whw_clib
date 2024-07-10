#include "basics.h"

double ** bb_matrix_alloc(int M, int N)
{
	double **A = (double **) calloc(M, sizeof(double *));
	for (int i = 0; i < M; i++) A[i] = (double *) calloc(N, sizeof(double));
	return A;
}

void bb_matrix_free(double **A, int M)
{
	for (int i = 0; i < M; i++) free(A[i]);
	free(A);
}

// Routines for variable length arrays of void pointers
void array_voidstar_init(array_voidstar *a, int size_)
{
	if (size_ > 0)
	{
		(*a).e = (void **) calloc(size_, sizeof(void *));
		(*a).mem = size_;
		for (int i = 0; i < (*a).mem; i++) (*a).e[i] = NULL;
	}
	else 
	{
		(*a).e = NULL;
		(*a).mem = 0;
	}
	(*a).len = 0;
}

void add_mem_array_voidstar(array_voidstar *a)
{
	if ((*a).mem > 0) {}
	else (*a).mem = 1;
	(*a).mem <<= 1;
	void **ne = (void **) calloc((*a).mem, sizeof(void *));
	for (int i = 0; i < (*a).len; i++)
	{
		ne[i] = (*a).e[i];
	}
	for (int i = (*a).len; i < (*a).mem; i++)
	{
		ne[i] = NULL;
	}
	if ((*a).e != NULL) free((*a).e);
	(*a).e = ne;
}

void add_mem_array_voidstar_until(array_voidstar *a, int i)
{
	if ((*a).mem > 0) {}
	else (*a).mem = 1;
	while ((*a).mem <= i)
	{
		(*a).mem <<= 1;
	}
	void **ne = (void **) calloc((*a).mem, sizeof(void *));
	for (int ii = 0; ii < (*a).len; ii++)
	{
		ne[ii] = (*a).e[ii];
	}
	for (int ii = (*a).len; ii < (*a).mem; ii++) ne[ii] = NULL;
	if ((*a).e != NULL) free((*a).e);
	(*a).e = ne;
}

void add2array_voidstar(array_voidstar *a, void *i)
{
	if ((*a).len < (*a).mem) {}
	else
	{
		add_mem_array_voidstar(a);
	}
	(*a).e[(*a).len] = i;
	(*a).len += 1;
}

void print_array_voidstar(array_voidstar a, aarray_char *fmt)
{
	for (int i = 0; i < a.len; i++)
	{
		printf((*fmt).e[i].e, a.e[i]);
	}
	printf("\n");
}

void remove_last_array_voidstar(array_voidstar *a, void (*free_elem)(void *))
{
	if ((*a).len > 0)
	{
		(*a).len -= 1;
		if ((*a).e[(*a).len] != NULL) 
		{
			if (free_elem != NULL) free_elem((*a).e[(*a).len]);
			free((*a).e[(*a).len]); // RESUME: check this!
			(*a).e[(*a).len] = NULL;
		}
		if ((*a).len < ((*a).mem >> 2)) contract_array_voidstar(a, free_elem);
	}
}

void remove_array_voidstar(array_voidstar *a, int n, void (*free_elem)(void *))
{
	if (n < (*a).len)
	{
		(*a).len -= 1;
		if (free_elem != NULL)
		{
			if ((*a).e[n] != NULL) 
			{
				free_elem((*a).e[n]);
				free((*a).e[n]); // RESUME: check consistency across other programs
			}
		}
		(*a).e[n] = (*a).e[(*a).len];
		(*a).e[(*a).len] = NULL;
		if ((*a).len < ((*a).mem >> 2)) contract_array_voidstar(a, free_elem);
	}
}

// Shrink array to reduce memory usage
void contract_array_voidstar(array_voidstar *a, void (*free_elem)(void *))
{
	if ((*a).mem > 0 && (*a).e != NULL)
	{
		int newsize = (*a).mem >> 1;
		void **data = (void **) calloc(newsize, sizeof(void *));
		int ulim = (*a).len < newsize ? (*a).len : newsize;
		for (int i = 0; i < ulim; i++) 
		{
			data[i] = (*a).e[i];
		}
		for (int i = ulim; i < newsize; i++)
		{
			data[i] = NULL;
		}
		if (free_elem != NULL)
		{
			for (int i = 0; i < (*a).mem; i++)
			{
				if ((*a).e[i] != NULL) 
				{
					free_elem((*a).e[i]);
					free((*a).e[i]); // RESUME: check this!
				}
			}
		}
		free((*a).e);
		(*a).e = data;
		(*a).mem = newsize;
	}
}

void free_elem_triv(void *) {}

void free_array_voidstar(array_voidstar *a, void (*free_elem)(void *))
{
	if ((*a).mem > 0)
	{
		if ((*a).e != NULL) 
		{
			if (free_elem != NULL)
			{
				for (int i = 0; i < (*a).len; i++) 
				{
					if ((*a).e[i] != NULL) 
					{
						free_elem((*a).e[i]);
						free((*a).e[i]);
					}
				}
			}
			free((*a).e);
			(*a).e = NULL;
		}
		(*a).len = 0;
		(*a).mem = 0;
	}
}

// NOTE: this function should be rarely used (and all pointers in the array should be freed independently.)
void reset_array_voidstar(array_voidstar *a)
{
	(*a).len = 0;
}

void transcribe_array_voidstar(array_voidstar *a, array_voidstar *b)
{
	add_mem_array_voidstar_until(a, (*b).len);
	(*a).len = (*b).len;
	for (int i = 0; i < (*b).len; i++)
	{
		(*a).e[i] = (*b).e[i];
	}
}

// methods for variable length arrays of integers
void array_int_init(array_int *a, int size_)
{
	if (size_ > 0)
	{
		(*a).e = (int *) calloc(size_, sizeof(int));
		(*a).mem = size_;
		(*a).len = 0;
	}
	else 
	{
		(*a).e = NULL;
		(*a).mem = 0;
		(*a).len = 0;
	}
}

void array_int_init_range(array_int *a, int i0, int i1)
{
	int len = i1 - i0;
	if (len > 0)
	{
		array_int_init(a, len);
		for (int i = i0; i < i1; i++)
		{
			add2array_int(a, i);
		}
	}
}

void transcribe_array_int(array_int *src, array_int *dest)
{
	array_int_init(dest, (*src).mem);
	(*dest).len = (*src).len;
	for (int i = 0; i < (*src).len; i++)
	{
		(*dest).e[i] = (*src).e[i];
	}
}

int array_int_search(array_int *a, int elem)
{
	for (int i = 0; i < (*a).len; i++)
	{
		if ((*a).e[i] == elem) return i;
	}
	return -1;
}

char array_int_contains(array_int *a, int elem)
{
	for (int i = 0; i < (*a).len; i++)
	{
		if ((*a).e[i] == elem) return 1;
	}
	return 0;
}

void add_mem_array_int(array_int *a)
{
	if ((*a).mem > 0) {}
	else (*a).mem = 1;
	(*a).mem <<= 1;
	int *ne = (int *) calloc((*a).mem, sizeof(int));
	for (int i = 0; i < (*a).len; i++)
	{
		ne[i] = (*a).e[i];
	}
	if ((*a).e != NULL) free((*a).e);
	(*a).e = ne;
}

void add_mem_array_int_until(array_int *a, int i)
{
	if ((*a).mem > i) return;
	if ((*a).mem > 0) {}
	else (*a).mem = 1;
	while ((*a).mem <= i)
	{
		(*a).mem <<= 1;
	}
	int *ne = (int *) calloc((*a).mem, sizeof(int));
	for (int ii = 0; ii < (*a).len; ii++)
	{
		ne[ii] = (*a).e[ii];
	}
	if ((*a).e != NULL) free((*a).e);
	(*a).e = ne;
}

void add2array_int(array_int *a, int i)
{
	if ((*a).len < (*a).mem) {}
	else
	{
		add_mem_array_int(a);
	}
	(*a).e[(*a).len] = i;
	(*a).len += 1;
}

void print_array_int(array_int a)
{
	for (int i = 0; i < a.len; i++)
	{
		printf("%d ", a.e[i]);
	}
	printf("\n");
}

void remove_last_array_int(array_int *a)
{
	if ((*a).len > 0)
	{
		(*a).len -= 1;
		if ((*a).len < ((*a).mem >> 2)) contract_array_int(a);
	}
}

void remove_array_int(array_int *a, int n)
{
	if (n < (*a).len && n > -1)
	{
		(*a).len -= 1;
		(*a).e[n] = (*a).e[(*a).len];
		if ((*a).len < ((*a).mem >> 2)) contract_array_int(a);
	}
}

void contract_array_int(array_int *a)
{
	if ((*a).mem > 0 && (*a).e != NULL) 
	{
		int newsize = (*a).mem >> 1;
		int *data = (int *) calloc(newsize, sizeof(int));
		int ulim = (*a).len < newsize ? (*a).len : newsize;
		for (int i = 0; i < ulim; i++) 
		{
			data[i] = (*a).e[i];
		}
		free((*a).e);
		(*a).e = data;
		(*a).mem = newsize;
	}
}

void fprintf_array_int(array_int *a, FILE *ofile)
{
	ofile = ofile != NULL ? ofile : stdout;
	for (int i = 0; i < (*a).len; i++) fprintf(ofile, "%d ", (*a).e[i]);
	fprintf(ofile, "\n");
}

void free_array_int(array_int *a)
{
	if ((*a).mem > 0)
	{
		if ((*a).e != NULL) free((*a).e);
		(*a).e = NULL;
		(*a).len = 0;
		(*a).mem = 0;
	}
}

void reset_array_int(array_int *a)
{
	(*a).len = 0;
}

void merge_array_int(array_int *a, int *e, int *c, int imin, int imax, int jmin, int jmax)
{
	if ((*a).e[imin] < (*a).e[jmin])
	{
		int il = imin;
		while ((*a).e[il] < (*a).e[jmin] && il < imax)
		{
			e[*c] = (*a).e[il];
			il += 1;
			*c += 1;
		}
		merge_array_int(a, e, c, il, imax, jmin, jmax);
	}
	else
	{
		int ir = jmin;
		while ((*a).e[ir] < (*a).e[imin] && ir < jmax)
		{
			e[*c] = (*a).e[ir];
			ir += 1;
			*c += 1;
		}
		merge_array_int(a, e, c, imin, imax, ir, jmax);
	}
}

// sort a variable length array along the index interval [i, j)
void merge_sort_array_int(array_int *a, int *e, int i, int j)
{
	int jmi = j - i;
	if (jmi > 2)
	{
		int mdpt = (i + j) >> 1;
		merge_sort_array_int(a, e, i, mdpt);
		merge_sort_array_int(a, e, mdpt, j);
		// merge the two ordered lists
		int c = i;
		merge_array_int(a, e, &c, i, mdpt, mdpt, j);
	}
	else if (jmi == 2)
	{
		int ip1 = i + 1;
		if ((*a).e[ip1] > (*a).e[i]) {}
		else
		{
			int tmp = (*a).e[ip1];
			(*a).e[ip1] = (*a).e[i];
			(*a).e[i] = tmp;
		}
	}
}

char is_sorted_array_int(array_int *a)
{
	if ((*a).len > 1)
	{
		char state_inc = 1;
		char state_dec = 1;
		int im1 = 0;
		for (int i = 1; i < (*a).len; i++)
		{
			if ((*a).e[i] > (*a).e[im1]) state_dec = 0;
			if ((*a).e[i] < (*a).e[im1]) state_inc = 0;
			if (state_inc || state_dec) {}
			else return 0;
			im1 = i;
		}
		if (state_inc && !state_dec) return SORTED_INCREASING;
		if (state_dec && !state_inc) return SORTED_DECREASING;
		if (state_dec && state_inc) return SORTED_CONST;
	}
	else return 1;
}

void sort_array_int(array_int *a)
{
	char status = is_sorted_array_int(a);
	if (status == 0)
	{
		int *e = (int *) calloc((*a).mem, sizeof(int));
		merge_sort_array_int(a, e, 0, (*a).len);
		free((*a).e);
		(*a).e = e;
	}
	else if (status == SORTED_DECREASING) // RESUME: #define this!
	{
		int hlength = (*a).len >> 1;
		int ii = (*a).len;
		for (int i = 0; i < hlength; i++)
		{
			ii -= 1;
			int aux = (*a).e[i];
			(*a).e[i] = (*a).e[ii];
			(*a).e[ii] = aux;
		}
	}
}

// methods for variable length arrays of variable length arrays
// aarray_int #beginning
void aarray_int_init(aarray_int *aa, int mem)
{
	aarray_int_init_precise(aa, mem, 0);
}

// An alternative to the last method that can be tuned to improve performance
void aarray_int_init_precise(aarray_int *aa, int mem1, int mem2)
{
	(*aa).mem = mem1;
	(*aa).e = (array_int *) calloc(mem1, sizeof(array_int));
	for (int i = 0; i < (*aa).mem; i++)
	{
		array_int_init(&(*aa).e[i], mem2);
	}
	(*aa).len = 0;
}

void transcribe_aarray_int(aarray_int *aa, aarray_int *aa_)
{
	(*aa_).mem = (*aa).mem;
	(*aa_).e = (array_int *) calloc((*aa).mem, sizeof(array_int));
	(*aa_).len = (*aa).len;
	for (int i = 0; i < (*aa).len; i++) 
	{
		transcribe_array_int(&((*aa).e[i]), &((*aa_).e[i]));
	}
}

void add_mem_aarray_int(aarray_int *aa)
{
	int init_mem = (*aa).mem;
	if ((*aa).mem > 0) {}
	else (*aa).mem = 1;
	(*aa).mem <<= 1;
	array_int *ne = (array_int *) calloc((*aa).mem, sizeof(array_int));
	int i = 0;
	while ((*aa).e[i].mem > 0 && i < init_mem)
	{
		ne[i] = (*aa).e[i];
		i += 1;
	}
	for (int ii = i; ii < (*aa).mem; ii++)
	{
		ne[ii].e = NULL;
		ne[ii].mem = 0;
		ne[ii].len = 0;
	}
	if ((*aa).e != NULL) free((*aa).e);
	(*aa).e = ne;
}

void contract_aarray_int(aarray_int *a)
{
	if ((*a).mem > 0 && (*a).e != NULL)
	{
		int newsize = (*a).mem >> 1;
		array_int *data = (array_int *) calloc(newsize, sizeof(array_int));
		int ulim = (*a).len < newsize ? (*a).len : newsize;
		for (int i = 0; i < ulim; i++) 
		{
			data[i] = (*a).e[i];
		}
		for (int i = ulim; i < newsize; i++)
		{
			data[i] = (*a).e[i];
		}
		free((*a).e);
		(*a).e = data;
		(*a).mem = newsize;
	}
}

void add_mem_aarray_int_until(aarray_int *aa, int i)
{
	if ((*aa).mem > i) return;
	int init_mem = (*aa).mem;
	if ((*aa).mem > 0) {}
	else (*aa).mem = 1;
	while ((*aa).mem <= i) (*aa).mem <<= 1;
	array_int *ne = (array_int *) calloc((*aa).mem, sizeof(array_int));
	int ii = 0; 
	while (ii < init_mem && (*aa).e[ii].mem > 0)
	{
		ne[ii] = (*aa).e[ii];
		ii += 1;
	}
	for (int i = ii; i < (*aa).mem; i++)
	{
		ne[i].e = NULL;
		ne[i].mem = 0;
		ne[i].len = 0;
	}
	if ((*aa).e != NULL) free((*aa).e);
	(*aa).e = ne;
}

void extend_aarray_int(aarray_int *aa)
{
	extend_aarray_int_precise(aa, INIT_A_MEM);
}

void extend_aarray_int_precise(aarray_int *aa, int init_mem)
{
	if ((*aa).len == (*aa).mem)
	{
		add_mem_aarray_int(aa);
	}
	if ((*aa).e[(*aa).len].e != NULL) 
	{
		(*aa).e[(*aa).len].len = 0;
	}
	else
	{
		array_int_init(&((*aa).e[(*aa).len]), init_mem);
	}
	(*aa).len += 1;

}

// RESUME: Check this! 001
void add2aarray_int(aarray_int *aa, array_int a)
{
	if (a.mem > 0) {}
	else 
	{
		printf("Error: attempting to add empty array_int to array of array_int (aarray_int.)\n");
		exit(EXIT_FAILURE);
	}
	if ((*aa).len == (*aa).mem)
	{
		add_mem_aarray_int(aa);
	}
	if ((*aa).e[(*aa).len].mem == 0) {}
	else 
	{
		// NOTE: alternatively, consider swapping the last element of aa with the next NULL element in 
		// 	memory, which might improve performance in some settings.
		free_array_int(&(*aa).e[(*aa).len]);
	}
	(*aa).e[(*aa).len] = a;
	(*aa).len += 1;
}

void add2aarray_int_elem(aarray_int *aa, int ai, int i)
{
	add_mem_aarray_int_until(aa, ai);
	add2array_int(&((*aa).e[ai]), i);
	(*aa).len = (*aa).len > ai? (*aa).len: ai + 1;
}



void reset_aarray_int_elem(int i, aarray_int *aa)
{
	reset_array_int(&((*aa).e[i]));
	free((*aa).e[i].e);
	(*aa).e[i].mem = 0;
	(*aa).e[i].e = NULL;
}

// RESUME: test this
void remove_aarray_int(aarray_int *aa, int i)
{
	if (i < (*aa).len && i > -1)
	{
		(*aa).len -= 1;
		// transcribe the last element of (*aa) onto the i-th element
		array_int aux = (*aa).e[i];
		(*aa).e[i] = (*aa).e[(*aa).len];
		(*aa).e[(*aa).len] = aux;
		reset_aarray_int_elem((*aa).len, aa);
		if ((*aa).len > (*aa).mem >> 2) {}
		else contract_aarray_int(aa);
	}
	else 
	{
		printf("Error: attempting to remove non-existent aarray_int element %d of %d\n", i, (*aa).len);
		exit(EXIT_FAILURE);
	}
}

void fprintf_aarray_int(aarray_int *aa, FILE *ofile)
{
  if (ofile != NULL)
    {
      for (int i = 0; i < (*aa).len; i++)
	{
	  fprintf(ofile, "%d ", (*aa).e[i].len);
	  for (int ii = 0; ii < (*aa).e[i].len; ii++) fprintf(ofile, "%d ", (*aa).e[i].e[ii]);
	  fprintf(ofile, "\n");
	}
    }
}

// Assumes aarray_int has been initialized
void load_aarray_int(aarray_int *aa, char *fname)
{
  FILE *ifile = fopen(fname, "r");
  if (ifile != NULL)
    {
      while (1)
	{
	  int n_elem;
	  int status = fscanf(ifile, "%d", &n_elem);
	  if (status != EOF)
	    {
	      int i = (*aa).len;
	      extend_aarray_int(aa);
	      for (int ni = 0; ni < n_elem; ni++)
		{
		  int ii;
		  fscanf(ifile, "%d", &ii);
		  add2array_int(&((*aa).e[i]), ii);
		}
	    }
	  else break;
	}
      fclose(ifile);
    }
}

void free_aarray_int(aarray_int *aa)
{
	if ((*aa).mem > 0)
	{
		int i = 0;
		while (i < (*aa).len)
		{
			free_array_int(&((*aa).e[i]));
			// free(&((*aa).e[i])); // RESUME: Check this! 000
			i += 1;
		}
		free((*aa).e);
	}
}

int aarray_int_max(aarray_int *a)
{
	int max = -RAND_MAX;
	for (int i = 0; i < (*a).len; i++)
	{
		for (int ni = 0; ni < (*a).e[i].len; ni++)
		{
			if ((*a).e[i].e[ni] < max) {}
			else max = (*a).e[i].e[ni];
		}
	}
	return max;
}

void print_aarray_int(aarray_int aa)
{
	for (int i = 0; i < aa.len; i++)
	{
		for (int ii = 0; ii < aa.e[i].len; ii++)
		{
			printf("%d ", aa.e[i].e[ii]);
		}
		printf("\n");
	}
}

// aarray_int #ending

// aarray_double #beginning
void aarray_double_init(aarray_double *aa, int mem)
{
	aarray_double_init_precise(aa, mem, 0);
}

void transcribe_aarray_double(aarray_double *src, aarray_double *dest)
{
	(*dest).mem = (*src).mem;
	(*dest).e = (array_double *) calloc((*dest).mem, sizeof(array_double));
	(*dest).len = (*src).len;
	for (int i = 0; i < (*src).len; i++)
	{
		transcribe_array_double(&((*src).e[i]), &((*dest).e[i]));
	}
}

// An alternative to the last method that can be tuned to improve performance
void aarray_double_init_precise(aarray_double *aa, int mem1, int mem2)
{
	(*aa).mem = mem1;
	(*aa).e = (array_double *) calloc(mem1, sizeof(array_double));
	for (int i = 0; i < (*aa).mem; i++)
	{
		array_double_init(&(*aa).e[i], mem2);
	}
	(*aa).len = 0;
}

void add_mem_aarray_double(aarray_double *aa)
{
	int init_mem = (*aa).mem;
	if ((*aa).mem > 0) {}
	else (*aa).mem = 1;
	(*aa).mem <<= 1;
	array_double *ne = (array_double *) calloc((*aa).mem, sizeof(array_double));
	int i = 0;
	while ((*aa).e[i].mem > 0 && i < init_mem)
	{
		ne[i] = (*aa).e[i];
		i += 1;
	}
	for (int ii = i; ii < (*aa).mem; ii++)
	{
		ne[ii].e = NULL;
		ne[ii].mem = 0;
		ne[ii].len = 0;
	}
	if ((*aa).e != NULL) free((*aa).e);
	(*aa).e = ne;
}

void contract_aarray_double(aarray_double *a)
{
	if ((*a).mem > 0 && (*a).e != NULL)
	{
		int newsize = (*a).mem >> 1;
		array_double *data = (array_double *) calloc(newsize, sizeof(array_double));
		int ulim = (*a).len < newsize ? (*a).len : newsize;
		for (int i = 0; i < ulim; i++) 
		{
			data[i] = (*a).e[i];
		}
		for (int i = ulim; i < newsize; i++)
		{
			data[i] = (*a).e[i];
		}
		free((*a).e);
		(*a).e = data;
		(*a).mem = newsize;
	}
}

void add_mem_aarray_double_until(aarray_double *aa, int i)
{
	if ((*aa).mem > i) return;
	int init_mem = (*aa).mem;
	if ((*aa).mem > 0) {}
	else (*aa).mem = 1;
	while ((*aa).mem <= i) (*aa).mem <<= 1;
	array_double *ne = (array_double *) calloc((*aa).mem, sizeof(array_double));
	int ii = 0; 
	while (ii < init_mem && (*aa).e[ii].mem > 0)
	{
		ne[ii] = (*aa).e[ii];
		ii += 1;
	}
	for (int i = ii; i < (*aa).mem; i++)
	{
		ne[i].e = NULL;
		ne[i].mem = 0;
		ne[i].len = 0;
	}
	if ((*aa).e != NULL) free((*aa).e);
	(*aa).e = ne;
}

void extend_aarray_double(aarray_double *aa)
{
	extend_aarray_double_precise(aa, INIT_A_MEM);
}

void extend_aarray_double_precise(aarray_double *aa, int init_mem)
{
	if ((*aa).len == (*aa).mem)
	{
		add_mem_aarray_double(aa);
	}
	if ((*aa).e[(*aa).len].e != NULL) 
	{
		(*aa).e[(*aa).len].len = 0;
	}
	else
	{
		array_double_init(&((*aa).e[(*aa).len]), init_mem);
	}
	(*aa).len += 1;

}

// RESUME: Check this! 001
void add2aarray_double(aarray_double *aa, array_double a)
{
	if (a.mem > 0) {}
	else 
	{
		printf("Error: attempting to add empty array_double to array of array_double (aarray_double.)\n");
		exit(EXIT_FAILURE);
	}
	if ((*aa).len == (*aa).mem)
	{
		add_mem_aarray_double(aa);
	}
	if ((*aa).e[(*aa).len].mem == 0) {}
	else 
	{
		// NOTE: alternatively, consider swapping the last element of aa with the next NULL element in 
		// 	memory, which might improve performance in some settings.
		free_array_double(&(*aa).e[(*aa).len]);
	}
	(*aa).e[(*aa).len] = a;
	(*aa).len += 1;
}

void reset_aarray_double_elem(int i, aarray_double *aa)
{
	reset_array_double(&((*aa).e[i]));
}

// RESUME: test this
void remove_aarray_double(aarray_double *aa, int i)
{
	if (i < (*aa).len && i > -1)
	{
		(*aa).len -= 1;
		// transcribe the last element of (*aa) onto the i-th element
		array_double aux = (*aa).e[i];
		(*aa).e[i] = (*aa).e[(*aa).len];
		(*aa).e[(*aa).len] = aux;
		reset_aarray_double_elem((*aa).len, aa);
		if ((*aa).len > (*aa).mem >> 2) {}
		else contract_aarray_double(aa);
	}
	else 
	{
		printf("Error: attempting to remove non-existent aarray_double element %d of %d\n", i, (*aa).len);
		exit(EXIT_FAILURE);
	}
}


void fprintf_aarray_double(aarray_double *aa, FILE *ofile)
{
    if (ofile != NULL)
    {
      for (int i = 0; i < (*aa).len; i++)
	{
	  fprintf(ofile, "%d ", (*aa).e[i].len);
	  for (int ii = 0; ii < (*aa).e[i].len; ii++) fprintf(ofile, "%g ", (*aa).e[i].e[ii]);
	  fprintf(ofile, "\n");
	}
    }
}

void load_aarray_double(aarray_double *aa, char *fname)
{
  FILE *ifile = fopen(fname, "r");
  if (ifile != NULL)
    {
      while (1)
	{
	  int n_elem;
	  int status = fscanf(ifile, "%d", &n_elem);
	  if (status != EOF)
	    {
	      int i = (*aa).len;
	      extend_aarray_double(aa);
	      for (int ni = 0; ni < n_elem; ni++)
		{
		  double ii;
		  fscanf(ifile, "%lg", &ii);
		  add2array_double(&((*aa).e[i]), (double) ii);
		}
	    }
	  else break;
	}
      fclose(ifile);
    }
}

void free_aarray_double(aarray_double *aa)
{
	if ((*aa).mem > 0)
	{
		int i = 0;
		while (i < (*aa).mem && (*aa).e[i].mem > 0)
		{
			free_array_double(&((*aa).e[i]));
			// free(&((*aa).e[i])); // RESUME: Check this! 000
			i += 1;
		}
		free((*aa).e);
	}
}

void print_aarray_double(aarray_double aa)
{
	for (int i = 0; i < aa.len; i++)
	{
		for (int ii = 0; ii < aa.e[i].len; ii++)
		{
			printf("%g ", aa.e[i].e[ii]);
		}
		printf("\n");
	}
}


// aarray_double #ending

// Methods for linked lists
linked_list linked_list_init()
{
	linked_list nll;
	nll.next = NULL;
	nll.data = NULL;
	return nll;
}

void add2linked_list(linked_list *ll, void *elem)
{
	linked_list *new = (linked_list *) calloc(1, sizeof(linked_list));
	(*new).data = (*ll).data;
	(*new).next = (*ll).next;
	(*ll).data = elem;
	(*ll).next = new;
}

void insert_linked_list(linked_list *ll, void *elem)
{
	linked_list *new = (linked_list *) calloc(1, sizeof(linked_list));
	(*new).data = elem;
	(*new).next = (*ll).next;
	(*ll).next = new;
}

void* pop_linked_list(linked_list *ll)
{
	void *popped = (*ll).data;
	linked_list *next_ll = (*ll).next;
	(*ll) = (*next_ll);
	free(next_ll);
	return popped;
}

void free_linked_list(linked_list *ll, void (*ff)(void *))
{
	while ((*ll).next != NULL) 
	{
		void *data = pop_linked_list(ll);
		if (ff != NULL) ff(data);
	}
}

dlinked_list dlinked_list_init()
{
	dlinked_list dll;
	dll.next = NULL;
	dll.prev = NULL;
	dll.data = NULL;
}

void add2dlinked_list(dlinked_list *dll, void *elem)
{
	dlinked_list *ndll = (dlinked_list *) calloc(1, sizeof(dlinked_list));
	(*ndll).next = (*dll).next;
	(*ndll).prev = dll;
	(*ndll).data = (*dll).data;
	(*dll).next = ndll;
	(*dll).data = elem;
}

void insert_dlinked_list(dlinked_list *dll, void *elem)
{
	dlinked_list *ndll = (dlinked_list *) calloc(1, sizeof(dlinked_list));
	(*ndll).next = (*dll).next;
	(*ndll).prev = dll;
	(*ndll).data = elem;
	(*dll).next = ndll;
}

void *pop_dlinked_list(dlinked_list *dll)
{
	void *elem = (*dll).data;
	(*dll).data = (*(*dll).next).data;
	(*dll).next = (*(*dll).next).next;
	(*(*(*dll).next).next).prev = dll;
	free((*dll).next);
}

void free_dlinked_list(dlinked_list *dll, void (*ff)(void *))
{
	while ((*dll).next != NULL)
	{
		void *data = pop_dlinked_list(dll);
		if (ff != NULL) ff(data);
	}
}

// methods for variable dimension box lists
void move_boxlist_elem(boxlist *bl, int elem, int new_box_index)
{
	int old_box_index = (*bl).addr.e[elem].e[0];
	if (new_box_index != old_box_index) {}
	else return;
	int eci = (*bl).addr.e[elem].e[1];
	remove_array_int(&(*bl).boxc.e[old_box_index], eci);
	// int lci = (*bl).boxc.e[old_box_index].len;
	int l_elem = (*bl).boxc.e[old_box_index].e[eci];
	(*bl).addr.e[l_elem].e[1] = eci;
	(*bl).addr.e[elem].e[0] = new_box_index;
	(*bl).addr.e[elem].e[1] = (*bl).boxc.e[new_box_index].len;
	add2array_int(&(*bl).boxc.e[new_box_index], elem);
}

void reset_boxlist(boxlist *bl)
{
	(*bl).addr.len = 0;
	for (int i = 0; i < (*bl).boxc.len; i++)
	{
		(*bl).boxc.e[i].len = 0;
	}
}

int flatten_vdim(int dim, int *m, int *box_index)
{
	int j = dim - 1;
	int fi = box_index[j];
	while (j > 0)
	{
		j -= 1;
		fi *= m[j];
		fi += box_index[j];
	}
	return fi;
}

void unflatten_box_index(int dim, int *m, int *flatbi, int *unflatbi)
{
	int last = dim - 1;
	int aux = (*flatbi);
	for (int i = 0; i < last; i++)
	{
		unflatbi[last] = aux / m[i];
		unflatbi[i] = aux % m[i];
		aux = unflatbi[last];
	}
}

void unflatten_box_index3D(int *m, int *flatbi, int *unflatbi)
{
	unflatbi[2] = (*flatbi) / m[0];
	unflatbi[0] = (*flatbi) % m[0];
	int aux = unflatbi[2];
	unflatbi[2] = aux / m[1];
	unflatbi[1] = aux % m[1];
}

void add2box(boxlist *bl, int elem, int *box_index)
{
	// Consider using a more parsimonious data structure for the boxlist address elements
	//	(because they all have equal length)
	array_int new_addr;
	array_int_init(&new_addr, 3);
	int fi = flatten_vdim((*bl).dim, (*bl).m, box_index);
	new_addr.len = 3;
	new_addr.e[0] = fi;
	new_addr.e[1] = (*bl).boxc.e[fi].len;
	new_addr.e[2] = elem;
	add2aarray_int_elem(&((*bl).boxc), fi, (*bl).addr.len);
	add2aarray_int(&((*bl).addr), new_addr); 
}

void add2box_preflattened(boxlist *bl, int elem, int fbox_index)
{
	array_int new_addr;
	array_int_init(&new_addr, 3);
	new_addr.len = 3;
	new_addr.e[0] = fbox_index;
	new_addr.e[1] = (*bl).boxc.e[fbox_index].len;
	new_addr.e[2] = elem;
	add2aarray_int_elem(&((*bl).boxc), fbox_index, (*bl).addr.len);
	add2aarray_int(&((*bl).addr), new_addr);
}

int boxlist_box_size(boxlist *bl, int *box_index)
{
	int fi = flatten_vdim((*bl).dim, (*bl).m, box_index);
	return (*bl).boxc.e[fi].len;
}

char check_boxlist_consistency(boxlist *bl)
{
	for (int i = 0; i < (*bl).addr.len; i++)
	{
		int fi = (*bl).addr.e[i].e[0];
		int ci = (*bl).addr.e[i].e[1];
		if ((*bl).boxc.e[fi].e[ci] == i) {}
		else
		{
			return 0;
		}
	}
	for (int i = 0; i < (*bl).boxc.len; i++)
	{
		for (int ci = 0; ci < (*bl).boxc.e[i].len; ci++)
		{
			int ii = (*bl).boxc.e[i].e[ci];
			int cii = (*bl).addr.e[ii].e[1];
			int fii = (*bl).addr.e[ii].e[0];
			if (cii == ci && fii == i) {}
			else
			{
				return 0;
			}
		}
	}
	return 1;
}

void remove_boxlist_elem_by_addr_pf(boxlist *bl, int fi, int content_index)
{
	if ((*bl).boxc.e[fi].len > content_index) {}
	else
	{
		printf("Warning: attempting to remove nonexistent entry %d from a boxlist containing %d elements\n", content_index, (*bl).boxc.e[fi].len);
		return;
	}
	int addri = (*bl).boxc.e[fi].e[content_index];
	// Remove (*bl).addr.e[addri] by copying data from the last entry
	remove_aarray_int(&(*bl).addr, addri);
	// 	(Update the address that appears in the box specified by addr.e[addri])
	if ((*bl).addr.len > 0)
	{
		int bi = (*bl).addr.e[addri].e[0];
		int ci = (*bl).addr.e[addri].e[1];
		int id_ = (*bl).addr.e[addri].e[2];
		(*bl).boxc.e[bi].e[ci] = addri;
	}
	// Remove the element from the associated box
	remove_array_int(&((*bl).boxc.e[fi]), content_index);
	if ((*bl).boxc.e[fi].len > 0 && content_index < (*bl).boxc.e[fi].len)
	{
		int last_addr = (*bl).boxc.e[fi].e[content_index];
		(*bl).addr.e[last_addr].e[1] = content_index;
	}
}

void remove_boxlist_elem(boxlist *bl, int i)
{
	array_int addr_i = (*bl).addr.e[i];
	int fi = addr_i.e[0];
	int ci = addr_i.e[1];
	remove_boxlist_elem_by_addr_pf(bl, fi, ci);
}

// this function has an effect on both addr and boxc
//	- remove an element from the box list
void remove_boxlist_elem_by_addr(boxlist *bl, int *box_index, int content_index)
{
	int fi = flatten_vdim((*bl).dim, (*bl).m, box_index);
	remove_boxlist_elem_by_addr_pf(bl, fi, content_index);
	/*if (check_boxlist_consistency(bl)) {}
	else
	{
		printf("Error: boxlist failed consistency check trying to remove element %d from box %d %d of length %d\n", content_index, box_index[0], box_index[1], (*bl).boxc.e[fi].len);
		exit(EXIT_FAILURE);
	}*/

}

void free_boxlist(boxlist *bl)
{
	free_aarray_int(&((*bl).boxc));
	free_aarray_int(&((*bl).addr));
	free((*bl).m);
}

void boxlist_init(boxlist *bl, int *m, int dim)
{
	(*bl).dim = dim;
	(*bl).m = (int *) calloc(dim, sizeof(int));
	for (int i = 0; i < dim; i++) (*bl).m[i] = m[i];
	if (dim > 0)
	{
		int Nboxes = (*bl).m[0];
		for (int i = 1; i < dim; i++)
		{
			Nboxes *= (*bl).m[i];
		}
		aarray_int_init_precise(&((*bl).addr), 1, 3);
		aarray_int_init_precise(&((*bl).boxc), Nboxes, 1);
		(*bl).boxc.len = Nboxes;
	}
}

void init_counter(int *counter, int dim)
{
	for (int i = 0; i < dim; i++) counter[i] = 0;
}

void advance_counter(int *counter, int *m, int dim)
{
	int i = 0;
	while (i < dim)
	{
		counter[i] += 1;
		if (counter[i] < m[i])
		{
			break;
		}
		else
		{
			counter[i] = 0;
			i += 1;
		}
	}
}

void print_boxlist_counts(boxlist *bl)
{
	int counter[(*bl).dim];
	init_counter(counter, (*bl).dim);
	for (int i = 0; i < (*bl).boxc.len; i++)
	{
		printf("%d ", (*bl).boxc.e[i].len);
		advance_counter(counter, (*bl).m, (*bl).dim);
		for (int ii = 0; ii < (*bl).dim; ii++)
		{
			if (counter[ii] == 0) printf("\n");
			else break;
		}
	}
}

void print_boxlist(boxlist *bl)
{
	int counter[(*bl).dim];
	init_counter(counter, (*bl).dim);
	for (int i = 0; i < (*bl).boxc.len; i++)
	{
		printf("Box (");
		for (int ii = 0; ii < (*bl).dim; ii++)
		{
			printf(" %d ", counter[i]);
		}
		printf("): ");
		for (int ii = 0; ii < (*bl).boxc.e[i].len; ii++)
		{
			printf("%d ", (*bl).boxc.e[i].e[ii]);
		}
		printf("\n");
		advance_counter(counter, (*bl).m, (*bl).dim);
	}
}

// Methods for 3D box lists
void move_boxlist3D_elem(boxlist3D *bl, int elem, int ijk)
{
	int eci = (*bl).addr.e[elem].e[1];
	int old_ijk = (*bl).addr.e[elem].e[0];
	remove_array_int(&(*bl).boxc.e[old_ijk], eci);
	int l_elem = (*bl).boxc.e[old_ijk].e[eci];
	(*bl).addr.e[l_elem].e[1] = eci;
	(*bl).addr.e[elem].e[0] = ijk;
	(*bl).addr.e[elem].e[1] = (*bl).boxc.e[ijk].len;
	add2array_int(&(*bl).boxc.e[ijk], eci);
}

void reset_boxlist3D(boxlist3D *bl)
{
	(*bl).addr.len = 0;
	for (int i = 0; i < (*bl).boxc.len; i++)
	{
		(*bl).boxc.e[i].len = 0;
	}
}

void add2box3D(boxlist3D *bl, int elem, int i, int j, int k)
{
	int ijk = i + (*bl).m[0] * (j + (*bl).m[1] * k);
	array_int new_addr;
	array_int_init(&new_addr, 3);
	new_addr.len = 3;
	new_addr.e[0] = ijk;
	new_addr.e[1] = (*bl).boxc.e[ijk].len;
	new_addr.e[2] = elem;
	add2aarray_int_elem(&((*bl).boxc), ijk, (*bl).addr.len);
	add2aarray_int(&((*bl).addr), new_addr);
}

void remove_boxlist3D_elem(boxlist3D *bl, int fi, int content_index)
{
	if ((*bl).boxc.e[fi].len > content_index) {}
	else
	{
		printf("Warning: attempting to remove nonexistent entry %d from a boxlist containing %d elements\n", content_index, (*bl).boxc.e[fi].len);
		return;
	}
	int addri = (*bl).boxc.e[fi].e[content_index];
	// Remove (*bl).addr.e[addri] by copying data from the last entry
	(*bl).addr.len -= 1;
	transcribe_array_int(&((*bl).addr.e[addri]), &((*bl).addr.e[(*bl).addr.len]));
	// 	(Update the address that appears in the box specified by addr.e[addri])
	int bi = (*bl).addr.e[addri].e[0];
	int ci = (*bl).addr.e[addri].e[1];
	int id_ = (*bl).addr.e[addri].e[2];
	(*bl).boxc.e[bi].e[ci] = addri;
	// Remove the element from the associated box
	remove_array_int(&((*bl).boxc.e[fi]), content_index);
	if ((*bl).boxc.e[fi].len > 0)
	{
		int last_addr = (*bl).boxc.e[fi].e[content_index];
		(*bl).addr.e[last_addr].e[1] = content_index;
	}
}

void remove_boxlist3D_elem_unflat(boxlist3D *bl, int i, int j, int k, int content_index)
{
	int fi = i + (*bl).m[0] * (j + (*bl).m[1] * k);
	remove_boxlist3D_elem(bl, fi, content_index);
}

void free_boxlist3D(boxlist3D *bl)
{
	free_aarray_int(&((*bl).addr));
	free_aarray_int(&((*bl).boxc));
}

boxlist3D boxlist3D_init(int m, int n, int o)
{
	boxlist3D bl;
	bl.m[0] = m;
	bl.m[1] = n;
	bl.m[2] = o;
	int Nboxes = m * n * o;
	aarray_int_init_precise(&(bl.addr), 1, 3);
	aarray_int_init_precise(&(bl.boxc), Nboxes, 1);
	bl.boxc.len = Nboxes;
	return bl;
}

// Methods for 2D box lists
void move_boxlist2D_elem(boxlist2D *bl, int elem, int ij)
{
	int eci = (*bl).addr.e[elem].e[1];
	int old_ij = (*bl).addr.e[elem].e[0];
	remove_array_int(&(*bl).boxc.e[old_ij], eci);
	int l_elem = (*bl).boxc.e[old_ij].e[eci];
	(*bl).addr.e[l_elem].e[1] = eci;
	(*bl).addr.e[elem].e[0] = ij;
	(*bl).addr.e[elem].e[1] = (*bl).boxc.e[ij].len;
	add2array_int(&(*bl).boxc.e[ij], eci);
}

void reset_boxlist2D(boxlist2D *bl)
{
	(*bl).addr.len = 0;
	for (int i = 0; i < (*bl).boxc.len; i++)
	{
		(*bl).boxc.e[i].len = 0;
	}
}

void add2box2D(boxlist2D *bl, int elem, int i, int j)
{
	int ij = i + (*bl).m[0] * j;
	array_int new_addr;
	array_int_init(&new_addr, 3);
	new_addr.len = 3;
	new_addr.e[0] = ij;
	new_addr.e[1] = (*bl).boxc.e[ij].len;
	new_addr.e[2] = elem;
	add2aarray_int_elem(&((*bl).boxc), ij, (*bl).addr.len);
	add2aarray_int(&((*bl).addr), new_addr);
}

void remove_boxlist2D_elem(boxlist2D *bl, int fi, int content_index)
{
	if ((*bl).boxc.e[fi].len > content_index) {}
	else
	{
		printf("Warning: attempting to remove nonexistent entry %d from a boxlist containing %d elements\n", content_index, (*bl).boxc.e[fi].len);
		return;
	}
	int addri = (*bl).boxc.e[fi].e[content_index];
	// Remove (*bl).addr.e[addri] by copying data from the last entry
	(*bl).addr.len -= 1;
	transcribe_array_int(&((*bl).addr.e[addri]), &((*bl).addr.e[(*bl).addr.len]));
	// 	(Update the address that appears in the box specified by addr.e[addri])
	int bi = (*bl).addr.e[addri].e[0];
	int ci = (*bl).addr.e[addri].e[1];
	int id_ = (*bl).addr.e[addri].e[2];
	(*bl).boxc.e[bi].e[ci] = addri;
	// Remove the element from the associated box
	remove_array_int(&((*bl).boxc.e[fi]), content_index);
	if ((*bl).boxc.e[fi].len > 0)
	{
		int last_addr = (*bl).boxc.e[fi].e[content_index];
		(*bl).addr.e[last_addr].e[1] = content_index;
	}
}

void remove_boxlist2D_elem_unflat(boxlist2D *bl, int i, int j, int content_index)
{
	int fi = i + (*bl).m[0] * j;
	remove_boxlist2D_elem(bl, fi, content_index);
}

void free_boxlist2D(boxlist2D *bl)
{
	free_aarray_int(&((*bl).addr));
	free_aarray_int(&((*bl).boxc));
}

boxlist2D boxlist2D_init(int m, int n)
{
	boxlist2D bl;
	bl.m[0] = m;
	bl.m[1] = n;
	int Nboxes = m * n;
	aarray_int_init_precise(&(bl.addr), 1, 3);
	aarray_int_init_precise(&(bl.boxc), Nboxes, 1);
	bl.boxc.len = Nboxes;
	return bl;
}



// Methods for neighbor lists
// RESUME: update all instances of nbrlists being initialized to be consistent with the new 
//		convention.
void transcribe_nbrlist(nbrlist *nb1, nbrlist *nb2)
{
	transcribe_aarray_int(&((*nb1).v), &((*nb2).v));
	transcribe_aarray_int(&((*nb1).i_of), &((*nb2).i_of));
}

int nbrlist_i_of(nbrlist *nbl, int i, int j)
{
	int i_of_ij = -1;
	for (int ni = 0; ni < (*nbl).v.e[j].len; ni++)
	{
		int ii = (*nbl).v.e[j].e[ni];
		if (ii == i) 
		{
			i_of_ij = ni;
			break;
		}
	}
	return i_of_ij;
}

char nbrlist_has_edge(nbrlist *nbl, int i, int j)
{
	if ((*nbl).v.e[i].len <= (*nbl).v.e[j].len) {}
	else return nbrlist_has_edge(nbl, j, i);
	for (int ni = 0; ni < (*nbl).v.e[i].len; ni++)
	{
		int ii = (*nbl).v.e[i].e[ni];
		if (ii == j) return 1;
	}
	return 0;
}

void nbrlist_init_precise(nbrlist *nbl, int mem)
{
	if (nbl == NULL)
	{
		nbl = (nbrlist *) calloc(1, sizeof(nbrlist));
	}
	aarray_int_init_precise(&(*nbl).v, mem, 1);
	aarray_int_init_precise(&(*nbl).i_of, mem, 1);
}

void nbrlist_init(nbrlist *nbl)
{
	nbrlist_init_precise(nbl, 1);
}

void prep_nbrlist(nbrlist *nbl)
{
	if ((*nbl).v.mem > (*nbl).v.len) {}
	else
	{
		add_mem_aarray_int(&((*nbl).v));
		add_mem_aarray_int(&((*nbl).i_of));
	}
}

void check_nbrlist(nbrlist *nbl)
{
	if ((*nbl).v.mem >= (*nbl).v.len) {}
	else
	{
		add_mem_aarray_int_until(&((*nbl).v), (*nbl).v.len);
		add_mem_aarray_int_until(&((*nbl).i_of), (*nbl).v.len);
	}

}

void identify_v_nbrlist(nbrlist *nbl, int i, int j)
{
	// Determine the neighbors of both i and j
	for (int nj = 0; nj < (*nbl).v.e[j].len; nj++)
	{
		int jj = (*nbl).v.e[j].e[nj];
		char adj_i = 0;
		for (int ni = 0; ni < (*nbl).v.e[i].len; ni++)
		{
			int ii = (*nbl).v.e[i].e[ni];
			if (ii == jj)
			{
				adj_i = 1;
				break;
			}
		}
		if (adj_i) {}
		else
		{
			add_edge_nbrlist(nbl, i, jj);
		}
	}
	remove_vertex_nbrlist(nbl, j);
}

void extend_nbrlist(nbrlist *nbl)
{
	prep_nbrlist(nbl);
	(*nbl).v.len += 1;
	(*nbl).i_of.len += 1;
}

void set_len_nbrlist(nbrlist *nbl, int len)
{
	(*nbl).v.len = len;
	(*nbl).i_of.len = len;
	check_nbrlist(nbl);
}

void extend_nbrlist_n(nbrlist *nbl, int N)
{
	int init_len = (*nbl).v.len;
	(*nbl).v.len += N;
	(*nbl).i_of.len += N;
	check_nbrlist(nbl);
	/*
	(*nbl).v.len -= N;
	for (int i = 0; i < N; i++)
	{
		array_int nbrs;
		array_int i_of_;
		array_int_init(&nbrs, 1);
		array_int_init(&i_of_, 1);
		add2aarray_int(&((*nbl).v), nbrs);
		add2aarray_int(&((*nbl).i_of), i_of_);

	}
	*/
}

void ensure_nbrlist_size_n(nbrlist *nbl, int N)
{
	if ((*nbl).v.len >= N) {}
	else
	{
		extend_nbrlist_n(nbl, N - (*nbl).v.len);
	}
}

int add_edge_nbrlist_safe(nbrlist *nbl, int vertex1, int vertex2)
{
	int vsl[2] = {vertex1, vertex2};
	if ((*nbl).v.e[vertex1].len < (*nbl).v.e[vertex2].len) {}
	else
	{
		vsl[0] = vertex2;
		vsl[1] = vertex1;
	}
	for (int i = 0; i < (*nbl).v.e[vsl[0]].len; i++)
	{
		if ((*nbl).v.e[vsl[0]].e[i] != vsl[1]) {}
		else
		{
			if (vsl[0] == vertex1) return i;
			else
			{
				return (*nbl).i_of.e[vsl[0]].e[i];
			}
		}
	}
	add_edge_nbrlist(nbl, vertex1, vertex2);
	return -1;
}

void add_edge_nbrlist(nbrlist *nbl, int vertex1, int vertex2)
{
	add2aarray_int_elem(&((*nbl).i_of), vertex1, (*nbl).v.e[vertex2].len);
	add2aarray_int_elem(&((*nbl).i_of), vertex2, (*nbl).v.e[vertex1].len);
	add2aarray_int_elem(&((*nbl).v), vertex1, vertex2);
	add2aarray_int_elem(&((*nbl).v), vertex2, vertex1);
}

void add_edge_nbrlist_1way(nbrlist *nbl, int vertex1, int vertex2)
{
	int back_edge = -1;
	for (int i = 0; i < (*nbl).v.e[vertex2].len; i++)
	{
		if ((*nbl).v.e[vertex2].e[i] == vertex1)
		{
			back_edge = i;
			(*nbl).i_of.e[vertex2].e[back_edge] = (*nbl).v.e[vertex1].len;
			break;
		}
	}
	add2aarray_int_elem(&(*nbl).i_of, vertex1, back_edge);
	add2aarray_int_elem(&(*nbl).v, vertex1, vertex2);

}

void remove_edge_nbrlist_1way(nbrlist *nbl, int vertex, int local_nbr_index)
{
	int nvertex = (*nbl).v.e[vertex].e[local_nbr_index];
	int nlocal_nbr_index = (*nbl).i_of.e[vertex].e[local_nbr_index];
	if (nlocal_nbr_index > -1)
	{
		(*nbl).i_of.e[nvertex].e[nlocal_nbr_index] = -1;
	}
	remove_array_int(&(*nbl).v.e[vertex], local_nbr_index);
	remove_array_int(&(*nbl).i_of.e[vertex], local_nbr_index);
	if ((*nbl).v.e[vertex].len > 0)
	{
		int _vertex = (*nbl).v.e[vertex].e[local_nbr_index];
		int _local_nbr_index = (*nbl).i_of.e[vertex].e[local_nbr_index];
		(*nbl).i_of.e[_vertex].e[_local_nbr_index] = local_nbr_index;
	}
}

void check_nbrlist_local_consistency(nbrlist *nbl, int v)
{
	for (int ni = 0; ni < (*nbl).v.e[v].len; ni++)
	{
		int ii = (*nbl).v.e[v].e[ni];
		int i_v_ii = (*nbl).i_of.e[v].e[ni];
		int i_ii_v = (*nbl).i_of.e[ii].e[i_v_ii];
		if (i_ii_v == ni) {}
		else
		{
			printf("Neighborlist failed consistency check at vertex %d (neighbor index = %d != %d)\n", v, ni, i_ii_v);
			exit(EXIT_FAILURE);
		}
	}
}

void remove_edge_nbrlist(nbrlist *nbl, int vertex, int local_nbr_index)
{
	if (local_nbr_index > -1 && vertex > -1) {}
	else
	{
		printf("Error (remove_edge_nbrlist): attempting to remove non-existent edge v.e[%d].e[%d]\n", vertex, local_nbr_index);
		exit(EXIT_FAILURE);
	}
	int vertex_ = (*nbl).v.e[vertex].e[local_nbr_index];
	int local_nbr_index_ = (*nbl).i_of.e[vertex].e[local_nbr_index];
	if (local_nbr_index_ > -1)
	{
		remove_array_int(&(*nbl).v.e[vertex_], local_nbr_index_);
		remove_array_int(&(*nbl).i_of.e[vertex_], local_nbr_index_);
		if ((*nbl).v.e[vertex_].len > local_nbr_index_)
		{
			int vertex__ = (*nbl).v.e[vertex_].e[local_nbr_index_];
			int local_nbr_index__ = (*nbl).i_of.e[vertex_].e[local_nbr_index_];
			(*nbl).i_of.e[vertex__].e[local_nbr_index__] = local_nbr_index_;
		}
	}
	remove_array_int(&(*nbl).v.e[vertex], local_nbr_index);
	remove_array_int(&(*nbl).i_of.e[vertex], local_nbr_index);
	if ((*nbl).v.e[vertex].len > local_nbr_index)
	{
		int _vertex = (*nbl).v.e[vertex].e[local_nbr_index];
		int _local_nbr_index = (*nbl).i_of.e[vertex].e[local_nbr_index];
		(*nbl).i_of.e[_vertex].e[_local_nbr_index] = local_nbr_index;
	}
}

void remove_vertex_nbrlist(nbrlist *nbl, int vertex)
{
	// Remove 'vertex' from the neighbor lists of its neighbors
	for (int i = 0; i < (*nbl).v.e[vertex].len; i++)
	{
		int vertex_ = (*nbl).v.e[vertex].e[i];
		int local_index_ = (*nbl).i_of.e[vertex].e[i];
		remove_edge_nbrlist_1way(nbl, vertex_, local_index_);
	}
	// Remove the vertex by transcribing the last entry in (*nbl).v and (*nbl).i_of:
	remove_aarray_int(&(*nbl).v, vertex);
	remove_aarray_int(&(*nbl).i_of, vertex);
	// Replace references to vertex (*nbl).v.len with references to 'vertex'
	for (int i = 0; i < (*nbl).v.e[vertex].len; i++)
	{
		int nv = (*nbl).v.e[vertex].e[i];
		int nv_i = (*nbl).i_of.e[vertex].e[i];
		(*nbl).v.e[nv].e[nv_i] = vertex;
	}
}

void remove_edges_vertex_nbrlist(nbrlist *nbl, int vertex)
{
	int i = (*nbl).v.e[vertex].len;
	while (i > 0)
	{
		i -= 1;
		remove_edge_nbrlist(nbl, vertex, i);
	} 
}

void remove_all_edges_nbrlist(nbrlist *nbl)
{
	for (int i = 0; i < (*nbl).v.len; i++)
	{
		(*nbl).v.e[i].len = 0;
		(*nbl).i_of.e[i].len = 0;
	}
}

void fprintf_nbrlist(nbrlist *nbl, FILE *ofile)
{
	if (ofile != NULL)
	{
		for (int i = 0; i < (*nbl).v.len; i++)
		{
			fprintf(ofile, "%d ", (*nbl).v.e[i].len);
			for (int ni = 0; ni < (*nbl).v.e[i].len; ni++)
			{
				fprintf(ofile, "%d ", (*nbl).v.e[i].e[ni]);
			}
			fprintf(ofile, "\n");
		}
	}
}

void write_nbrlist(nbrlist *nbl, char *ofname)
{
	FILE *ofile = fopen(ofname, "w");
	if (ofile != NULL)
	{
		for (int i = 0; i < (*nbl).v.len; i++)
		{
			fprintf(ofile, "%d ", (*nbl).v.e[i].len);
			for (int ni = 0; ni < (*nbl).v.e[i].len; ni++)
			{
				fprintf(ofile, "%d ", (*nbl).v.e[i].e[ni]);
			}
			fprintf(ofile, "\n");
		}
		fclose(ofile);
	}
}

void free_nbrlist(nbrlist *nbl)
{
	free_aarray_int(&((*nbl).v));
	free_aarray_int(&((*nbl).i_of));
}

int N_vertices_nbrlist(nbrlist *nbl)
{
	return (*nbl).v.len;
}

int N_nbors_nbrlist(nbrlist *nbl, int v_index)
{
	return (*nbl).v.e[v_index].len;
}

// Methods for edge_wtd_graph

void edge_wtd_graph_init(edge_wtd_graph *g, int size)
{
	nbrlist_init_precise(&((*g).top), size);
	aarray_int_init(&((*g).edge_wts), size);
}

void free_edge_wtd_graph(edge_wtd_graph *g)
{
	free_aarray_int(&((*g).edge_wts));
	free_nbrlist(&((*g).top));
}

void transcribe_edge_wtd_graph(edge_wtd_graph *src, edge_wtd_graph *dest)
{
	transcribe_nbrlist(&((*src).top), &((*dest).top));
	transcribe_aarray_int(&((*src).edge_wts), &((*dest).edge_wts));
}

void fprintf_edge_wtd_graph(edge_wtd_graph *g, char *prefix)
{
	char buf[256];
	sprintf(buf, "%s.top", prefix);
	FILE *ofile = fopen(buf, "w");
	if (ofile != NULL)
	{
		fprintf_nbrlist(&((*g).top), ofile);
		fclose(ofile);
	}
	sprintf(buf, "%s.wts", prefix);
	ofile = fopen(buf, "w");
	if (ofile != NULL)
	{
		fprintf_aarray_int(&((*g).edge_wts), ofile);
		fclose(ofile);
	}
}

int add_edge_edge_wtd_graph_safe(edge_wtd_graph *g, int i, int j, int w)
{
	for (int ni = 0; ni < (*g).top.v.e[i].len; ni++)
	{
		int ii = (*g).top.v.e[i].e[ni];
		if (ii == j)
		{
			return ni;
		}
	}
	add_edge_edge_wtd_graph(g, i, j, w);
	return -1;
}

void add_edge_edge_wtd_graph(edge_wtd_graph *g, int i, int j, int w)
{
	add_edge_nbrlist(&((*g).top), i, j);
	add2array_int(&((*g).edge_wts.e[i]), w);
	add2array_int(&((*g).edge_wts.e[j]), w);
}

void extend_edge_wtd_graph(edge_wtd_graph *g)
{
	extend_nbrlist(&((*g).top));
	extend_aarray_int(&((*g).edge_wts));
}

void remove_edges_vertex_edge_wtd_graph(edge_wtd_graph *g, int i)
{
	for (int ni = 0; ni < (*g).top.v.e[i].len; ni++)
	{
		int ii = (*g).top.v.e[i].e[ni];
		int i_i_ii = (*g).top.i_of.e[i].e[ni];
		remove_array_int(&((*g).edge_wts.e[ii]), i_i_ii);
	}
	remove_edges_vertex_nbrlist(&((*g).top), i);
	(*g).edge_wts.e[i].len = 0;
}

void remove_vertex_edge_wtd_graph(edge_wtd_graph *g, int i)
{
	for (int ni = 0; ni < (*g).top.v.e[i].len; ni++)
	{
		int ii = (*g).top.v.e[i].e[ni];
		int i_i_ii = (*g).top.i_of.e[i].e[ni];
		remove_array_int(&((*g).edge_wts.e[ii]), i_i_ii);
	}
	remove_vertex_nbrlist(&((*g).top), i);
	remove_aarray_int(&((*g).edge_wts), i);
}

void edge_wtd_graph_merge(edge_wtd_graph *g, int i, int j)
{
	for (int ni = 0; ni < (*g).top.v.e[i].len; ni++)
	{
		int ii = (*g).top.v.e[i].e[ni];
		if (ii != j) {}
		else continue;
		int nj = add_edge_edge_wtd_graph_safe(g, j, ii, (*g).edge_wts.e[i].e[ni]);
		if (nj > -1)
		{
			int i_j_ii = (*g).top.i_of.e[j].e[nj];
			(*g).edge_wts.e[j].e[nj] += (*g).edge_wts.e[i].e[ni];
			(*g).edge_wts.e[ii].e[i_j_ii] += (*g).edge_wts.e[i].e[ni];
		}
	}
	remove_vertex_edge_wtd_graph(g, i);
}

// Methods for contracted neighbor lists RESUME: test these! (contr_nbrlist)
// 	
void contr_nbrlist_init(contr_nbrlist *cnb, nbrlist *top, aarray_int *wts)
{
	(*cnb).prec = top;
	//edge_wtd_graph_init(&((*cnb).top), (*top).v.len);
	transcribe_nbrlist(top, &((*cnb).top.top));
	if (wts != NULL) transcribe_aarray_int(wts, &((*cnb).top.edge_wts));
	else
	{
		// Initialize edge_wts to have unit weight for each edge
		aarray_int_init(&((*cnb).top.edge_wts), (*top).v.len);
		for (int i = 0; i < (*top).v.len; i++)
		{
			array_int wts_i;
			array_int_init(&wts_i, (*top).v.e[i].len);
			wts_i.len = (*top).v.e[i].len;
			for (int ni = 0; ni < (*top).v.e[i].len; ni++) wts_i.e[ni] = 1;
			add2aarray_int(&((*cnb).top.edge_wts), wts_i);
		}
	}
	aarray_int_init(&((*cnb).fibers), (*top).v.len);
	array_int_init(&((*cnb).fiber_addr), (*top).v.len);
	array_int_init(&((*cnb).map), (*top).v.len);
	for (int i = 0; i < (*top).v.len; i++)
	{
		add2array_int(&((*cnb).map), i);
		add2array_int(&((*cnb).fiber_addr), 0);
		extend_aarray_int(&((*cnb).fibers));
		add2array_int(&((*cnb).fibers.e[i]), i);
	}
}

void free_contr_nbrlist(contr_nbrlist *cnb)
{
	free_edge_wtd_graph(&((*cnb).top));
	free_aarray_int(&((*cnb).fibers));
	free_array_int(&((*cnb).fiber_addr));
	free_array_int(&((*cnb).map));
}

void transcribe_contr_nbrlist(contr_nbrlist *src, contr_nbrlist *dest)
{
	transcribe_edge_wtd_graph(&((*src).top), &((*dest).top));
	transcribe_aarray_int(&((*src).fibers), &((*dest).fibers));
	(*dest).prec = (*src).prec;
	transcribe_array_int(&((*src).map), &((*dest).map));
	transcribe_array_int(&((*src).fiber_addr), &((*dest).fiber_addr));
}

char contr_nbrlist_fiber_corresp(contr_nbrlist *cnb1, contr_nbrlist *cnb2, array_int *pvs)
{
  //	printf("contr_nbrlist_fiber_corresp:\n");
	if (pvs != NULL) {}
	else return contr_nbrlist_vertex_corresp(cnb1, cnb2);
	// Check that the correspondence associated with the pvs vertices is well defined
	int corresp[(*cnb1).top.top.v.len];
	char match = 1;
	for (int ci = 0; ci < (*cnb1).top.top.v.len; ci++) corresp[ci] = -1;
	for (int pi = 0; pi < (*pvs).len; pi++)
	{
		int i = (*pvs).e[pi];
		int ci1 = (*cnb1).map.e[i];
		int ci2 = (*cnb2).map.e[i];
		if (corresp[ci1] == -1)
		{
			if ((*cnb1).fibers.e[ci1].len == (*cnb2).fibers.e[ci2].len) {}
			else 
			{
				match = 0;
				break;
			}
			corresp[ci1] = ci2;
		}
		else
		{
			if (corresp[ci1] == ci2) {}
			else 
			{
				match = 0;
				break;
			}
		}
	}
	//	printf("(done)\n");
	return match;
}

char contr_nbrlist_vertex_corresp(contr_nbrlist *cnb1, contr_nbrlist *cnb2)
{
	printf("contr_nbrlist_vertex_corresp:\n");
	char match = 1;
	// This presumes that (*cnb1).prec is equivalent to (*cnb2).prec
	if ((*(*cnb1).prec).v.len == (*(*cnb2).prec).v.len) {}
	else match = 0;
	if ((*cnb1).top.top.v.len == (*cnb2).top.top.v.len) {}
	else match = 0;
	if (match)
	{
		int corresp[(*cnb1).top.top.v.len];
		for (int ci = 0; ci < (*cnb1).top.top.v.len; ci++) 
		{
			corresp[ci] = -1;
			int ri = (*cnb1).fibers.e[ci].e[0];
			int cii = (*cnb2).map.e[ri];
			if ((*cnb1).fibers.e[ci].len == (*cnb2).fibers.e[cii].len) {}
			else
			{
				match = 0;
				break;
			}
		}
		if (match)
		{
			for (int i = 0; i < (*(*cnb1).prec).v.len; i++) 
			{
				int ci1 = (*cnb1).map.e[i];
				int ci2 = (*cnb2).map.e[i];
				if (corresp[ci1] == -1) corresp[ci1] = ci2;
				else if (corresp[ci1] == ci2) {}
				else
				{
					match = 0;
					break;
				}
			}
		}
	}
	printf("(done)\n");
	return match;
}

// RESUME: Test this!
void contr_nbrlist_expand(contr_nbrlist *cnb, int ci, aarray_int *wts)
{
	array_int old_fiber = (*cnb).fibers.e[ci];
	array_int new_fiber;
	array_int_init(&new_fiber, 1);
	new_fiber.len = 1;
	new_fiber.e[0] = old_fiber.e[0];
	(*cnb).fibers.e[ci] = new_fiber;
	int init_len = (*cnb).fibers.len;
	for (int fi = 1; fi < old_fiber.len; fi++)
	{
		int i = old_fiber.e[fi];
		extend_edge_wtd_graph(&((*cnb).top));
		array_int fiber_i;
		array_int_init(&fiber_i, 1);
		add2array_int(&fiber_i, i);
		(*cnb).map.e[i] = (*cnb).fibers.len;
		(*cnb).fiber_addr.e[i] = 0;
		add2aarray_int(&((*cnb).fibers), fiber_i);
	}
	for (int fi = 1; fi < old_fiber.len; fi++)
	{
		int j = old_fiber.e[fi];
		int cj = (*cnb).map.e[j];
		for (int ni = 0; ni < (*(*cnb).prec).v.e[j].len; ni++)
		{
			int ii = (*(*cnb).prec).v.e[j].e[ni];
			int cii = (*cnb).map.e[ii];
			int wt = wts != NULL ? (*wts).e[j].e[ni] : 1;
			if (cii >= init_len)
			{
				if (cii < cj) add_edge_edge_wtd_graph(&((*cnb).top), cj, cii, wt);
			}
			else if (cii != ci)
			{
				int njii = add_edge_edge_wtd_graph_safe(&((*cnb).top), cj, cii, wt);
				if (njii > -1)
				{
					int i_j_ii = (*cnb).top.top.i_of.e[j].e[njii];
					(*cnb).top.edge_wts.e[j].e[njii] += wt;
					(*cnb).top.edge_wts.e[ii].e[i_j_ii] += wt;
				}
			}
		}
	}
	//remove_aarray_int(&((*cnb).fibers), ci);
	free_array_int(&old_fiber);
	/*if ((*cnb).fibers.len > ci)
	{
		for (int fi = 0; fi < (*cnb).fibers.e[ci].len; fi++) (*cnb).map.e[(*cnb).fibers.e[ci].e[fi]] = ci;
	}*/

	//remove_vertex_edge_wtd_graph(&((*cnb).top), ci);
	remove_edges_vertex_edge_wtd_graph(&((*cnb).top), ci);
	int ri = (*cnb).fibers.e[ci].e[0];
	for (int ni = 0; ni < (*(*cnb).prec).v.e[ri].len; ni++)
	{
		int ii = (*(*cnb).prec).v.e[ri].e[ni];
		int cii = (*cnb).map.e[ii];
		int wt = wts != NULL ? (*wts).e[ri].e[ni] : 1;
		if (cii >= init_len) add_edge_edge_wtd_graph(&((*cnb).top), ci, cii, wt);
		else
		{
			int nii = add_edge_edge_wtd_graph_safe(&((*cnb).top), ci, cii, wt);
			if (nii > -1)
			{
				int i_i_ii = (*cnb).top.top.i_of.e[ci].e[nii];
				(*cnb).top.edge_wts.e[ci].e[nii] += wt;
				(*cnb).top.edge_wts.e[cii].e[i_i_ii] += wt;
			}
		}
	}
}

void contr_nbrlist_merge_cluster(contr_nbrlist *cnb, int ci, int cj)
{
	for (int fi = 0; fi < (*cnb).fibers.e[ci].len; fi++)
	{
		(*cnb).fiber_addr.e[(*cnb).fibers.e[ci].e[fi]] = (*cnb).fibers.e[cj].len;
		add2array_int(&((*cnb).fibers.e[cj]), (*cnb).fibers.e[ci].e[fi]);
		(*cnb).map.e[(*cnb).fibers.e[ci].e[fi]] = cj;
	}
	remove_aarray_int(&((*cnb).fibers), ci);
	if (ci < (*cnb).fibers.len)
	{
		for (int fi = 0; fi < (*cnb).fibers.e[ci].len; fi++)
		{
			(*cnb).map.e[(*cnb).fibers.e[ci].e[fi]] = ci;
		}
	}
	edge_wtd_graph_merge(&((*cnb).top), ci, cj);
}

void contr_nbrlist_merge(contr_nbrlist *cnb, int i, int j)
{
	int ci = (*cnb).map.e[i];
	int cj = (*cnb).map.e[j];
	contr_nbrlist_merge_cluster(cnb, ci, cj);
}

// Methods for dir_graph
//
void dir_graph_init_precise(dir_graph *dgr, int mem)
{
	aarray_int_init_precise(&((*dgr).out), mem, 1);
	aarray_int_init_precise(&((*dgr).in), mem, 1);
	aarray_int_init_precise(&((*dgr).i_of_io), mem, 1);
	aarray_int_init_precise(&((*dgr).i_of_oi), mem, 1);
}

void dir_graph_init(dir_graph *dgr)
{
	dir_graph_init_precise(dgr, 1);
}

void prep_dir_graph(dir_graph *dgr)
{
	if ((*dgr).out.len < (*dgr).out.mem)
       	{
		if ((*dgr).out.len > ((*dgr).out.mem >> 2)) {}
		else
		{
			contract_aarray_int(&((*dgr).out));
			contract_aarray_int(&((*dgr).in));
			contract_aarray_int(&((*dgr).i_of_io));
			contract_aarray_int(&((*dgr).i_of_oi));
		}
	}
	else
	{
		add_mem_aarray_int(&((*dgr).out));
		add_mem_aarray_int(&((*dgr).in));
		add_mem_aarray_int(&((*dgr).i_of_io));
		add_mem_aarray_int(&((*dgr).i_of_oi));
	}
}

void check_dir_graph(dir_graph *dgr)
{
	if ((*dgr).out.len <= (*dgr).out.mem) {}
	else
	{
		add_mem_aarray_int_until(&((*dgr).out), (*dgr).out.len);
		add_mem_aarray_int_until(&((*dgr).in), (*dgr).out.len);
		add_mem_aarray_int_until(&((*dgr).i_of_io), (*dgr).out.len);
		add_mem_aarray_int_until(&((*dgr).i_of_oi), (*dgr).out.len);
	}
}

void extend_dir_graph(dir_graph *dgr)
{
	prep_dir_graph(dgr);
	(*dgr).out.len += 1;
	(*dgr).in.len += 1;
	(*dgr).i_of_io.len += 1;
	(*dgr).i_of_oi.len += 1;
}

void extend_dir_graph_n(dir_graph *dgr, int N)
{
	(*dgr).out.len += N;
	(*dgr).in.len += N;
	(*dgr).i_of_io.len += N;
	(*dgr).i_of_oi.len += N;
	check_dir_graph(dgr);
}

void ensure_dir_graph_size_n(dir_graph *dgr, int N)
{
	if ((*dgr).out.len >= N) {}
	else extend_dir_graph_n(dgr, N - (*dgr).out.len);
}

void set_len_dir_graph(dir_graph *dgr, int len)
{
	(*dgr).out.len = len;
	(*dgr).in.len = len;
	(*dgr).i_of_io.len = len;
	(*dgr).i_of_oi.len = len;
	check_dir_graph(dgr);
}

void add_edge_dir_graph(dir_graph *dgr, int vertex1, int vertex2)
{
	add2array_int(&((*dgr).i_of_io.e[vertex2]), (*dgr).out.e[vertex1].len);
	add2array_int(&((*dgr).i_of_oi.e[vertex1]), (*dgr).in.e[vertex2].len);
	add2array_int(&((*dgr).out.e[vertex1]), vertex2);
	add2array_int(&((*dgr).in.e[vertex2]), vertex1);
}

void add_edge_dir_graph_safe(dir_graph *dgr, int vertex1, int vertex2)
{
	// Check that vertex2 does not already appear in the list of outgoing edges from vertex 1
	char new_edge = 1;
	for (int ii = 0; ii < (*dgr).out.e[vertex1].len; ii++)
	{
		if (vertex2 != (*dgr).out.e[vertex1].e[ii]) {}
		else 
		{
			new_edge = 0;
			break;
		}
	}
	if (new_edge) add_edge_dir_graph(dgr, vertex1, vertex2);
}

void remove_out_dir_graph_exp(aarray_int *out, aarray_int *in, aarray_int *i_of_oi, aarray_int *i_of_io, int vertex, int local_nbr_index)
{
	// Update: (*dgr).out.e[vertex], (*dgr).in.e[vertex2], (*dgr).i_of_oi.e[vertex], (*dgr).i_of_io.e[vertex_1], (*dgr).i_of_oi.e[vertex_2]
	int vertex2 = (*out).e[vertex].e[local_nbr_index];
	int i_of_12 = (*i_of_oi).e[vertex].e[local_nbr_index];
	remove_array_int(&((*out).e[vertex]), local_nbr_index);
	remove_array_int(&((*in).e[vertex2]), i_of_12);
	remove_array_int(&((*i_of_oi).e[vertex]), local_nbr_index);
	remove_array_int(&((*i_of_io).e[vertex2]), i_of_12);
	int vertex_ = (*out).e[vertex].e[local_nbr_index];
	int i_of_1_ = (*i_of_oi).e[vertex].e[local_nbr_index];
	(*i_of_io).e[vertex_].e[i_of_1_] = local_nbr_index;
	int vertex2_ = (*in).e[vertex2].e[i_of_12];
	int i_of_22_ = (*i_of_io).e[vertex2].e[i_of_12];
	(*i_of_oi).e[vertex2_].e[i_of_22_] = i_of_12;

}

void remove_out_dir_graph(dir_graph *dgr, int vertex, int local_nbr_index)
{
	remove_out_dir_graph_exp(&((*dgr).out), &((*dgr).in), &((*dgr).i_of_oi), &((*dgr).i_of_io), vertex, local_nbr_index);
	// Update: (*dgr).out.e[vertex], (*dgr).in.e[vertex2], (*dgr).i_of_oi.e[vertex], (*dgr).i_of_io.e[vertex_1], (*dgr).i_of_oi.e[vertex_2]
	/*
	 * int vertex2 = (*dgr).out.e[vertex].e[local_nbr_index];
	int i_of_12 = (*dgr).i_of_oi.e[vertex].e[local_nbr_index];
	remove_array_int(&((*dgr).out.e[vertex]), local_nbr_index);
	remove_array_int(&((*dgr).in.e[vertex2]), i_of_12);
	remove_array_int(&((*dgr).i_of_oi.e[vertex]), local_nbr_index);
	remove_array_int(&((*dgr).i_of_io.e[vertex2]), i_of_12);
	int vertex_ = (*dgr).out.e[vertex].e[local_nbr_index];
	int i_of_1_ = (*dgr).i_of_oi.e[vertex].e[local_nbr_index];
	(*dgr).i_of_io.e[vertex_].e[i_of_i_] = local_nbr_index;
	int vertex2_ = (*dgr).in.e[vertex2].e[i_of_12];
	int i_of_22_ = (*dgr).i_of_io.e[vertex2].e[i_of_12];
	(*dgr).i_of_oi.e[vertex2_].e[i_of_22_] = i_of_12;
	*/
}

void remove_in_dir_graph(dir_graph *dgr, int vertex, int local_nbr_index)
{
	remove_out_dir_graph_exp(&((*dgr).in), &((*dgr).out), &((*dgr).i_of_io), &((*dgr).i_of_oi), vertex, local_nbr_index);
}

void remove_all_edges_dir_graph(dir_graph *dgr)
{
	for (int i = 0; i < (*dgr).out.len; i++)
	{
		(*dgr).out.e[i].len = 0;
		(*dgr).in.e[i].len = 0;
		(*dgr).i_of_io.e[i].len = 0;
		(*dgr).i_of_oi.e[i].len = 0;
	}	
}

void remove_edges_vertex_dir_graph(dir_graph *dgr, int vertex)
{
	int ii_out = (*dgr).out.e[vertex].len;
	while (ii_out > 0)
	{
		ii_out -= 1;
		remove_out_dir_graph(dgr, vertex, ii_out);
	}
	int ii_in = (*dgr).in.e[vertex].len;
	while (ii_in > 0)
	{
		ii_in -= 1;
		remove_in_dir_graph(dgr, vertex, ii_in);
	}
}

void remove_vertex_dir_graph(dir_graph *dgr, int vertex)
{
	remove_edges_vertex_dir_graph(dgr, vertex);
	remove_aarray_int(&((*dgr).out), vertex);
	remove_aarray_int(&((*dgr).in), vertex);
	remove_aarray_int(&((*dgr).i_of_oi), vertex);
	remove_aarray_int(&((*dgr).i_of_io), vertex);
	for (int ii = 0; ii < (*dgr).out.e[vertex].len; ii++)
	{
		int j_ = (*dgr).out.e[vertex].e[ii];
		int i_of_ij_ = (*dgr).i_of_oi.e[vertex].e[ii];
		(*dgr).in.e[j_].e[i_of_ij_] = vertex;

	}
	for (int ii = 0; ii < (*dgr).in.e[vertex].len; ii++)
	{
		int _j = (*dgr).in.e[vertex].e[ii];
		int i_of_i_j = (*dgr).i_of_io.e[vertex].e[ii];
		(*dgr).out.e[_j].e[i_of_i_j] = vertex;
	}
}

void fprintf_dir_graph(dir_graph *dgr, FILE *ofile)
{
	if (ofile != NULL)
	{
		for (int i = 0; i < (*dgr).out.len; i++)
		{
			fprintf(ofile, "%d %d ", (*dgr).out.e[i].len, (*dgr).in.e[i].len);
			for (int ii = 0; ii < (*dgr).out.e[i].len; ii++)
			{
				fprintf(ofile, "%d ", (*dgr).out.e[i].e[ii]);
			}
			for (int ii = 0; ii < (*dgr).in.e[i].len; ii++)
			{
				fprintf(ofile, "%d ", (*dgr).in.e[i].e[ii]);
			}
			fprintf(ofile, "\n");
		}

	}
}

void write_dir_graph(dir_graph *dgr, char *ofname)
{
	FILE *ofile = fopen(ofname, "w");
	if (ofile != NULL)
	{
		for (int i = 0; i < (*dgr).out.len; i++)
		{
			fprintf(ofile, "%d %d ", (*dgr).out.e[i].len, (*dgr).in.e[i].len);
			for (int ii = 0; ii < (*dgr).out.e[i].len; ii++)
			{
				fprintf(ofile, "%d ", (*dgr).out.e[i].e[ii]);
			}
			for (int ii = 0; ii < (*dgr).in.e[i].len; ii++)
			{
				fprintf(ofile, "%d ", (*dgr).in.e[i].e[ii]);
			}
			fprintf(ofile, "\n");
		}
		fclose(ofile);
	}
}

void free_dir_graph(dir_graph *dgr)
{
	free_aarray_int(&((*dgr).out));
	free_aarray_int(&((*dgr).in));
	free_aarray_int(&((*dgr).i_of_oi));
	free_aarray_int(&((*dgr).i_of_io));
}
// end (methods for dir_graph)

// Methods for variable length arrays of doubles
// RESUME: make sure that all instances of array_double_init are consistent with the new implementation
void array_double_init(array_double *a, int size_)
{
	if (size_ > 0)
	{
		(*a).e = (double *) calloc(size_, sizeof(double));
	}
	else 
	{
		(*a).e = NULL;
	}
	(*a).mem = size_;
	(*a).len = 0;
}

void transcribe_array_double(array_double *src, array_double *dest)
{
	array_double_init(dest, (*src).mem);
	(*dest).len = (*src).len;
	for (int i = 0; i < (*src).len; i++)
	{
		(*dest).e[i] = (*src).e[i];
	}
}

void add_mem_array_double(array_double *a)
{
	if ((*a).mem > 0) {}
	else (*a).mem = 1;
	(*a).mem <<= 1;
	double *ne = (double *) calloc((*a).mem, sizeof(double));
	for (int i = 0; i < (*a).len; i++)
	{
		ne[i] = (*a).e[i];
	}
	if ((*a).e != NULL) free((*a).e);
	(*a).e = ne;
}

void add_mem_array_double_until(array_double *a, int i)
{
	if ((*a).mem > i) return;
	if ((*a).mem > 0) {}
	else (*a).mem = 1;
	while ((*a).mem <= i)
	{
		(*a).mem <<= 1;
	}
	double *ne = (double *) calloc((*a).mem, sizeof(double));
	for (int ii = 0; ii < (*a).len; ii++)
	{
		ne[ii] = (*a).e[ii];
	}
	if ((*a).e != NULL) free((*a).e);
	(*a).e = ne;
}

void sum_array_double(array_double *a, array_double *b)
{
	int min_len = (*a).len <= (*b).len ? (*a).len : (*b).len;
	for (int i = 0; i < min_len; i++) (*a).e[i] += (*b).e[i];
}

void scale_array_double(array_double *a, double s)
{
	for (int i = 0; i < (*a).len; i++) (*a).e[i] *= s;
}

void add2array_double(array_double *a, double i)
{
	if ((*a).len < (*a).mem) {}
	else
	{
		add_mem_array_double(a);
	}
	(*a).e[(*a).len] = i;
	(*a).len += 1;
}

void print_array_double(array_double a)
{
	for (int i = 0; i < a.len; i++)
	{
		printf("%g ", a.e[i]);
	}
	printf("\n");
}

void remove_last_array_double(array_double *a)
{
	if ((*a).len > 0)
	{
		(*a).len -= 1;
		if ((*a).len < ((*a).mem >> 2)) contract_array_double(a);
	}
}

void remove_array_double(array_double *a, int n)
{
	if (n < (*a).len)
	{
		(*a).len -= 1;
		(*a).e[n] = (*a).e[(*a).len];
		if ((*a).len < ((*a).mem >> 2)) contract_array_double(a);
	}
}

void contract_array_double(array_double *a)
{
	if ((*a).mem > 0) 
	{
		int newsize = (*a).mem >> 1;
		double *data = (double *) calloc(newsize, sizeof(double));
		int ulim = (*a).len < newsize ? (*a).len : newsize;
		for (int i = 0; i < ulim; i++) 
		{
			data[i] = (*a).e[i];
		}
		free((*a).e);
		(*a).e = data;
		(*a).mem = newsize;
	}
}

void fprintf_array_double(array_double *a, FILE *ofile)
{
	ofile = ofile != NULL ? ofile : stdout;
	for (int i = 0; i < (*a).len; i++) fprintf(ofile, "%g ", (*a).e[i]);
	fprintf(ofile, "\n");
}

void free_array_double(array_double *a)
{
	if ((*a).mem > 0)
	{
		if ((*a).e != NULL) free((*a).e);
		(*a).e = NULL;
		(*a).len = 0;
		(*a).mem = 0;
	}
}

void reset_array_double(array_double *a)
{
	(*a).len = 0;
}

void merge_array_double(array_double *a, double *e, int *c, int imin, int imax, int jmin, int jmax)
{
	if ((*a).e[imin] < (*a).e[jmin])
	{
		int il = imin;
		while ((*a).e[il] < (*a).e[jmin] && il < imax)
		{
			e[*c] = (*a).e[il];
			il += 1;
			*c += 1;
		}
		merge_array_double(a, e, c, il, imax, jmin, jmax);
	}
	else
	{
		int ir = jmin;
		while ((*a).e[ir] < (*a).e[imin] && ir < jmax)
		{
			e[*c] = (*a).e[ir];
			ir += 1;
			*c += 1;
		}
		merge_array_double(a, e, c, imin, imax, ir, jmax);
	}
}

// sort a variable length array along the index interval [i, j)
void merge_sort_array_double(array_double *a, double *e, int i, int j)
{
	int jmi = j - i;
	if (jmi > 2)
	{
		int mdpt = (i + j) / 2;
		merge_sort_array_double(a, e, i, mdpt);
		merge_sort_array_double(a, e, mdpt, j);
		// merge the two ordered lists
		int c = i;
		merge_array_double(a, e, &c, i, mdpt, mdpt, j);
	}
	else if (jmi == 2)
	{
		int ip1 = i + 1;
		if ((*a).e[ip1] > (*a).e[i]) {}
		else
		{
			int tmp = (*a).e[ip1];
			(*a).e[ip1] = (*a).e[i];
			(*a).e[i] = tmp;
		}
	}
}

void sort_array_double(array_double *a)
{
	double *e = (double *) calloc((*a).mem, sizeof(double));
	merge_sort_array_double(a, e, 0, (*a).len);
	free((*a).e);
	(*a).e = e;
}

// Methods for variable length arrays of chars
void array_char_init(array_char *a, int size_)
{
	if (size_ > 0)
	{
		(*a).e = (char *) calloc(size_, sizeof(char));
		(*a).mem = size_;
		(*a).len = 0;
	}
	else 
	{
		(*a).e = NULL;
		(*a).mem = 0;
		(*a).len = 0;
	}
}

char array_char_local_match(array_char *a, int i, char *b, int b_len)
{
	char *ae = &((*a).e[i]);
	for (int ii = 0; ii < b_len; ii++)
	{
		if (ae[ii] == b[ii]) {}
		else return 0;
	}
	return 1;
}

char array_char_contains_substring(array_char *a, char *b, int b_len)
{
	if ((*a).len >= b_len)
	{
		int ulim = (*a).len - b_len;
		for (int i = 0; i <= ulim; i++)
		{
			if (array_char_local_match(a, i, b, b_len)) return 1;
		}
	}
	return 0;
}

void transcribe_array_char(array_char *src, array_char *dest)
{
	array_char_init(dest, (*src).mem);
	(*dest).len = (*src).len;
	for (int i = 0; i < (*src).len; i++)
	{
		(*dest).e[i] = (*src).e[i];
	}
}

int array_char_search(array_char *a, char elem)
{
	for (int i = 0; i < (*a).len; i++)
	{
		if ((*a).e[i] == elem) return i;
	}
	return -1;
}

char array_char_contains(array_char *a, char elem)
{
	for (int i = 0; i < (*a).len; i++)
	{
		if ((*a).e[i] == elem) return 1;
	}
	return 0;
}


void add_mem_array_char(array_char *a)
{
	if ((*a).mem > 0) {}
	else (*a).mem = 1;
	(*a).mem <<= 1;
	char *ne = (char *) calloc((*a).mem, sizeof(char));
	for (int i = 0; i < (*a).len; i++)
	{
		ne[i] = (*a).e[i];
	}
	if ((*a).e != NULL) free((*a).e);
	(*a).e = ne;
}

void add_mem_array_char_until(array_char *a, int i)
{
	if ((*a).mem > i) return;
	if ((*a).mem > 0) {}
	else (*a).mem = 1;
	while ((*a).mem <= i)
	{
		(*a).mem <<= 1;
	}
	char *ne = (char *) calloc((*a).mem, sizeof(char));
	for (int ii = 0; ii < (*a).len; ii++)
	{
		ne[ii] = (*a).e[ii];
	}
	if ((*a).e != NULL) free((*a).e);
	(*a).e = ne;
}

void append_array_char(array_char *a, char *buf, int len)
{
	if (len > -1) {}
	else
	{
		len = str_len(buf);
	}
	add_mem_array_char_until(a, (*a).len + len);
	for (int i = 0; i < len; i++)
	{
		(*a).e[(*a).len] = buf[i];
		(*a).len += 1;
	}
}

void add2array_char(array_char *a, char i)
{
	if ((*a).len < (*a).mem) {}
	else
	{
		add_mem_array_char(a);
	}
	(*a).e[(*a).len] = i;
	(*a).len += 1;
}

void print_array_char(array_char a)
{
	for (int i = 0; i < a.len; i++)
	{
		printf("%c ", a.e[i]);
	}
	printf("\n");
}

void remove_last_array_char(array_char *a)
{
	if ((*a).len > 0)
	{
		(*a).len -= 1;
		if ((*a).len < ((*a).mem >> 2)) contract_array_char(a);
	}
}

void remove_array_char(array_char *a, int n)
{
	if (n < (*a).len)
	{
		(*a).len -= 1;
		(*a).e[n] = (*a).e[(*a).len];
		if ((*a).len < ((*a).mem >> 2)) contract_array_char(a);
	}
}

void contract_array_char(array_char *a)
{
	if ((*a).mem > 0 && (*a).e != NULL) 
	{
		int newsize = (*a).mem >> 1;
		char *data = (char *) calloc(newsize, sizeof(char));
		int ulim = (*a).len < newsize ? (*a).len : newsize;
		for (int i = 0; i < ulim; i++) 
		{
			data[i] = (*a).e[i];
		}
		free((*a).e);
		(*a).e = data;
		(*a).mem = newsize;
	}
}

void fprintf_array_char(array_char *a, FILE *ofile)
{
	ofile = ofile != NULL ? ofile : stdout;
	for (int i = 0; i < (*a).len; i++) fprintf(ofile, "%d ", (*a).e[i]);
	fprintf(ofile, "\n");
}

void free_array_char(array_char *a)
{
	if ((*a).mem > 0)
	{
		if ((*a).e != NULL) free((*a).e);
		(*a).len = 0;
		(*a).e = NULL;
		(*a).mem = 0;
	}
}

void reset_array_char(array_char *a)
{
	(*a).len = 0;
}

void merge_array_char(array_char *a, char *e, int *c, int imin, int imax, int jmin, int jmax)
{
	if ((*a).e[imin] < (*a).e[jmin])
	{
		int il = imin;
		while ((*a).e[il] < (*a).e[jmin] && il < imax)
		{
			e[*c] = (*a).e[il];
			il += 1;
			*c += 1;
		}
		merge_array_char(a, e, c, il, imax, jmin, jmax);
	}
	else
	{
		int ir = jmin;
		while ((*a).e[ir] < (*a).e[imin] && ir < jmax)
		{
			e[*c] = (*a).e[ir];
			ir += 1;
			*c += 1;
		}
		merge_array_char(a, e, c, imin, imax, ir, jmax);
	}
}

// sort a variable length array along the index interval [i, j)
void merge_sort_array_char(array_char *a, char *e, int i, int j)
{
	int jmi = j - i;
	if (jmi > 2)
	{
		int mdpt = (i + j) / 2;
		merge_sort_array_char(a, e, i, mdpt);
		merge_sort_array_char(a, e, mdpt, j);
		// merge the two ordered lists
		int c = i;
		merge_array_char(a, e, &c, i, mdpt, mdpt, j);
	}
	else if (jmi == 2)
	{
		int ip1 = i + 1;
		if ((*a).e[ip1] > (*a).e[i]) {}
		else
		{
			int tmp = (*a).e[ip1];
			(*a).e[ip1] = (*a).e[i];
			(*a).e[i] = tmp;
		}
	}
}

void sort_array_char(array_char *a)
{
	if ((*a).e != NULL) {}
	else return;
	char *e = (char *) calloc((*a).mem, sizeof(char));
	merge_sort_array_char(a, e, 0, (*a).len);
	free((*a).e);
	(*a).e = e;
}

// aarray_char #beginning
void aarray_char_init(aarray_char *aa, int mem)
{
	aarray_char_init_precise(aa, mem, 0);
}

void transcribe_aarray_char(aarray_char *src, aarray_char *dest)
{
	(*dest).mem = (*src).mem;
	(*dest).e = (array_char *) calloc((*dest).mem, sizeof(array_char));
	(*dest).len = (*src).len;
	for (int i = 0; i < (*src).len; i++)
	{
		transcribe_array_char(&((*src).e[i]), &((*dest).e[i]));
	}
}

// An alternative to the last method that can be tuned to improve performance
void aarray_char_init_precise(aarray_char *aa, int mem1, int mem2)
{
	(*aa).mem = mem1;
	(*aa).e = (array_char *) calloc(mem1, sizeof(array_char));
	for (int i = 0; i < (*aa).mem; i++)
	{
		array_char_init(&(*aa).e[i], mem2);
	}
	(*aa).len = 0;
}

void add_mem_aarray_char(aarray_char *aa)
{
	int init_mem = (*aa).mem;
	if ((*aa).mem > 0) {}
	else (*aa).mem = 1;
	(*aa).mem <<= 1;
	array_char *ne = (array_char *) calloc((*aa).mem, sizeof(array_char));
	int i = 0;
	while ((*aa).e[i].mem > 0 && i < init_mem)
	{
		ne[i] = (*aa).e[i];
		i += 1;
	}
	for (int ii = i; ii < (*aa).mem; ii++)
	{
		ne[ii].e = NULL;
		ne[ii].mem = 0;
		ne[ii].len = 0;
	}
	if ((*aa).e != NULL) free((*aa).e);
	(*aa).e = ne;
}

void contract_aarray_char(aarray_char *a)
{
	if ((*a).mem > 0 && (*a).e != NULL)
	{
		int newsize = (*a).mem >> 1;
		array_char *data = (array_char *) calloc(newsize, sizeof(array_char));
		int ulim = (*a).len < newsize ? (*a).len : newsize;
		for (int i = 0; i < ulim; i++) 
		{
			data[i] = (*a).e[i];
		}
		for (int i = ulim; i < newsize; i++)
		{
			data[i] = (*a).e[i];
		}
		free((*a).e);
		(*a).e = data;
		(*a).mem = newsize;
	}
}

void add_mem_aarray_char_until(aarray_char *aa, int i)
{
	if ((*aa).mem > i) return;
	int init_mem = (*aa).mem;
	if ((*aa).mem > 0) {}
	else (*aa).mem = 1;
	while ((*aa).mem <= i) (*aa).mem <<= 1;
	array_char *ne = (array_char *) calloc((*aa).mem, sizeof(array_char));
	int ii = 0; 
	while (ii < init_mem && (*aa).e[ii].mem > 0)
	{
		ne[ii] = (*aa).e[ii];
		ii += 1;
	}
	for (int i = ii; i < (*aa).mem; i++)
	{
		ne[i].e = NULL;
		ne[i].mem = 0;
		ne[i].len = 0;
	}
	if ((*aa).e != NULL) free((*aa).e);
	(*aa).e = ne;
}

void extend_aarray_char(aarray_char *aa)
{
	extend_aarray_char_precise(aa, INIT_A_MEM);
}

void extend_aarray_char_precise(aarray_char *aa, int init_mem)
{
	if ((*aa).len == (*aa).mem)
	{
		add_mem_aarray_char(aa);
	}
	if ((*aa).e[(*aa).len].e != NULL) 
	{
		(*aa).e[(*aa).len].len = 0;
	}
	else
	{
		array_char_init(&((*aa).e[(*aa).len]), init_mem);
	}
	(*aa).len += 1;

}

// RESUME: Check this! 001
void add2aarray_char(aarray_char *aa, array_char a)
{
	if (a.mem > 0) {}
	else 
	{
		printf("Error: attempting to add empty array_char to array of array_char (aarray_char.)\n");
		exit(EXIT_FAILURE);
	}
	if ((*aa).len == (*aa).mem)
	{
		add_mem_aarray_char(aa);
	}
	if ((*aa).e[(*aa).len].mem == 0) {}
	else 
	{
		// NOTE: alternatively, consider swapping the last element of aa with the next NULL element in 
		// 	memory, which might improve performance in some settings.
		free_array_char(&(*aa).e[(*aa).len]);
	}
	(*aa).e[(*aa).len] = a;
	(*aa).len += 1;
}

void reset_aarray_char_elem(int i, aarray_char *aa)
{
	reset_array_char(&((*aa).e[i]));
}

// RESUME: test this
void remove_aarray_char(aarray_char *aa, int i)
{
	if (i < (*aa).len && i > -1)
	{
		(*aa).len -= 1;
		// transcribe the last element of (*aa) onto the i-th element
		array_char aux = (*aa).e[i];
		(*aa).e[i] = (*aa).e[(*aa).len];
		(*aa).e[(*aa).len] = aux;
		reset_aarray_char_elem((*aa).len, aa);
		if ((*aa).len > (*aa).mem >> 2) {}
		else contract_aarray_char(aa);
	}
	else 
	{
		printf("Error: attempting to remove non-existent aarray_char element %d of %d\n", i, (*aa).len);
		exit(EXIT_FAILURE);
	}
}

void fprintf_aarray_char(aarray_char *aa, FILE *ofile)
{
    if (ofile != NULL)
    {
      for (int i = 0; i < (*aa).len; i++)
	{
	  fprintf(ofile, "%d ", (*aa).e[i].len);
	  for (int ii = 0; ii < (*aa).e[i].len; ii++) fprintf(ofile, "%d ", (*aa).e[i].e[ii]);
	  fprintf(ofile, "\n");
	}
    }
}

void load_aarray_char(aarray_char *aa, char *fname)
{
  FILE *ifile = fopen(fname, "r");
  if (ifile != NULL)
    {
      while (1)
	{
	  int n_elem;
	  int status = fscanf(ifile, "%d", &n_elem);
	  if (status != EOF)
	    {
	      int i = (*aa).len;
	      extend_aarray_char(aa);
	      for (int ni = 0; ni < n_elem; ni++)
		{
		  int ii;
		  fscanf(ifile, "%d", &ii);
		  add2array_char(&((*aa).e[i]), (char) ii);
		}
	    }
	  else break;
	}
      fclose(ifile);
    }
}

void free_aarray_char(aarray_char *aa)
{
	if ((*aa).mem > 0)
	{
		int i = 0;
		while (i < (*aa).mem && (*aa).e[i].mem > 0)
		{
			free_array_char(&((*aa).e[i]));
			// free(&((*aa).e[i])); // RESUME: Check this! 000
			i += 1;
		}
		free((*aa).e);
	}
}

void print_aarray_char(aarray_char aa)
{
	for (int i = 0; i < aa.len; i++)
	{
		for (int ii = 0; ii < aa.e[i].len; ii++)
		{
			printf("%d ", aa.e[i].e[ii]);
		}
		printf("\n");
	}
}

// aarray_char #ending

char one_bit(char *a, int *b)
{
	return ((*a) >> (*b)) & 1;
}

void remove_121_int(array_int *inj, array_int *left_inv, int i)
{
	int imi = (*inj).e[i];
	(*left_inv).e[imi] = -1;
	remove_array_int(inj, i);
	imi = (*inj).e[i];
	(*left_inv).e[imi] = i;
}

void add_121_int(array_int *inj, array_int *left_inv, int i)
{
	if (i < (*left_inv).len)
	{
		if ((*left_inv).e[i] == -1)
		{
			(*left_inv).e[i] = (*inj).len;
		}
		else
		{
			printf("Error (add_121_int): proposed map %d -> %d conflicts with existing map %d -> %d (broken injectivity)\n", (*inj).len, i, (*left_inv).e[i], i);
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		add_mem_array_int_until(left_inv, i);
		for (int ii = (*left_inv).len; ii < i; ii++)
		{
			(*left_inv).e[ii] = -1;
		}
		(*left_inv).e[i] = (*inj).len - 1;
	}
	add2array_int(inj, i);
}

void toupper_string(char *a)
{
	int pos = 0;
	while (a[pos] != '\0')
	{
		a[pos] = toupper(a[pos]);
		pos += 1;
	}
}

void tolower_string(char *a)
{
	int pos = 0;
	while (a[pos] != '\0')
	{
		a[pos] = tolower(a[pos]);
		pos += 1;
	}
}

void load_array_template(FILE *infile, void *aptr, char *prompt, char *mode)
{
	if (infile != NULL)
	{
		char prompt_[256];
		sprintf(prompt_, "%s", prompt);
		toupper_string(prompt_);
		int prompt_len = str_len(prompt_);
		char mode_[256];
		sprintf(mode_, "%s", mode);
		tolower_string(mode_);
		char buf[256];
		while (fscanf(infile, "%s\n", buf) != EOF)
		{
			toupper_string(buf);
			if (strcmp(buf, prompt_) == 0)
			{
				long int current_pos;
				current_pos = ftell(infile);
				if (aptr != NULL)
				{
					if (strcmp(mode_, "double") == 0)
					{
						double field;
						array_double *a = (array_double *) aptr;
						while (fscanf(infile, "%lg", &field) == 1)
						{
							add2array_double(a, field);
						}
					}
					else if (strcmp(mode_, "int") == 0)
					{
						int field;
						array_int *a = (array_int *) aptr;
						while (fscanf(infile, "%d", &field) == 1)
						{
							add2array_int(a, field);
						}
					}
					else if (strcmp(mode_, "char") == 0)
					{
						// Determine where the array data terminates
						long int end_pos;
						int n_args = 1;
						while (n_args != EOF)
						{
							end_pos = ftell(infile);
							n_args = fscanf(infile, "%s", buf);
							toupper_string(buf);
							if (strcmp(buf, "END") == 0)
							{
								break;
							}
						}
						fseek(infile, current_pos, SEEK_SET);
						char field;
						array_char *a = (array_char *) aptr;
						while (current_pos != end_pos)
						{
							add2array_char(a, fgetc(infile));
							current_pos = ftell(infile);
						}
					}
				}
				else
				{
					printf("Attempting to load data into non-existing array_double (prompt: %s)\n", prompt);
				}
				return;
			}
		}
	}
}

void load_array_double(FILE *infile, array_double *a, char *prompt)
{
	load_array_template(infile, (void *) a, prompt, "double");
}

void load_array_int(FILE *infile, array_int *a, char *prompt)
{
	load_array_template(infile, (void *) a, prompt, "int");
}

void load_array_char(FILE *infile, array_char *a, char *prompt)
{
	load_array_template(infile, (void *) a, prompt, "char");
}

// Operations between arrays
void elementwise_product_array_double(array_double *a1, array_double *a2, array_double *a3)
{
	int minlen = ((*a1).len >= (*a2).len? (*a1).len: (*a2).len);
	if ((*a3).mem >= minlen) {}
	else
	{
		if ((*a3).mem > 0) 
		{
			free((*a3).e);
		}
		else (*a3).mem = 16;
		do
		{
			(*a3).mem <<= 1;
		} 
		while ((*a3).mem < minlen);
		(*a3).e = (double *) calloc((*a3).mem, sizeof(double));
	}
	(*a3).len = minlen;
	for (int i = 0; i < minlen; i++)
	{
		(*a3).e[i] = (*a1).e[i] * (*a2).e[i];
	}
}

void elementwise_product_array_int(array_int *a1, array_int *a2, array_int *a3)
{
	int minlen = ((*a1).len >= (*a2).len? (*a1).len: (*a2).len);
	if ((*a3).mem >= minlen) {}
	else
	{
		if ((*a3).mem > 0) 
		{
			free((*a3).e);
		}
		else (*a3).mem = 16;
		do
		{
			(*a3).mem <<= 1;
		} 
		while ((*a3).mem < minlen);
		(*a3).e = (int *) calloc((*a3).mem, sizeof(int));
	}
	(*a3).len = minlen;
	for (int i = 0; i < minlen; i++)
	{
		(*a3).e[i] = (*a1).e[i] * (*a2).e[i];
	}
}

// Operations involving strings/character arrays
void array_char_split(array_char *a, char *div, aarray_char *diva) // RESUME: TEST THIS
{
	int div_len = str_len(div);
	if (diva != NULL) {}
	else 
	{
		printf("Sorry! array_char_split doesn't work with null aarray_char instances yet. Try passing the address of an uninitialized instance, or an allocated pointer.\n");
	}
	if ((*diva).mem > 0)
	{
		free_aarray_char(diva);
	}
	(*diva).e = (array_char *) calloc(1, sizeof(array_char));
	array_char_init(&(*diva).e[0], 1);
	(*diva).len = 1;
	(*diva).mem = 1;
	int next = 0;
	int i = 0;
	int count = 0;
	while (i < (*a).len)
	{
		if (div[next] != (*a).e[i])
		{
			// Copy the matching substring into the last element of 'diva'
			for (int ii = 0; ii < next; ii++)
			{
				add2array_char(&((*diva).e[count]), div[ii]);
			}
			// Reset the read index
			next = 0;
			// Add the last character of 'a' into the last element of 'diva'
			add2array_char(&((*diva).e[count]), (*a).e[i]);
		}
		else
		{
			// Check if the buffer is complete
			next += 1;
			if (next == div_len)
			{
				// Add a new element to diva
				count = (*diva).len;
				extend_aarray_char(diva);
				next = 0;
			}
		}
		i += 1;
	}
	return;
}

void array_char_init_from_string(array_char *a, char *str) // RESUME: TEST THIS
{
	if (a != NULL) {}
	else
	{
		printf("Sorry! array_char_init_from_string doesn't work with null pointers yet. (Make sure to pass either the address of a possibly uninitialized instance or an allocated pointer.)\n");
		exit(EXIT_FAILURE);
	}
	int strlen = 0;
	while (str[strlen] != '\0') strlen += 1;
	if ((*a).mem > 0) free((*a).e);
	int mem = 1;
	while (mem < strlen) mem <<= 1;
	(*a).mem = mem;
	(*a).e = (char *) calloc((*a).mem, sizeof(char));
	(*a).len = strlen;
	sprintf((*a).e, "%s", str);
}

void concatenate_array_char(array_char *a, array_char *b) // RESUME: TEST THIS
{
	int new_len = (*a).len + (*b).len;
	if ((*a).mem >= new_len) {}
	else
	{
		add_mem_array_char_until(a, new_len);
	}
	for (int i = 0; i < (*b).len; i++)
	{
		(*a).e[(*a).len] = (*b).e[i];
		(*a).len += 1;
	}
}

int str_len(char *a)
{
	int len = 0;
	while (a[len] != '\0') len += 1;
	return len;
}

void merge_aarray_int(aarray_int *a, array_int *e, int *c, int imin, int imax, int jmin, int jmax, char (*order)(array_int *, array_int *))
{
	if (order(&(*a).e[imin], &(*a).e[jmin]) == -1)
	{
		int il = imin;
		while (order(&(*a).e[il], &(*a).e[jmin]) == -1 && il < imax)
		{
			e[*c] = (*a).e[il];
			il += 1;
			*c += 1;
		}
		merge_aarray_int(a, e, c, il, imax, jmin, jmax, order);
	}
	else
	{
		int ir = jmin;
		while (order(&(*a).e[ir], &(*a).e[imin]) == -1 && ir < jmax)
		{
			e[*c] = (*a).e[ir];
			ir += 1;
			*c += 1;
		}
		merge_aarray_int(a, e, c, imin, imax, ir, jmax, order);
	}
}

void merge_sort_aarray_int(aarray_int *a, array_int *e, int i, int j, char (*order)(array_int *, array_int *))
{
	int jmi = j - i;
	if (jmi > 2)
	{
		int mdpt = (i + j) / 2;
		merge_sort_aarray_int(a, e, i, mdpt, order);
		merge_sort_aarray_int(a, e, mdpt, j, order);
		// merge the two ordered lists
		int c = i;
		merge_aarray_int(a, e, &c, i, mdpt, mdpt, j, order);
	}
	else if (jmi == 2)
	{
		int ip1 = i + 1;
		if (order(&(*a).e[ip1], &(*a).e[i]) == 1) {}
		else
		{
			array_int tmp = (*a).e[ip1];
			(*a).e[ip1] = (*a).e[i];
			(*a).e[i] = tmp;
		}
	}
}

void sort_aarray_int(aarray_int *a, char (*order)(array_int *, array_int *))
{
	array_int *e = (array_int *) calloc((*a).mem, sizeof(array_int));
	merge_sort_aarray_int(a, e, 0, (*a).len, order);
	free((*a).e);
	(*a).e = e;
}

void merge_aarray_double(aarray_double *a, array_double *e, int *c, int imin, int imax, int jmin, int jmax, char (*order)(array_double *, array_double *))
{
	if (order(&(*a).e[imin], &(*a).e[jmin]) == -1)
	{
		int il = imin;
		while (order(&(*a).e[il], &(*a).e[jmin]) == -1 && il < imax)
		{
			e[*c] = (*a).e[il];
			il += 1;
			*c += 1;
		}
		merge_aarray_double(a, e, c, il, imax, jmin, jmax, order);
	}
	else
	{
		int ir = jmin;
		while (order(&(*a).e[ir], &(*a).e[imin]) == -1 && ir < jmax)
		{
			e[*c] = (*a).e[ir];
			ir += 1;
			*c += 1;
		}
		merge_aarray_double(a, e, c, imin, imax, ir, jmax, order);
	}
}

void merge_sort_aarray_double(aarray_double *a, array_double *e, int i, int j, char (*order)(array_double *, array_double *))
{
	int jmi = j - i;
	if (jmi > 2)
	{
		int mdpt = (i + j) / 2;
		merge_sort_aarray_double(a, e, i, mdpt, order);
		merge_sort_aarray_double(a, e, mdpt, j, order);
		// merge the two ordered lists
		int c = i;
		merge_aarray_double(a, e, &c, i, mdpt, mdpt, j, order);
	}
	else if (jmi == 2)
	{
		int ip1 = i + 1;
		if (order(&(*a).e[ip1], &(*a).e[i]) == 1) {}
		else
		{
			array_double tmp = (*a).e[ip1];
			(*a).e[ip1] = (*a).e[i];
			(*a).e[i] = tmp;
		}
	}
}

void sort_aarray_double(aarray_double *a, char (*order)(array_double *, array_double *))
{
	array_double *e = (array_double *) calloc((*a).mem, sizeof(array_double));
	merge_sort_aarray_double(a, e, 0, (*a).len, order);
	free((*a).e);
	(*a).e = e;
}

void merge_aarray_char(aarray_char *a, array_char *e, int *c, int imin, int imax, int jmin, int jmax, char (*order)(array_char *, array_char *))
{
	if (order(&(*a).e[imin], &(*a).e[jmin]) == -1)
	{
		int il = imin;
		while (order(&(*a).e[il], &(*a).e[jmin]) == -1 && il < imax)
		{
			e[*c] = (*a).e[il];
			il += 1;
			*c += 1;
		}
		merge_aarray_char(a, e, c, il, imax, jmin, jmax, order);
	}
	else
	{
		int ir = jmin;
		while (order(&(*a).e[ir], &(*a).e[imin]) == -1 && ir < jmax)
		{
			e[*c] = (*a).e[ir];
			ir += 1;
			*c += 1;
		}
		merge_aarray_char(a, e, c, imin, imax, ir, jmax, order);
	}
}

void merge_sort_aarray_char(aarray_char *a, array_char *e, int i, int j, char (*order)(array_char *, array_char *))
{
	int jmi = j - i;
	if (jmi > 2)
	{
		int mdpt = (i + j) / 2;
		merge_sort_aarray_char(a, e, i, mdpt, order);
		merge_sort_aarray_char(a, e, mdpt, j, order);
		// merge the two ordered lists
		int c = i;
		merge_aarray_char(a, e, &c, i, mdpt, mdpt, j, order);
	}
	else if (jmi == 2)
	{
		int ip1 = i + 1;
		if (order(&(*a).e[ip1], &(*a).e[i]) == 1) {}
		else
		{
			array_char tmp = (*a).e[ip1];
			(*a).e[ip1] = (*a).e[i];
			(*a).e[i] = tmp;
		}
	}
}

void sort_aarray_char(aarray_char *a, char (*order)(array_char *, array_char *))
{
	array_char *e = (array_char *) calloc((*a).mem, sizeof(array_char));
	merge_sort_aarray_char(a, e, 0, (*a).len, order);
	free((*a).e);
	(*a).e = e;
}

void merge_array_voidstar(array_voidstar *a, void **e, int *c, int *imin, int imax, int *jmin, int jmax, char (*order)(void *, void *))
{
	if ((*imin) < imax && (*jmin) < jmax)
	{
		if (order((*a).e[*imin], (*a).e[*jmin]) != 1)
		{
			e[(*c)] = (*a).e[(*imin)];
			(*c) += 1;
			(*imin) += 1;
		}
		else
		{
			e[(*c)] = (*a).e[(*jmin)];
			(*c) += 1;
			(*jmin) += 1;
		}
		merge_array_voidstar(a, e, c, imin, imax, jmin, jmax, order);
	}
	else 
	{
		int llim, ulim;
		if ((*imin) < imax)
		{
			llim = (*imin);
			ulim = imax;
		}
		else if ((*jmin) < jmax)
		{
			llim = (*jmin);
			ulim = imax;
		}
		else return;
		for (int i = llim; i < ulim; i++)
		{
			e[(*c)] = (*a).e[i];
			(*c) += 1;
		}
	}
}

void merge_array_voidstar_permutation(array_voidstar *a, int *pi, int *buf, int i0, int i1, int mdpt, char (*order)(void *, void *))
{
	if (order((*a).e[pi[mdpt]], (*a).e[pi[mdpt - 1]]) != -1) 
	{
		return;
	}
	else if (order((*a).e[pi[i0]], (*a).e[pi[i1 - 1]]) != -1)
	{
		// case mdpt - i0 > i1 - mdpt:
		// p[i0 + 1: mdpt] -> p[mdpt : i1] -> p[i0 : mdpt - 1] 
		// p[i0] -> p[mdpt - 1]
		// case mdpt - i0 == i1 - mdpt:
		// p[i0 : mdpt] <-> p[mdpt : i1]
		int mdptmi0 = mdpt - i0;
		int i1mmdpt = i1 - mdpt;
		int *seg_a = &(buf[i0]);
		int *seg_b = &pi[mdpt];
		for (int i = 0; i < i1mmdpt; i++)
		{
			//int ipi0 = i + i0;
			seg_a[i] = seg_b[i];
			//buf[ipi0] = (*p).e[i + mdpt];
		}
		seg_a = &(buf[i1 - mdptmi0]);
		seg_b = &pi[i0];
		for (int i = 0; i < mdptmi0; i++)
		{
			seg_a[i] = seg_b[i];
		}
		for (int i = i0; i < i1; i++) pi[i] = buf[i];
		return;
	}
	int term = mdpt;
	int c0 = i0;
	int c1 = mdpt;
	int c = i0;
	while ((c0 < mdpt) && (c1 < i1))
	{
		if (order((*a).e[pi[c0]], (*a).e[pi[c1]]) != 1)
		{
			buf[c] = pi[c0];
			c0 += 1;
			c += 1;
		}
		else
		{
			buf[c] = pi[c1];
			c1 += 1;
			c += 1;
		}
	}
	if (c0 < mdpt)
	{
		while (c < i1)
		{
			buf[c] = pi[c0];
			c += 1;
			c0 += 1;
		}
	}
	if (c1 < i1)
	{
		if (c == c1) {}
		else
		{
			printf("Something weird happened in merge_array_int_permutation! left half of array folded in, but c = %d != c1 = %d!\n", c, c1);
			exit(EXIT_FAILURE);
		}
		while (c < i1)
		{
			buf[c] = pi[c1];
			c += 1;
			c1 += 1;
		}
	}
	for (int i = i0; i < i1; i++)
	{
		pi[i] = buf[i];
	}

}

void merge_sort_array_voidstar_permutation(array_voidstar *a, int *pi, int *buf, int i0, int i1, char (*order)(void *, void *))
{
	int i1mi0 = i1 - i0;
	if (i1mi0 > 2)
	{
		int mdpt = (i0 + i1 + 1) >> 1; 
		merge_sort_array_voidstar_permutation(a, pi, buf, i0, mdpt, order);
		merge_sort_array_voidstar_permutation(a, pi, buf, mdpt, i1, order);
		merge_array_voidstar_permutation(a, pi, buf, i0, i1, mdpt, order);
	}
	else
	{
		if (i1mi0 < 2) return;
		else
		{
			int i0p1 = i0 + 1;
			if (order((*a).e[pi[i0]], (*a).e[pi[i0p1]]) != 1) 
			{
				buf[i0] = pi[i0];
				buf[i0p1] = pi[i0p1];
			}
			else
			{
				buf[i0] = pi[i0p1];
				buf[i0p1] = pi[i0];
			}
			pi[i0] = buf[i0];
			pi[i0p1] = buf[i0p1];
		}
	}

}

void sort_array_voidstar_permutation(array_voidstar *a, int *pi, char (*order)(void *, void *))
{
	int aux_pi[(*a).len];
	merge_sort_array_voidstar_permutation(a, pi, &(aux_pi[0]), 0, (*a).len, order);
}

// sort a variable length array along the index interval [i, j)
void merge_sort_array_voidstar(array_voidstar *a, void **e, int i, int j, char (*order)(void *, void *))
{
	int jmi = j - i;
	if (jmi > 2)
	{
		int mdpt = (i + j) >> 1;
		merge_sort_array_voidstar(a, e, i, mdpt, order);
		merge_sort_array_voidstar(a, e, mdpt, j, order);
		// merge the two ordered lists
		int c = i;
		int imin = i;
		int jmin = mdpt;
		merge_array_voidstar(a, e, &c, &imin, mdpt, &jmin, j, order);
		for (int ii = i; ii < j; ii++)
		{
			(*a).e[ii] = e[ii];
		}
	}
	else if (jmi == 2)
	{
		int ip1 = i + 1;
		if (order((*a).e[ip1], (*a).e[i]) != -1) 
		{
			e[ip1] = (*a).e[ip1];
			e[i] = (*a).e[i];
		}
		else
		{
			e[ip1] = (*a).e[i];
			e[i] = (*a).e[ip1];
			(*a).e[ip1] = e[ip1];
			(*a).e[i] = e[i];
		}
	}
}

void sort_array_voidstar(array_voidstar *a, char (*order)(void *, void *))
{
	void **e = (void **) calloc((*a).len, sizeof(void *));
	for (int i = 0; i < (*a).len; i++) e[i] = (*a).e[i];
	merge_sort_array_voidstar(a, e, 0, (*a).len, order);
	for (int i = 0; i < (*a).len; i++)
	{
		(*a).e[i] = e[i];
	}
	free((*a).e);
	(*a).e = e;
}

// Note: this method assumes that the array provided has already been sorted
char search_ordered_array_voidstar(array_voidstar *a, void *e, char (*order)(void *, void *))
{
	if ((*a).len > 0) {}
	else return 0;
	int ubnd = (*a).len;
	int lbnd = 0;
	while (ubnd != lbnd)
	{
		if (order((*a).e[lbnd], e) != 0 && order((*a).e[ubnd], e) != 0) {}
		else return 1;
		int guess = (lbnd + ubnd) / 2;
		if (order((*a).e[guess], e) == -1)
		{
			lbnd = guess;
		}
		else if (order((*a).e[guess], e) == 1)
		{
			ubnd = guess;
		}
		else
		{
			return 1;
		}
	}
	return 0;
}

// Transcriptions of voidstar arrays
// Note: it seems silly to have to repeat the body of these programs for each implementation.
//		(There's probably a more efficient way of expressing this code: like a 'c'-based emulation of 
//		interfaces or object oriented principles with void pointers or similar.)
void transcribe_array_int2voidstar(array_int *a, array_voidstar *b)
{
	if ((*b).mem < (*a).len)
	{
		add_mem_array_voidstar_until(b, (*a).len);
	}
	for (int i = 0; i < (*a).len; i++)
	{
		(*b).e[i] = (void *) &(*a).e[i];
	}
	(*b).len = (*a).len;
}

void transcribe_array_double2voidstar(array_double *a, array_voidstar *b)
{
	if ((*b).mem < (*a).len)
	{
		add_mem_array_voidstar_until(b, (*a).len);
	}
	for (int i = 0; i < (*a).len; i++)
	{
		(*b).e[i] = (void *) &(*a).e[i];
	}
	(*b).len = (*a).len;
}

void transcribe_array_char2voidstar(array_char *a, array_voidstar *b)
{
	if ((*b).mem < (*a).len)
	{
		add_mem_array_voidstar_until(b, (*a).len);
	}
	for (int i = 0; i < (*a).len; i++)
	{
		(*b).e[i] = (void *) &(*a).e[i];
	}
	(*b).len = (*a).len;
}

void transcribe_aaray_int2voidstar(aarray_int *a, array_voidstar *b)
{
	if ((*b).mem < (*a).len)
	{
		add_mem_array_voidstar_until(b, (*a).len);
	}
	for (int i = 0; i < (*a).len; i++)
	{
		(*b).e[i] = (void *) &(*a).e[i];
	}
	(*b).len = (*a).len;
}

void transcribe_aarray_char2voidstar(aarray_char *a, array_voidstar *b)
{
	if ((*b).mem < (*a).len)
	{
		add_mem_array_voidstar_until(b, (*a).len);
	}
	for (int i = 0; i < (*a).len; i++)
	{
		(*b).e[i] = (void *) &(*a).e[i];
	}
	(*b).len = (*a).len;
}

char order_lexical_array_int(array_int *e1, array_int *e2) 
{
	for (int i = 0; i < (*e1).len; i++)
	{
		if (i < (*e2).len) {}
		else
		{
			return 1;
		}
		if ((*e1).e[i] > (*e2).e[i]) return 1;
		else if ((*e1).e[i] < (*e2).e[i]) return -1;
	}
	return 0;
}

char order_lexical_array_double(array_double *e1, array_double *e2) 
{
	for (int i = 0; i < (*e1).len; i++)
	{
		if (i < (*e2).len) {}
		else
		{
			return 1;
		}
		if ((*e1).e[i] > (*e2).e[i]) return 1;
		else if ((*e1).e[i] < (*e2).e[i]) return -1;
	}
	return 0;
}

char order_lexical_array_char(array_char *e1, array_char *e2) 
{
	for (int i = 0; i < (*e1).len; i++)
	{
		if (i < (*e2).len) {}
		else
		{
			return 1;
		}
		if ((*e1).e[i] > (*e2).e[i]) return 1;
		else if ((*e1).e[i] < (*e2).e[i]) return -1;
	}
	return 0;
}

char order_induced_voidstar_lexical_array_double(void *a, void *b)
{
	return order_lexical_array_double(a, b);
}

char order_induced_voidstar_lexical_array_int(void *a, void *b)
{
	return order_lexical_array_int(a, b);
}

char order_induced_voidstar_lexical_array_char(void *a, void *b)
{
	return order_lexical_array_char(a, b);
}

char has_substring(char *str, char *sub)
{
	int len_str = str_len(str);
	int len_sub = str_len(sub);
	if (len_sub < len_str)
	{
		for (int i = 0; i < len_str - len_sub; i++)
		{
			char match = 1;
			for (int ii = 0; ii < len_sub; ii++)
			{
				if (str[i + ii] == sub[ii]) {}
				else
				{
					match = 0;
					break;
				}
			}
			if (!match) {}
			else return 1;
		}
	}
	return 0;
}

void array_double_diff(const double *e1, const double *e2, double *e3, int len)
{
	for (int i = 0; i < len; i++)
	{
		e3[i] = e1[i] - e2[i];
	}
}

double array_double_dot(const double *e1, const double *e2, int len)
{
	double dp = 0;
	for (int i = 0; i < len; i++) dp += e1[i] * e2[i];
	return dp;
}

double array_double_norm(double *e1, int len)
{
	double normsq = array_double_dot(e1, e1, len);
	return sqrt(normsq);
}

void array_bit_int_init_zero(array_bit *abit, int length)
{
	array_int *data = (array_int *) calloc(1, sizeof(array_int));
       	int N_ints = (length >> LG_BLOCK_SIZE_INT) + 1;
	array_int_init(data, N_ints);
	(*data).len = N_ints;
	for (int i = 0; i < N_ints; i++) (*data).e[i] = 0;
	(*abit).data = (void *) data;
	(*abit).len = length;
	(*abit).mem = ((*data).mem << LG_BLOCK_SIZE_INT);
}

// Methods for arrays of bits
void array_bit_int_init(array_bit *abit, int length)
{
	array_int *data = (array_int *) calloc(1, sizeof(array_int));
	int n_blocks = (length >> LG_BLOCK_SIZE_INT) + 1;
	array_int_init(data, n_blocks);
	(*data).len = n_blocks;
	(*abit).data = (void *) data;
	(*abit).len = 0;
	(*abit).mem = ((*data).mem << LG_BLOCK_SIZE_INT);
}

void add_mem_array_bit_int(array_bit *abit)
{
	add_mem_array_int((array_int *) (*abit).data);
	(*abit).mem <<= 1;
}

void add2array_bit_int(array_bit *abit, char bit)
{
	if ((*abit).len != (*abit).mem) {}
	else add_mem_array_bit_int(abit);
	int addr = ((*abit).len >> array_bit_int_lg_block_size);
	int bit_addr = (*abit).len & array_bit_int_addr_mask;
	array_int *data = (array_int *) (*abit).data;
	if (bit == 0) (*data).e[addr] &= array_bit_int_cmasks[bit_addr];
	else (*data).e[addr] |= array_bit_int_masks[bit_addr];
	(*abit).len += 1;
}

void array_bit_int_set(array_bit *abit, int bit_pos, char bit_val)
{
	if ((*abit).mem > bit_pos) {}
	else
	{
		while ((*abit).mem <= bit_pos) add_mem_array_bit_int(abit);
	}
	int addr = (bit_pos >> array_bit_int_lg_block_size);
	int bit_addr = bit_pos & array_bit_int_addr_mask;
	array_int *data = (array_int *) (*abit).data;
	if (bit_val == 0) (*data).e[addr] &= array_bit_int_cmasks[bit_addr];
	else (*data).e[addr] |= array_bit_int_masks[bit_addr];
}

char array_bit_int_get(array_bit *abit, int bit_pos)
{
	int addr = (bit_pos >> array_bit_int_lg_block_size);
	int bit_addr = bit_pos & array_bit_int_addr_mask;
	array_int *data = (array_int *) (*abit).data;
	return ((*data).e[addr] & array_bit_int_masks[bit_addr]) != 0;
}

void array_bit_int_ll(array_bit *abit)
{
	// This is a bit involved...
	array_int *data = (array_int *) (*abit).data;
	int i = (*data).len;
	int im1 = (*data).len - 1;
	while (im1 > 0)
	{
		(*data).e[im1] <<= 1;
		i = im1;
		im1 -= 1;
		// Load the last bit from (*data).e[im1] into the register of (*data).e[i]
		(*data).e[i] |= ((*data).e[im1] & array_bit_int_masks[31]) != 0;
	}
	(*data).e[0] <<= 1;
}

void array_bit_int_rr(array_bit *abit)
{
	array_int *data = (array_int *) (*abit).data;
	int im1 = 0;
	for (int i = 1; i < (*data).len; i++)
	{
		(*data).e[im1] >>= 1;
		if (((*data).e[i] & 1) == 0) {}
		else
		{
			(*data).e[im1] |= array_bit_int_masks[31];
		}
		im1 = i;
	}
}

void array_bit_int_and(array_bit *abit1, array_bit *abit2, array_bit *abit3)
{
	int len = (*abit1).len > (*abit2).len? (*abit2).len: (*abit1).len;
	if ((*abit3).mem >= len) {}
	else while ((*abit3).mem < len) add_mem_array_bit_int(abit3);
	array_int *data1 = (array_int *) (*abit1).data;
	array_int *data2 = (array_int *) (*abit2).data;
	array_int *data3 = (array_int *) (*abit3).data;
	for (int i = 0; i < len; i++)
	{
		(*data3).e[i] = (*data1).e[i] & (*data2).e[i];
	}
	(*abit3).len = len;
}

void array_bit_int_or(array_bit *abit1, array_bit *abit2, array_bit *abit3)
{
	int len = (*abit1).len > (*abit2).len? (*abit2).len: (*abit1).len;
	if ((*abit3).mem >= len) {}
	else while ((*abit3).mem < len) add_mem_array_bit_int(abit3);
	array_int *data1 = (array_int *) (*abit1).data;
	array_int *data2 = (array_int *) (*abit2).data;
	array_int *data3 = (array_int *) (*abit3).data;
	for (int i = 0; i < len; i++)
	{
		(*data3).e[i] = (*data1).e[i] | (*data2).e[i];
	}
	(*abit3).len = len;

}

void free_array_bit_int(array_bit *abit)
{
	array_int *data = (array_int *) (*abit).data;
	free_array_int(data);
	free(data);
}

void display_array_bit_int(array_bit *abit)
{
	array_int *data = (array_int *) (*abit).data;
	int i = (*abit).len;
	do
	{
		i -= 1;
		char val = array_bit_int_get(abit, i);
		printf("%d", val);
	} while (i > 0);
	printf("\n");
}

char aarray_int_contains_set(aarray_int *aa, int *set, int len)
{
	for (int i = 0; i < (*aa).len; i++)
	{
		if ((*aa).e[i].len == len)
		{
			char equal = 1;
			for (int ii = 0; ii < len; ii++)
			{
				char match = 0;
				for (int iii = 0; iii < len; iii++)
				{
					if (set[ii] != (*aa).e[i].e[iii]) {}
					else 
					{
						match = 1;
						break;
					}
				}
				if (!match) 
				{
					equal = 0;
					break;
				}
			}
			if (!equal) {}
			else return 1;
		}
	}
	return 0;
}

// RESUME: Test this!
int aarray_int_sorted_place(aarray_int *aa, int *seq, int len, char (*order_gt)(void *, void *))
{
	int llim = 0;
	int ulim = (*aa).len - 1;
	int diff = (*aa).len;
	int mid = (llim + ulim) / 2;
	void *oseq[2];
	oseq[0] = (void *) seq;
	oseq[1] = (void *) &len;
	void *lseq[2];
	void *useq[2];
	void *mseq[2];
	char changed = 1;
	while (changed)
	{
		lseq[0] = (void *) (*aa).e[llim].e;
		lseq[1] = (void *) &((*aa).e[llim].len);
		useq[0] = (void *) (*aa).e[ulim].e;
		useq[1] = (void *) &((*aa).e[ulim].len);
		mseq[0] = (void *) (*aa).e[mid].e;
		mseq[1] = (void *) &((*aa).e[mid].len);
		char mid_comp = order_gt((void *) mseq, (void *) oseq);
		if (mid_comp == 1)
		{
			llim = mid;
		}
		else if (mid_comp == -1)
		{
			ulim = mid;
		}
		else
		{
			return mid;
		}
		int nmid = (llim + ulim) / 2;
		if (nmid != mid)
		{
			mid = nmid;
		}
		else
		{
			return llim;
		}
	}
}

char aarray_int_sorted_contains(aarray_int *aa, int *seq, int len, char (*order)(void *, void *))
{
	int pos = aarray_int_sorted_place(aa, seq, len, order);
	if ((*aa).e[pos].len == len) {}
	else
	{
		return 0;
	}
	for (int i = 0; i < len; i++)
	{
		if ((*aa).e[pos].e[i] == seq[i]) {}
		else
		{
			return 0;
		}
	}
	return 1;
}

char aarray_int_contains(aarray_int *aa, int *seq, int len)
{
	for (int i = 0; i < (*aa).len; i++)
	{
		if ((*aa).e[i].len != len) {}
		else
		{
			char match = 1;
			for (int ii = 0; ii < len; ii++)
			{
				if ((*aa).e[i].e[ii] == seq[ii]) {}
				else 
				{
					match = 0;
					break;
				}
			}
			if (!match) {}
			else
			{
				return 1;
			}
		}
	}
	return 0;
}

int array_int_min(int *a, int len)
{
	int min_a = a[0];
	for (int i = 1; i < len; i++)
	{
		if (a[i] >= min_a) {}
		else
		{
			min_a = a[i];
		}
	}
	return min_a;
}

int array_int_max(int *a, int len)
{
	int max_a = a[0];
	for (int i = 1; i < len; i++)
	{
		if (a[i] <= max_a) {}
		else
		{
			max_a = a[i];
		}
	}
	return max_a;
}

int *parse_int(void *a)
{
	return ((int *) a);
}
char *parse_char(void *a)
{
	return ((char *) a);
}
double *parse_double(void *a)
{
	return ((double *) a);
}

array_int *parse_array_int(void *a)
{
	return ((array_int *) a);	
}
array_char *parse_array_char(void *a)
{
	return ((array_char *) a);
}
array_double *parse_array_double(void *a)
{
	return ((array_double *) a);
}
array_voidstar *parse_array_voidstar(void *a)
{
	return ((array_voidstar *) a);
}
aarray_int *parse_aarray_int(void *a)
{
	return ((aarray_int *) a);
}
aarray_char *parse_aarray_char(void *a)
{
	return ((aarray_char *) a);
}
aarray_double *parse_aarray_double(void *a)
{
	return ((aarray_double *) a);
}
nbrlist *parse_nbrlist(void *a)
{
	return ((nbrlist *) a);
}

double array_double_min(double *a, int len)
{
	double min_a = a[0];
	for (int i = 1; i < len; i++)
	{
		if (a[i] >= min_a) {}
		else
		{
			min_a = a[i];
		}
	}
	return min_a;
}

double array_double_max(double *a, int len)
{
	double max_a = a[0];
	for (int i = 1; i < len; i++)
	{
		if (a[i] <= max_a) {}
		else
		{
			max_a = a[i];
		}
	}
	return max_a;
}

char array_char_min(char *a, int len)
{
	char min_a = a[0];
	for (int i = 1; i < len; i++)
	{
		if (a[i] >= min_a) {}
		else
		{
			min_a = a[i];
		}
	}
	return min_a;
}

char array_char_max(char *a, int len)
{
	char max_a = a[0];
	for (int i = 1; i < len; i++)
	{
		if (a[i] <= max_a) {}
		else
		{
			max_a = a[i];
		}
	}
	return max_a;
}

// This might be easier with 'cyclable' arrays
void merge_array_int_permutation(array_int *a, array_int *p, int *buf, int i0, int i1, int mdpt)
{
	if ((*a).e[(*p).e[mdpt]] >= (*a).e[(*p).e[mdpt - 1]]) 
	{
		return;
	}
	else if ((*a).e[(*p).e[i0]] >= (*a).e[(*p).e[i1 - 1]])
	{
		// case mdpt - i0 > i1 - mdpt:
		// p[i0 + 1: mdpt] -> p[mdpt : i1] -> p[i0 : mdpt - 1] 
		// p[i0] -> p[mdpt - 1]
		// case mdpt - i0 == i1 - mdpt:
		// p[i0 : mdpt] <-> p[mdpt : i1]
		int mdptmi0 = mdpt - i0;
		int i1mmdpt = i1 - mdpt;
		int *seg_a = &(buf[i0]);
		int *seg_b = &((*p).e[mdpt]);
		for (int i = 0; i < i1mmdpt; i++)
		{
			//int ipi0 = i + i0;
			seg_a[i] = seg_b[i];
			//buf[ipi0] = (*p).e[i + mdpt];
		}
		seg_a = &(buf[i1 - mdptmi0]);
		seg_b = &((*p).e[i0]);
		for (int i = 0; i < mdptmi0; i++)
		{
			seg_a[i] = seg_b[i];
		}
		for (int i = i0; i < i1; i++) (*p).e[i] = buf[i];
		return;
	}
	int term = mdpt;
	int c0 = i0;
	int c1 = mdpt;
	int c = i0;
	while ((c0 < mdpt) && (c1 < i1))
	{
		
		if ((*a).e[(*p).e[c0]] <= (*a).e[(*p).e[c1]])
		{
			buf[c] = (*p).e[c0];
			c0 += 1;
			c += 1;
		}
		else
		{
			buf[c] = (*p).e[c1];
			c1 += 1;
			c += 1;
		}
	}
	if (c0 < mdpt)
	{
		while (c < i1)
		{
			buf[c] = (*p).e[c0];
			c += 1;
			c0 += 1;
		}
	}
	if (c1 < i1)
	{
		if (c == c1) {}
		else
		{
			printf("Something weird happened in merge_array_int_permutation! left half of array folded in, but c = %d != c1 = %d!\n", c, c1);
			exit(EXIT_FAILURE);
		}
		while (c < i1)
		{
			buf[c] = (*p).e[c1];
			c += 1;
			c1 += 1;
		}
	}
	for (int i = i0; i < i1; i++)
	{
		(*p).e[i] = buf[i];
	}
}

void merge_array_double_permutation(array_double *a, array_int *p, int *buf, int i0, int i1, int mdpt)
{
	if (i0 < mdpt && mdpt < i1) {}
	else return;
	if ((*a).e[(*p).e[mdpt]] >= (*a).e[(*p).e[mdpt - 1]]) 
	{
		return;
	}
	else if ((*a).e[(*p).e[i0]] >= (*a).e[(*p).e[i1 - 1]])
	{
		// case mdpt - i0 > i1 - mdpt:
		// p[i0 + 1: mdpt] -> p[mdpt : i1] -> p[i0 : mdpt - 1] 
		// p[i0] -> p[mdpt - 1]
		// case mdpt - i0 == i1 - mdpt:
		// p[i0 : mdpt] <-> p[mdpt : i1]
		int mdptmi0 = mdpt - i0;
		int i1mmdpt = i1 - mdpt;
		int *seg_a = &(buf[i0]);
		int *seg_b = &((*p).e[mdpt]);
		for (int i = 0; i < i1mmdpt; i++)
		{
			//int ipi0 = i + i0;
			seg_a[i] = seg_b[i];
			//buf[ipi0] = (*p).e[i + mdpt];
		}
		seg_a = &(buf[i1 - mdptmi0]);
		seg_b = &((*p).e[i0]);
		for (int i = 0; i < mdptmi0; i++)
		{
			seg_a[i] = seg_b[i];
			// buf[i + i1mmdpt] = (*p).e[i + i0];
		}
		for (int i = i0; i < i1; i++) (*p).e[i] = buf[i];
		return;
	}
	int term = mdpt;
	int c0 = i0;
	int c1 = mdpt;
	int c = i0;
	while ((c0 < mdpt) && (c1 < i1))
	{
		if ((*a).e[(*p).e[c0]] <= (*a).e[(*p).e[c1]])
		{
			buf[c] = (*p).e[c0];
			c0 += 1;
			c += 1;
		}
		else
		{
			buf[c] = (*p).e[c1];
			c1 += 1;
			c += 1;
		}
	}
	if (c0 < mdpt)
	{
		while (c < i1)
		{
			buf[c] = (*p).e[c0];
			c += 1;
			c0 += 1;
		}
	}
	if (c1 < i1)
	{
		if (c == c1) {}
		else
		{
			printf("Something weird happened in merge_array_int_permutation! left half of array folded in, but c = %d != c1 = %d!\n", c, c1);
			exit(EXIT_FAILURE);
		}
		while (c < i1)
		{
			buf[c] = (*p).e[c1];
			c += 1;
			c1 += 1;
		}
	}
	for (int i = i0; i < i1; i++)
	{
		(*p).e[i] = buf[i];
	}

}

void merge_sort_array_double_permutation(array_double *a, array_int *p, int *buf, int i0, int i1)
{
	int i1mi0 = i1 - i0;
	if (i1mi0 > 2)
	{
		int mdpt = (i0 + i1 + 1) >> 1; 
		merge_sort_array_double_permutation(a, p, buf, i0, mdpt);
		merge_sort_array_double_permutation(a, p, buf, mdpt, i1);
		merge_array_double_permutation(a, p, buf, i0, i1, mdpt);
	}
	else
	{
		if (i1mi0 < 2) return;
		else
		{
			int i0p1 = i0 + 1;
			if ((*a).e[(*p).e[i0]] < (*a).e[(*p).e[i0p1]]) 
			{
				buf[i0] = (*p).e[i0];
				buf[i0p1] = (*p).e[i0p1];
			}
			else
			{
				buf[i0] = (*p).e[i0p1];
				buf[i0p1] = (*p).e[i0];
			}
			(*p).e[i0] = buf[i0];
			(*p).e[i0p1] = buf[i0p1];
		}
	}
}

void sort_array_double_permutation(array_double *a, array_int *p)
{
	if ((*a).len == (*p).len) {}
	else
	{
		printf("Error (sort_array_double_permutation): permutation must be pre-allocated\n");
		exit(EXIT_FAILURE);
	}
	int *buf = (int *) calloc((*p).len, sizeof(int));
	merge_sort_array_double_permutation(a, p, buf, 0, (*a).len);
	free(buf);
}

void merge_sort_array_int_permutation(array_int *a, array_int *p, int *buf, int i0, int i1)
{
	int i1mi0 = i1 - i0;
	if (i1mi0 > 2)
	{
		int mdpt = (i0 + i1 + 1) >> 1; 
		merge_sort_array_int_permutation(a, p, buf, i0, mdpt);
		merge_sort_array_int_permutation(a, p, buf, mdpt, i1);
		merge_array_int_permutation(a, p, buf, i0, i1, mdpt);
	}
	else
	{
		if (i1mi0 < 2) return;
		else
		{
			int i0p1 = i0 + 1;
			if ((*a).e[(*p).e[i0]] <= (*a).e[(*p).e[i0p1]]) 
			{
				buf[i0] = (*p).e[i0];
				buf[i0p1] = (*p).e[i0p1];
			}
			else
			{
				buf[i0] = (*p).e[i0p1];
				buf[i0p1] = (*p).e[i0];
			}
			(*p).e[i0] = buf[i0];
			(*p).e[i0p1] = buf[i0p1];
		}
	}
}

void sort_array_int_permutation(array_int *a, array_int *p)
{
	if ((*a).len == (*p).len) {}
	else
	{
		array_int_init(p, (*a).len);
		for (int i = 0; i < (*a).len; i++) (*p).e[i] = i;
		(*p).len = (*a).len;
	}
	int *buf = (int *) calloc((*p).len, sizeof(int));
	merge_sort_array_int_permutation(a, p, buf, 0, (*a).len);
	free(buf);
}

void load_nbrlist(nbrlist *nbl, char *ifname)
{
	FILE *ifile = fopen(ifname, "r");
	if (ifile != NULL)
	{
		while (1)
		{
			int n_nbrs;
			int status = fscanf(ifile, "%d", &n_nbrs);
			if (status != EOF)
			{
				int i = (*nbl).v.len;
				extend_nbrlist(nbl);
				for (int ni = 0; ni < n_nbrs; ni++)
				{
					int ii;
					fscanf(ifile, "%d", &ii);
					if (i > ii) add_edge_nbrlist(nbl, ii, i);
				}
			}
			else break;
		}
		fclose(ifile);
	}
}

void hash_table_int_init(hash_table_int *ht, int min_data_size, basics_Uint64 min_elem_size)
{
	(*ht).largest_bin = 0;
	(*ht).lg_data_len = 0;
	int data_len = 1;
	while (data_len < min_data_size)
	{
		data_len <<= 1;
		(*ht).lg_data_len += 1;
	}
	array_voidstar_init(&((*ht).data), data_len);
	(*ht).data.len = data_len;
	if ((*ht).data.len != (*ht).data.mem) 
	{
		printf("Error (hash_table_int_init): something weird happened! c.f. %d vs. %d a;sldjf;lsjkdf;jaw\n", (*ht).data.len, (*ht).data.mem);
		exit(EXIT_FAILURE);
	}
	(*ht).size_mask = (*ht).data.len - 1;
	int lg_elem_size = (*ht).lg_data_len + 2;
	lg_elem_size = lg_elem_size < 65 ? lg_elem_size : 64;
	basics_Uint64 elem_size = 1 << lg_elem_size;
	while (elem_size < min_elem_size)
	{
		elem_size <<= 1;
		lg_elem_size += 1;
	}
	(*ht).elem_mask = elem_size - 1;
	//(*ht).elem_mask = elem_size - 1;
	(*ht).gen_a = rand() & (*ht).elem_mask;
	//(*ht).gen_b = rand() & ((*ht).elem_mask);
	(*ht).gen_b = rand() & (*ht).elem_mask;
	(*ht).gen_b |= 1;
	(*ht).gen_a |= 1;
	array_int_init(&((*ht).addr_0), 0);
	array_int_init(&((*ht).addr_1), 0);
	array_int_init(&((*ht).elem), 0);
	array_voidstar_init(&((*ht).values), 0);
	(*ht).n_bits_ignored = lg_elem_size - (*ht).lg_data_len;
}

// RESUME
void free_hash_table_int(hash_table_int *ht, void (*free_value_func)(void *))
{
	free_array_int(&((*ht).addr_0));
	free_array_int(&((*ht).addr_1));
	free_array_int(&((*ht).elem));
	if (free_value_func != NULL)
	{
		for (int i = 0; i < (*ht).values.len; i++) free_value_func((*ht).values.e[i]);
	}
	free_array_voidstar(&((*ht).values), NULL);
	for (int i = 0; i < (*ht).data.len; i++)
	{
		if ((*ht).data.e[i] != NULL)
		{
			free_array_int((array_int *) (*ht).data.e[i]);
			free((*ht).data.e[i]);
		}
	}
	free_array_voidstar(&((*ht).data), NULL);
}

// Resume: test this!
void resize_hash_table_int(hash_table_int *ht, int new_min_size)
{
	hash_table_int htp;
	int nmsp1 = new_min_size + 1;
	if ((*ht).elem_mask <= nmsp1)
	{
		(*ht).elem_mask += 1;
		while ((*ht).elem_mask <= nmsp1)
		{
			(*ht).elem_mask <<= 1;
		}
		(*ht).elem_mask -= 1;
	}
	hash_table_int_init(&htp, new_min_size, (*ht).elem_mask);
	for (int i = 0; i < (*ht).elem.len; i++)
	{
		add2hash_table_int(&htp, (*ht).elem.e[i], (*ht).values.e[i]);
	}
	free_hash_table_int(ht, NULL);
	(*ht) = htp;
}

void transcribe_hash_table_int(hash_table_int *src, hash_table_int *dest)
{
	array_voidstar_init(&((*dest).data), (*src).data.len);
	(*dest).data.len = (*src).data.len;
	for (int i = 0; i < (*src).data.len; i++)
	{
		if ((*src).data.e[i] != NULL)
		{
			array_int *sei = (array_int *) (*src).data.e[i];
			array_int *dei = (array_int *) calloc(1, sizeof(array_int));
			transcribe_array_int(sei, dei);
			(*dest).data.e[i] = dei;
		}
	}
	transcribe_array_int(&((*src).addr_0), &((*dest).addr_0));
	transcribe_array_int(&((*src).addr_1), &((*dest).addr_1));
	transcribe_array_int(&((*src).elem), &((*dest).elem));
	transcribe_array_voidstar(&((*src).values), &((*dest).values));
	(*dest).n_bits_ignored = (*src).n_bits_ignored;
	(*dest).elem_mask = (*src).elem_mask;
	(*dest).lg_elem_size = (*src).lg_elem_size;
	(*dest).size_mask = (*src).size_mask;
	(*dest).gen_a = (*src).gen_a;
	(*dest).gen_b = (*src).gen_b;
	(*dest).largest_bin = (*src).largest_bin;
}

int hash_table_int_map(hash_table_int *ht, int n)
{
	return (((*ht).gen_a * n + (*ht).gen_b) & (*ht).elem_mask) >> (*ht).n_bits_ignored;
}

char query_hash_table_int(hash_table_int *ht, int n, int *elem_addr, void **val)
{
	int addr_0 = hash_table_int_map(ht, n);
	if ((*ht).data.e[addr_0] != NULL) {}
	else return 0;
	// Check if 'n' already exists in the hash table
	array_int *ehn = (array_int *) (*ht).data.e[addr_0];
	char present = 0;
	for (int i = 0; i < (*ehn).len; i++)
	{
		if ((*ht).elem.e[(*ehn).e[i]] != n) {}
		else
		{
			(*elem_addr) = (*ehn).e[i];
			present = 1;
			if (val != NULL) (*val) = (*ht).values.e[(*ehn).e[i]];
			break;
		}
	}
	return present;
}

// RESUME: include actual nontrivial data in both hash table implementations
void add2hash_table_int(hash_table_int *ht, int n, void *data)
{
	if ((*ht).elem.len < ((*ht).data.len << 1)) {}
	else
	{
		// resize the hash table
		resize_hash_table_int(ht, ((*ht).data.len << 1));
	}
	int hash_n = hash_table_int_map(ht, n);
	if ((*ht).data.e[hash_n] != NULL) {}
	else
	{
		(*ht).data.e[hash_n] = (array_int *) calloc(1, sizeof(array_int));
		array_int_init((array_int *) (*ht).data.e[hash_n], 1);
	}
	// Check if 'n' already exists in the hash table
	array_int *ehn = (array_int *) (*ht).data.e[hash_n];
	char present = 0;
	for (int i = 0; i < (*ehn).len; i++)
	{
		if ((*ht).elem.e[(*ehn).e[i]] != n) {}
		else
		{
			present = 1;
			break;
		}
	}
	if (present) {}
	else
	{
		add2array_int(&((*ht).addr_0), hash_n);
		add2array_int(&((*ht).addr_1), (*ehn).len);
		add2array_int(ehn, (*ht).elem.len);
		(*ht).largest_bin = (*ehn).len < (*ht).largest_bin ? (*ht).largest_bin : (*ehn).len;
		add2array_int(&((*ht).elem), n);
		add2array_voidstar(&((*ht).values), data);
	}
}

void remove_hash_table_int_exp(hash_table_int *ht, int elem_addr, void (*free_value_func)(void *))
{
	int addr_0 = (*ht).addr_0.e[elem_addr];
	array_int *ehn = (array_int *) (*ht).data.e[addr_0];
	int addr_1 = (*ht).addr_1.e[elem_addr];
	remove_array_int(ehn, addr_1);
	if ((*ehn).len > addr_1)
	{
		int index_p = (*ehn).e[addr_1];
		int n_p = (*ht).elem.e[index_p];
		(*ht).addr_1.e[index_p] = addr_1;
	}
	else if ((*ehn).len > 0) {}
	else
	{
		free_array_int(ehn);
		free((*ht).data.e[addr_0]);
		(*ht).data.e[addr_0] = NULL;
	}
	remove_array_int(&((*ht).elem), elem_addr);
	remove_array_voidstar(&((*ht).values), elem_addr, free_value_func);
	remove_array_int(&((*ht).addr_0), elem_addr);
	remove_array_int(&((*ht).addr_1), elem_addr);
	if ((*ht).elem.len > elem_addr)
	{
		int n_p = (*ht).elem.e[elem_addr];
		int hash_n_p = (*ht).addr_0.e[elem_addr];
		array_int *ehn_p = (array_int *) (*ht).data.e[hash_n_p];
		(*ehn_p).e[(*ht).addr_1.e[elem_addr]] = elem_addr;
	}
	if ((*ht).elem.len < ((*ht).data.len >> 2)) resize_hash_table_int(ht, ((*ht).data.len >> 1) - 1);
}

char remove_hash_table_int(hash_table_int *ht, int n, void (*free_value_func)(void *))
{
	void *val_addr;
	int elem_addr;
	char present = query_hash_table_int(ht, n, &elem_addr, &val_addr);
	if (present) remove_hash_table_int_exp(ht, elem_addr, free_value_func); 
	return present;
}

void transcribe_hash_table_int_str(hash_table_int_str *src, hash_table_int_str *dest, char copy_mode)
{
	(*dest).key_data_src = copy_mode == BASICS_H_OTHER ? (*src).key_data_src : &((*dest).elem);
	(*dest).largest_bin = (*src).largest_bin;
	(*dest).lg_data_len = (*src).lg_data_len;
	array_voidstar_init(&((*dest).data), (*src).data.len);
	array_voidstar_init(&((*dest).values), (*src).values.len);
	(*dest).values.len = (*src).values.len;
	(*dest).data.len = (*src).data.len;
	for (int i = 0; i < (*src).data.len; i++)
	{
		if ((*src).data.e[i] != NULL)
		{
			(*dest).data.e[i] = (array_int *) calloc(1, sizeof(array_int));
			transcribe_array_int((array_int *) (*src).data.e[i], (array_int *) (*dest).data.e[i]);
		}
	}
	(*dest).n_bits_ignored = (*src).n_bits_ignored;
	(*dest).elem_mask = (*src).elem_mask;
	(*dest).size_mask = (*src).size_mask;
	(*dest).gen_a = (*src).gen_a;
	(*dest).gen_b = (*src).gen_b;
	(*dest).str_gen_a = (*src).str_gen_a;
	(*dest).str_gen_b = (*src).str_gen_b;
	transcribe_array_int(&((*src).addr_0), &((*dest).addr_0));
	transcribe_array_int(&((*src).addr_1), &((*dest).addr_1));
	transcribe_array_int(&((*src).lens), &((*dest).lens));
	transcribe_array_voidstar(&((*src).values), &((*dest).values));
	if (copy_mode == BASICS_H_OTHER) transcribe_array_voidstar(&((*src).elem), &((*dest).elem));
	else
	{
		array_voidstar_init(&((*dest).elem), (*src).elem.len);
		(*dest).elem.len = (*src).elem.len;
		for (int i = 0; i < (*src).elem.len; i++)
		{
			int *ei = (int *) (*src).elem.e[i];
			int *ei_ = (int *) calloc((*src).lens.e[i], sizeof(int));
			for (int ii = 0; ii < (*src).lens.e[i]; ii++)
			{
				ei_[ii] = ei[ii];
			}
			(*dest).elem.e[i] = (void *) ei_;
		}
	}
}

int int_str_ht_size(hash_table_int_str *ht)
{
	return (*ht).data.len;
}

void hash_table_int_str_init(hash_table_int_str *ht, int min_ht_size, basics_Uint64 max_int_size, char src_mode)
{
	(*ht).key_data_src = &((*ht).elem);
	(*ht).largest_bin = 0;
	(*ht).lg_data_len = 0;
	int data_len = 1;
	while (data_len < min_ht_size)
	{
		data_len <<= 1;
		(*ht).lg_data_len += 1;
	}
	array_voidstar_init(&((*ht).data), data_len);
	(*ht).data.len = data_len;
	(*ht).n_bits_ignored = (*ht).lg_data_len + 2;
	(*ht).elem_mask = 1 << (*ht).n_bits_ignored;
	while ((*ht).elem_mask < max_int_size)
	{
		(*ht).n_bits_ignored += 1;
		(*ht).elem_mask <<= 1;
	}
	(*ht).n_bits_ignored -= (*ht).lg_data_len;
	(*ht).elem_mask -= 1;
	(*ht).size_mask = data_len - 1;
	(*ht).gen_a = rand() & (*ht).elem_mask;
	(*ht).gen_b = rand() & (((*ht).elem_mask << 1) + 1);
	(*ht).str_gen_a = rand() & (*ht).elem_mask;
	(*ht).str_gen_b = rand() & (((*ht).elem_mask << 1) + 1);
	array_int_init(&((*ht).addr_0), 0);
	array_int_init(&((*ht).addr_1), 0);
	array_int_init(&((*ht).lens), 0);
	array_voidstar_init(&((*ht).elem), 0);
	array_voidstar_init(&((*ht).values), 0);
}

void free_hash_table_int_str(hash_table_int_str *ht, void (*free_value_func)(void *))
{
	if ((*ht).key_data_src == &((*ht).elem))
	{
		for (int i = 0; i < (*ht).elem.len; i++)
		{
			free((int *) (*ht).elem.e[i]);
		}
	}
	free_array_voidstar(&((*ht).elem), NULL);
	if (free_value_func != NULL)
	{
		for (int i = 0; i < (*ht).values.len; i++) free_value_func((*ht).values.e[i]);
	}
	free_array_voidstar(&((*ht).values), NULL);
	free_array_int(&((*ht).addr_0));
	free_array_int(&((*ht).addr_1));
	free_array_int(&((*ht).lens));
	for (int i = 0; i < (*ht).data.len; i++)
	{
		if ((*ht).data.e[i] != NULL)
		{
			array_int *ei = (array_int *) (*ht).data.e[i];
			free_array_int(ei);
			free(ei);
		}
	}
	free_array_voidstar(&((*ht).data), NULL);
}

char hash_table_int_str_self_srcd(hash_table_int_str *ht)
{
	return (*ht).key_data_src == &((*ht).elem);
}

void hash_table_int_str_transfer_exp(hash_table_int_str *ht0, int elem_addr, hash_table_int_str *ht1)
{
	void *src0 = (*ht0).key_data_src;
	void *src1 = (*ht1).key_data_src;
	(*ht0).key_data_src = NULL;
	(*ht1).key_data_src = NULL;
	add2hash_table_int_str(ht1, (int *) (*ht0).elem.e[elem_addr], (*ht0).lens.e[elem_addr], (*ht0).values.e[elem_addr]);
	remove_hash_table_int_str_exp(ht0, elem_addr, NULL);
	(*ht0).key_data_src = src0;
	(*ht1).key_data_src = src1;
}

void hash_table_int_str_transfer(hash_table_int_str *ht0, int *key, int key_len, hash_table_int_str *ht1)
{
	void *val;
	int elem_addr;
	char present = query_hash_table_int_str(ht0, key, key_len, &elem_addr, &val); // NOTE: it might be simpler and possibly more efficient to set the 'elem' index instead of the table indices (addr_0/addr_1).
	if (present) hash_table_int_str_transfer_exp(ht0, elem_addr, ht1); // RESUME: change this to refer to elem_addr
}

// NOTE: consider storing actual keys in a separate, aggregated/global hash table/set
// and storing a simple reference in local instances. The advantage of this approach 
// is that it would save memory with a relatively small cost in performance. Also,
// performance losses could be offset to some extent (especially with larger strings)
// by eliminating the need to allocate largish arrays.
void resize_hash_table_int_str(hash_table_int_str *ht, int new_min_size)
{
	hash_table_int_str htp;
	char copy_mode = (*ht).key_data_src == &((*ht).elem) ? BASICS_H_SELF : BASICS_H_OTHER;
	hash_table_int_str_init(&htp, new_min_size, (*ht).elem_mask, BASICS_H_OTHER);
	htp.key_data_src = (*ht).key_data_src;
	//if (copy_mode == BASICS_H_OTHER) htp.key_data_src = (*ht).key_data_src;
	for (int i = 0; i < (*ht).elem.len; i++)
	{
		int *ei = (int *) (*ht).elem.e[i];
		int len_ei = (*ht).lens.e[i];
		add2hash_table_int_str(&htp, ei, len_ei, (*ht).values.e[i]);
	}
	if (copy_mode == BASICS_H_SELF) 
	{
		(*ht).key_data_src = NULL;
	}
	free_hash_table_int_str(ht, NULL);
	(*ht) = htp;
	if (copy_mode == BASICS_H_SELF)
	{
		(*ht).key_data_src = &((*ht).elem);
	}
}

void add2hash_table_int_str(hash_table_int_str *ht, int *n, int len, void *value)
{
	if ((*ht).elem.len < ((*ht).data.len << 1)) {}
	else
	{
		resize_hash_table_int_str(ht, (*ht).elem.len); // this should increase the size of the data field to the next power of two greater than data.len
	}
	int hash_n = hash_table_int_str_map(ht, n, len);
	char present = 0;
	if ((*ht).data.e[hash_n] != NULL) 
	{
		array_int *ehn = (array_int *) (*ht).data.e[hash_n];
		for (int i = 0; i < (*ehn).len; i++)
		{
			int index = (*ehn).e[i];
			int *si = (int *) (*ht).elem.e[index];
			int si_len = (*ht).lens.e[index];
			if (si_len != len) {}
			else
			{
				if (arr_int_comp(si, n, len) != 0) {}
				else
				{
					present = 1;
					break;
				}
			}
		}
	}
	else
	{
		(*ht).data.e[hash_n] = (array_int *) calloc(1, sizeof(array_int));
		array_int_init((array_int *) (*ht).data.e[hash_n], 1);
	}
	if (present) {}
	else
	{
		array_int *ehn = (array_int *) (*ht).data.e[hash_n];
		add2array_int(&((*ht).addr_0), hash_n);
		add2array_int(&((*ht).addr_1), (*ehn).len);
		add2array_int(ehn, (*ht).elem.len);
		(*ht).largest_bin = (*ht).largest_bin > (*ehn).len ? (*ht).largest_bin : (*ehn).len;
		if ((*ht).key_data_src == &((*ht).elem))
		{
			int *n_ = (int *) calloc(len, sizeof(int));
			for (int i = 0; i < len; i++) n_[i] = n[i];
			add2array_voidstar(&((*ht).elem), n_);
		}
		else 
		{
			add2array_voidstar(&((*ht).elem), n);
		}
		add2array_voidstar(&((*ht).values), value);
		add2array_int(&((*ht).lens), len);
	}
}

void remove_hash_table_int_str_exp(hash_table_int_str *ht, int index, void (*free_value_func)(void *))
{
	int addr_0 = (*ht).addr_0.e[index];
	int addr_1 = (*ht).addr_1.e[index];
	array_int *ehn = (array_int *) (*ht).data.e[addr_0];
	char check_largest_bin = (*ehn).len == (*ht).largest_bin;
	remove_array_int(ehn, addr_1);
	if ((*ehn).len > addr_1)
	{
		int index_p = (*ehn).e[addr_1];
		(*ht).addr_1.e[index_p] = addr_1;
	}
	else if ((*ehn).len > 0) {}
	else
	{
		free_array_int(ehn);
		free((*ht).data.e[addr_0]);
		(*ht).data.e[addr_0] = NULL;
	}
	remove_array_int(&((*ht).addr_1), index);
	remove_array_int(&((*ht).addr_0), index);
	remove_array_int(&((*ht).lens), index);
	if ((*ht).key_data_src == &((*ht).elem))
	{
		free((int *) (*ht).elem.e[index]);
	}
	remove_array_voidstar(&((*ht).elem), index, NULL);
	remove_array_voidstar(&((*ht).values), index, free_value_func);
	if ((*ht).addr_0.len > index)
	{
		int hash_n_p = (*ht).addr_0.e[index];
		int addr_1_p = (*ht).addr_1.e[index];
		array_int *ehn_p = (array_int *) (*ht).data.e[hash_n_p];
		(*ehn_p).e[addr_1_p] = index;
	}
}

char remove_hash_table_int_str(hash_table_int_str *ht, int *n, int len, void (*free_value_func)(void *))
{
	void *val_addr;
	int elem_addr;
	int present = query_hash_table_int_str(ht, n, len, &elem_addr, &val_addr);
	if (present)
	{
		remove_hash_table_int_str_exp(ht, elem_addr, free_value_func); // RESUME: update remove_hash_table_int_str to accept 'elem_addr'
		if ((*ht).elem.len < ((*ht).data.len >> 2)) resize_hash_table_int_str(ht, ((*ht).data.len >> 2));
	}
	return present;
}

int arr_int_comp(int *a, int *b, int len)
{
	return len > 0 ? (a[0] > b[0] ? 1 : (b[0] > a[0] ? -1 : arr_int_comp(&a[1], &b[1], len - 1))) : 0;
}

char query_hash_table_int_str(hash_table_int_str *ht, int *n, int len, int *elem_addr, void **val)
{
	int addr0 = hash_table_int_str_map(ht, n, len);
	//(*addr0) = hash_table_int_str_map(ht, n, len);
	//if ((*ht).data.e[(*addr0)] != NULL)
	if ((*ht).data.e[addr0] != NULL)
	{
		//array_int *ehn = (array_int *) (*ht).data.e[(*addr0)];
		array_int *ehn = (array_int *) (*ht).data.e[addr0];
		for (int i = 0; i < (*ehn).len; i++)
		{
			int index = (*ehn).e[i];
			if ((*ht).lens.e[index] == len)
			{
				if (arr_int_comp(n, (int *) (*ht).elem.e[index], len) == 0) 
				{
					//(*addr1) = i;
					(*elem_addr) = (*ehn).e[i];
					if (val != NULL)
					{
						(*val) = (*ht).values.e[index];
					}
					return 1;
				}
			}
		}
	}
	return 0;
}

int hash_table_int_str_map(hash_table_int_str *ht, int *n, int len)
{
	int a = (*ht).str_gen_b;
	for (int i = 0; i < len; i++)
	{
		a = ((*ht).str_gen_a * a + n[i]) & (*ht).elem_mask;
	}
	return ((a * (*ht).gen_a + (*ht).gen_b) & (*ht).elem_mask) >> (*ht).n_bits_ignored;
}

// Union find methods
void union_find_init(union_find *uf, int size)
{
	array_int_init(&((*uf).membership), size);
	array_int_init(&((*uf).cluster_addr), size);
	array_int_init(&((*uf).cluster_size), size);
	aarray_int_init(&((*uf).clusters), size);
	(*uf).membership.len = size;
	(*uf).cluster_size.len = size;
	(*uf).cluster_addr.len = size;
	for (int i = 0; i < size; i++)
	{
		(*uf).membership.e[i] = i;
		extend_aarray_int(&((*uf).clusters));
		(*uf).cluster_addr.e[i] = -1;
		(*uf).cluster_size.e[i] = 1;
	}
}

void free_union_find(union_find *uf)
{
	free_array_int(&((*uf).cluster_addr));
	free_array_int(&((*uf).cluster_size));
	free_array_int(&((*uf).membership));
	free_aarray_int(&((*uf).clusters));
}

void union_find_merge_cluster(union_find *uf, int ci, int cj)
{
	if ((*uf).cluster_addr.e[ci] == -1 && (*uf).cluster_addr.e[cj] == -1) {}
	else
	{
		printf("Error (union_find_merge_cluster): attempting to merge non-terminal points %d->%d->%d, %d->%d->%d\n", (*uf).cluster_addr.e[ci], ci, (*uf).membership.e[ci], (*uf).cluster_addr.e[cj], cj, (*uf).membership.e[cj]);
		exit(EXIT_FAILURE);
	}
	if (ci != cj) {}
	else return;
	(*uf).membership.e[ci] = cj;
	(*uf).cluster_addr.e[ci] = (*uf).clusters.e[cj].len;
	add2array_int(&((*uf).clusters.e[cj]), ci);
}

void union_find_union(union_find *uf, int i, int j)
{
	int ci, cj;
	int depth_ci = union_find_find(uf, i, &ci);
	int depth_cj = union_find_find(uf, j, &cj);
       	if (ci == cj) {}
	else
	{
		int c_upper, c_lower;
		if ((*uf).cluster_size.e[cj] <= (*uf).cluster_size.e[ci]) 
		{
			c_upper = ci;
			c_lower = cj;
		}
		else 
		{
			c_upper = cj;
			c_lower = ci;
		}
		union_find_merge_cluster(uf, c_lower, c_upper);
		(*uf).cluster_size.e[c_upper] += (*uf).cluster_size.e[c_lower];
	}
}

int union_find_find(union_find *uf, int i, int *ci)
{
	(*ci) = (*uf).membership.e[i];
	if ((*ci) == i) return 0;
	int addr_i = (*uf).cluster_addr.e[i];
	if (addr_i > -1) {}
	else
	{
		printf("Something weird happened! a;lsdja;lfa\n");
		exit(EXIT_FAILURE);
	}
	remove_array_int(&((*uf).clusters.e[(*ci)]), addr_i);
	if ((*uf).clusters.e[(*ci)].len > addr_i)
	{
		int j = (*uf).clusters.e[(*ci)].e[addr_i];
		if ((*uf).membership.e[j] != j) {}
		else
		{
			printf("Something weird happened! a;lsjkl;jjoawdk\n");
			exit(EXIT_FAILURE);
		}
		(*uf).cluster_addr.e[j] = addr_i;
	}
	int depth = 1 + union_find_find(uf, (*ci), ci);
	// change the membership of 'i' to (*ci), and add 'i' to the associated cluster
	(*uf).cluster_addr.e[i] = (*uf).clusters.e[(*ci)].len;
	add2array_int(&((*uf).clusters.e[(*ci)]), i);
	return depth;
}

// RESUME: test these!
int union_find_order_root(union_find *uf, int j)
{
	int c = 1;
	for (int i = 0; i < (*uf).clusters.e[j].len; i++)
	{
		c += union_find_order_root(uf, (*uf).clusters.e[j].e[i]);
	}
	return c;
}

int union_find_order(union_find *uf, int i)
{
	int ci;
	union_find_find(uf, i, &ci);
	return (*uf).cluster_size.e[ci];
	//return union_find_order_root(uf, ci);
}

// 'Bare-bones' routines

void union_find_bb_init(array_int *uf, int size)
{
	array_int_init(uf, size);
	(*uf).len = size;
	for (int i = 0; i < size; i++) (*uf).e[i] = i;
}

void union_find_bb_union(array_int *uf, int i, int j)
{
	int ri, rj;
	int di = union_find_bb_find(uf, i, &ri);
	int dj = union_find_bb_find(uf, j, &rj);
	if (ri != rj)
	{
		if (di <= dj) (*uf).e[ri] = (*uf).e[i] = rj;
		else (*uf).e[rj] = (*uf).e[j] = ri;
	}
}

int union_find_bb_find(array_int *uf, int i, int *ri)
{
	if ((*uf).e[i] != i)
	{
		int depth = union_find_bb_find(uf, (*uf).e[i], ri) + 1;
		(*uf).e[i] = (*ri);
		return depth;
	}
	else 
	{
		(*ri) = i;
		return 0;
	}
}

