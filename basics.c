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
		(*a).len = 0;
	}
	else 
	{
		(*a).e = NULL;
	}
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

void remove_array_voidstar(array_voidstar *a, int n, void (*free_elem)(void *))
{
	if (n < (*a).len)
	{
		(*a).len -= 1;
		if ((*a).e[n] != NULL) 
		{
			if (free_elem != NULL) free_elem((*a).e[n]);
			free((*a).e[n]); // RESUME: check this!
		}
		(*a).e[n] = (*a).e[(*a).len];
		(*a).e[(*a).len] = NULL;
		if ((*a).len < (*a).mem >> 2) contract_array_voidstar(a, free_elem);
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


void free_array_voidstar(array_voidstar *a, void (*free_elem)(void *))
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
	}
}

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

void remove_array_int(array_int *a, int n)
{
	if (n < (*a).len && n > -1)
	{
		(*a).len -= 1;
		(*a).e[n] = (*a).e[(*a).len];
		if ((*a).len < (*a).mem >> 2) contract_array_int(a);
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


void free_array_int(array_int *a)
{
	if ((*a).e != NULL) free((*a).e);
}

void reset_array_int(array_int *a)
{
	(*a).len = 0;
}

void transcribe_array_int(array_int *a, array_int *b)
{
	add_mem_array_int_until(a, (*b).len);
	(*a).len = (*b).len;
	for (int i = 0; i < (*b).len; i++)
	{
		(*a).e[i] = (*b).e[i];
	}
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
		int mdpt = (i + j) / 2;
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

void sort_array_int(array_int *a)
{
	int *e = (int *) calloc((*a).mem, sizeof(int));
	merge_sort_array_int(a, e, 0, (*a).len);
	free((*a).e);
	(*a).e = e;
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

void free_aarray_int(aarray_int *aa)
{
	int i = 0;
	while (i < (*aa).mem && (*aa).e[i].mem > 0)
	{
		free_array_int(&((*aa).e[i]));
		// free(&((*aa).e[i])); // RESUME: Check this! 000
		i += 1;
	}
	free((*aa).e);
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

void free_aarray_double(aarray_double *aa)
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
	return nll;
}

void add2linked_list(linked_list *ll, void *elem)
{
	if ((*ll).next != NULL) {}
	else 
	{
		(*ll).next = elem;
		return;
	}
	linked_list *new = (linked_list *) calloc(1, sizeof(linked_list));
	(*new).data = (*ll).data;
	(*new).next = (*ll).next;
	(*ll).data = elem;
	(*ll).next = new;
}

void* pop_linked_list(linked_list **ll)
{
	void *popped = (**ll).data;
	linked_list *next_ll = (**ll).next;
	free(*ll);
	*ll = next_ll;
	return popped;
}

void free_linked_list(linked_list **ll)
{
	while ((**ll).next != NULL) pop_linked_list(ll);
}


// 

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
	//	(as each has a known length)
	array_int new_addr;
	array_int_init(&new_addr, 3);
	int fi = flatten_vdim((*bl).dim, (*bl).m, box_index);
	new_addr.len = 3;
	new_addr.e[0] = fi;
	new_addr.e[1] = (*bl).boxc.e[fi].len;
	new_addr.e[2] = elem;
	add2aarray_int_elem(&((*bl).boxc), fi, (*bl).addr.len);
	add2aarray_int(&((*bl).addr), new_addr); // WARNING: POTENTIAL MEMORY PROBLEM (try to test this)
	// RESUME: Consider either revising the basic aarray structure so that additional arrays can be 
	//		incorporated directly without transcribing their contents, or clearing new_addr at the end of this.
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

// this function has an effect on both addr and boxc
//	- remove an element from the box list
void remove_boxlist_elem(boxlist *bl, int *box_index, int content_index)
{
	int fi = flatten_vdim((*bl).dim, (*bl).m, box_index);
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

void free_boxlist(boxlist *bl)
{
	free_aarray_int(&((*bl).boxc));
	free_aarray_int(&((*bl).addr));
	free((*bl).m);
}

boxlist boxlist_init(int dim, int *m)
{
	boxlist bl;
	bl.dim = dim;
	bl.m = (int *) calloc(dim, sizeof(int));
	for (int i = 0; i < dim; i++) bl.m[i] = m[i];
	if (dim > 0)
	{
		int Nboxes = bl.m[0];
		for (int i = 1; i < dim; i++)
		{
			Nboxes *= bl.m[i];
		}
		aarray_int_init_precise(&(bl.addr), 1, 3);
		aarray_int_init_precise(&(bl.boxc), Nboxes, 1);
		bl.boxc.len = Nboxes;
	}
	return bl;
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

void extend_nbrlist(nbrlist *nbl)
{
	prep_nbrlist(nbl);
	(*nbl).v.len += 1;
	(*nbl).i_of.len = (*nbl).v.len;
}

void set_len_nbrlist(nbrlist *nbl, int len)
{
	(*nbl).v.len = len;
	(*nbl).i_of.len = len;
	check_nbrlist(nbl);
}

void extend_nbrlist_n(nbrlist *nbl, int N)
{
	(*nbl).v.len += N;
	(*nbl).i_of.len = (*nbl).v.len;
	check_nbrlist(nbl);
}

void ensure_nbrlist_size_n(nbrlist *nbl, int N)
{
	if ((*nbl).v.len >= N) {}
	else
	{
		extend_nbrlist_n(nbl, N - (*nbl).v.len);
	}
}

void add_edge_nbrlist_safe(nbrlist *nbl, int vertex1, int vertex2)
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
			return;
		}
	}
	add_edge_nbrlist(nbl, vertex1, vertex2);
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
	add2aarray_int_elem(&(*nbl).i_of, vertex1, -1);
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
}

void remove_edge_nbrlist(nbrlist *nbl, int vertex, int local_nbr_index)
{
	int vertex_ = (*nbl).v.e[vertex].e[local_nbr_index];
	int local_nbr_index_ = (*nbl).i_of.e[vertex].e[local_nbr_index];
	if (local_nbr_index_ > -1)
	{
		remove_array_int(&(*nbl).v.e[vertex_], local_nbr_index_);
		remove_array_int(&(*nbl).i_of.e[vertex_], local_nbr_index_);
	}
	remove_array_int(&(*nbl).v.e[vertex], local_nbr_index);
	remove_array_int(&(*nbl).i_of.e[vertex], local_nbr_index);
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

void remove_edges_nbrlist_vertex(nbrlist *nbl, int vertex)
{
	int i = (*nbl).v.e[vertex].len;
	do
	{
		i -= 1;
		remove_edge_nbrlist(nbl, vertex, i);
	} 
	while (i > 0);
}

void remove_all_edges_nbrlist(nbrlist *nbl)
{
	for (int i = 0; i < (*nbl).v.len; i++)
	{
		(*nbl).v.e[i].len = 0;
		(*nbl).i_of.e[i].len = 0;
	}
}

void write_nbrlist(nbrlist *nbl, char *ofname)
{
	FILE *ofile = fopen(ofname, "w");
	if (ofile != NULL)
	{
		for (int i = 0; i < (*nbl).v.len; i++)
		{
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

void remove_array_double(array_double *a, int n)
{
	if (n < (*a).len)
	{
		(*a).len -= 1;
		(*a).e[n] = (*a).e[(*a).len];
		if ((*a).len < (*a).mem >> 2) contract_array_double(a);
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


void free_array_double(array_double *a)
{
	if ((*a).e != NULL) free((*a).e);
}

void reset_array_double(array_double *a)
{
	(*a).len = 0;
}

void transcribe_array_double(array_double *a, array_double *b)
{
	add_mem_array_double_until(a, (*b).len);
	(*a).len = (*b).len;
	for (int i = 0; i < (*b).len; i++)
	{
		(*a).e[i] = (*b).e[i];
	}
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

void remove_array_char(array_char *a, int n)
{
	if (n < (*a).len)
	{
		(*a).len -= 1;
		(*a).e[n] = (*a).e[(*a).len];
		if ((*a).len < (*a).mem >> 2) contract_array_char(a);
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

void free_array_char(array_char *a)
{
	if ((*a).e != NULL) free((*a).e);
}

void reset_array_char(array_char *a)
{
	(*a).len = 0;
}

void transcribe_array_char(array_char *a, array_char *b)
{
	add_mem_array_char_until(a, (*b).len);
	(*a).len = (*b).len;
	for (int i = 0; i < (*b).len; i++)
	{
		(*a).e[i] = (*b).e[i];
	}
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

void free_aarray_char(aarray_char *aa)
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

void merge_array_voidstar(array_voidstar *a, void **e, int *c, int imin, int imax, int jmin, int jmax, char (*order)(void *, void *))
{
	if (order((*a).e[imin], (*a).e[jmin]) == -1)
	{
		int il = imin;
		while (order((*a).e[il], (*a).e[jmin]) == -1 && il < imax)
		{
			e[*c] = (*a).e[il];
			il += 1;
			*c += 1;
		}
		merge_array_voidstar(a, e, c, il, imax, jmin, jmax, order);
	}
	else
	{
		int ir = jmin;
		while (order((*a).e[ir], (*a).e[imin]) == -1 && ir < jmax)
		{
			e[*c] = (*a).e[ir];
			ir += 1;
			*c += 1;
		}
		merge_array_voidstar(a, e, c, imin, imax, ir, jmax, order);
	}
}

// sort a variable length array along the index interval [i, j)
void merge_sort_array_voidstar(array_voidstar *a, void **e, int i, int j, char (*order)(void *, void *))
{
	int jmi = j - i;
	if (jmi > 2)
	{
		int mdpt = (i + j) / 2;
		merge_sort_array_voidstar(a, e, i, mdpt, order);
		merge_sort_array_voidstar(a, e, mdpt, j, order);
		// merge the two ordered lists
		int c = i;
		merge_array_voidstar(a, e, &c, i, mdpt, mdpt, j, order);
	}
	else if (jmi == 2)
	{
		int ip1 = i + 1;
		if (order((*a).e[ip1], (*a).e[i]) == 1) {}
		else
		{
			void *tmp = (*a).e[ip1];
			(*a).e[ip1] = (*a).e[i];
			(*a).e[i] = tmp;
		}
	}
}

void sort_array_voidstar(array_voidstar *a, char (*order)(void *, void *))
{
	void **e = (void **) calloc((*a).mem, sizeof(void *));
	merge_sort_array_voidstar(a, e, 0, (*a).len, order);
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

void array_double_diff(double *e1, double *e2, double *e3, int len)
{
	for (int i = 0; i < len; i++)
	{
		e3[i] = e1[i] - e2[i];
	}
}

double array_double_dot(double *e1, double *e2, int len)
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

// Methods for arrays of bits
void array_bit_int_init(array_bit *abit, int length)
{
	array_int *data = (array_int *) calloc(1, sizeof(array_int));
	int block_size = 8 * sizeof(int); 
	array_int_init(data, length / block_size + 1);
	(*abit).data = (void *) data;
	(*abit).block_size = block_size;
	(*abit).len = 0;
	(*abit).mem = block_size * (*data).mem;
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
	int addr = (*abit).len / (*abit).block_size;
	int bit_addr = (*abit).len % (*abit).block_size;
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
	int addr = bit_pos / (*abit).block_size;
	int bit_addr = bit_pos % (*abit).block_size;
	array_int *data = (array_int *) (*abit).data;
	if (bit_val == 0) (*data).e[addr] &= array_bit_int_cmasks[bit_addr];
	else (*data).e[addr] |= array_bit_int_masks[bit_addr];
}

char array_bit_int_get(array_bit *abit, int bit_pos)
{
	int addr = bit_pos / (*abit).block_size;
	int bit_addr = bit_pos % (*abit).block_size;
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
