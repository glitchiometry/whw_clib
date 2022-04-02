#include "vector.h"

#define TRANSCRIBE_AB for (int aux_i = 0; aux_i < len; aux_i++) a[aux_i] = b[aux_i];


vector_real vector_real_init(int dim)
{
	vector_real vr;
	vr.dim = dim;
	vr.e = (double *) calloc(dim, sizeof(double));
	return vr;
}

void vectors_real_init(vectors_real *vrs, int mem, int dim)
{
	(*vrs).mem = mem;
	(*vrs).dim = dim;
	(*vrs).e = (double ***) calloc(mem, sizeof(double **));
	for (int i = 0; i < mem; i++)
	{
		(*vrs).e[i] = (double **) calloc(1, sizeof(double *));
		(*(*vrs).e[i]) = (double *) calloc(dim, sizeof(double));
	}
	(*vrs).len = 0;
}

double vectors_real_get(vectors_real *vrs, int elem, int coord)
{
	return (*(*vrs).e[elem])[coord];
}

void vectors_real_set(vectors_real *vrs, int elem, int coord, double val)
{
	(*(*vrs).e[elem])[coord] = val;
}

void free_vector_real(vector_real *vr)
{
	if ((*vr).e != NULL) free((*vr).e);
}

void free_vectors_real(vectors_real *vrs)
{
	for (int i = 0; i < (*vrs).mem; i++)
	{
		free(*(*vrs).e[i]);
		free((*vrs).e[i]);
	}
	free((*vrs).e);
}

void reset_mem_vectors_real(vectors_real *vrs, int new_mem)
{
	//printf("Attempting to reset memory\n");
	if (new_mem > (*vrs).len) {}
	else
	{
		printf("Warning: attempting to store a collection of vectors with insufficient memory\n");
	}
	double *** ne = (double ***) calloc(new_mem, sizeof(double **));
	for (int i = 0; i < new_mem; i++)
	{
		ne[i] = (double **) calloc(1, sizeof(double *));
	}
	for (int i = (*vrs).mem; i < new_mem; i++)
	{
		ne[i][0] = (double *) calloc((*vrs).dim, sizeof(double));
	}
	//printf("Attempting to transcribe vectors\n");
	for (int i = 0; i < (*vrs).mem; i++)
	{
		ne[i][0] = (*vrs).e[i][0];
		/*for (int ii = 0; ii < (*vrs).dim; ii++)
		{
			ne[i][0][ii] = (*vrs).e[i][0][ii];
		}*/
	}
	//printf("Freeing original memory\n");
	for (int i = 0; i < (*vrs).mem; i++)
	{
		free((*vrs).e[i]);
	}
	free((*vrs).e);
	//printf("Replacing old reference with new\n");
	(*vrs).e = ne;
	(*vrs).mem = new_mem;
}

void prep_vectors_real(vectors_real *vrs)
{
	if ((*vrs).len < (*vrs).mem) {}
	else
	{
		add_mem_vectors_real(vrs);
	}
}

void add_mem_vectors_real(vectors_real *vrs)
{
	reset_mem_vectors_real(vrs, (*vrs).mem << 1);
}

void add_mem_vectors_real_until(vectors_real *vrs, int lim)
{
	int new_mem = (*vrs).mem;
	while (new_mem < lim) new_mem <<= 1;
	reset_mem_vectors_real(vrs, new_mem);
}

void add2vectors_real(vectors_real *vrs, double *pt)
{
	prep_vectors_real(vrs);
	for (int i = 0; i < (*vrs).dim; i++)
	{
		(*(*vrs).e[(*vrs).len])[i] = pt[i];
	}
	(*vrs).len += 1;
}

void add2vectors_real_zeros(vectors_real *vrs)
{
	prep_vectors_real(vrs);
	for (int i = 0; i < (*vrs).dim; i++) (*vrs).e[(*vrs).len][0][i] = 0;
	(*vrs).len += 1;
}

void remove_vectors_real(vectors_real *vrs, int vertex)
{
	(*vrs).len -= 1;
	double *aux_ptr = *(*vrs).e[vertex];
	*(*vrs).e[vertex] = *(*vrs).e[(*vrs).len];
	*(*vrs).e[(*vrs).len] = aux_ptr;
}

void array_vector_real_init(array_vector_real *avr, int mem)
{
	(*avr).mem = mem;
	(*avr).e = (vector_real *) calloc(mem, sizeof(vector_real));
	for (int i = 0; i < (*avr).mem; i++)
	{
		(*avr).e[i].e = NULL;
		(*avr).e[i].dim = 0;
	}
	(*avr).len = 0;
}

void free_array_vector_real(array_vector_real *avr)
{
	for (int i = 0; i < (*avr).len; i++)
	{
		free_vector_real(&(*avr).e[i]);
	}
	if ((*avr).e != NULL) free((*avr).e);
}

void reset_mem_array_vector_real(array_vector_real *avr, int new_mem)
{
	vector_real *ne = (vector_real *) calloc(new_mem, sizeof(vector_real));
	for (int i = 0; i < (*avr).len; i++)
	{
		ne[i] = (*avr).e[i];
		/*ne[i] = vector_real_init((*avr).e[i].dim);
		for (int ii = 0; ii < (*avr).e[i].dim; ii++)
		{
			ne[i].e[ii] = (*avr).e[i].e[ii];
		}
		free_vector_real(&((*avr).e[i]));
		*/
	}
	for (int i = (*avr).len; i < new_mem; i++)
	{
		ne[i].e = NULL;
		ne[i].dim = 0;
	}
	free((*avr).e);
	(*avr).e = ne;
	(*avr).mem = new_mem;
}

void add_mem_array_vector_real(array_vector_real *avr)
{
	(*avr).mem <<= 1;
	reset_mem_array_vector_real(avr, (*avr).mem);
}

void add_mem_array_vector_real_until(array_vector_real *avr, int lim)
{
	while ((*avr).mem < lim) (*avr).mem <<= 1;
	reset_mem_array_vector_real(avr, (*avr).mem);
}

void extend_array_vector_real(array_vector_real *avr, int dim)
{
	if ((*avr).mem == (*avr).len)
	{
		add_mem_array_vector_real(avr);
	}
	(*avr).e[(*avr).len] = vector_real_init(dim);
	(*avr).len += 1;
}

void add2array_vector_real(array_vector_real *avr, vector_real elem)
{
	if ((*avr).mem > (*avr).len) {}
	else
	{
		add_mem_array_vector_real(avr);	
	}
	(*avr).e[(*avr).len] = elem;
	(*avr).len += 1;
}

void aarray_vector_real_init(aarray_vector_real *aavr, int mem)
{
	(*aavr).mem = mem;
	(*aavr).e = (array_vector_real *) calloc(mem, sizeof(array_vector_real));
	for (int i = 0; i < mem; i++)
	{
		(*aavr).e[i].e = NULL;
		(*aavr).e[i].mem = 0;
		(*aavr).e[i].len = 0;
	}
	(*aavr).len = 0;
}

void free_aarray_vector_real(aarray_vector_real *aavr)
{
	for (int i = 0; i < (*aavr).mem; i++)
	{
		free_array_vector_real(&(*aavr).e[i]);
	}
}

void add2aarray_vector_real(aarray_vector_real *aavr, array_vector_real elem)
{
	prep_aarray_vector_real(aavr);
	(*aavr).e[(*aavr).len] = elem;
	(*aavr).len += 1;
}

void reset_mem_aarray_vector_real(aarray_vector_real *aavr)
{
	if ((*aavr).mem >= (*aavr).len) {}
	else
	{
		printf("Error (reset_mem_aarray_vector_real): attempting to store aarray_vector_real with insufficient memory\n");
		exit(EXIT_FAILURE);
	}
	array_vector_real *ne = (array_vector_real *) calloc((*aavr).mem, sizeof(array_vector_real));
	for (int i = 0; i < (*aavr).len; i++)
	{
		ne[i] = (*aavr).e[i];
	}
	for (int i = (*aavr).len; i < (*aavr).mem; i++)
	{
		ne[i].e = NULL;
		ne[i].mem = 0;
		ne[i].len = 0;
	}
	if ((*aavr).e != NULL) free((*aavr).e);
	(*aavr).e = ne;
}

void add_mem_aarray_vector_real(aarray_vector_real *aavr)
{
	if ((*aavr).mem > 0) {}
	else (*aavr).mem = 1;
	(*aavr).mem <<= 1;
	reset_mem_aarray_vector_real(aavr);
}

void add_mem_aarray_vector_real_until(aarray_vector_real *aavr, int len)
{
	if ((*aavr).mem > 0) {}
	else (*aavr).mem = 1;
	while ((*aavr).mem < len) 
	{
		(*aavr).mem <<= 1;
	}
	reset_mem_aarray_vector_real(aavr);
}

void prep_aarray_vector_real(aarray_vector_real *aavr)
{
	if ((*aavr).mem >= (*aavr).len) {}
	else add_mem_aarray_vector_real(aavr);
}

void extend_aarray_vector_real(aarray_vector_real *aavr)
{
	prep_aarray_vector_real(aavr);
	array_vector_real_init(&((*aavr).e[(*aavr).len]), 1);
	(*aavr).len += 1;
}

void xform_p_vector_real(vector_real *v1, vector_real *v3)
{
	if ((*v1).dim == (*v3).dim)
	{
		for (int i = 0; i < (*v1).dim; i++)
		{
			(*v3).e[i] = (*v1).e[i] + (*v3).e[i];
		}
	}
	else
	{
		printf("Error: attempting to add vectors with incompatible dimensions\n");
	}
}

vector_real plus_vector_real(vector_real *v1, vector_real *v2)
{
	if ((*v1).dim == (*v2).dim)
	{
		vector_real v3 = vector_real_init((*v1).dim);
		for (int i = 0; i < (*v1).dim; i++)
		{
			v3.e[i] = (*v1).e[i] + (*v2).e[i];
		}
		return v3;
	}
	else
	{
		printf("Error: attempting to add vectors with unequal lengths\n");
		exit(EXIT_FAILURE);
	}
}

vector_real minus_vector_real(vector_real *v1, vector_real *v2)
{
	if ((*v1).dim == (*v2).dim)
	{
		vector_real v3 = vector_real_init((*v1).dim);
		for (int i = 0; i < (*v1).dim; i++)
		{
			v3.e[i] = (*v1).e[i] - (*v2).e[i];
		}
		return v3;
	}
	else
	{
		printf("Error: attempting to add vectors with unequal lengths\n");
		exit(EXIT_FAILURE);
	}
}

vector_real multiply_vector_real(vector_real *v, double s)
{
	vector_real vp = vector_real_init((*v).dim);
	for (int i = 0; i < (*v).dim; i++)
	{
		vp.e[i] = (*v).e[i] * s;
	}
	return vp;
}

double vector_real_dot(vector_real *v1, vector_real *v2)
{
	if ((*v1).dim == (*v2).dim)
	{
		double dp = 0;
		for (int i = 0; i < (*v1).dim; i++)
		{
			dp += (*v1).e[i] * (*v2).e[i];
		}
		return dp;
	}
	else
	{
		printf("Error: attempting to compute dot product between vectors with different lengths\n");
		exit(EXIT_FAILURE);
	}
}

double vectors_real_dot(vectors_real *vs, vector_real *v2, int i)
{
	if ((*vs).dim == (*v2).dim)
	{
		double dp = 0;
		for (int ii = 0; ii < (*vs).dim; ii++)
		{
			dp += (*(*vs).e[i])[ii] * (*v2).e[ii];
		}
		return dp;
	}
	else
	{
		printf("Error: attempting to compute dot product between vectors with different lengths\n");
		exit(EXIT_FAILURE);
	}
}

double vectors_reals_dot(vectors_real *v1s, vectors_real *v2s, int i, int j)
{
	if ((*v1s).dim == (*v2s).dim)
	{
		double dp = 0;
		for (int ii = 0; ii < (*v1s).dim; ii++)
		{
			dp += (*(*v1s).e[i])[ii] * (*(*v2s).e[j])[ii];
		}
		return dp;
	}
	else
	{
		printf("Error: attempting to compute dot product between vectors with different lengths\n");
		exit(EXIT_FAILURE);
	}
}

vector_int vector_int_init(int dim)
{
	vector_int vr;
	vr.dim = dim;
	vr.e = (int *) calloc(dim, sizeof(int));
	return vr;
}

void vectors_int_init(vectors_int *vrs, int mem, int dim)
{
	(*vrs).mem = mem;
	(*vrs).dim = dim;
	(*vrs).e = (int ***) calloc(mem, sizeof(int **));
	for (int i = 0; i < mem; i++)
	{
		(*vrs).e[i] = (int **) calloc(1, sizeof(int *));
		*(*vrs).e[i] = (int *) calloc(dim, sizeof(int));
	}
	(*vrs).len = 0;
}

int vectors_int_get(vectors_int *vrs, int elem, int coord)
{
	return (*(*vrs).e[elem])[coord];
}

void vectors_int_set(vectors_int *vrs, int elem, int coord, int val)
{
	(*(*vrs).e[elem])[coord] = val;
}

void array_vector_int_init(array_vector_int *avr, int mem)
{
	(*avr).mem = mem;
	(*avr).e = (vector_int *) calloc(mem, sizeof(vector_int));
	for (int i = 0; i < mem; i++)
	{
		(*avr).e[i].e = NULL;
		(*avr).e[i].dim = 0;
	}
	(*avr).len = 0;
}

void free_array_vector_int(array_vector_int *avr)
{
	for (int i = 0; i < (*avr).len; i++)
	{
		if ((*avr).e[i].dim > 0)
		{
			free_vector_int(&(*avr).e[i]);
		}
	}
	if ((*avr).e != NULL) free((*avr).e);
}

void free_vector_int(vector_int *vr)
{
	if ((*vr).e != NULL) free((*vr).e);
}

void free_vectors_int(vectors_int *vrs)
{
	for (int i = 0; i < (*vrs).mem; i++)
	{
		free(*(*vrs).e[i]);
		free((*vrs).e[i]);
	}
	free((*vrs).e);
}

void reset_mem_vectors_int(vectors_int *vrs)
{
	if ((*vrs).mem > (*vrs).len) {}
	else
	{
		printf("Warning: attempting to store a collection of vectors with insufficient memory\n");
	}
	int ***ne = (int ***) calloc((*vrs).mem, sizeof(int **));
	for (int i = 0; i < (*vrs).mem; i++)
	{
		ne[i] = (int **) calloc(1, sizeof(int *));
		*ne[i] = (int *) calloc((*vrs).dim, sizeof(int));
	}
	for (int i = 0; i < (*vrs).len; i++)
	{
		for (int ii = 0; ii < (*vrs).dim; i++)
		{
			(*ne[i])[ii] = (*(*vrs).e[i])[ii];
		}
	}
	for (int i = 0; i < (*vrs).mem; i++)
	{
		free(*(*vrs).e[i]);
		free((*vrs).e[i]);
	}
	free((*vrs).e);
	(*vrs).e = ne;
}

void prep_vectors_int(vectors_int *vrs)
{
	if ((*vrs).len < (*vrs).mem) {}
	else
	{
		add_mem_vectors_int(vrs);
	}
}

void add_mem_vectors_int(vectors_int *vrs)
{
	(*vrs).mem <<= 1;
	reset_mem_vectors_int(vrs);
}

void add_mem_vectors_int_until(vectors_int *vrs, int lim)
{
	while ((*vrs).mem < lim) (*vrs).mem <<= 1;
	reset_mem_vectors_int(vrs);
}

void add2vectors_int(vectors_int *vrs, int *pt)
{
	prep_vectors_int(vrs);
	for (int i = 0; i < (*vrs).dim; i++)
	{
		(*(*vrs).e[(*vrs).len])[i] = pt[i];
	}
	(*vrs).len += 1;
}

void add2vectors_int_zeros(vectors_int *vrs)
{
	prep_vectors_int(vrs);
	for (int i = 0; i < (*vrs).dim; i++) (*(*vrs).e[(*vrs).len])[i] = 0;
	(*vrs).len += 1;
}

void remove_vectors_int(vectors_int *vrs, int vertex)
{
	(*vrs).len -= 1;
	int *aux_ptr = *(*vrs).e[vertex];
	*(*vrs).e[vertex] = *(*vrs).e[(*vrs).len];
	*(*vrs).e[(*vrs).len] = aux_ptr;
	/*for (int i = 0; i < (*vrs).dim; i++)
	{
		(*vrs).e[vertex][i] = (*vrs).e[(*vrs).len][i];
	}*/
}

void reset_mem_array_vector_int(array_vector_int *avr)
{
	vector_int *ne = (vector_int *) calloc((*avr).mem, sizeof(vector_int));
	for (int i = 0; i < (*avr).len; i++)
	{
		// ne[i] = vector_int_init((*avr).e[i].dim);
		ne[i] = (*avr).e[i];
		/*for (int ii = 0; ii < (*avr).e[i].dim; ii++)
		{
			ne[i].e[ii] = (*avr).e[i].e[ii];
		}
		free_vector_int(&((*avr).e[i]));*/
	}
	for (int i = (*avr).len; i < (*avr).mem; i++)
	{
		ne[i].e = NULL;
		ne[i].dim = 0;
	}
	free((*avr).e);
	(*avr).e = ne;
}

void add_mem_array_vector_int(array_vector_int *avr)
{
	(*avr).mem <<= 1;
	reset_mem_array_vector_int(avr);
}

void add_mem_array_vector_int_until(array_vector_int *avr, int lim)
{
	while ((*avr).mem < lim) (*avr).mem <<= 1;
	reset_mem_array_vector_int(avr);
}

void extend_array_vector_int(array_vector_int *avr, int dim)
{
	if ((*avr).mem == (*avr).len)
	{
		add_mem_array_vector_int(avr);
	}
	(*avr).e[(*avr).len] = vector_int_init(dim);
	(*avr).len += 1;
}

void add2array_vector_int(array_vector_int *avr, vector_int elem)
{
	if ((*avr).mem > (*avr).len) {}
	else
	{
		add_mem_array_vector_int(avr);	
	}
	(*avr).e[(*avr).len] = elem;
	(*avr).len += 1;
}

void aarray_vector_int_init(aarray_vector_int *aavr, int mem)
{
	(*aavr).mem = mem;
	(*aavr).e = (array_vector_int *) calloc(mem, sizeof(array_vector_int));
	(*aavr).len = 0;
}

void free_aarray_vector_int(aarray_vector_int *aavr)
{
	for (int i = 0; i < (*aavr).len; i++)
	{
		free_array_vector_int(&(*aavr).e[i]);
	}
}

void add2aarray_vector_int(aarray_vector_int *aavr, array_vector_int elem)
{
	prep_aarray_vector_int(aavr);
	(*aavr).e[(*aavr).len] = elem;
	(*aavr).len += 1;
}

void reset_mem_aarray_vector_int(aarray_vector_int *aavr)
{
	if ((*aavr).mem >= (*aavr).len) {}
	else
	{
		printf("Error (reset_mem_aarray_vector_int): attempting to store aarray_vector_int with insufficient memory\n");
		exit(EXIT_FAILURE);
	}
	array_vector_int *ne = (array_vector_int *) calloc((*aavr).mem, sizeof(array_vector_int));
	for (int i = 0; i < (*aavr).len; i++)
	{
		ne[i] = (*aavr).e[i];
	}
	for (int i = (*aavr).len; i < (*aavr).mem; i++)
	{
		ne[i].e = NULL;
		ne[i].mem = 0;
		ne[i].len = 0;
	}
	if ((*aavr).e != NULL) free((*aavr).e);
	(*aavr).e = ne;
}

void add_mem_aarray_vector_int(aarray_vector_int *aavr)
{
	if ((*aavr).mem > 0) {}
	else (*aavr).mem = 1;
	(*aavr).mem <<= 1;
	reset_mem_aarray_vector_int(aavr);
}

void add_mem_aarray_vector_int_until(aarray_vector_int *aavr, int len)
{
	if ((*aavr).mem > 0) {}
	else (*aavr).mem = 1;
	while ((*aavr).mem < len) 
	{
		(*aavr).mem <<= 1;
	}
	reset_mem_aarray_vector_int(aavr);
}

void prep_aarray_vector_int(aarray_vector_int *aavr)
{
	if ((*aavr).mem >= (*aavr).len) {}
	else add_mem_aarray_vector_int(aavr);
}

void extend_aarray_vector_int(aarray_vector_int *aavr)
{
	prep_aarray_vector_int(aavr);
	array_vector_int_init(&((*aavr).e[(*aavr).len]), 1);
	(*aavr).len += 1;
}

array_vectors_real array_vectors_real_init(int mem)
{
	array_vectors_real avr;
	avr.mem = mem;
	avr.e = (vectors_real *) calloc(mem, sizeof(vectors_real));
	avr.len = 0;
	return avr;
}

void reset_mem_array_vectors_real(array_vectors_real *avr)
{
	vectors_real *ne = (vectors_real *) calloc((*avr).mem, sizeof(vectors_real));
	for (int i = 0; i < (*avr).len; i++)
	{
		ne[i] = (*avr).e[i];
		/*vectors_real_init(&ne[i], (*avr).e[i].mem, (*avr).e[i].dim);
		for (int ii = 0; ii < (*avr).e[i].len; ii++)
		{
			for (int iii = 0; iii < (*avr).e[i].dim; iii++)
			{
				ne[i].e[ii][iii] = (*avr).e[i].e[ii][iii];
			}
		}
		free_vectors_real(&(*avr).e[i]);*/
	}
	free((*avr).e);
	(*avr).e = ne;
}

void add_mem_array_vectors_real(array_vectors_real *avr)
{
	(*avr).mem <<= 1;
	reset_mem_array_vectors_real(avr);
}

void add_mem_array_vectors_real_until(array_vectors_real *avr, int lim)
{
	while ((*avr).mem < lim) (*avr).mem <<= 1;
	reset_mem_array_vectors_real(avr);
}

void extend_array_vectors_real(array_vectors_real *avr, int mem, int dim)
{
	if ((*avr).len < (*avr).mem) {}
	else
	{
		add_mem_array_vectors_real(avr);
	}
	vectors_real_init(&(*avr).e[(*avr).len], mem, dim);
	(*avr).len += 1;
}

void free_array_vectors_real(array_vectors_real *avr)
{
	for (int i = 0; i < (*avr).len; i++)
	{
		free_vectors_real(&(*avr).e[i]);
	}
	free((*avr).e);
}

void xform_p_vector_int(vector_int *v1, vector_int *v3)
{
	if ((*v1).dim == (*v3).dim)
	{
		for (int i = 0; i < (*v1).dim; i++)
		{
			(*v3).e[i] = (*v1).e[i] + (*v3).e[i];
		}
	}
	else
	{
		printf("Error: attempting to add vectors with incompatible dimensions\n");
	}
}

vector_int plus_vector_int(vector_int *v1, vector_int *v2)
{
	if ((*v1).dim == (*v2).dim)
	{
		vector_int v3 = vector_int_init((*v1).dim);
		for (int i = 0; i < (*v1).dim; i++)
		{
			v3.e[i] = (*v1).e[i] + (*v2).e[i];
		}
		return v3;
	}
	else
	{
		printf("Error: attempting to add vectors with unequal lengths\n");
		exit(EXIT_FAILURE);
	}
}

vector_int minus_vector_int(vector_int *v1, vector_int *v2)
{
	if ((*v1).dim == (*v2).dim)
	{
		vector_int v3 = vector_int_init((*v1).dim);
		for (int i = 0; i < (*v1).dim; i++)
		{
			v3.e[i] = (*v1).e[i] - (*v2).e[i];
		}
		return v3;
	}
	else
	{
		printf("Error: attempting to add vectors with unequal lengths\n");
		exit(EXIT_FAILURE);
	}
}

vector_int multiply_vector_int(vector_int *v, int s)
{
	vector_int vp = vector_int_init((*v).dim);
	for (int i = 0; i < (*v).dim; i++)
	{
		vp.e[i] = (*v).e[i] * s;
	}
	return vp;
}

int vector_int_dot(vector_int *v1, vector_int *v2)
{
	if ((*v1).dim == (*v2).dim)
	{
		int dp = 0;
		for (int i = 0; i < (*v1).dim; i++)
		{
			dp += (*v1).e[i] * (*v2).e[i];
		}
		return dp;
	}
	else
	{
		printf("Error: attempting to compute dot product between vectors with different lengths\n");
		exit(EXIT_FAILURE);
	}
}

int vectors_int_dot(vectors_int *vs, vector_int *v2, int i)
{
	if ((*vs).dim == (*v2).dim)
	{
		int dp = 0;
		for (int ii = 0; ii < (*vs).dim; ii++)
		{
			dp += (*(*vs).e[i])[ii] * (*v2).e[ii];
		}
		return dp;
	}
	else
	{
		printf("Error: attempting to compute dot product between vectors with different lengths\n");
		exit(EXIT_FAILURE);
	}
}

int vectors_ints_dot(vectors_int *v1s, vectors_int *v2s, int i, int j)
{
	if ((*v1s).dim == (*v2s).dim)
	{
		int dp = 0;
		for (int ii = 0; ii < (*v1s).dim; ii++)
		{
			dp += (*(*v1s).e[i])[ii] * (*(*v2s).e[j])[ii];
		}
		return dp;
	}
	else
	{
		printf("Error: attempting to compute dot product between vectors with different lengths\n");
		exit(EXIT_FAILURE);
	}
}

array_vectors_int array_vectors_int_init(int mem)
{
	array_vectors_int avr;
	avr.mem = mem;
	avr.e = (vectors_int *) calloc(mem, sizeof(vectors_int));
	avr.len = 0;
	return avr;
}

void reset_mem_array_vectors_int(array_vectors_int *avr)
{
	vectors_int *ne = (vectors_int *) calloc((*avr).mem, sizeof(vectors_int));
	for (int i = 0; i < (*avr).len; i++)
	{
		ne[i] = (*avr).e[i];
		/*vectors_int_init(&ne[i], (*avr).e[i].mem, (*avr).e[i].dim);
		for (int ii = 0; ii < (*avr).e[i].len; ii++)
		{
			for (int iii = 0; iii < (*avr).e[i].dim; iii++)
			{
				ne[i].e[ii][iii] = (*avr).e[i].e[ii][iii];
			}
		}
		free_vectors_int(&(*avr).e[i]); */
	}
	free((*avr).e);
	(*avr).e = ne;
}

void add_mem_array_vectors_int(array_vectors_int *avr)
{
	(*avr).mem <<= 1;
	reset_mem_array_vectors_int(avr);
}

void add_mem_array_vectors_int_until(array_vectors_int *avr, int lim)
{
	while ((*avr).mem < lim) (*avr).mem <<= 1;
	reset_mem_array_vectors_int(avr);
}

void extend_array_vectors_int(array_vectors_int *avr, int mem, int dim)
{
	if ((*avr).len < (*avr).mem) {}
	else
	{
		add_mem_array_vectors_int(avr);
	}
	vectors_int_init(&(*avr).e[(*avr).len], mem, dim);
	(*avr).len += 1;
}

void free_array_vectors_int(array_vectors_int *avr)
{
	for (int i = 0; i < (*avr).len; i++)
	{
		free_vectors_int(&(*avr).e[i]);
	}
	free((*avr).e);
}

void set_equal_double(double *a, double *b, int len)
{
	TRANSCRIBE_AB
}

void set_equal_int(int *a, int *b, int len)
{
	TRANSCRIBE_AB
}

vector_char vector_char_init(int dim)
{
	vector_char vr;
	vr.dim = dim;
	vr.e = (char *) calloc(dim, sizeof(char));
	return vr;
}

void vectors_char_init(vectors_char *vrs, int mem, int dim)
{
	(*vrs).mem = mem;
	(*vrs).dim = dim;
	(*vrs).e = (char ***) calloc(mem, sizeof(char **));
	for (int i = 0; i < mem; i++)
	{
		(*vrs).e[i] = (char **) calloc(1, sizeof(char *));
		*(*vrs).e[i] = (char *) calloc(dim, sizeof(char));
	}
	(*vrs).len = 0;
}

char vectors_char_get(vectors_char *vrs, int elem, int coord)
{
	return (*(*vrs).e[elem])[coord];
}

void vectors_char_set(vectors_char *vrs, int elem, int coord, char val)
{
	(*(*vrs).e[elem])[coord] = val;
}

void array_vector_char_init(array_vector_char *avr, int mem)
{
	(*avr).mem = mem;
	(*avr).e = (vector_char *) calloc(mem, sizeof(vector_char));
	for (int i = 0; i < mem; i++)
	{
		(*avr).e[i].e = NULL;
		(*avr).e[i].dim = 0;
	}
	(*avr).len = 0;
}

void free_array_vector_char(array_vector_char *avr)
{
	for (int i = 0; i < (*avr).len; i++)
	{
		free_vector_char(&(*avr).e[i]);
	}
	if ((*avr).e != NULL) free((*avr).e);
}


void free_vector_char(vector_char *vr)
{
	if ((*vr).e != NULL) free((*vr).e);
}

void free_vectors_char(vectors_char *vrs)
{
	for (int i = 0; i < (*vrs).mem; i++)
	{
		free(*(*vrs).e[i]);
		free((*vrs).e[i]);
	}
	free((*vrs).e);
}

void reset_mem_vectors_char(vectors_char *vrs)
{
	if ((*vrs).mem > (*vrs).len) {}
	else
	{
		printf("Warning: attempting to store a collection of vectors with insufficient memory\n");
	}
	char ***ne = (char ***) calloc((*vrs).mem, sizeof(char **));
	for (int i = 0; i < (*vrs).mem; i++)
	{
		ne[i] = (char **) calloc(1, sizeof(char *));
		*ne[i] = (char *) calloc((*vrs).dim, sizeof(char));
	}
	for (int i = 0; i < (*vrs).len; i++)
	{
		for (int ii = 0; ii < (*vrs).dim; i++)
		{
			(*ne[i])[ii] = (*(*vrs).e[i])[ii];
		}
	}
	for (int i = 0; i < (*vrs).mem; i++)
	{
		free(*(*vrs).e[i]);
		free((*vrs).e[i]);
	}
	free((*vrs).e);
	(*vrs).e = ne;
}

void prep_vectors_char(vectors_char *vrs)
{
	if ((*vrs).len < (*vrs).mem) {}
	else
	{
		add_mem_vectors_char(vrs);
	}
}

void add_mem_vectors_char(vectors_char *vrs)
{
	(*vrs).mem <<= 1;
	reset_mem_vectors_char(vrs);
}

void add_mem_vectors_char_until(vectors_char *vrs, int lim)
{
	while ((*vrs).mem < lim) (*vrs).mem <<= 1;
	reset_mem_vectors_char(vrs);
}

void add2vectors_char(vectors_char *vrs, char *pt)
{
	prep_vectors_char(vrs);
	for (int i = 0; i < (*vrs).dim; i++)
	{
		(*(*vrs).e[(*vrs).len])[i] = pt[i];
	}
	(*vrs).len += 1;
}


void add2vectors_char_zeros(vectors_char *vrs)
{
	prep_vectors_char(vrs);
	for (int i = 0; i < (*vrs).dim; i++) (*vrs).e[(*vrs).len][i] = 0;
	(*vrs).len += 1;
}

void remove_vectors_char(vectors_char *vrs, int vertex)
{
	(*vrs).len -= 1;
	char *aux_ptr = *(*vrs).e[vertex];
	*(*vrs).e[vertex] = *(*vrs).e[(*vrs).len];
	*(*vrs).e[(*vrs).len] = aux_ptr;
	/*for (int i = 0; i < (*vrs).dim; i++)
	{
		(*vrs).e[vertex][i] = (*vrs).e[(*vrs).len][i];
	}*/
}

void reset_mem_array_vector_char(array_vector_char *avr)
{
	vector_char *ne = (vector_char *) calloc((*avr).mem, sizeof(vector_char));
	for (int i = 0; i < (*avr).len; i++)
	{
		ne[i] = (*avr).e[i];
		/*ne[i] = vector_char_init((*avr).e[i].dim);
		for (int ii = 0; ii < (*avr).e[i].dim; ii++)
		{
			ne[i].e[ii] = (*avr).e[i].e[ii];
		}
		free_vector_char(&((*avr).e[i]));*/
	}
	for (int i = (*avr).len; i < (*avr).mem; i++)
	{
		ne[i].e = NULL;
		ne[i].dim = 0;
	}
	free((*avr).e);
	(*avr).e = ne;
}

void add_mem_array_vector_char(array_vector_char *avr)
{
	(*avr).mem <<= 1;
	reset_mem_array_vector_char(avr);
}

void add_mem_array_vector_char_until(array_vector_char *avr, int lim)
{
	while ((*avr).mem < lim) (*avr).mem <<= 1;
	reset_mem_array_vector_char(avr);
}

void extend_array_vector_char(array_vector_char *avr, int dim)
{
	if ((*avr).mem == (*avr).len)
	{
		add_mem_array_vector_char(avr);
	}
	(*avr).e[(*avr).len] = vector_char_init(dim);
	(*avr).len += 1;
}

void add2array_vector_char(array_vector_char *avr, vector_char elem)
{
	if ((*avr).mem > (*avr).len) {}
	else
	{
		add_mem_array_vector_char(avr);	
	}
	(*avr).e[(*avr).len] = elem;
	(*avr).len += 1;
}

void aarray_vector_char_init(aarray_vector_char *aavr, int mem)
{
	(*aavr).mem = mem;
	(*aavr).e = (array_vector_char *) calloc(mem, sizeof(array_vector_char));
	(*aavr).len = 0;
}

void free_aarray_vector_char(aarray_vector_char *aavr)
{
	for (int i = 0; i < (*aavr).len; i++)
	{
		free_array_vector_char(&(*aavr).e[i]);
	}
}

void add2aarray_vector_char(aarray_vector_char *aavr, array_vector_char elem)
{
	prep_aarray_vector_char(aavr);
	(*aavr).e[(*aavr).len] = elem;
	(*aavr).len += 1;
}

void reset_mem_aarray_vector_char(aarray_vector_char *aavr)
{
	if ((*aavr).mem >= (*aavr).len) {}
	else
	{
		printf("Error (reset_mem_aarray_vector_char): attempting to store aarray_vector_char with insufficient memory\n");
		exit(EXIT_FAILURE);
	}
	array_vector_char *ne = (array_vector_char *) calloc((*aavr).mem, sizeof(array_vector_char));
	for (int i = 0; i < (*aavr).len; i++)
	{
		ne[i] = (*aavr).e[i];
	}
	for (int i = (*aavr).len; i < (*aavr).mem; i++)
	{
		ne[i].e = NULL;
		ne[i].mem = 0;
		ne[i].len = 0;
	}
	if ((*aavr).e != NULL) free((*aavr).e);
	(*aavr).e = ne;
}

void add_mem_aarray_vector_char(aarray_vector_char *aavr)
{
	if ((*aavr).mem > 0) {}
	else (*aavr).mem = 1;
	(*aavr).mem <<= 1;
	reset_mem_aarray_vector_char(aavr);
}

void add_mem_aarray_vector_char_until(aarray_vector_char *aavr, int len)
{
	if ((*aavr).mem > 0) {}
	else (*aavr).mem = 1;
	while ((*aavr).mem < len) 
	{
		(*aavr).mem <<= 1;
	}
	reset_mem_aarray_vector_char(aavr);
}

void prep_aarray_vector_char(aarray_vector_char *aavr)
{
	if ((*aavr).mem >= (*aavr).len) {}
	else add_mem_aarray_vector_char(aavr);
}

void extend_aarray_vector_char(aarray_vector_char *aavr)
{
	prep_aarray_vector_char(aavr);
	array_vector_char_init(&((*aavr).e[(*aavr).len]), 1);
	(*aavr).len += 1;
}

void xform_p_vector_char(vector_char *v1, vector_char *v3)
{
	if ((*v1).dim == (*v3).dim)
	{
		for (int i = 0; i < (*v1).dim; i++)
		{
			(*v3).e[i] = (*v1).e[i] + (*v3).e[i];
		}
	}
	else
	{
		printf("Error: attempting to add vectors with incompatible dimensions\n");
	}
}

vector_char plus_vector_char(vector_char *v1, vector_char *v2)
{
	if ((*v1).dim == (*v2).dim)
	{
		vector_char v3 = vector_char_init((*v1).dim);
		for (int i = 0; i < (*v1).dim; i++)
		{
			v3.e[i] = (*v1).e[i] + (*v2).e[i];
		}
		return v3;
	}
	else
	{
		printf("Error: attempting to add vectors with unequal lengths\n");
		exit(EXIT_FAILURE);
	}
}

vector_char minus_vector_char(vector_char *v1, vector_char *v2)
{
	if ((*v1).dim == (*v2).dim)
	{
		vector_char v3 = vector_char_init((*v1).dim);
		for (int i = 0; i < (*v1).dim; i++)
		{
			v3.e[i] = (*v1).e[i] - (*v2).e[i];
		}
		return v3;
	}
	else
	{
		printf("Error: attempting to add vectors with unequal lengths\n");
		exit(EXIT_FAILURE);
	}
}

vector_char multiply_vector_char(vector_char *v, char s)
{
	vector_char vp = vector_char_init((*v).dim);
	for (int i = 0; i < (*v).dim; i++)
	{
		vp.e[i] = (*v).e[i] * s;
	}
	return vp;
}

char vector_char_dot(vector_char *v1, vector_char *v2)
{
	if ((*v1).dim == (*v2).dim)
	{
		char dp = 0;
		for (int i = 0; i < (*v1).dim; i++)
		{
			dp += (*v1).e[i] * (*v2).e[i];
		}
		return dp;
	}
	else
	{
		printf("Error: attempting to compute dot product between vectors with different lengths\n");
		exit(EXIT_FAILURE);
	}
}

char vectors_char_dot(vectors_char *vs, vector_char *v2, int i)
{
	if ((*vs).dim == (*v2).dim)
	{
		char dp = 0;
		for (int ii = 0; ii < (*vs).dim; ii++)
		{
			dp += (*(*vs).e[i])[ii] * (*v2).e[ii];
		}
		return dp;
	}
	else
	{
		printf("Error: attempting to compute dot product between vectors with different lengths\n");
		exit(EXIT_FAILURE);
	}
}

char vectors_chars_dot(vectors_char *v1s, vectors_char *v2s, int i, int j)
{
	if ((*v1s).dim == (*v2s).dim)
	{
		char dp = 0;
		for (int ii = 0; ii < (*v1s).dim; ii++)
		{
			dp += (*(*v1s).e[i])[ii] * (*(*v2s).e[j])[ii];
		}
		return dp;
	}
	else
	{
		printf("Error: attempting to compute dot product between vectors with different lengths\n");
		exit(EXIT_FAILURE);
	}
}

array_vectors_char array_vectors_char_init(int mem)
{
	array_vectors_char avr;
	avr.mem = mem;
	avr.e = (vectors_char *) calloc(mem, sizeof(vectors_char));
	avr.len = 0;
	return avr;
}

void reset_mem_array_vectors_char(array_vectors_char *avr)
{
	vectors_char *ne = (vectors_char *) calloc((*avr).mem, sizeof(vectors_char));
	for (int i = 0; i < (*avr).len; i++)
	{
		ne[i] = (*avr).e[i];
		/*vectors_char_init(&ne[i], (*avr).e[i].mem, (*avr).e[i].dim);
		for (int ii = 0; ii < (*avr).e[i].len; ii++)
		{
			for (int iii = 0; iii < (*avr).e[i].dim; iii++)
			{
				ne[i].e[ii][iii] = (*avr).e[i].e[ii][iii];
			}
		}
		free_vectors_char(&(*avr).e[i]); */
	}
	free((*avr).e);
	(*avr).e = ne;
}

void add_mem_array_vectors_char(array_vectors_char *avr)
{
	(*avr).mem <<= 1;
	reset_mem_array_vectors_char(avr);
}

void add_mem_array_vectors_char_until(array_vectors_char *avr, int lim)
{
	while ((*avr).mem < lim) (*avr).mem <<= 1;
	reset_mem_array_vectors_char(avr);
}

void extend_array_vectors_char(array_vectors_char *avr, int mem, int dim)
{
	if ((*avr).len < (*avr).mem) {}
	else
	{
		add_mem_array_vectors_char(avr);
	}
	vectors_char_init(&(*avr).e[(*avr).len], mem, dim);
	(*avr).len += 1;
}

void free_array_vectors_char(array_vectors_char *avr)
{
	for (int i = 0; i < (*avr).len; i++)
	{
		free_vectors_char(&(*avr).e[i]);
	}
	free((*avr).e);
}

void set_equal_char(char *a, char *b, int len)
{
	TRANSCRIBE_AB
}

double vector_real_norm(vector_real *vr)
{
	double len = 0;
	for (int i = 0; i < (*vr).dim; i++)
	{
		len += (*vr).e[i] * (*vr).e[i];
	}
	return sqrt(len);
}
double vector_int_norm(vector_int *vi)
{
	double len = 0;
	for (int i = 0; i < (*vi).dim; i++)
	{
		len += (double) (*vi).e[i] * (*vi).e[i];
	}
	return sqrt(len);
}
