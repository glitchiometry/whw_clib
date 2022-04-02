#include "matrix.h"

matrix_real matrix_real_init(int m, int n)
{
	matrix_real vr;
	vr.m = m;
	vr.n = n;
	vr.e = (double **) calloc(m, sizeof(double *));
	for (int i = 0; i < m; i++)
	{
		vr.e[i] = (double *) calloc(n, sizeof(double));
	}
	return vr;
}

matrices_real matrices_real_init(int mem, int m, int n)
{
	matrices_real vrs;
	vrs.mem = mem;
	vrs.m = m;
	vrs.n = n;
	vrs.e = (double ***) calloc(mem, sizeof(double **));
	for (int i = 0; i < mem; i++)
	{
		vrs.e[i] = (double **) calloc(m, sizeof(double *));
		for (int ii = 0; ii < m; ii++)
		{
			vrs.e[i][ii] = (double *) calloc(n, sizeof(double));
		}
	}
	vrs.len = 0;
	return vrs;
}

array_matrix_real array_matrix_real_init(int mem)
{
	array_matrix_real avr;
	avr.mem = mem;
	avr.e = (matrix_real *) calloc(mem, sizeof(matrix_real));
	avr.len = 0;
	return avr;
}

void free_matrix_real(matrix_real *vr)
{
	for (int i = 0; i < (*vr).m; i++) free((*vr).e[i]);
	free((*vr).e);
}

void free_matrices_real(matrices_real *vrs)
{
	for (int i = 0; i < (*vrs).mem; i++)
	{
		for (int ii = 0; ii < (*vrs).m; ii++) free((*vrs).e[i][ii]);
		free((*vrs).e[i]);
	}
	free((*vrs).e);
}

void reset_mem_matrices_real(matrices_real *vrs)
{
	if ((*vrs).mem > (*vrs).len) {}
	else
	{
		printf("Warning: attempting to store a collection of matrices with insufficient memory\n");
	}
	double *** ne = (double ***) calloc((*vrs).mem, sizeof(double **));
	for (int i = 0; i < (*vrs).mem; i++)
	{
		ne[i] = (double **) calloc((*vrs).m, sizeof(double *));
		for (int ii = 0; ii < (*vrs).m; ii++)
		{
			ne[i][ii] = (double *) calloc((*vrs).n, sizeof(double));
		}
	}
	for (int i = 0; i < (*vrs).len; i++)
	{
		for (int ii = 0; ii < (*vrs).m; ii++)
		{
			for (int iii = 0; iii < (*vrs).n; iii++) ne[i][ii][iii] = (*vrs).e[i][ii][iii];
		}
	}
	for (int i = 0; i < (*vrs).mem; i++)
	{
		for (int ii = 0; ii < (*vrs).m; ii++) free((*vrs).e[i][ii]);
		free((*vrs).e[i]);
	}
	free((*vrs).e);
	(*vrs).e = ne;
}

void prep_matrices_real(matrices_real *vrs)
{
	if ((*vrs).len < (*vrs).mem) {}
	else
	{
		add_mem_matrices_real(vrs);
	}
}

void add_mem_matrices_real(matrices_real *vrs)
{
	(*vrs).mem <<= 1;
	reset_mem_matrices_real(vrs);
}

void add_mem_matrices_real_until(matrices_real *vrs, int lim)
{
	while ((*vrs).mem < lim) (*vrs).mem <<= 1;
	reset_mem_matrices_real(vrs);
}

void add2matrices_real(matrices_real *vrs, double **pt)
{
	prep_matrices_real(vrs);
	for (int i = 0; i < (*vrs).m; i++)
	{
		for (int ii = 0; ii < (*vrs).n; ii++) (*vrs).e[(*vrs).len][i][ii] = pt[i][ii];
	}
	(*vrs).len += 1;
}

void add2matrices_real_zeros(matrices_real *vrs)
{
	prep_matrices_real(vrs);
	for (int i = 0; i < (*vrs).m; i++) {
		for (int ii = 0; ii < (*vrs).n; ii++) (*vrs).e[(*vrs).len][i][ii] = 0;
	}
	(*vrs).len += 1;
}

void remove_matrices_real(matrices_real *vrs, int vertex)
{
	(*vrs).len -= 1;
	for (int i = 0; i < (*vrs).m; i++)
	{
		for (int ii = 0; ii < (*vrs).n; ii++) (*vrs).e[vertex][i][ii] = (*vrs).e[(*vrs).len][i][ii];
	}
}

void reset_mem_array_matrix_real(array_matrix_real *avr)
{
	matrix_real *ne = (matrix_real *) calloc((*avr).mem, sizeof(matrix_real));
	for (int i = 0; i < (*avr).len; i++)
	{
		ne[i] = matrix_real_init((*avr).e[i].m, (*avr).e[i].n);
		for (int ii = 0; ii < (*avr).e[i].m; ii++)
		{
			for (int iii = 0; iii < (*avr).e[i].n; iii++) ne[i].e[ii][iii] = (*avr).e[i].e[ii][iii];
		}
		free_matrix_real(&((*avr).e[i]));
	}
	free((*avr).e);
	(*avr).e = ne;
}

void add_mem_array_matrix_real(array_matrix_real *avr)
{
	(*avr).mem <<= 1;
	reset_mem_array_matrix_real(avr);
}

void add_mem_array_matrix_real_until(array_matrix_real *avr, int lim)
{
	while ((*avr).mem < lim) (*avr).mem <<= 1;
	reset_mem_array_matrix_real(avr);
}

void extend_array_matrix_real(array_matrix_real *avr, int m, int n)
{
	if ((*avr).mem == (*avr).len)
	{
		add_mem_array_matrix_real(avr);
	}
	(*avr).e[(*avr).len] = matrix_real_init(m, n);
	(*avr).len += 1;
}



matrix_int matrix_int_init(int m, int n)
{
	matrix_int vr;
	vr.m = m;
	vr.n = n;
	vr.e = (int **) calloc(m, sizeof(int *));
	for (int i = 0; i < m; i++)
	{
		vr.e[i] = (int *) calloc(n, sizeof(int));
	}
	return vr;
}

matrices_int matrices_int_init(int mem, int m, int n)
{
	matrices_int vrs;
	vrs.mem = mem;
	vrs.m = m;
	vrs.n = n;
	vrs.e = (int ***) calloc(mem, sizeof(int **));
	for (int i = 0; i < mem; i++)
	{
		vrs.e[i] = (int **) calloc(m, sizeof(int *));
		for (int ii = 0; ii < m; ii++) vrs.e[i][ii] = (int *) calloc(n, sizeof(int));
	}
	vrs.len = 0;
	return vrs;
}

array_matrix_int array_matrix_int_init(int mem)
{
	array_matrix_int avr;
	avr.mem = mem;
	avr.e = (matrix_int *) calloc(mem, sizeof(matrix_int));
	avr.len = 0;
	return avr;
}

void free_matrix_int(matrix_int *vr)
{
	for (int i = 0; i < (*vr).m; i++) free((*vr).e[i]);
	free((*vr).e);
}

void free_matrices_int(matrices_int *vrs)
{
	for (int i = 0; i < (*vrs).mem; i++)
	{
		for (int ii = 0; ii < (*vrs).m; ii++) free((*vrs).e[i][ii]);
		free((*vrs).e[i]);
	}
	free((*vrs).e);
}

void reset_mem_matrices_int(matrices_int *vrs)
{
	if ((*vrs).mem > (*vrs).len) {}
	else
	{
		printf("Warning: attempting to store a collection of matrices with insufficient memory\n");
	}
	int *** ne = (int ***) calloc((*vrs).mem, sizeof(int **));
	for (int i = 0; i < (*vrs).mem; i++)
	{
		ne[i] = (int **) calloc((*vrs).m, sizeof(int *));
		for (int ii = 0; ii < (*vrs).m; ii++) ne[i][ii] = (int *) calloc((*vrs).n, sizeof(int));
	}
	for (int i = 0; i < (*vrs).len; i++)
	{
		for (int ii = 0; ii < (*vrs).m; ii++)
		{
			for (int iii = 0; iii < (*vrs).n; iii++) ne[i][ii][iii] = (*vrs).e[i][ii][iii];
		}
	}
	for (int i = 0; i < (*vrs).mem; i++)
	{
		for (int ii = 0; ii < (*vrs).m; ii++) free((*vrs).e[i][ii]);
		free((*vrs).e[i]);
	}
	free((*vrs).e);
	(*vrs).e = ne;
}

void prep_matrices_int(matrices_int *vrs)
{
	if ((*vrs).len < (*vrs).mem) {}
	else
	{
		add_mem_matrices_int(vrs);
	}
}

void add_mem_matrices_int(matrices_int *vrs)
{
	(*vrs).mem <<= 1;
	reset_mem_matrices_int(vrs);
}

void add_mem_matrices_int_until(matrices_int *vrs, int lim)
{
	while ((*vrs).mem < lim) (*vrs).mem <<= 1;
	reset_mem_matrices_int(vrs);
}

void add2matrices_int(matrices_int *vrs, int **pt)
{
	prep_matrices_int(vrs);
	for (int i = 0; i < (*vrs).m; i++)
	{
		for (int ii = 0; ii < (*vrs).n; ii++) (*vrs).e[(*vrs).len][i][ii] = pt[i][ii];
	}
	(*vrs).len += 1;
}

void remove_matrices_int(matrices_int *vrs, int vertex)
{
	(*vrs).len -= 1;
	for (int i = 0; i < (*vrs).m; i++)
	{
		for (int ii = 0; ii < (*vrs).n; ii++) (*vrs).e[vertex][i][ii] = (*vrs).e[(*vrs).len][i][ii];
	}
}

void reset_mem_array_matrix_int(array_matrix_int *avr)
{
	matrix_int *ne = (matrix_int *) calloc((*avr).mem, sizeof(matrix_int));
	for (int i = 0; i < (*avr).len; i++)
	{
		ne[i] = matrix_int_init((*avr).e[i].m, (*avr).e[i].n);
		for (int ii = 0; ii < (*avr).e[i].m; ii++)
		{
			for (int iii = 0; iii < (*avr).e[i].n; iii++) ne[i].e[ii][iii] = (*avr).e[i].e[ii][iii];
		}
		free_matrix_int(&((*avr).e[i]));
	}
	free((*avr).e);
	(*avr).e = ne;
}

void add_mem_array_matrix_int(array_matrix_int *avr)
{
	(*avr).mem <<= 1;
	reset_mem_array_matrix_int(avr);
}

void add_mem_array_matrix_int_until(array_matrix_int *avr, int lim)
{
	while ((*avr).mem < lim) (*avr).mem <<= 1;
	reset_mem_array_matrix_int(avr);
}

void extend_array_matrix_int(array_matrix_int *avr, int m, int n)
{
	if ((*avr).mem == (*avr).len)
	{
		add_mem_array_matrix_int(avr);
	}
	(*avr).e[(*avr).len] = matrix_int_init(m, n);
	(*avr).len += 1;
}

void free_array_matrix_real(array_matrix_real *avr)
{
	for (int i = 0; i < (*avr).len; i++)
	{
		free_matrix_real(&(*avr).e[i]);
	}
}

array_matrices_real array_matrices_real_init(int mem)
{
	array_matrices_real avr;
	avr.mem = mem;
	avr.e = (matrices_real *) calloc(mem, sizeof(matrices_real));
	avr.len = 0;
	return avr;
}

void reset_mem_array_matrices_real(array_matrices_real *avr)
{
	matrices_real *ne = (matrices_real *) calloc((*avr).mem, sizeof(matrices_real));
	for (int i = 0; i < (*avr).len; i++)
	{
		ne[i] = matrices_real_init((*avr).e[i].mem, (*avr).e[i].m, (*avr).e[i].n);
		for (int ii = 0; ii < (*avr).e[i].len; ii++)
		{
			for (int iii = 0; iii < (*avr).e[i].m; iii++)
			{
				for (int iv = 0; iv < (*avr).e[i].n; iv++) ne[i].e[ii][iii][iv] = (*avr).e[i].e[ii][iii][iv];
			}
		}
		free_matrices_real(&(*avr).e[i]);
	}
	free(&(*avr).e);
	(*avr).e = ne;
}

void add_mem_array_matrices_real(array_matrices_real *avr)
{
	(*avr).mem <<= 1;
	reset_mem_array_matrices_real(avr);
}

void add_mem_array_matrices_real_until(array_matrices_real *avr, int lim)
{
	while ((*avr).mem < lim) (*avr).mem <<= 1;
	reset_mem_array_matrices_real(avr);
}

void extend_array_matrices_real(array_matrices_real *avr, int mem, int m, int n)
{
	if ((*avr).len < (*avr).mem) {}
	else
	{
		add_mem_array_matrices_real(avr);
	}
	(*avr).e[(*avr).len] = matrices_real_init(mem, m, n);
	(*avr).len += 1;
}

void free_array_matrices_real(array_matrices_real *avr)
{
	for (int i = 0; i < (*avr).len; i++)
	{
		free_matrices_real(&(*avr).e[i]);
	}
	free(&(*avr).e);
}

array_matrices_int array_matrices_int_init(int mem)
{
	array_matrices_int avr;
	avr.mem = mem;
	avr.e = (matrices_int *) calloc(mem, sizeof(matrices_int));
	avr.len = 0;
	return avr;
}

void reset_mem_array_matrices_int(array_matrices_int *avr)
{
	matrices_int *ne = (matrices_int *) calloc((*avr).mem, sizeof(matrices_int));
	for (int i = 0; i < (*avr).len; i++)
	{
		ne[i] = matrices_int_init((*avr).e[i].mem, (*avr).e[i].m, (*avr).e[i].n);
		for (int ii = 0; ii < (*avr).e[i].len; ii++)
		{
			for (int iii = 0; iii < (*avr).e[i].m; iii++)
			{
				for (int iv = 0; iv < (*avr).e[i].n; iv++) ne[i].e[ii][iii][iv] = (*avr).e[i].e[ii][iii][iv];
			}
		}
		free_matrices_int(&(*avr).e[i]);
	}
	free(&(*avr).e);
	(*avr).e = ne;
}

void add_mem_array_matrices_int(array_matrices_int *avr)
{
	(*avr).mem <<= 1;
	reset_mem_array_matrices_int(avr);
}

void add_mem_array_matrices_int_until(array_matrices_int *avr, int lim)
{
	while ((*avr).mem < lim) (*avr).mem <<= 1;
	reset_mem_array_matrices_int(avr);
}

void extend_array_matrices_int(array_matrices_int *avr, int mem, int m, int n)
{
	if ((*avr).len < (*avr).mem) {}
	else
	{
		add_mem_array_matrices_int(avr);
	}
	(*avr).e[(*avr).len] = matrices_int_init(mem, m, n);
	(*avr).len += 1;
}

void free_array_matrices_int(array_matrices_int *avr)
{
	for (int i = 0; i < (*avr).len; i++)
	{
		free_matrices_int(&(*avr).e[i]);
	}
	free(&(*avr).e);
}
