/*
Sources: 
	* Algorithms 4th Edition (Sedgewick & Wayne): 
		* Data structures: resizable arrays, heaps, queues, stacks, hash functions, etc.
*/

#ifndef BASICS_H
#define BASICS_H
#include "stdio.h"
#include "stdlib.h"
#include "ctype.h"
#include "string.h"
#include "math.h"
#define INIT_A_MEM 8

double ** bb_matrix_alloc(int M, int N);
void bb_matrix_free(double **A, int M);

// A variable length array of void pointers
typedef struct
{
	void **e;
	int mem;
	int len;
} array_voidstar;

typedef struct
{
	void *data;
	int block_size;
	int len;
	int mem;
} array_bit;

static const int array_bit_int_masks[] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304, 8388608, 16777216, 33554432, 67108864, 134217728, 268435456, 536870912, 1073741824, -2147483648};

static const int array_bit_int_cmasks[] = {-2, -3, -5, -9, -17, -33, -65, -129, -257, -513, -1025, -2049, -4097, -8193, -16385, -32769, -65537, -131073, -262145, -524289, -1048577, -2097153, -4194305, -8388609, -16777217, -33554433, -67108865, -134217729, -268435457, -536870913, -1073741825, 2147483647}; 

void array_bit_int_init(array_bit *abit, int length);
void add_mem_array_bit_int(array_bit *abit);
void add2array_bit_int(array_bit *abit, char bit);
void array_bit_int_set(array_bit *abit, int bit_pos, char bit_val);
char array_bit_int_get(array_bit *abit, int bit_pos);
void array_bit_int_ll(array_bit *abit);
void array_bit_int_rr(array_bit *abit);
void array_bit_int_and(array_bit *abit1, array_bit *abit2, array_bit *abit3);
void array_bit_int_or(array_bit *abit1, array_bit *abit2, array_bit *abit3);
void free_array_bit_int(array_bit *abit);

void display_array_bit_int(array_bit *abit);

// A variable length array of integers
typedef struct 
{
	int *e;
	int mem;
	int len;
} array_int;

// A variable length array of variable length arrays
typedef struct
{
	array_int *e;
	int len;
	int mem;
} aarray_int;

// A variable length array of doubles
typedef struct 
{
	double *e;
	int mem;
	int len;
} array_double;

// A variable length array of variable length arrays
typedef struct
{
	array_double *e;
	int len;
	int mem;
} aarray_double;

// A variable length array of chars
typedef struct 
{
	char *e;
	int mem;
	int len;
} array_char;

// A variable length array of variable length arrays
typedef struct
{
	array_char *e;
	int len;
	int mem;
} aarray_char;

// A more general linked list:
typedef struct gll
{
	void *data;
	struct gll *next;
} linked_list;

// A box list with a variable number of dimensions
typedef struct 
{
	// the dimension
	int dim;
	// the number of boxes along each dimension
	int *m;
	// the content of each box
	aarray_int boxc;
	// the address(es) of each element
	aarray_int addr;
} boxlist;

// A box list with three dimensions
typedef struct
{
	int m[3];
	aarray_int boxc;
	aarray_int addr;
} boxlist3D;

// A box list with two dimensions
typedef struct
{
	int m[2];
	aarray_int boxc;
	aarray_int addr;
} boxlist2D;


// A neighbor list or graph
typedef struct
{
	// a list of neighbors for each vertex
	aarray_int v;
	// index of each vertex w.r.t. its neighbors
	aarray_int i_of;
} nbrlist; 

typedef struct generic_iter
{
	int state;
	char (*advance)(struct generic_iter *, void *);
} generic_iterator;

// methods for variable length arrays
void array_voidstar_init(array_voidstar *a, int size_);
void add_mem_array_voidstar(array_voidstar *a);
void add_mem_array_voidstar_until(array_voidstar *a, int i);
void add2array_voidstar(array_voidstar *a, void *i);
void print_array_voidstar(array_voidstar a, aarray_char *format);
void remove_array_voidstar(array_voidstar *a, int n, void (*free_elem)(void *));
void contract_array_voidstar(array_voidstar *a, void (*free_elem)(void *));
void free_array_voidstar(array_voidstar *a, void (*free_elem)(void *));
void reset_array_voidstar(array_voidstar *a);
void transcribe_array_voidstar(array_voidstar *a, array_voidstar *b);
void merge_array_voidstar(array_voidstar *a, void **e, int *c, int imin, int imax, int jmin, int jmax, char (*order)(void *, void *));
void merge_sort_array_voidstar(array_voidstar *a, void **e, int i, int j, char (*order)(void *, void *));
void sort_array_voidstar(array_voidstar *a, char (*order)(void *, void *));

void array_int_init(array_int *a, int size_);
char array_int_contains(array_int *a, int elem);
int array_int_search(array_int *a, int elem);
void add_mem_array_int(array_int *a);
void add_mem_array_int_until(array_int *a, int i);
void add2array_int(array_int *a, int i);
void print_array_int(array_int a);
void remove_array_int(array_int *a, int n);
void contract_array_int(array_int *a);
void free_array_int(array_int *a);
void reset_array_int(array_int *a);
void transcribe_array_int(array_int *a, array_int *b);
void merge_array_int(array_int *a, int *e, int *c, int imin, int imax, int jmin, int jmax);
void merge_sort_array_int(array_int *a, int *e, int i, int j);
void sort_array_int(array_int *a);

void array_double_init(array_double *a, int size_);
void add_mem_array_double(array_double *a);
void add_mem_array_double_until(array_double *a, int i);
void add2array_double(array_double *a, double i);
void print_array_double(array_double a);
void remove_array_double(array_double *a, int n);
void contract_array_double(array_double *a);
void free_array_double(array_double *a);
void reset_array_double(array_double *a);
void transcribe_array_double(array_double *a, array_double *b);
void merge_array_double(array_double *a, double *e, int *c, int imin, int imax, int jmin, int jmax);
void merge_sort_array_double(array_double *a, double *e, int i, int j);
void sort_array_double(array_double *a);

void array_char_init(array_char *a, int size_);
char array_char_contains(array_char *a, char elem);
int array_char_search(array_char *a, char elem);
void add_mem_array_char(array_char *a);
void add_mem_array_char_until(array_char *a, int i);
void add2array_char(array_char *a, char i);
void print_array_char(array_char a);
void remove_array_char(array_char *a, int n);
void contract_array_char(array_char *a);
void free_array_char(array_char *a);
void reset_array_char(array_char *a);
void transcribe_array_char(array_char *a, array_char *b);
void merge_array_char(array_char *a, char *e, int *c, int imin, int imax, int jmin, int jmax);
void merge_sort_array_char(array_char *a, char *e, int i, int j);
void sort_array_char(array_char *a);

// methods for variable length arrays of variable length arrays
void aarray_int_init(aarray_int *aa, int mem);
void aarray_int_init_precise(aarray_int *aa, int mem1, int mem2);
void add_mem_aarray_int(aarray_int *aa);
void add_mem_aarray_int_until(aarray_int *aa, int i);
void extend_aarray_int(aarray_int *aa);
void add2aarray_int(aarray_int *aa, array_int a);
void add2aarray_int_elem(aarray_int *aa, int ai, int i);
void reset_aarray_int_elem(int i, aarray_int *aa);
void remove_aarray_int(aarray_int *aa, int i);
void print_aarray_int(aarray_int aa);
void free_aarray_int(aarray_int *aa);
void print_aarray_int(aarray_int aa);
void merge_aarray_int(aarray_int *a, array_int *e, int *c, int imin, int imax, int jmin, int jmax, char (*order)(array_int *, array_int *));
void merge_sort_aarray_int(aarray_int *a, array_int *e, int i, int j, char (*order)(array_int *, array_int *));
void sort_aarray_int(aarray_int *a, char (*order)(array_int *, array_int *));

void aarray_double_init(aarray_double *aa, int mem);
void aarray_double_init_precise(aarray_double *aa, int mem1, int mem2);
void add_mem_aarray_double(aarray_double *aa);
void add_mem_aarray_double_until(aarray_double *aa, int i);
void extend_aarray_double(aarray_double *aa);
void add2aarray_double(aarray_double *aa, array_double a);
void add2aarray_double_elem(aarray_double *aa, int ai, double i);
void reset_aarray_double_elem(int i, aarray_double *aa);
void remove_aarray_double(aarray_double *aa, int i);
void print_aarray_double(aarray_double aa);
void free_aarray_double(aarray_double *aa);
void print_aarray_double(aarray_double aa);
void merge_aarray_double(aarray_double *a, array_double *e, int *c, int imin, int imax, int jmin, int jmax, char (*order)(array_double *, array_double *));
void merge_sort_aarray_double(aarray_double *a, array_double *e, int i, int j, char (*order)(array_double *, array_double *));
void sort_aarray_double(aarray_double *a, char (*order)(array_double *, array_double *));

void aarray_char_init(aarray_char *aa, int mem);
void aarray_char_init_precise(aarray_char *aa, int mem1, int mem2);
void add_mem_aarray_char(aarray_char *aa);
void add_mem_aarray_char_until(aarray_char *aa, int i);
void extend_aarray_char(aarray_char *aa);
void add2aarray_char(aarray_char *aa, array_char a);
void add2aarray_char_elem(aarray_char *aa, int ai, char i);
void reset_aarray_char_elem(int i, aarray_char *aa);
void remove_aarray_char(aarray_char *aa, int i);
void print_aarray_char(aarray_char aa);
void free_aarray_char(aarray_char *aa);
void print_aarray_char(aarray_char aa);
void merge_aarray_char(aarray_char *a, array_char *e, int *c, int imin, int imax, int jmin, int jmax, char (*order)(array_char *, array_char *));
void merge_sort_aarray_char(aarray_char *a, array_char *e, int i, int j, char (*order)(array_char *, array_char *));
void sort_aarray_char(aarray_char *a, char (*order)(array_char *, array_char *));

// Methods for linked lists
linked_list linked_list_init();
void add2linked_list(linked_list *ll, void *elem);
void* pop_linked_list(linked_list **ll);
void free_linked_list(linked_list **ll);

// methods for variable dimension box lists
int flatten_vdim(int dim, int *m, int *box_index); 
void unflatten_box_index(int dim, int *m, int *flatbi, int *unflatbi);

void add2box(boxlist *bl, int elem, int *box_index);
void add2box_preflattened(boxlist *bl, int elem, int fbox_index);
void remove_boxlist_elem(boxlist *bl, int *box_index, int content_index);
void free_boxlist(boxlist *bl);
boxlist boxlist_init(int dim, int *m); // RESUME: consider replacing this with a void version
void move_boxlist_elem(boxlist *bl, int elem, int box_index);
void reset_boxlist(boxlist *bl);
void print_boxlist_counts(boxlist *bl);
void print_boxlist(boxlist *bl);

// methods for 3D box lists
void add2box3D(boxlist3D *bl, int elem, int i, int j, int k);
void remove_boxlist3D_elem(boxlist3D *bl, int fi, int content_index);
void remove_boxlist3D_elem_unflat(boxlist3D *bl, int i, int j, int k, int content_index);
void free_boxlist3D(boxlist3D *bl);
boxlist3D boxlist3D_init(int m, int n, int o);
void move_boxlist3D_elem(boxlist3D *bl, int elem, int ijk); 
void reset_boxlist3D(boxlist3D *bl);

// methods for 2D box lists
void add2box2D(boxlist2D *bl, int elem, int i, int j);
void remove_boxlist2D_elem(boxlist2D *bl, int fi, int content_index);
void remove_boxlist2D_elem_unflat(boxlist2D *bl, int i, int j, int content_index);
void free_boxlist2D(boxlist2D *bl);
boxlist2D boxlist2D_init(int m, int n);
void move_boxlist2D_elem(boxlist2D *bl, int elem, int ij); 
void reset_boxlist2D(boxlist2D *bl);


// methods for neighbor lists
// RESUME: convert versions of initializers with non-void return values to void routines acting 
//		on pointers.
void nbrlist_init_precise(nbrlist *nbl, int mem);
void nbrlist_init(nbrlist *nbl);
void prep_nbrlist(nbrlist *nbl);
void check_nbrlist(nbrlist *nbl);
void extend_nbrlist(nbrlist *nbl);
void extend_nbrlist_n(nbrlist *nbl, int N);
void ensure_nbrlist_size_n(nbrlist *nbl, int N);
void set_len_nbrlist(nbrlist *nbl, int len);
void add_edge_nbrlist(nbrlist *nbl, int vertex1, int vertex2);
void add_edge_nbrlist_safe(nbrlist *nbl, int vertex1, int vertex2);
void remove_edge_nbrlist_1way(nbrlist *nbl, int vertex, int local_nbr_index);
void remove_edge_nbrlist(nbrlist *nbl, int vertex, int local_nbr_index);
void remove_all_edges_nbrlist(nbrlist *nbl);
void remove_edges_vertex_nbrlist(nbrlist *nbl, int vertex);
void remove_vertex_nbrlist(nbrlist *nbl, int vertex);
void write_nbrlist(nbrlist *nbl, char *ofname);
void free_nbrlist(nbrlist *nbl);
int N_vertices_nbrlist(nbrlist *nbl);
int N_nbors_nbrlist(nbrlist *nbl, int v_index);

// misc (to be sorted)
void unflatten_box_index3D(int *m, int *flatbi, int *unflatbi);

char one_bit(char *a, int *b);
void remove_121_int(array_int *inj, array_int *left_inverse, int i);
void add_121_int(array_int *inj, array_int *left_inv, int i);
void toupper_string(char *a);
void tolower_string(char *a);
int str_len(char *a);

// methods to load arrays
void load_array_template(FILE *infile, void *aptr, char *prompt, char *mode);
void load_array_double(FILE *infile, array_double *a, char *prompt);
void load_array_int(FILE *infile, array_int *a, char *prompt);
void load_array_char(FILE *infile, array_char *a, char *prompt);

// operations between arrays
void elementwise_product_array_double(array_double *a1, array_double *a2, array_double *a3);
void elementwise_product_array_int(array_int *a1, array_int *a2, array_int *a3);

// operations on character arrays
void array_char_split(array_char *a, char *div, aarray_char *diva);
void array_char_init_from_string(array_char *a, char *str);
void concatenate_array_char(array_char *a, array_char *b);

// binary searches on ordered arrays
char search_ordered_array_voidstar(array_voidstar *a, void *e, char (*order)(void *, void *));

// transcriptions of arrays to void-star arrays
void transcribe_array_int2voidstar(array_int *a, array_voidstar *b);
void transcribe_array_double2voidstar(array_double *, array_voidstar *b);
void transcribe_array_char2voidstar(array_char *, array_voidstar *b);
void transcribe_aarray_int2voidstar(aarray_int *a, array_voidstar *b);
void transcribe_aarray_char2voidstar(aarray_char *a, array_voidstar *b);
void transcribe_aarray_double2voidstar(aarray_double *a, array_voidstar *b);

// orderings
char order_lexical_array_int(array_int *a, array_int *b);
char order_lexical_array_double(array_double *a, array_double *b);
char order_lexical_array_char(array_char *a, array_char *b);
char order_induced_voidstar_lexical_array_double(void *a, void *b);
char order_induced_voidstar_lexical_array_int(void *a, void *b);
char order_induced_voidstar_lexical_array_char(void *a, void *b);

void array_double_diff(double *e1, double *e2, double *e3, int len);
double array_double_dot(double *e1, double *e2, int len);
double array_double_norm(double *ad, int len);

char aarray_int_contains_set(aarray_int *aa, int *set, int len);
int aarray_int_sorted_place(aarray_int *aa, int *seq, int len, char (*order_gt)(void *, void *));
char aarray_int_sorted_contains(aarray_int *aa, int *seq, int len, char (*order_gt)(void *, void *));
char aarray_int_contains(aarray_int *aa, int *seq, int len);

int array_int_min(int *a, int len);
int array_int_max(int *a, int len);

int *parse_int(void *a);
double *parse_double(void *a);
char *parse_char(void *a);

array_int *parse_array_int(void *a);
array_char *parse_array_char(void *a);
array_double *parse_array_double(void *a);
array_voidstar *parse_array_voidstar(void *a);
aarray_int *parse_aarray_int(void *a);
aarray_char *parse_aarray_char(void *a);
aarray_double *parse_aarray_double(void *a);
nbrlist *parse_nbrlist(void *a);

double array_double_min(double *a, int len);
double array_double_max(double *a, int len);
char array_char_min(char *a, int len);
char array_char_max(char *a, int len);



#endif
