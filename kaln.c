/* The MIT License

   Copyright (c) 2003-2006, 2008, 2009, by Heng Li <lh3lh3@gmail.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include "kaln.h"

#define FROM_M 0
#define FROM_I 1
#define FROM_D 2

typedef struct {
	int i, j;
	unsigned char ctype;
} path_t;

int aln_sm_blosum62[] = {
/*	 A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  *  X */
	 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-4, 0,
	-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-4,-1,
	-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,-4,-1,
	-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,-4,-1,
	 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-4,-2,
	-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,-4,-1,
	-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,-4,-1,
	 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-4,-1,
	-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,-4,-1,
	-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-4,-1,
	-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-1,
	-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,-4,-1,
	-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-4,-1,
	-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-4,-1,
	-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-4,-2,
	 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2,-4, 0,
	 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-4, 0,
	-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-2,
	-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-4,-1,
	 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-4,-1,
	-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1,-4,
	 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-4,-1
};

int aln_sm_blast[] = {
	1, -3, -3, -3, -2,
	-3, 1, -3, -3, -2,
	-3, -3, 1, -3, -2,
	-3, -3, -3, 1, -2,
	-2, -2, -2, -2, -2
};

ka_param_t ka_param_blast = {  5,  2,   5, 2, aln_sm_blast, 5, 50 };
ka_param_t ka_param_aa2aa = { 10,  2,  10, 2, aln_sm_blosum62, 22, 50 };

static uint32_t *ka_path2cigar32(const path_t *path, int path_len, int *n_cigar)
{
	int i, n;
	uint32_t *cigar;
	unsigned char last_type;

	if (path_len == 0 || path == 0) {
		*n_cigar = 0;
		return 0;
	}

	last_type = path->ctype;
	for (i = n = 1; i < path_len; ++i) {
		if (last_type != path[i].ctype) ++n;
		last_type = path[i].ctype;
	}
	*n_cigar = n;
	cigar = (uint32_t*)calloc(*n_cigar, 4);

	cigar[0] = 1u << 4 | path[path_len-1].ctype;
	last_type = path[path_len-1].ctype;
	for (i = path_len - 2, n = 0; i >= 0; --i) {
		if (path[i].ctype == last_type) cigar[n] += 1u << 4;
		else {
			cigar[++n] = 1u << 4 | path[i].ctype;
			last_type = path[i].ctype;
		}
	}

	return cigar;
}

/***************************/
/* START OF common_align.c */
/***************************/

#define SET_INF(s) (s).M = (s).I = (s).D = MINOR_INF;

#define set_M(MM, cur, p, sc)							\
{														\
	if ((p)->M >= (p)->I) {								\
		if ((p)->M >= (p)->D) {							\
			(MM) = (p)->M + (sc); (cur)->Mt = FROM_M;	\
		} else {										\
			(MM) = (p)->D + (sc); (cur)->Mt = FROM_D;	\
		}												\
	} else {											\
		if ((p)->I > (p)->D) {							\
			(MM) = (p)->I + (sc); (cur)->Mt = FROM_I;	\
		} else {										\
			(MM) = (p)->D + (sc); (cur)->Mt = FROM_D;	\
		}												\
	}													\
}
#define set_I(II, cur, p)								\
{														\
	if ((p)->M - gap_open > (p)->I) {					\
		(cur)->It = FROM_M;								\
		(II) = (p)->M - gap_open - gap_ext;				\
	} else {											\
		(cur)->It = FROM_I;								\
		(II) = (p)->I - gap_ext;						\
	}													\
}
#define set_end_I(II, cur, p)							\
{														\
	if (gap_end_ext >= 0) {								\
		if ((p)->M - gap_end_open > (p)->I) {			\
			(cur)->It = FROM_M;							\
			(II) = (p)->M - gap_end_open - gap_end_ext;	\
		} else {										\
			(cur)->It = FROM_I;							\
			(II) = (p)->I - gap_end_ext;				\
		}												\
	} else set_I(II, cur, p);							\
}
#define set_D(DD, cur, p)								\
{														\
	if ((p)->M - gap_open > (p)->D) {					\
		(cur)->Dt = FROM_M;								\
		(DD) = (p)->M - gap_open - gap_ext;				\
	} else {											\
		(cur)->Dt = FROM_D;								\
		(DD) = (p)->D - gap_ext;						\
	}													\
}
#define set_end_D(DD, cur, p)							\
{														\
	if (gap_end_ext >= 0) {								\
		if ((p)->M - gap_end_open > (p)->D) {			\
			(cur)->Dt = FROM_M;							\
			(DD) = (p)->M - gap_end_open - gap_end_ext;	\
		} else {										\
			(cur)->Dt = FROM_D;							\
			(DD) = (p)->D - gap_end_ext;				\
		}												\
	} else set_D(DD, cur, p);							\
}

typedef struct {
	uint8_t Mt:3, It:2, Dt:2;
} dpcell_t;

typedef struct {
	int M, I, D;
} dpscore_t;

/***************************
 * banded global alignment *
 ***************************/
uint32_t *ka_global_core(uint8_t *seq1, int len1, uint8_t *seq2, int len2, const ka_param_t *ap, int *_score, int *n_cigar)
{
	int i, j;
	dpcell_t **dpcell, *q;
	dpscore_t *curr, *last, *s;
	int b1, b2, tmp_end;
	int *mat, end, max = 0;
	uint8_t type, ctype;
	uint32_t *cigar = 0;

	int gap_open, gap_ext, gap_end_open, gap_end_ext, b;
	int *score_matrix, N_MATRIX_ROW;

	/* initialize some align-related parameters. just for compatibility */
	gap_open = ap->gap_open;
	gap_ext = ap->gap_ext;
	gap_end_open = ap->gap_end_open;
	gap_end_ext = ap->gap_end_ext;
	b = ap->band_width;
	score_matrix = ap->matrix;
	N_MATRIX_ROW = ap->row;

	*n_cigar = 0;
	if (len1 == 0 || len2 == 0) return 0;

	/* calculate b1 and b2 */
	if (len1 > len2) {
		b1 = len1 - len2 + b;
		b2 = b;
	} else {
		b1 = b;
		b2 = len2 - len1 + b;
	}
	if (b1 > len1) b1 = len1;
	if (b2 > len2) b2 = len2;
	--seq1; --seq2;

	/* allocate memory */
	end = (b1 + b2 <= len1)? (b1 + b2 + 1) : (len1 + 1);
	dpcell = (dpcell_t**)malloc(sizeof(dpcell_t*) * (len2 + 1));
	for (j = 0; j <= len2; ++j)
		dpcell[j] = (dpcell_t*)malloc(sizeof(dpcell_t) * end);
	for (j = b2 + 1; j <= len2; ++j)
		dpcell[j] -= j - b2;
	curr = (dpscore_t*)malloc(sizeof(dpscore_t) * (len1 + 1));
	last = (dpscore_t*)malloc(sizeof(dpscore_t) * (len1 + 1));
	
	/* set first row */
	SET_INF(*curr); curr->M = 0;
	for (i = 1, s = curr + 1; i < b1; ++i, ++s) {
		SET_INF(*s);
		set_end_D(s->D, dpcell[0] + i, s - 1);
	}
	s = curr; curr = last; last = s;

	/* core dynamic programming, part 1 */
	tmp_end = (b2 < len2)? b2 : len2 - 1;
	for (j = 1; j <= tmp_end; ++j) {
		q = dpcell[j]; s = curr; SET_INF(*s);
		set_end_I(s->I, q, last);
		end = (j + b1 <= len1 + 1)? (j + b1 - 1) : len1;
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		++s; ++q;
		for (i = 1; i != end; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]); /* this will change s->M ! */
			set_I(s->I, q, last + i);
			set_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + i - 1, mat[seq1[i]]);
		set_D(s->D, q, s - 1);
		if (j + b1 - 1 > len1) { /* bug fixed, 040227 */
			set_end_I(s->I, q, last + i);
		} else s->I = MINOR_INF;
		s = curr; curr = last; last = s;
	}
	/* last row for part 1, use set_end_D() instead of set_D() */
	if (j == len2 && b2 != len2 - 1) {
		q = dpcell[j]; s = curr; SET_INF(*s);
		set_end_I(s->I, q, last);
		end = (j + b1 <= len1 + 1)? (j + b1 - 1) : len1;
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		++s; ++q;
		for (i = 1; i != end; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]); /* this will change s->M ! */
			set_I(s->I, q, last + i);
			set_end_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + i - 1, mat[seq1[i]]);
		set_end_D(s->D, q, s - 1);
		if (j + b1 - 1 > len1) { /* bug fixed, 040227 */
			set_end_I(s->I, q, last + i);
		} else s->I = MINOR_INF;
		s = curr; curr = last; last = s;
		++j;
	}

	/* core dynamic programming, part 2 */
	for (; j <= len2 - b2 + 1; ++j) {
		SET_INF(curr[j - b2]);
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		end = j + b1 - 1;
		for (i = j - b2 + 1, q = dpcell[j] + i, s = curr + i; i != end; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]);
			set_I(s->I, q, last + i);
			set_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + i - 1, mat[seq1[i]]);
		set_D(s->D, q, s - 1);
		s->I = MINOR_INF;
		s = curr; curr = last; last = s;
	}

	/* core dynamic programming, part 3 */
	for (; j < len2; ++j) {
		SET_INF(curr[j - b2]);
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		for (i = j - b2 + 1, q = dpcell[j] + i, s = curr + i; i < len1; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]);
			set_I(s->I, q, last + i);
			set_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + len1 - 1, mat[seq1[i]]);
		set_end_I(s->I, q, last + i);
		set_D(s->D, q, s - 1);
		s = curr; curr = last; last = s;
	}
	/* last row */
	if (j == len2) {
		SET_INF(curr[j - b2]);
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		for (i = j - b2 + 1, q = dpcell[j] + i, s = curr + i; i < len1; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]);
			set_I(s->I, q, last + i);
			set_end_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + len1 - 1, mat[seq1[i]]);
		set_end_I(s->I, q, last + i);
		set_end_D(s->D, q, s - 1);
		s = curr; curr = last; last = s;
	}

	*_score = last[len1].M;
	if (n_cigar) { /* backtrace */
		path_t *p, *path = (path_t*)malloc(sizeof(path_t) * (len1 + len2 + 2));
		i = len1; j = len2;
		q = dpcell[j] + i;
		s = last + len1;
		max = s->M; type = q->Mt; ctype = FROM_M;
		if (s->I > max) { max = s->I; type = q->It; ctype = FROM_I; }
		if (s->D > max) { max = s->D; type = q->Dt; ctype = FROM_D; }

		p = path;
		p->ctype = ctype; p->i = i; p->j = j; /* bug fixed 040408 */
		++p;
		do {
			switch (ctype) {
			case FROM_M: --i; --j; break;
			case FROM_I: --j; break;
			case FROM_D: --i; break;
			}
			q = dpcell[j] + i;
			ctype = type;
			switch (type) {
			case FROM_M: type = q->Mt; break;
			case FROM_I: type = q->It; break;
			case FROM_D: type = q->Dt; break;
			}
			p->ctype = ctype; p->i = i; p->j = j;
			++p;
		} while (i || j);
		cigar = ka_path2cigar32(path, p - path - 1, n_cigar);
		free(path);
	}

	/* free memory */
	for (j = b2 + 1; j <= len2; ++j)
		dpcell[j] += j - b2;
	for (j = 0; j <= len2; ++j)
		free(dpcell[j]);
	free(dpcell);
	free(curr); free(last);

	return cigar;
}

#define get_k(b, i, kk) ((i) < (b)? (kk)/3 : (kk)/3+(i)-(b))
#define set_u(u, b, i, k) { int x=(i)-(b); x=x>0?x:0; (u)=((k)-x+1)*3; }

ka_probpar_t ka_probpar_def = { 0.0001, 0.1, 10 };

int ka_prob_extend(uint8_t *_ref, int l_ref, uint8_t *_query, int l_query, float *_qual,
				   const ka_probpar_t *c, int *state, uint8_t *q)
{
	double **f, **b, m[9], f_sink, sI, sM, pf, pb;
	float *qual;
	uint8_t *ref, *query;
	int bw, bw2, i, k;
	/*** initialization ***/
	ref = _ref - 1;
	query = _query - 1;
	qual = _qual - 1;
	bw = c->bw;
	bw2 = c->bw * 2 + 1;
	f = calloc(l_query+1, sizeof(void*));
	b = calloc(l_query+1, sizeof(void*));
	for (i = 0; i <= l_query; ++i) {
		f[i] = calloc(bw2 * 3 + 6, sizeof(double));
		b[i] = calloc(bw2 * 3 + 6, sizeof(double));
	}
	// initialize m
	sM = sI = 1. / (2 * l_query + 1);
	m[0*3+0] = (1 - c->d - c->d) * (1 - sM); m[0*3+1] = m[0*3+2] = c->d * (1 - sM);
	m[1*3+0] = (1 - c->e) * (1 - sI); m[1*3+1] = c->e * (1 - sI); m[1*3+2] = 0.;
	m[2*3+0] = 1 - c->e; m[2*3+1] = 0.; m[2*3+2] = c->e;
	/*** forward ***/
	// f[0]
	set_u(k, bw, 0, 0);
	f[0][k] = 1.;
	// f[1..l_query]; core loop
	for (i = 1; i <= l_query; ++i) {
		double *fi = f[i], *fi1 = f[i-1];
		int beg = 0, end = l_ref, x;
		x = i - bw; beg = beg > x? beg : x;
		x = i + bw; end = end < x? end : x;
		for (k = beg; k <= end; ++k) {
			int u, v11, v01, v10;
			double e;
			e = (ref[k] == 4 || query[i] == 4)? 1. : ref[k] == query[i]? 1. - qual[i] : qual[i]; // FIXME: can be better
			set_u(u, bw, i, k); set_u(v11, bw, i-1, k-1); set_u(v10, bw, i-1, k); set_u(v01, bw, i, k-1);
			fi[u+0] = e * (m[0] * fi1[v11+0] + m[3] * fi1[v11+1] + m[6] * fi1[v11+2]);
			fi[u+1] = .25 * (m[1] * fi1[v10+0] + m[4] * fi1[v10+1]);
			fi[u+2] = m[2] * fi[v01+0] + m[8] * fi[v01+2];
		}
	}
	// sink
	for (k = 0, f_sink = 0.; k <= l_ref; ++k) {
		int u;
		set_u(u, bw, l_query, k);
		if (u < 3 || u >= bw2*3+3) continue;
		f_sink += f[l_query][u+0] * sM + f[l_query][u+1] * sI;
	}
	pf = f_sink;
	/*** backward ***/
	// sink
	for (k = 0; k <= l_ref; ++k) {
		int u;
		double *bi = b[l_query];
		set_u(u, bw, l_query, k);
		if (u < 3 || u >= bw2*3+3) continue;
		bi[u+0] = sM; bi[u+1] = sI;
	}
	// b[l_query-1..1]; core loop
	for (i = l_query - 1; i > 0; --i) {
		int beg = 0, end = l_ref, x;
		double *bi = b[i], *bi1 = b[i+1];
		x = i - bw; beg = beg > x? beg : x;
		x = i + bw; end = end < x? end : x;
		for (k = end; k >= beg; --k) {
			int u, v11, v01, v10;
			double e;
			e = k >= l_ref? 0 : (ref[k+1] == 4 || query[i+1] == 4)? 1. : ref[k+1] == query[i+1]? 1. - qual[i+1] : qual[i+1];
			set_u(u, bw, i, k); set_u(v11, bw, i+1, k+1); set_u(v10, bw, i+1, k); set_u(v01, bw, i, k+1);
			bi[u+0] = e * m[0] * bi1[v11+0] + .25 * m[1] * bi1[v10+1] + m[2] * bi[v01+2];
			bi[u+1] = e * m[3] * bi1[v11+0] + .25 * m[4] * bi1[v10+1];
			bi[u+2] = e * m[6] * bi1[v11+0] + m[8] * bi[v01+2];
		}
	}
	{ // the "zero" coordinate (i==k==0)
		double *bi = b[0], *bi1 = b[1];
		int u, v11, v01, v10;
		double e = (ref[1] == 4 || query[1] == 4)? 1. : ref[1] == query[1]? 1. - qual[1] : qual[1];
		set_u(u, bw, 0, 0); set_u(v11, bw, 1, 1); set_u(v10, bw, 1, 0); set_u(v01, bw, 0, 1);
		bi[u+0] = e * m[0] * bi1[v11+0] + .25 * m[1] * bi1[v10+1] + m[2] * bi[v01+2];
	}
	set_u(k, bw, 0, 0);
	pb = b[0][k];
	{ // pf should equal pb as they are both the likelihood of data but calculated in different ways
		double z = fabs(pb/pf - 1.);
		if (z > 1e-7) fprintf(stderr, "[%s] %lg!=%lg\n", __func__, pf, pb);
	}
	/*** MAP ***/
	return 0;
}

#ifdef _MAIN
int main()
{
	int l_ref = 5, l_query = 4;
	uint8_t *ref = (uint8_t*)"\0\1\3\3\1";
	uint8_t *query = (uint8_t*)"\0\1\3\1";
//	uint8_t *query = "\1\3\3\1";
	static float qual[4] = {.2, .01, .01, .01};
	ka_prob_extend(ref, l_ref, query, l_query, qual, &ka_probpar_def, 0, 0);
	return 0;
}
#endif
