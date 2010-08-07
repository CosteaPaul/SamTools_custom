#ifndef BCF_PROB1_H
#define BCF_PROB1_H

#include "bcf.h"

struct __bcf_p1aux_t;
typedef struct __bcf_p1aux_t bcf_p1aux_t;

typedef struct {
	int rank0;
	double f_em, f_exp, f_flat, p_ref, g[3];
} bcf_p1rst_t;

#define MC_PTYPE_FULL  1
#define MC_PTYPE_COND2 2
#define MC_PTYPE_FLAT  3

#ifdef __cplusplus
extern "C" {
#endif

	bcf_p1aux_t *bcf_p1_init(int n);
	void bcf_p1_init_prior(bcf_p1aux_t *ma, int type, double theta);
	void bcf_p1_destroy(bcf_p1aux_t *ma);
	int bcf_p1_cal(bcf1_t *b, bcf_p1aux_t *ma, bcf_p1rst_t *rst);
	int bcf_p1_call_gt(const bcf_p1aux_t *ma, double f0, int k);
	void bcf_p1_dump_afs(bcf_p1aux_t *ma);

#ifdef __cplusplus
}
#endif

#endif
