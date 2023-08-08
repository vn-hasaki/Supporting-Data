#include "gmxpre.h"
#include "config.h"
#include <assert.h>

#include "gromacs/legacyheaders/force.h"
#include "gromacs/legacyheaders/gmx_omp_nthreads.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_consts.h"
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_common.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/utility/smalloc.h"
#include <slave.h>

#define UNROLLI    NBNXN_CPU_CLUSTER_I_SIZE
#define UNROLLJ    NBNXN_CPU_CLUSTER_I_SIZE

/* We could use nbat->xstride and nbat->fstride, but macros might be faster */
#define X_STRIDE   3
#define F_STRIDE   3
/* Local i-atom buffer strides */
#define XI_STRIDE  3
#define FI_STRIDE  3

extern nbnxn_pairlist_t     *nbl_s;
extern nbnxn_atomdata_t     *nbat_s;
extern nbnxn_atomdata_output_t *out;
extern real                    *fshift_p;
extern interaction_const_t  *ic_s;
extern rvec                 *shift_vec_s;

__thread_local int inx, section, res, st, ed, num;

void
nbnxn_kernel_ElecQSTab_VdwLJ_VgrpF_ref_slave()
{
    real *f = out->f;
    real *Vvdw = out->Vvdw;
    real *Vc = out->Vc;
    const nbnxn_ci_t   *nbln;
    const nbnxn_cj_t   *l_cj;
    const int          *type;
    const real         *q;
    const real         *shiftvec;
    const real         *x;
    const real         *nbfp;
    real                rcut2;
#ifdef VDW_CUTOFF_CHECK
    real                rvdw2;
#endif
    int                 ntype2;
    real                facel;
    real               *nbfp_i;
    int                 n, ci, ci_sh;
    int                 ish, ishf;
    gmx_bool            do_LJ, half_LJ, do_coul, do_self;
    int                 cjind0, cjind1, cjind;
    int                 ip, jp;

    real                xi[UNROLLI*XI_STRIDE];
    real                fi[UNROLLI*FI_STRIDE];
    real                qi[UNROLLI];

#ifdef CALC_ENERGIES
#ifndef ENERGY_GROUPS

    real       Vvdw_ci, Vc_ci;
#else
    int        egp_mask;
    int        egp_sh_i[UNROLLI];
#endif
#endif
#ifdef LJ_POT_SWITCH
    real       swV3, swV4, swV5;
    real       swF2, swF3, swF4;
#endif
#ifdef LJ_EWALD
    real        lje_coeff2, lje_coeff6_6, lje_vc;
    const real *ljc;
#endif

#ifdef CALC_COUL_RF
    real       k_rf2;
#ifdef CALC_ENERGIES
    real       k_rf, c_rf;
#endif
#endif
#ifdef CALC_COUL_TAB
    real       tabscale;
#ifdef CALC_ENERGIES
    real       halfsp;
#endif
#ifndef GMX_DOUBLE
    const real *tab_coul_FDV0;
#else
    const real *tab_coul_F;
    const real *tab_coul_V;
#endif
#endif

    int ninner;

#ifdef COUNT_PAIRS
    int npair = 0;
#endif

#ifdef LJ_POT_SWITCH
    swV3 = ic_s->vdw_switch.c3;
    swV4 = ic_s->vdw_switch.c4;
    swV5 = ic_s->vdw_switch.c5;
    swF2 = 3*ic_s->vdw_switch.c3;
    swF3 = 4*ic_s->vdw_switch.c4;
    swF4 = 5*ic_s->vdw_switch.c5;
#endif

#ifdef LJ_EWALD
    lje_coeff2   = ic_s->ewaldcoeff_lj*ic_s->ewaldcoeff_lj;
    lje_coeff6_6 = lje_coeff2*lje_coeff2*lje_coeff2/6.0;
    lje_vc       = ic_s->sh_lj_ewald;

    ljc          = nbat_s->nbfp_comb;
#endif

#ifdef CALC_COUL_RF
    k_rf2 = 2*ic_s->k_rf;
#ifdef CALC_ENERGIES
    k_rf = ic_s->k_rf;
    c_rf = ic_s->c_rf;
#endif
#endif
#ifdef CALC_COUL_TAB
    tabscale = ic_s->tabq_scale;
#ifdef CALC_ENERGIES
    halfsp = 0.5/ic_s->tabq_scale;
#endif

#ifndef GMX_DOUBLE
    tab_coul_FDV0 = ic_s->tabq_coul_FDV0;
#else
    tab_coul_F    = ic_s->tabq_coul_F;
    tab_coul_V    = ic_s->tabq_coul_V;
#endif
#endif

#ifdef ENERGY_GROUPS
    egp_mask = (1<<nbat_s->neg_2log) - 1;
#endif


    rcut2               = ic_s->rcoulomb*ic_s->rcoulomb;
#ifdef VDW_CUTOFF_CHECK
    rvdw2               = ic_s->rvdw*ic_s->rvdw;
#endif

    ntype2              = nbat_s->ntype*2;
    nbfp                = nbat_s->nbfp;
    q                   = nbat_s->q;
    type                = nbat_s->type;
    facel               = ic_s->epsfac;
    shiftvec            = shift_vec_s[0];
    x                   = nbat_s->x;

    l_cj = nbl_s->cj;

    ninner = 0;

    inx = _PEN;
//#ifndef CALC_ENERGIES
    //printf("---------------------------nbnxn_kernel_ElecQSTab_VdwLJ_VgrpF_ref_slave slave id is %d---------------------------------\n", inx);
//#endif
    section = nbl_s->nci / 64;
	res     = nbl_s->nci % 64;
	st = (inx < res) ? (inx * (section + 1))
								  : (inx *  section + res);
	ed = (inx < res) ? (st       + (section + 1))
								  : (st       +  section);
    if (section == 0)
    {
        st = inx;
        ed = inx + 1;
        if (inx >= res)
        {
            return;
        }
        
    }
    //printf("---------------------------st is %d, ed is %d---------------------------------\n", st, ed);
    for (num = st; num < ed; num++)
    {
        int i, d;

        nbln = &nbl_s->ci[num];

        ish              = (nbln->shift & NBNXN_CI_SHIFT);
        /* x, f and fshift are assumed to be stored with stride 3 */
        ishf             = ish*DIM;
        cjind0           = nbln->cj_ind_start;
        cjind1           = nbln->cj_ind_end;
        /* Currently only works super-cells equal to sub-cells */
        ci               = nbln->ci;
        ci_sh            = (ish == CENTRAL ? ci : -1);
        if (inx == 64)
            printf("---------------------------block 1---------------------------------\n");
        /* We have 5 LJ/C combinations, but use only three inner loops,
         * as the other combinations are unlikely and/or not much faster:
         * inner half-LJ + C for half-LJ + C / no-LJ + C
         * inner LJ + C      for full-LJ + C
         * inner LJ          for full-LJ + no-C / half-LJ + no-C
         */
        do_LJ   = (nbln->shift & NBNXN_CI_DO_LJ(0));
        do_coul = (nbln->shift & NBNXN_CI_DO_COUL(0));
        half_LJ = ((nbln->shift & NBNXN_CI_HALF_LJ(0)) || !do_LJ) && do_coul;
#ifdef LJ_EWALD
        do_self = TRUE;
#else
        do_self = do_coul;
#endif
        if (inx == 64)
            printf("---------------------------block 2---------------------------------\n");
#ifdef CALC_ENERGIES
#ifndef ENERGY_GROUPS
        Vvdw_ci = 0;
        Vc_ci   = 0;
#else
        for (i = 0; i < UNROLLI; i++)
        {
            egp_sh_i[i] = ((nbat_s->energrp[ci]>>(i*nbat_s->neg_2log)) & egp_mask)*nbat_s->nenergrp;
        }
#endif
#endif
        if (inx == 64)
            printf("---------------------------block 3---------------------------------\n");
        for (i = 0; i < UNROLLI; i++)
        {
            for (d = 0; d < DIM; d++)
            {
                xi[i*XI_STRIDE+d] = x[(ci*UNROLLI+i)*X_STRIDE+d] + shiftvec[ishf+d];
                fi[i*FI_STRIDE+d] = 0;
            }

            qi[i] = facel*q[ci*UNROLLI+i];
        }
        if (inx == 64)
            printf("---------------------------block 4---------------------------------\n");
#ifdef CALC_ENERGIES
        if (do_self)
        {
            real Vc_sub_self;

#ifdef CALC_COUL_RF
            Vc_sub_self = 0.5*c_rf;
#endif
#ifdef CALC_COUL_TAB
#ifdef GMX_DOUBLE
            Vc_sub_self = 0.5*tab_coul_V[0];
#else
            Vc_sub_self = 0.5*tab_coul_FDV0[2];
#endif
#endif

            if (l_cj[nbln->cj_ind_start].cj == ci_sh)
            {
                for (i = 0; i < UNROLLI; i++)
                {
                    int egp_ind;
#ifdef ENERGY_GROUPS
                    egp_ind = egp_sh_i[i] + ((nbat_s->energrp[ci]>>(i*nbat_s->neg_2log)) & egp_mask);
#else
                    egp_ind = 0;
#endif
                    /* Coulomb self interaction */
                    Vc[egp_ind]   -= qi[i]*q[ci*UNROLLI+i]*Vc_sub_self;

#ifdef LJ_EWALD
                    /* LJ Ewald self interaction */
                    Vvdw[egp_ind] += 0.5*nbat_s->nbfp[nbat_s->type[ci*UNROLLI+i]*(nbat_s->ntype + 1)*2]/6*lje_coeff6_6;
#endif
                }
            }
        }
#endif  /* CALC_ENERGIES */
        if (inx == 64)
            printf("---------------------------block 5---------------------------------\n");
        cjind = cjind0;
        while (cjind < cjind1 && nbl_s->cj[cjind].excl != 0xffff)
        {
#define CHECK_EXCLS
            if (half_LJ)
            {
#define CALC_COULOMB
#define HALF_LJ
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_inner.h"
#undef HALF_LJ
#undef CALC_COULOMB
            }
            else if (do_coul)
            {
#define CALC_COULOMB
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_inner.h"
#undef CALC_COULOMB
            }
            else
            {
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_inner.h"
            }
#undef CHECK_EXCLS
            cjind++;
        }

        for (; (cjind < cjind1); cjind++)
        {
            if (half_LJ)
            {
#define CALC_COULOMB
#define HALF_LJ
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_inner.h"
#undef HALF_LJ
#undef CALC_COULOMB
            }
            else if (do_coul)
            {
#define CALC_COULOMB
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_inner.h"
#undef CALC_COULOMB
            }
            else
            {
#include "gromacs/mdlib/nbnxn_kernels/nbnxn_kernel_ref_inner.h"
            }
        }
        ninner += cjind1 - cjind0;

        /* Add accumulated i-forces to the force array */
        for (i = 0; i < UNROLLI; i++)
        {
            for (d = 0; d < DIM; d++)
            {
                f[(ci*UNROLLI+i)*F_STRIDE+d] += fi[i*FI_STRIDE+d];
            }
        }
#ifdef CALC_SHIFTFORCES
        if (fshift_p != NULL)
        {
            /* Add i forces to shifted force list */
            for (i = 0; i < UNROLLI; i++)
            {
                for (d = 0; d < DIM; d++)
                {
                    fshift_p[ishf+d] += fi[i*FI_STRIDE+d];
                }
            }
        }
#endif

#ifdef CALC_ENERGIES
#ifndef ENERGY_GROUPS
        *Vvdw += Vvdw_ci;
        *Vc   += Vc_ci;
#endif
#endif
    }

#ifdef COUNT_PAIRS
    printf("atom pairs %d\n", npair);
#endif
}
