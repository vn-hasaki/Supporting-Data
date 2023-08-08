#include <slave.h>
#include <simd.h>
#define CPE
#undef rpcc
#define LWPF_UNIT U(LIST)
#define LWPF_KERNELS _K(ALL) K(PRE) K(PUTBACK) K(L_CX) K(L_CY) K(T1) K(T2) K(MKCLUS) K(SUBC) K(MISS) K(NEXT)
#include "lwpf.h"
#define lwpf_start(x)
#define lwpf_stop(x)


// ======================== Header from gmx =======================
#include "gromacs/mdlib/nbnxn_internal.h"
#include "make_pairlist_sw64.h"
// #define NOINLINE

int nbl_old_nci;
int nbl_old_ncj;
int nbl_nci_host[NCPE];
int nbl_ncj_host[NCPE];
struct nbnxn_ci nbl_ci_host[MAX_NCI_64SLAVE]; // __attribute__((aligned(128)));
struct nbnxn_cj nbl_cj_host[MAX_NCJ_64SLAVE]; // __attribute__((aligned(128)));

static __thread_local int nbl_nci, l_nbl_nci;
static __thread_local int nbl_ncj, l_nbl_ncj;
static __thread_local struct nbnxn_ci *nbl_ci_host_ptr;
static __thread_local struct nbnxn_cj *nbl_cj_host_ptr;


#define REG_PUTR(var, dst) asm volatile ("putr %0,%1\n"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1\n"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile("getr %0\n":"=r"(var))
#define REG_GETC(var) asm volatile("getc %0\n":"=r"(var))

#define sqr(x) ((x)*(x))

#define max(a,b) ((a > b)? a:b)
#define min(a,b) ((a < b)? a:b)

#define pe_syn() {asm volatile("memb");}
#define pe_get(src, dst, size) { \
    volatile unsigned long get_reply; \
    get_reply = 0; \
    athread_get(PE_MODE, (src), (dst), (size), &get_reply, 0, 0, 0); \
    while(get_reply != 1); \
    pe_syn(); \
}
#define pe_put(src, dst, size) { \
    volatile unsigned long put_reply; \
    put_reply = 0; \
    athread_put(PE_MODE, (src), (dst), (size), &put_reply, 0, 0); \
    while(put_reply != 1); \
    pe_syn(); \
}
#define bcast_get(src, dst, size) { \
    volatile unsigned long get_reply; \
    if (NCPE == 64) { \
        athread_syn(ARRAY_SCOPE, 0xffff); \
        if (_MYID == 0) { \
            get_reply = 0; \
            athread_get(BCAST_MODE, (src), (dst), (size), &get_reply, 0xff, 0, 0); \
            while(get_reply != 1); \
        } \
        athread_syn(ARRAY_SCOPE, 0xffff); \
    } \
    else { \
        pe_get(src, dst, size); \
    } \
}


// ======================== Funcs only for MPE =======================

#if 0
#define OVER_ALLOC_FAC 1.19
#define over_alloc_large(n) (int)(OVER_ALLOC_FAC*(n) + 1000)
#define over_alloc_small(n) (int)(OVER_ALLOC_FAC*(n) + 8000)
#define srenew(ptr, nelem) \
    (ptr) = save_realloc(#ptr, __FILE__, __LINE__, (ptr), (nelem), sizeof(*(ptr)))

void *save_realloc(const char *name, const char *file, int line,
                   void *ptr, size_t nelem, size_t elsize);

static void nb_realloc_ci(nbnxn_pairlist_t *nbl, int n)
{
    nbl->ci_nalloc = over_alloc_small(n);
    nbnxn_realloc_void((void **)&nbl->ci,
                       nbl->nci*sizeof(*nbl->ci),
                       nbl->ci_nalloc*sizeof(*nbl->ci),
                       nbl->alloc, nbl->free);
}

static void check_subcell_list_space_simple(nbnxn_pairlist_t *nbl,
                                            int               ncell)
{
    int cj_max;

    cj_max = nbl->ncj + ncell;

    if (cj_max > nbl->cj_nalloc)
    {
        nbl->cj_nalloc = over_alloc_small(cj_max);
        nbnxn_realloc_void((void **)&nbl->cj,
                           nbl->ncj*sizeof(*nbl->cj),
                           nbl->cj_nalloc*sizeof(*nbl->cj),
                           nbl->alloc, nbl->free);
    }
}

#endif

// =========================================================================================

#define P_CPE_SIZE  16
#define T_CPE_SIZE 8
#define TU_SIZE 2
#ifdef NOINLINE
static __attribute__((__noinline__)) void* miss_get( int *c_t, void *cache_data, void *data, int data_size, int id)
#else
static void* miss_get( int *c_t, void *cache_data, void *data, int data_size, int id)
#endif
{
    lwpf_start(MISS);
    int t_s = T_CPE_SIZE;
    int p_s = P_CPE_SIZE;
    int tu_size = TU_SIZE;
    int p_ps = id / t_s;
    int cache_p_ps = p_ps % p_s;
    int cache_t_ps = id - p_ps * t_s;
    char *C_D = (char*)(cache_data);
    char *D = (char*)(data);
    int s = data_size / sizeof(char);
    int ps = cache_p_ps * (tu_size + 1);
    int get_cache_ps = c_t[ps + tu_size];
    int i,j,k;
    for(i = 0;i < tu_size;i++)
    {
        if(c_t[ps + i] == p_ps)
        {
            lwpf_stop(MISS);
            return (void*)(&C_D[(((cache_p_ps * tu_size + i) * t_s) + cache_t_ps) * s]);
        }
    }

    {
        volatile unsigned long get_reply, put_reply;
        get_reply = 0;
        athread_get(PE_MODE, D + (p_ps * t_s) * s, C_D + ((cache_p_ps * tu_size + get_cache_ps ) * t_s) * s, 
                                    data_size * t_s, (void*)&get_reply,0,0,0);
        while(get_reply != 1);
        c_t[ps + get_cache_ps] = p_ps;
    }
    c_t[ps + tu_size] = (get_cache_ps + 1) % tu_size;
    lwpf_stop(MISS);
    return (void*)(&C_D[(((cache_p_ps * tu_size + get_cache_ps) * t_s) + cache_t_ps) * s]);
}


static inline void copy_rvec(const rvec a, rvec b)
{
    b[XX] = a[XX];
    b[YY] = a[YY];
    b[ZZ] = a[ZZ];
}

static inline void copy_mat(gmx_cxx_const matrix a, matrix b)
{
    copy_rvec(a[XX], b[XX]);
    copy_rvec(a[YY], b[YY]);
    copy_rvec(a[ZZ], b[ZZ]);
}

static void sort_cj_excl(nbnxn_cj_t *cj, int ncj,
                         nbnxn_list_work_t *work)
{
    int jnew, j;

    /* Make a list of the j-cells involving exclusions */
    jnew = 0;
    for (j = 0; j < ncj; j++)
    {
        if (cj[j].excl != NBNXN_INTERACTION_MASK_ALL)
        {
            work->cj[jnew++] = cj[j];
        }
    }
    /* Check if there are exclusions at all or not just the first entry */
    if (!((jnew == 0) ||
          (jnew == 1 && cj[0].excl != NBNXN_INTERACTION_MASK_ALL)))
    {
        for (j = 0; j < ncj; j++)
        {
            if (cj[j].excl == NBNXN_INTERACTION_MASK_ALL)
            {
                work->cj[jnew++] = cj[j];
            }
        }
        for (j = 0; j < ncj; j++)
        {
            cj[j] = work->cj[j];
        }
    }
}

static inline void set_icell_bb_simple(const nbnxn_bb_t *bb, int ci,
                                           real shx, real shy, real shz,
                                           nbnxn_bb_t *bb_ci)
{
    bb_ci->lower[BB_X] = bb[ci].lower[BB_X] + shx;
    bb_ci->lower[BB_Y] = bb[ci].lower[BB_Y] + shy;
    bb_ci->lower[BB_Z] = bb[ci].lower[BB_Z] + shz;
    bb_ci->upper[BB_X] = bb[ci].upper[BB_X] + shx;
    bb_ci->upper[BB_Y] = bb[ci].upper[BB_Y] + shy;
    bb_ci->upper[BB_Z] = bb[ci].upper[BB_Z] + shz;
}

static unsigned int get_imask(gmx_bool rdiag, int ci, int cj)
{
    return (rdiag && ci == cj ? NBNXN_INTERACTION_MASK_DIAG : NBNXN_INTERACTION_MASK_ALL);
}

#ifdef NOINLINE
static __attribute__((__noinline__)) float subc_bb_dist2(int si, const nbnxn_bb_t *bb_i_ci,
                           int csj, const nbnxn_bb_t *bb_j_all)
#else
static float subc_bb_dist2(int si, const nbnxn_bb_t *bb_i_ci,
                           int csj, const nbnxn_bb_t *bb_j_all)
#endif
{
    lwpf_start(SUBC);
    const nbnxn_bb_t *bb_i, *bb_j;
    float             d2;
    float             dl, dh, dm, dm0;

    bb_i = bb_i_ci  +  si;
    bb_j = bb_j_all + csj;

    d2 = 0;

    dl  = bb_i->lower[BB_X] - bb_j->upper[BB_X];
    dh  = bb_j->lower[BB_X] - bb_i->upper[BB_X];
    dm  = max(dl, dh);
    dm0 = max(dm, 0);
    d2 += dm0*dm0;

    dl  = bb_i->lower[BB_Y] - bb_j->upper[BB_Y];
    dh  = bb_j->lower[BB_Y] - bb_i->upper[BB_Y];
    dm  = max(dl, dh);
    dm0 = max(dm, 0);
    d2 += dm0*dm0;

    dl  = bb_i->lower[BB_Z] - bb_j->upper[BB_Z];
    dh  = bb_j->lower[BB_Z] - bb_i->upper[BB_Z];
    dm  = max(dl, dh);
    dm0 = max(dm, 0);
    d2 += dm0*dm0;

    lwpf_stop(SUBC);
    return d2;
}

static void get_cell_range(real b0, real b1,
                           int nc, real c0, real s, real invs,
                           real d2, real r2, int *cf, int *cl)
{
    *cf = max((int)((b0 - c0)*invs), 0);

    while (*cf > 0 && d2 + sqr((b0 - c0) - (*cf-1+1)*s) < r2)
    {
        (*cf)--;
    }

    *cl = min((int)((b1 - c0)*invs), nc-1);
    while (*cl < nc-1 && d2 + sqr((*cl+1)*s - (b1 - c0)) < r2)
    {
        (*cl)++;
    }
}

static void icell_set_x_simple(int ci,
                               real shx, real shy, real shz,
                               int gmx_unused na_c,
                               int stride, const real *x,
                               nbnxn_list_work_t *work)
{
    int  ia, i;

    ia = ci*NBNXN_CPU_CLUSTER_I_SIZE;

    for (i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
    {
        work->x_ci[i*STRIDE_XYZ+XX] = x[(ia+i)*stride+XX] + shx;
        work->x_ci[i*STRIDE_XYZ+YY] = x[(ia+i)*stride+YY] + shy;
        work->x_ci[i*STRIDE_XYZ+ZZ] = x[(ia+i)*stride+ZZ] + shz;
    }
}

#ifdef NOINLINE
static __attribute__((__noinline__)) void make_cluster_list_simple(const nbnxn_grid_t *gridj,
                                     nbnxn_pairlist_t *nbl,
                                     int ci, int cjf, int cjl,
                                     gmx_bool remove_sub_diag,
                                     const real *x_j,
                                     real rl2, float rbb2,
                                     int *ndistc,
                                     real *x_ci_cache,
                                     int *x_ci_cache_ps,
                                     nbnxn_bb_t *bb_cache,
                                     int *bb_cache_ps       
                                     )
#else
static void make_cluster_list_simple(const nbnxn_grid_t *gridj,
                                     nbnxn_pairlist_t *nbl,
                                     int ci, int cjf, int cjl,
                                     gmx_bool remove_sub_diag,
                                     const real *x_j,
                                     real rl2, float rbb2,
                                     int *ndistc,
                                     real *x_ci_cache,
                                     int *x_ci_cache_ps,
                                     nbnxn_bb_t *bb_cache,
                                     int *bb_cache_ps       
                                     )
#endif
{
    lwpf_start(MKCLUS);
    const nbnxn_list_work_t *work;

    const nbnxn_bb_t        *bb_ci;
    const nbnxn_bb_t        *bb;
    const real              *x_ci;

    gmx_bool                 InRange;
    real                     d2;
    int                      cjf_gl, cjl_gl, cj;
    work = nbl->work;
    bb_ci = nbl->work->bb_ci;
    x_ci  = nbl->work->x_ci;

    real *lc_x_j;

    InRange = FALSE;
    while (!InRange && cjf <= cjl)
    {
        bb = miss_get( bb_cache_ps, bb_cache, gridj->bb, sizeof(nbnxn_bb_t), cjf);
        d2       = subc_bb_dist2(0, bb_ci, 0, bb);
        //d2       = subc_bb_dist2(0, bb_ci, cjf, gridj->bb);
        *ndistc += 2;

        /* Check if the distance is within the distance where
         * we use only the bounding box distance rbb,
         * or within the cut-off and there is at least one atom pair
         * within the cut-off.
         */


        if (d2 < rbb2)
        {
            InRange = TRUE;
        }
        else if (d2 < rl2)
        {
            int i, j;

            cjf_gl = gridj->cell0 + cjf;
            lc_x_j = miss_get( x_ci_cache_ps, x_ci_cache, x_j, sizeof(real) * NBNXN_CPU_CLUSTER_I_SIZE * STRIDE_XYZ, cjf_gl);
            for (i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE && !InRange; i++)
            {
                for (j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
                {
                    InRange = InRange ||
                        (sqr(x_ci[i*STRIDE_XYZ+XX] - lc_x_j[j*STRIDE_XYZ+XX]) +
                         sqr(x_ci[i*STRIDE_XYZ+YY] - lc_x_j[j*STRIDE_XYZ+YY]) +
                         sqr(x_ci[i*STRIDE_XYZ+ZZ] - lc_x_j[j*STRIDE_XYZ+ZZ]) < rl2);
                }
                //for (j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
                //{
                //    InRange = InRange ||
                //        (sqr(x_ci[i*STRIDE_XYZ+XX] - x_j[(cjf_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+XX]) +
                //         sqr(x_ci[i*STRIDE_XYZ+YY] - x_j[(cjf_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+YY]) +
                //         sqr(x_ci[i*STRIDE_XYZ+ZZ] - x_j[(cjf_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+ZZ]) < rl2);
                //}
            }
            *ndistc += NBNXN_CPU_CLUSTER_I_SIZE*NBNXN_CPU_CLUSTER_I_SIZE;
        }
        if (!InRange)
        {
            cjf++;
        }
    }
    if (!InRange)
    {
        lwpf_stop(MKCLUS);
        return;
    }

    InRange = FALSE;
    while (!InRange && cjl > cjf)
    {
        
        bb = miss_get( bb_cache_ps, bb_cache, gridj->bb, sizeof(nbnxn_bb_t), cjl);
        d2       = subc_bb_dist2(0, bb_ci, 0, bb);
        //d2       = subc_bb_dist2(0, bb_ci, cjl, gridj->bb);
        *ndistc += 2;

        /* Check if the distance is within the distance where
         * we use only the bounding box distance rbb,
         * or within the cut-off and there is at least one atom pair
         * within the cut-off.
         */



        if (d2 < rbb2)
        {
            InRange = TRUE;
        }
        else if (d2 < rl2)
        {
            int i, j;

            cjl_gl = gridj->cell0 + cjl;
            lc_x_j = miss_get( x_ci_cache_ps, x_ci_cache, x_j, sizeof(real) * NBNXN_CPU_CLUSTER_I_SIZE * STRIDE_XYZ, cjl_gl);
            for (i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE && !InRange; i++)
            {
                for (j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
                {
                    InRange = InRange ||
                        (sqr(x_ci[i*STRIDE_XYZ+XX] - lc_x_j[j*STRIDE_XYZ+XX]) +
                         sqr(x_ci[i*STRIDE_XYZ+YY] - lc_x_j[j*STRIDE_XYZ+YY]) +
                         sqr(x_ci[i*STRIDE_XYZ+ZZ] - lc_x_j[j*STRIDE_XYZ+ZZ]) < rl2);
                }
                //for (j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
                //{
                //    InRange = InRange ||
                //        (sqr(x_ci[i*STRIDE_XYZ+XX] - x_j[(cjl_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+XX]) +
                //         sqr(x_ci[i*STRIDE_XYZ+YY] - x_j[(cjl_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+YY]) +
                //         sqr(x_ci[i*STRIDE_XYZ+ZZ] - x_j[(cjl_gl*NBNXN_CPU_CLUSTER_I_SIZE+j)*STRIDE_XYZ+ZZ]) < rl2);
                //}
            }
            *ndistc += NBNXN_CPU_CLUSTER_I_SIZE*NBNXN_CPU_CLUSTER_I_SIZE;
        }
        if (!InRange)
        {
            cjl--;
        }
    }

    if (cjf <= cjl)
    {
        for (cj = cjf; cj <= cjl; cj++)
        {
            /* Store cj and the interaction mask */
            nbl->cj[l_nbl_ncj].cj   = gridj->cell0 + cj;
            nbl->cj[l_nbl_ncj].excl = get_imask(remove_sub_diag, ci, cj);
            l_nbl_ncj++;
        }
        /* Increase the closing index in i super-cell list */
        nbl->ci[l_nbl_nci].cj_ind_end = nbl_ncj + l_nbl_ncj;
    }
    lwpf_stop(MKCLUS);
}

static void set_ci_top_excls(const nbnxn_search_t nbs,
                             nbnxn_pairlist_t    *nbl,
                             gmx_bool             diagRemoved,
                             int                  na_ci_2log,
                             int                  na_cj_2log,
                             const nbnxn_ci_t    *nbl_ci,
                             const t_blocka      *excl)
{
    const int    *cell;
    int           ci;
    int           cj_ind_first, cj_ind_last;
    int           cj_first, cj_last;
    int           ndirect;
    int           i, ai, aj, si, eind, ge, se;
    int           found, cj_ind_0, cj_ind_1, cj_ind_m;
    int           cj_m;
    gmx_bool      Found_si;
    int           si_ind;
    nbnxn_excl_t *nbl_excl;
    int           inner_i, inner_e;

    cell = nbs->cell;

    if (nbl_ci->cj_ind_end == nbl_ci->cj_ind_start)
    {
        /* Empty list */
        return;
    }

    ci = nbl_ci->ci;

    cj_ind_first = nbl_ci->cj_ind_start;
    cj_ind_last  = nbl_ncj + l_nbl_ncj - 1;

    cj_first = nbl->cj[cj_ind_first - nbl_ncj].cj;
    cj_last  = nbl->cj[cj_ind_last - nbl_ncj].cj;

    /* Determine how many contiguous j-cells we have starting
     * from the first i-cell. This number can be used to directly
     * calculate j-cell indices for excluded atoms.
     */
    ndirect = 0;
    if (na_ci_2log == na_cj_2log)
    {
        while (cj_ind_first + ndirect <= cj_ind_last &&
               nbl->cj[cj_ind_first+ndirect - nbl_ncj].cj == ci + ndirect) // ****
        {
            ndirect++;
        }
    }

    /* Loop over the atoms in the i super-cell */
    for (i = 0; i < nbl->na_sc; i++)
    {
        ai = nbs->a[ci*nbl->na_sc+i];
        if (ai >= 0)
        {
            si  = (i>>na_ci_2log);

            /* Loop over the topology-based exclusions for this i-atom */
            for (eind = excl->index[ai]; eind < excl->index[ai+1]; eind++)
            {
                aj = excl->a[eind];

                if (aj == ai)
                {
                    /* The self exclusion are already set, save some time */
                    continue;
                }

                ge = cell[aj];

                /* Without shifts we only calculate interactions j>i
                 * for one-way pair-lists.
                 */
                if (diagRemoved && ge <= ci*nbl->na_sc + i)
                {
                    continue;
                }

                se = (ge >> na_cj_2log);

                /* Could the cluster se be in our list? */
                if (se >= cj_first && se <= cj_last)
                {
                    if (se < cj_first + ndirect)
                    {
                        /* We can calculate cj_ind directly from se */
                        found = cj_ind_first + se - cj_first;
                    }
                    else
                    {
                        /* Search for se using bisection */
                        found    = -1;
                        cj_ind_0 = cj_ind_first + ndirect;
                        cj_ind_1 = cj_ind_last + 1;
                        while (found == -1 && cj_ind_0 < cj_ind_1)
                        {
                            cj_ind_m = (cj_ind_0 + cj_ind_1)>>1;

                            cj_m = nbl->cj[cj_ind_m - nbl_ncj].cj;

                            if (se == cj_m)
                            {
                                found = cj_ind_m;
                            }
                            else if (se < cj_m)
                            {
                                cj_ind_1 = cj_ind_m;
                            }
                            else
                            {
                                cj_ind_0 = cj_ind_m + 1;
                            }
                        }
                    }

                    if (found >= 0)
                    {
                        inner_i = i  - (si << na_ci_2log);
                        inner_e = ge - (se << na_cj_2log);

                        nbl->cj[found - nbl_ncj].excl &= ~(1U<<((inner_i<<na_cj_2log) + inner_e));
                    }
                }
            }
        }
    }
}

static void close_ci_entry_simple(nbnxn_pairlist_t *nbl)
{
    int jlen;

    /* All content of the new ci entry have already been filled correctly,
     * we only need to increase the count here (for non empty lists).
     */
    jlen = nbl->ci[l_nbl_nci].cj_ind_end - nbl->ci[l_nbl_nci].cj_ind_start;
    if (jlen > 0)
    {
        sort_cj_excl(nbl->cj + (nbl->ci[l_nbl_nci].cj_ind_start - nbl_ncj), jlen, nbl->work);

        /* The counts below are used for non-bonded pair/flop counts
         * and should therefore match the available kernel setups.
         */
        if (!(nbl->ci[l_nbl_nci].shift & NBNXN_CI_DO_COUL(0)))
        {
            nbl->work->ncj_noq += jlen;
        }
        else if ((nbl->ci[l_nbl_nci].shift & NBNXN_CI_HALF_LJ(0)) ||
                 !(nbl->ci[l_nbl_nci].shift & NBNXN_CI_DO_LJ(0)))
        {
            nbl->work->ncj_hlj += jlen;
        }

        l_nbl_nci++;
    }

    lwpf_start(PUTBACK);
    if (l_nbl_nci == MAX_LNCI) {
        pe_put(nbl->ci, nbl_ci_host_ptr + _MYID * MAX_NCI_1SLAVE + nbl_nci, sizeof(struct nbnxn_ci) * l_nbl_nci);
        nbl_nci += l_nbl_nci;
        l_nbl_nci = 0;
    }
    lwpf_stop(PUTBACK);
}

static void new_ci_entry(nbnxn_pairlist_t *nbl, int ci, int shift, int flags)
{
    nbl->ci[l_nbl_nci].ci            = ci;
    nbl->ci[l_nbl_nci].shift         = shift;
    /* Store the interaction flags along with the shift */
    nbl->ci[l_nbl_nci].shift        |= flags;
    nbl->ci[l_nbl_nci].cj_ind_start  = l_nbl_ncj + nbl_ncj;
    nbl->ci[l_nbl_nci].cj_ind_end    = l_nbl_ncj + nbl_ncj;
}

#ifdef NOINLINE
static gmx_bool __attribute__((__noinline__)) next_ci(const nbnxn_grid_t *grid,
                        int conv,
                        int nth, int ci_block,
                        int *ci_x, int *ci_y,
                        int *ci_b, int *ci)
#else
static gmx_bool next_ci(const nbnxn_grid_t *grid,
                        int conv,
                        int nth, int ci_block,
                        int *ci_x, int *ci_y,
                        int *ci_b, int *ci)
#endif
{
    lwpf_start(NEXT);
    (*ci_b)++;
    (*ci)++;

    if (*ci_b == ci_block)
    {
        /* Jump to the next block assigned to this task */
        *ci   += (nth - 1)*ci_block;
        *ci_b  = 0;
    }

    if (*ci >= grid->nc*conv)
    {
        lwpf_stop(NEXT);
        return FALSE;
    }

    while (*ci >= grid->cxy_ind[*ci_x*grid->ncy + *ci_y + 1]*conv)
    {
        *ci_y += 1;
        if (*ci_y == grid->ncy)
        {
            *ci_x += 1;
            *ci_y  = 0;
        }
    }
    lwpf_stop(NEXT);

    return TRUE;
}

void nbnxn_make_pairlist_part_slave(make_pairlist_ref *params)
{
    if (_MYID >= NCPE)
        return;

    lwpf_start(ALL);
    lwpf_start(PRE);

    make_pairlist_ref l_params;
    bcast_get(params, &l_params, sizeof(make_pairlist_ref));
    params = &l_params;

    const nbnxn_grid_t *gridi;
    const nbnxn_grid_t *gridj;
    #define MAX_CXY_IND (2000)
    int l_gridi_cxy_ind[MAX_CXY_IND];
    int l_gridj_cxy_ind[MAX_CXY_IND];
    nbnxn_grid_t l_gridi, l_gridj;
    bcast_get(params->gridi, &l_gridi, sizeof(nbnxn_grid_t));
    gridi = &l_gridi;
    // bcast_get(gridi->cxy_ind, l_gridi_cxy_ind, sizeof(int) * gridi->cxy_nalloc);
    // if (l_gridi.cxy_nalloc > MAX_CXY_IND)
    //     exit(-1);
    // if (l_gridj.cxy_nalloc > MAX_CXY_IND)
    //     exit(-1);
    if (params->gridi == params->gridj) {
        gridj = gridi;
    }
    else {
        bcast_get(params->gridj, &l_gridj, sizeof(nbnxn_grid_t));
        gridj = &l_gridj;
        // bcast_get(gridj->cxy_ind, l_gridj_cxy_ind, sizeof(int) * gridj->cxy_nalloc);
    }
    #undef MAX_CXY_IND

    const int *flags_i = params->flags_i; // gld ?

    const nbnxn_bb_t *bb_i = params->bb_i;
    struct nbnxn_bb bb_i_ci;

    struct {
        int xstride;
        real *x;
    }l_nbat, *nbat;
    l_nbat.xstride = params->nbat->xstride;
    l_nbat.x       = params->nbat->x      ; // cache
    nbat = &l_nbat;

    nbnxn_pairlist_t l_nbl;
    nbnxn_list_work_t l_work;
    real l_work_x_ci[NBNXN_CPU_CLUSTER_I_SIZE * DIM];
    struct nbnxn_bb l_work_bb_ci;
    struct nbnxn_cj l_work_cj[MAX_1NCJ];
    struct nbnxn_ci l_ci[MAX_LNCI];
    struct nbnxn_cj l_cj[MAX_1NCJ];
    bcast_get(params->nbl, &l_nbl, sizeof(nbnxn_pairlist_t));
    nbnxn_pairlist_t *nbl = &l_nbl;
    nbl->work = &l_work;
    nbl->work->bb_ci = &l_work_bb_ci;
    nbl->work->x_ci = l_work_x_ci;
    nbl->work->cj = (struct nbnxn_cj*) l_work_cj;
    nbl->ci = l_ci;
    nbl->cj = l_cj;
    nbl_ci_host_ptr = nbl_ci_host;
    nbl_cj_host_ptr = nbl_cj_host;


    const t_blocka *excl = params->excl;
    int conv_i = params->conv_i, nth = params->nth, th = params->th;
    int ci_block = params->ci_block, cell0_i = params->cell0_i, nb_kernel_type = params->nb_kernel_type;
    int gridj_flag_shift = params->gridj_flag_shift, gridi_flag_shift = params->gridi_flag_shift;
    int  na_cj_2log = params->na_cj_2log, nsubpair_max = params->nsubpair_max, progBal = params->progBal;
    float rbb2 = params->rbb2;

    const float *bbcz_i = params->bbcz_i;

    const float *bbcz_j = params->bbcz_j; // cache

    float nsubpair_tot_est = params->nsubpair_tot_est;
    real rl_fep2 = params->rl_fep2;
    real rl2 = params->rl2;

    struct nbnxn_search l_nbs;
    bcast_get(params->nbs, &l_nbs, sizeof(struct nbnxn_search));
    const  nbnxn_search_t nbs = &l_nbs;

    ivec shp;
    // bcast_get(&(params->shp), &shp, sizeof(real) * DIM);
    shp[XX] = params->shp[XX];
    shp[YY] = params->shp[YY];
    shp[ZZ] = params->shp[ZZ];

    struct mklist_ind_t my_ind_start;
    pe_get(params->ind_start + _MYID, &my_ind_start, sizeof(struct mklist_ind_t));

    int ci_start   = my_ind_start.ci;
    int ci_x_start = my_ind_start.ci_x;
    int ci_y_start = my_ind_start.ci_y;
    int nloops     = my_ind_start.n;

    int mpi_rank = params->mpi_rank;

    int ci_b = -1;
    int ci   = ci_start;
    int ci_x = ci_x_start;
    int ci_y = ci_y_start;
    real              d2cx, d2z, d2z_cx, d2z_cy, d2zx, d2zxy, d2xy;
    int               cxf, cxl, cyf, cyf_x, cyl;
    real              bx0, bx1, by0, by1, bz0, bz1;
    int               c0, c1, cs, cf, cl;
    int               ci_xy, cj;
    int               tx, ty, tz;
    int               shift;
    real              shx, shy, shz;
    real              bz1_frac;
    int               cx, cy;
    int               ndistc = 0;
    int               ncpcheck = 0;

    float bbcz_j_cache[P_CPE_SIZE * T_CPE_SIZE * NNBSBB_D * TU_SIZE];
    int bbcz_j_cache_ps[P_CPE_SIZE * T_CPE_SIZE][TU_SIZE + 1];
    int i, j;
    for(i = 0; i <  P_CPE_SIZE * T_CPE_SIZE;i++)
    {
        for(j = 0;j < TU_SIZE;j++)
            bbcz_j_cache_ps[i][j] = -1;
        bbcz_j_cache_ps[i][TU_SIZE] = 0;
    }


    real x_ci_cache[NBNXN_CPU_CLUSTER_I_SIZE * STRIDE_XYZ * P_CPE_SIZE * T_CPE_SIZE * TU_SIZE];
    int x_ci_cache_ps[P_CPE_SIZE * T_CPE_SIZE][TU_SIZE + 1];
    for(i = 0; i <  P_CPE_SIZE * T_CPE_SIZE;i++)
    {
        for(j = 0;j < TU_SIZE;j++)
            x_ci_cache_ps[i][j] = -1;
        x_ci_cache_ps[i][TU_SIZE] = 0;
    }


    nbnxn_bb_t bb_cache[P_CPE_SIZE * T_CPE_SIZE * NNBSBB_D * TU_SIZE];
    int bb_cache_ps[P_CPE_SIZE * T_CPE_SIZE][TU_SIZE + 1];
    for(i = 0; i <  P_CPE_SIZE * T_CPE_SIZE;i++)
    {
        for(j = 0;j < TU_SIZE;j++)
            bb_cache_ps[i][j] = -1;
        bb_cache_ps[i][TU_SIZE] = 0;
    }

    int cxy_indj[100];

    matrix            box;
    copy_mat(nbs->box, box);

    if (_MYID == 0) {
        nbl_old_nci = nbl->nci;
        nbl_old_ncj = nbl->ncj;
    }
    nbl_nci = 0;
    nbl_ncj = 0;
    l_nbl_nci = 0;
    l_nbl_ncj = 0;

    lwpf_stop(PRE);

    for (; nloops; next_ci(gridi, conv_i, nth, ci_block, &ci_x, &ci_y, &ci_b, &ci))
    {
        if (flags_i[ci] == 0)
        {
            continue;
        }
        --nloops;

        d2cx = 0;
        pe_get(bb_i + ci, &bb_i_ci, sizeof(struct nbnxn_bb));
        if (gridj != gridi && shp[XX] == 0)
        {
            bx1 = bb_i_ci.upper[BB_X];

            if (bx1 < gridj->c0[XX])
            {
                d2cx = sqr(gridj->c0[XX] - bx1);

                if (d2cx >= rl2)
                {
                    continue;
                }
            }
        }

        ci_xy = ci_x*gridi->ncy + ci_y;

        /* Loop over shift vectors in three dimensions */
        for (tz = -shp[ZZ]; tz <= shp[ZZ]; tz++)
        {
            shz = tz*box[ZZ][ZZ];

            bz0 = bbcz_i[ci*NNBSBB_D  ] + shz;
            bz1 = bbcz_i[ci*NNBSBB_D+1] + shz;

            if (tz == 0)
            {
                d2z = 0;
            }
            else if (tz < 0)
            {
                d2z = sqr(bz1);
            }
            else
            {
                d2z = sqr(bz0 - box[ZZ][ZZ]);
            }

            d2z_cx = d2z + d2cx;

            if (d2z_cx >= rl2)
            {
                continue;
            }

            bz1_frac =
                bz1/((real)(gridi->cxy_ind[ci_xy+1] - gridi->cxy_ind[ci_xy]));
            if (bz1_frac < 0)
            {
                bz1_frac = 0;
            }
            /* The check with bz1_frac close to or larger than 1 comes later */

            for (ty = -shp[YY]; ty <= shp[YY]; ty++)
            {
                shy = ty*box[YY][YY] + tz*box[ZZ][YY];

                by0 = bb_i_ci.lower[BB_Y] + shy;
                by1 = bb_i_ci.upper[BB_Y] + shy;

                get_cell_range(by0, by1,
                               gridj->ncy, gridj->c0[YY], gridj->sy, gridj->inv_sy,
                               d2z_cx, rl2,
                               &cyf, &cyl);

                if (cyf > cyl)
                {
                    continue;
                }

                d2z_cy = d2z;
                if (by1 < gridj->c0[YY])
                {
                    d2z_cy += sqr(gridj->c0[YY] - by1);
                }
                else if (by0 > gridj->c1[YY])
                {
                    d2z_cy += sqr(by0 - gridj->c1[YY]);
                }

                for (tx = -shp[XX]; tx <= shp[XX]; tx++)
                {
                    shift = XYZ2IS(tx, ty, tz);

                    if (gridi == gridj && shift > CENTRAL)
                    {
                        continue;
                    }

                    shx = tx*box[XX][XX] + ty*box[YY][XX] + tz*box[ZZ][XX];

                    bx0 = bb_i_ci.lower[BB_X] + shx;
                    bx1 = bb_i_ci.upper[BB_X] + shx;

                    get_cell_range(bx0, bx1,
                                   gridj->ncx, gridj->c0[XX], gridj->sx, gridj->inv_sx,
                                   d2z_cy, rl2,
                                   &cxf, &cxl);

                    if (cxf > cxl)
                    {
                        continue;
                    }

                    new_ci_entry(nbl, cell0_i+ci, shift, flags_i[ci]);

                    if (shift == CENTRAL && gridi == gridj &&
                        cxf < ci_x)
                    {
                        /* Leave the pairs with i > j.
                         * x is the major index, so skip half of it.
                         */
                        cxf = ci_x;
                    }

                    set_icell_bb_simple(&bb_i_ci, 0, shx, shy, shz,
                            nbl->work->bb_ci);

                    icell_set_x_simple(cell0_i+ci, shx, shy, shz,
                                     gridi->na_c, nbat->xstride, nbat->x,
                                     nbl->work);

                    lwpf_start(L_CX);
                    for (cx = cxf; cx <= cxl; cx++)
                    {
                        d2zx = d2z;
                        if (gridj->c0[XX] + cx*gridj->sx > bx1)
                        {
                            d2zx += sqr(gridj->c0[XX] + cx*gridj->sx - bx1);
                        }
                        else if (gridj->c0[XX] + (cx+1)*gridj->sx < bx0)
                        {
                            d2zx += sqr(gridj->c0[XX] + (cx+1)*gridj->sx - bx0);
                        }

                        if (gridi == gridj &&
                            cx == 0 && shift == CENTRAL && cyf < ci_y)

                        {
                            /* Leave the pairs with i > j.
                             * Skip half of y when i and j have the same x.
                             */
                            cyf_x = ci_y;
                        }
                        else
                        {
                            cyf_x = cyf;
                        }

                        lwpf_start(L_CY);
                        int jsize = cyl - cyf_x + 1 + 1;
                        pe_get( &gridj->cxy_ind[cx*gridj->ncy + cyf_x], cxy_indj, sizeof(real) * jsize);
                        for (cy = cyf_x; cy <= cyl; cy++)
                        {
                            
                            c0 = cxy_indj[cy     - cyf_x];
                            c1 = cxy_indj[cy + 1 - cyf_x];
                            //c0 = gridj->cxy_ind[cx*gridj->ncy+cy];
                            //c1 = gridj->cxy_ind[cx*gridj->ncy+cy+1];

                            if (gridi == gridj &&
                                shift == CENTRAL && c0 < ci)
                            {
                                c0 = ci;
                            }

                            d2zxy = d2zx;
                            if (gridj->c0[YY] + cy*gridj->sy > by1)
                            {
                                d2zxy += sqr(gridj->c0[YY] + cy*gridj->sy - by1);
                            }
                            else if (gridj->c0[YY] + (cy+1)*gridj->sy < by0)
                            {
                                d2zxy += sqr(gridj->c0[YY] + (cy+1)*gridj->sy - by0);
                            }
                            if (c1 > c0 && d2zxy < rl2)
                            {
                                cs = c0 + (int)(bz1_frac*(c1 - c0));
                                if (cs >= c1)
                                {
                                    cs = c1 - 1;
                                }

                                d2xy = d2zxy - d2z;

                                /* Find the lowest cell that can possibly
                                 * be within range.
                                 */
                                cf = cs;

                                float *bbcz = miss_get( bbcz_j_cache_ps, bbcz_j_cache, bbcz_j, sizeof(float) * NNBSBB_D, cf);
                                lwpf_start(T1);
                                while (cf > c0 &&
                                       (bbcz[1] >= bz0 ||
                                        d2xy + sqr(bbcz[1] - bz0) < rl2))
                                {
                                    cf--;
                                    lwpf_start(T2);
                                    bbcz = miss_get( bbcz_j_cache_ps, bbcz_j_cache, bbcz_j, sizeof(float) * NNBSBB_D, cf);
                                    lwpf_stop(T2);
                                }
                                lwpf_stop(T1);


                                    
                                //while (cf > c0 &&
                                //       (bbcz_j[cf*NNBSBB_D+1] >= bz0 ||
                                //        d2xy + sqr(bbcz_j[cf*NNBSBB_D+1] - bz0) < rl2))
                                //{
                                //    cf--;
                                //}

                                /* Find the highest cell that can possibly
                                 * be within range.
                                 */
                                cl = cs;
                                

                                bbcz = miss_get( bbcz_j_cache_ps, bbcz_j_cache, bbcz_j, sizeof(float) * NNBSBB_D, cl);
                                lwpf_start(T1);
                                while (cl < c1-1 &&
                                       (bbcz[0] <= bz1 ||
                                        d2xy + sqr(bbcz[0] - bz1) < rl2))
                                {
                                    cl++;
                                    lwpf_start(T2);
                                    bbcz = miss_get( bbcz_j_cache_ps, bbcz_j_cache, bbcz_j, sizeof(float) * NNBSBB_D, cl);
                                    lwpf_stop(T2);
                                }
                                lwpf_stop(T1);


                                
                                
                                //while (cl < c1-1 &&
                                //       (bbcz_j[cl*NNBSBB_D] <= bz1 ||
                                //        d2xy + sqr(bbcz_j[cl*NNBSBB_D] - bz1) < rl2))
                                //{
                                //    cl++;
                                //}

                                if (gridi == gridj)
                                {
                                    /* We want each atom/cell pair only once,
                                     * only use cj >= ci.
                                     */
                                    if (shift == CENTRAL)
                                    {
                                        cf = max(cf, ci);
                                    }
                                }

                                if (cf <= cl)
                                {

                                    make_cluster_list_simple(gridj,
                                            nbl, ci, cf, cl,
                                            (gridi == gridj && shift == CENTRAL),
                                            nbat->x,
                                            rl2, rbb2,
                                            &ndistc,
                                            x_ci_cache,
                                            x_ci_cache_ps,
                                            bb_cache,
                                            bb_cache_ps
                                            );

                                    ncpcheck += cl - cf + 1;

                                }
                            }
                        }
                        lwpf_stop(L_CY);
                    }
                    lwpf_stop(L_CX);

                    /* Set the exclusions for this ci list */
                        set_ci_top_excls(nbs,
                                         nbl,
                                         shift == CENTRAL && gridi == gridj,
                                         gridj->na_c_2log,
                                         na_cj_2log,
                                         &(nbl->ci[l_nbl_nci]),
                                         excl);

                    /* Close this ci list */
                        close_ci_entry_simple(nbl);

                        lwpf_start(PUTBACK);
                        if (l_nbl_ncj > 0) {
                            pe_put(nbl->cj, nbl_cj_host_ptr + _MYID * MAX_NCJ_1SLAVE + nbl_ncj, sizeof(struct nbnxn_cj) * l_nbl_ncj);
                            nbl_ncj += l_nbl_ncj;
                            l_nbl_ncj = 0;
                        }
                        lwpf_stop(PUTBACK);

                }
            }
        }
    }

    /* flush left list */
    lwpf_start(PUTBACK);
    if (l_nbl_nci > 0) {
        pe_put(nbl->ci, nbl_ci_host_ptr + _MYID * MAX_NCI_1SLAVE + nbl_nci, sizeof(struct nbnxn_ci) * l_nbl_nci);
        nbl_nci += l_nbl_nci;
        l_nbl_nci = 0;
    }
    if (l_nbl_ncj > 0) {
        pe_put(nbl->cj, nbl_cj_host_ptr + _MYID * MAX_NCJ_1SLAVE + nbl_ncj, sizeof(struct nbnxn_cj) * l_nbl_ncj);
        nbl_ncj += l_nbl_ncj;
        l_nbl_ncj = 0;
    }
    lwpf_stop(PUTBACK);
    
    nbl_nci_host[_MYID] = nbl_nci;
    nbl_ncj_host[_MYID] = nbl_ncj;

    lwpf_stop(ALL);
}


void mklist_gather_cicj_slave(mklist_gather_ref *params)
{
    if (_MYID >= NCPE)
        return;

    mklist_gather_ref l_params;
    bcast_get(params, &l_params, sizeof(mklist_gather_ref));
    params = &l_params;

    int i, j;
    int nci_off[NCPE], ncj_off[NCPE];

    nci_off[0] = 0;
    ncj_off[0] = 0;
    for (i = 1; i < NCPE; ++i) {
        nci_off[i] = nci_off[i-1] + params->nci[i-1];
        ncj_off[i] = ncj_off[i-1] + params->ncj[i-1];
    }

    int nbl_nci_base = params->old_nci + nci_off[_MYID];
    int nbl_ncj_base = params->old_ncj + ncj_off[_MYID];

    struct nbnxn_ci ci64[MAX_LNCI];
    struct nbnxn_cj cj64[MAX_1NCJ];

    for (i = 0; i < params->nci[_MYID]; i += MAX_LNCI) {
        int i_size = (i + MAX_LNCI < params->nci[_MYID] ? MAX_LNCI : params->nci[_MYID] - i);
        pe_get(params->ci64 + _MYID * MAX_NCI_1SLAVE + i, ci64, sizeof(struct nbnxn_ci) * i_size);
        for (j = 0; j < i_size; ++j) {
            ci64[j].cj_ind_start += nbl_ncj_base;
            ci64[j].cj_ind_end   += nbl_ncj_base;
        }
        pe_put(ci64, params->ci + nbl_nci_base + i, sizeof(struct nbnxn_ci) * i_size);
    }

    for (j = 0; j < params->ncj[_MYID]; j += MAX_1NCJ) {
        int j_size = (j + MAX_1NCJ < params->ncj[_MYID] ? MAX_1NCJ : params->ncj[_MYID] - j);
        pe_get(params->cj64 + _MYID * MAX_NCJ_1SLAVE + j, cj64, sizeof(struct nbnxn_cj) * j_size);
        pe_put(cj64, params->cj + nbl_ncj_base + j, sizeof(struct nbnxn_cj) * j_size);
    }
}

