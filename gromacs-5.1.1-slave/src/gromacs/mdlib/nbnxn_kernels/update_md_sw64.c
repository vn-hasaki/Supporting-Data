#include <slave.h>
#include <simd.h>

#include "gromacs/mdlib/nbnxn_internal.h"
#include "make_pairlist_sw64.h"

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

static inline void rvec_sub(const rvec a, const rvec b, rvec c)
{
    real x, y, z;

    x = a[XX]-b[XX];
    y = a[YY]-b[YY];
    z = a[ZZ]-b[ZZ];

    c[XX] = x;
    c[YY] = y;
    c[ZZ] = z;
}

static inline real iprod(const rvec a, const rvec b)
{
    return (a[XX]*b[XX]+a[YY]*b[YY]+a[ZZ]*b[ZZ]);
}


#undef NCPE
#define NCPE (64)
#define N_DO_UPDATE_MD (512)
typedef struct do_update_md_ref{
    unsigned short eptVSite, eptShell;
    rvec *f, *v, *xprime, *x;
    real *invmass;
    real tcstat0_lambda;
    matrix M;
    unsigned short *ptype;
    double dt;
    int start, nrend;
    int nstpcouple;

    rvec gstat0_u;
    rvec accel0;
    ivec nFreeze0;
    
}do_update_md_ref;

void do_update_md_slave(do_update_md_ref *params)
{
    if (_MYID >= NCPE)
        return;

    do_update_md_ref l_params;
    bcast_get(params, &l_params, sizeof(do_update_md_ref));
    params = &l_params;

    double dt = params->dt;
    unsigned short eptVSite = params->eptVSite;
    unsigned short eptShell = params->eptShell;
    matrix M;
    memcpy(M, params->M, sizeof(matrix));
    int start = params->start;
    int nrend = params->nrend;
    int nstpcouple = params->nstpcouple;
    real tcstat0_lambda = params->tcstat0_lambda;

    rvec gstat0_u;
    memcpy(gstat0_u, params->gstat0_u, sizeof(rvec));
    rvec accel0;
    memcpy(accel0, params->accel0, sizeof(rvec));
    rvec nFreeze0;
    memcpy(nFreeze0, params->nFreeze0, sizeof(ivec));

    rvec l_f[N_DO_UPDATE_MD];
    rvec l_v[N_DO_UPDATE_MD];
    rvec l_xprime[N_DO_UPDATE_MD];
    rvec l_x[N_DO_UPDATE_MD];
    real l_invmass[N_DO_UPDATE_MD];
    unsigned short l_ptype[N_DO_UPDATE_MD + 2];
    unsigned short *ptype;
    rvec *f       = (rvec*) l_f;
    rvec *v       = (rvec*) l_v;
    rvec *xprime  = (rvec*) l_xprime;
    rvec *x       = (rvec*) l_x     ;
    real *invmass = l_invmass;

    double imass;
    rvec vrel;
    real vxi = 0, lg, u, vn, vnrel;

    int n, ni, d;
    int my_start = start + (nrend - start + NCPE - 1) / NCPE * _MYID;
    int my_end   = start + (nrend - start + NCPE - 1) / NCPE * (_MYID + 1);
    if (my_end > nrend)
        my_end = nrend;

    for (ni = my_start; ni < my_end; ni += N_DO_UPDATE_MD)
    {
        int ni_size = ni + N_DO_UPDATE_MD < my_end ? N_DO_UPDATE_MD : my_end - ni;
        pe_get(params->invmass + ni, invmass, sizeof(real) * ni_size);
        pe_get(params->x + ni, x, sizeof(rvec) * ni_size);
        pe_get(params->v + ni, v, sizeof(rvec) * ni_size);
        pe_get(params->f + ni, f, sizeof(rvec) * ni_size);
        {
            long ptype_addr = (long) (params->ptype + ni);
            if (ptype_addr & 2) {
                ptype_addr -= 2;
                ptype = l_ptype + 1;
            }
            else {
                ptype = l_ptype;
            }
            pe_get((void *) ptype_addr, l_ptype, sizeof(short) * ((ni_size + 2) / 2 * 2));
        }

        for (n = 0; n < ni_size; ++n) {
            imass = invmass[n];
            rvec_sub(v[n], gstat0_u, vrel);
            lg = tcstat0_lambda;

            for (d = 0; d < DIM; d++)
            {
                if ((ptype[n] != eptVSite) && (ptype[n] != eptShell) && !nFreeze0[d])
                {
                    vnrel = (lg*vrel[d] + dt*(imass*f[n][d] - 0.5*vxi*vrel[d]
                                - nstpcouple*iprod(M[d], vrel)))/(1 + 0.5*vxi*dt);
                    /* do not scale the mean velocities u */
                    vn             = gstat0_u[d] + accel0[d]*dt + vnrel;
                    v[n][d]        = vn;
                    xprime[n][d]   = x[n][d]+vn*dt;
                }
                else
                {
                    v[n][d]        = 0.0;
                    xprime[n][d]   = x[n][d];
                }
            }
        }

        pe_put(v, params->v + ni, sizeof(rvec) * ni_size);
        pe_put(xprime, params->xprime + ni, sizeof(rvec) * ni_size);

    }
}
