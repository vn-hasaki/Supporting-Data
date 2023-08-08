#include <slave.h>
#include <stdio.h>
#include <stdlib.h>
#include <simd.h>
#define CPE
#define LWPF_UNIT U(CASE2)
#define LWPF_KERNELS K(ALL) K(PRE) K(GET_CACHE) K(PUT_CACHE)  K(CACU)  K(COMPUTE1)  K(COMPUTE2) K(COMPUTE3) K(CACU1)  K(CACU2)
#undef rpcc
#include "param.h"
#include "cpe_print.h"
#include "lwpf2.h"
#include "/home/export/online3/swmore/opensource/cal/cal.h"
#define lwpf_start(x)
#define lwpf_stop(x)

#ifdef atom_cache_tail
#undef atom_cache_tail
#undef C_T
#endif
#define atom_cache_tail 3
#define C_T atom_cache_tail


#define REG_PUTR(var, dst) asm volatile ("putr %0,%1\n"::"r"(var),"r"(dst))
#define REG_PUTC(var, dst) asm volatile ("putc %0,%1\n"::"r"(var),"r"(dst))
#define REG_GETR(var) asm volatile("getr %0\n":"=r"(var))
#define REG_GETC(var) asm volatile("getc %0\n":"=r"(var))


#define test_id  1
#define TABQ_SIZE 1073
#define type_size 41
#define T_S  type_size
#define TYPE_CACHE_SIZE 3
#define T_C_S  TYPE_CACHE_SIZE

extern long f_cache_bmap_host[CPE_SIZE][F_CACHE_BMAP_SIZE];

#define swapABCD(in0, in1, in2, in3, ot0, ot1, ot2, ot3) {\
    doublev4 o0 = simd_vshff(in1,in0,68 ); \
    doublev4 o1 = simd_vshff(in1,in0,238); \
    doublev4 o2 = simd_vshff(in3,in2,68 ); \
    doublev4 o3 = simd_vshff(in3,in2,238); \
    ot0 = simd_vshff(o2,o0,136); \
    ot1 = simd_vshff(o2,o0,221); \
    ot2 = simd_vshff(o3,o1,136); \
    ot3 = simd_vshff(o3,o1,221); \
}
#define get_f_cache_bmap(bmap, x) ((bmap)[(x) >> 6] & (1lu << (((x) & 63))))
#define set_f_cache_bmap(bmap, x) (bmap)[(x) >> 6] |= (1lu << (((x) & 63)))
#define reset_f_cache(x) { \
    int i; \
    realv4 zero_v4; \
    zero_v4 = 0; \
    real *f_cache_x = (real*) (f_cache + x); \
    for (i = 0; i < (1<<C_T) * (NBNXN_CPU_CLUSTER_I_SIZE) * (DIM); i += 4) { \
        simd_store(zero_v4, f_cache_x + i); \
    } \
}
#define reset_realv4_f(p) { \
    int i; \
    realv4 zero_v4; \
    zero_v4 = 0; \
    real *p_real = (real*) (p); \
    simd_store(zero_v4, p_real + 0); \
    simd_store(zero_v4, p_real + 4); \
    simd_store(zero_v4, p_real + 8); \
}
#define add_realv4_f(a, b) { \
    int i; \
    realv4 a0v4, a4v4, a8v4; \
    realv4 b0v4, b4v4, b8v4; \
    real *a_real = (real*) (a); \
    real *b_real = (real*) (b); \
    simd_load(a0v4, a_real + 0); \
    simd_load(b0v4, b_real + 0); \
    simd_load(a4v4, a_real + 4); \
    simd_load(b4v4, b_real + 4); \
    simd_load(a8v4, a_real + 8); \
    simd_load(b8v4, b_real + 8); \
    a0v4 += b0v4; \
    a4v4 += b4v4; \
    a8v4 += b8v4; \
    simd_store(a0v4, a_real + 0); \
    simd_store(a4v4, a_real + 4); \
    simd_store(a8v4, a_real + 8); \
}

__thread_local volatile int rank;
__thread_local static realv4 one_fv4, skipmask1[4];


__thread_local long slave_count = 0;
#if 0
static void compute1_3_test( atom *a, atom *b, real *fi, real *f, const real *nbfp, real rcut2,
                                  int ci_sh, int ci, int cj, int ntype2,
                                  unsigned int excl)
{
    int i; 
    for (i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
    {
        int ai;
        int type_i_off;
        int j;
    
        ai = ci*NBNXN_CPU_CLUSTER_I_SIZE + i;
    
        type_i_off = a->type[i] * ntype2;
    
        for (j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
        {
            int aj;
            real dx, dy, dz;
            real rsq, rinv;
            real rinvsq, rinvsix;
            real c6, c12;
            real FrLJ6 = 0, FrLJ12 = 0, frLJ = 0, VLJ = 0;
    
            real qq;
            real fcoul;
            real fscal;
            real fx, fy, fz;
            real skipmask;
            int interact;

            lwpf_start(CACU1)
            interact = ((excl>>(i*NBNXN_CPU_CLUSTER_I_SIZE + j)) & 1);
            skipmask = !(cj == ci_sh && j <= i);
            aj = cj*NBNXN_CPU_CLUSTER_I_SIZE + j;
            lwpf_stop(CACU1)
            dx = a->x[i][XX] - b->x[j][XX] + a->shiftvec[XX];
            dy = a->x[i][YY] - b->x[j][YY] + a->shiftvec[YY];
            dz = a->x[i][ZZ] - b->x[j][ZZ] + a->shiftvec[ZZ];
            rsq = dx*dx + dy*dy + dz*dz;
            skipmask = (rsq >= rcut2) ? 0 : skipmask;
            rsq += (1 - interact)*NBNXN_AVOID_SING_R2_INC;
            rinv = 1.0/ sqrt(rsq);
            rinv = rinv * skipmask;
            rinvsq = rinv*rinv;
            {
                c6  = nbfp[type_i_off+b->type[j]*2 ];
                c12 = nbfp[type_i_off+b->type[j]*2+1];
                rinvsix = interact*rinvsq*rinvsq*rinvsq;
                FrLJ6 = c6*rinvsix;
                FrLJ12 = c12*rinvsix*rinvsix;
                frLJ = FrLJ12 - FrLJ6;
            }
                
            fscal = frLJ*rinvsq;
            fx = fscal*dx;
            fy = fscal*dy;
            fz = fscal*dz;
            fi[i*3 +XX] += fx;
            fi[i*3 +YY] += fy;
            fi[i*3 +ZZ] += fz;
    
            f[j*3 +XX] -= fx;
            f[j*3 +YY] -= fy;
            f[j*3 +ZZ] -= fz;
        }
    }
}
#endif

#if 0
static void compute1_2_test( atom *a, atom *b, real *fi, real *f,const real *nbfp,
                                  real rcut2, real k_rf2,
                                  int ci_sh, int ci, int cj, int ntype2,
                                  unsigned int excl)
{
    int i; 
    for (i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
    {
        int ai;
        int type_i_off;
        int j;
    
        ai = ci*NBNXN_CPU_CLUSTER_I_SIZE + i;
    
        type_i_off = a->type[i] * ntype2;
    
        for (j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
        {
            int aj;
            real dx, dy, dz;
            real rsq, rinv;
            real rinvsq, rinvsix;
            real c6, c12;
            real FrLJ6 = 0, FrLJ12 = 0, frLJ = 0, VLJ = 0;
    
            real qq;
            real fcoul;
            real fscal;
            real fx, fy, fz;
            real skipmask;
            int interact;
            interact = ((excl>>(i*NBNXN_CPU_CLUSTER_I_SIZE + j)) & 1);
            skipmask = !(cj == ci_sh && j <= i);
            aj = cj*NBNXN_CPU_CLUSTER_I_SIZE + j;
            dx = a->x[i][XX] - b->x[j][XX] + a->shiftvec[XX];
            dy = a->x[i][YY] - b->x[j][YY] + a->shiftvec[YY];
            dz = a->x[i][ZZ] - b->x[j][ZZ] + a->shiftvec[ZZ];
            rsq = dx*dx + dy*dy + dz*dz;
            skipmask = (rsq >= rcut2) ? 0 : skipmask;
            rsq += (1 - interact)*NBNXN_AVOID_SING_R2_INC;
            rinv = 1.0 / sqrt(rsq);
            rinv = rinv * skipmask;
            rinvsq = rinv*rinv;
            {
                c6  = nbfp[type_i_off+b->type[j]*2 ];
                c12 = nbfp[type_i_off+b->type[j]*2+1];
                rinvsix = interact*rinvsq*rinvsq*rinvsq;
                FrLJ6 = c6*rinvsix;
                FrLJ12 = c12*rinvsix*rinvsix;
                frLJ = FrLJ12 - FrLJ6;
            }
            qq = skipmask * a->q[i] * b->q[j];
            fcoul = qq*(interact*rinv*rinvsq - k_rf2);
            {
                fscal = frLJ*rinvsq + fcoul;
            }
            fx = fscal*dx;
            fy = fscal*dy;
            fz = fscal*dz;
            fi[i*3 +XX] += fx;
            fi[i*3 +YY] += fy;
            fi[i*3 +ZZ] += fz;
    
            f[j*3 +XX] -= fx;
            f[j*3 +YY] -= fy;
            f[j*3 +ZZ] -= fz;
        }
    }
}
#endif




static void printf_v4(realv4 *a)
{
    real c[4];
    simd_storeu(*a, c);
    cpe_printf("1 = %.f 2 =  %.f 3 = %.f 4 = %.f\n", c[0], c[1], c[2], c[3]);

}

static void compute1_3_test_vec( atom *a, atom *b, real *fi, real *f,const real *nbfp,
                                  real rcut2, int ci_sh, int ci, int cj, int ntype2,
                                  unsigned int excl)
{
    int i;
    int j;
    lwpf_start(COMPUTE3)

    real test[16 * 3], c6_array[4], c12_array[4];
    realv4 x1[4],x2[4], shift[4], xq;
    realv4 sk, inter, k_rf2_v4;
    real skipmask_array[4], interact_array[4];
    realv4 fj[4];

    realv4 zero_v4 = 0.0;
    simd_load(xq, a->q);
    shift[XX] = a->shiftvec[XX]; 
    shift[YY] = a->shiftvec[YY];
    shift[ZZ] = a->shiftvec[ZZ];
    simd_load(x2[0], &(a->x[0][0]));
    simd_load(x2[1], &(a->x[1][1]));
    simd_load(x2[2], &(a->x[2][2]));
    x1[0] = simd_vshff(x2[2], x2[1], (2 << 6) + (1 << 4) + (3 << 2) + 2);
    x1[2] = simd_vshff(x2[1], x2[0], (1 << 6) + (0 << 4) + (2 << 2) + 1);
    x1[1] = simd_vshff(x1[0], x1[2], (3 << 6) + (1 << 4) + (2 << 2) + 0);
    x1[0] = simd_vshff(x1[0], x2[0], (2 << 6) + (0 << 4) + (3 << 2) + 0);
    x1[2] = simd_vshff(x2[2], x1[2], (3 << 6) + (0 << 4) + (3 << 2) + 1);



    for(i = 0;i < 16 * 3;i++)
        test[i] = 0;
    
    
    for (j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
    {
        
        realv4 dx, dy, dz;
        realv4 rsq, rinv;
        realv4 rinvsq, rinvsix;
        realv4 c6, c12;
        realv4 FrLJ6 = 0, FrLJ12 = 0, frLJ = 0, VLJ = 0;
    
        realv4 qq;
        realv4 fcoul;
        realv4 fscal;
        realv4 fx, fy, fz;
        realv4 skipmask, interact;
        //int interact;
        //interact = ((excl>>(i*NBNXN_CPU_CLUSTER_I_SIZE + j)) & 1);
        //skipmask = !(cj == ci_sh && j <= i);

        for(i = 0;i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
        {
            int type_i_off;
            type_i_off = a->type[i] * ntype2;
    
            c6_array[i]       = nbfp[type_i_off+b->type[j]*2 ];
            c12_array[i]      = nbfp[type_i_off+b->type[j]*2+1];
        }
        simd_load(c6, c6_array);
        simd_load(c12, c12_array);

        if (cj == ci_sh) {
            skipmask = skipmask1[j];
        }
        else {
            skipmask = one_fv4;
        }
        if (excl == -1) {
            #pragma freq always
            interact = one_fv4;
        }
        else {
            for(i = 0;i < NBNXN_CPU_CLUSTER_I_SIZE; i++) {
                int inter = ((excl>>(i*NBNXN_CPU_CLUSTER_I_SIZE + j)) & 1);
                interact_array[i] = inter;
            }
            simd_load(interact, interact_array);
        }
        
        x2[XX] = b->x[j][XX]; 
        x2[YY] = b->x[j][YY];
        x2[ZZ] = b->x[j][ZZ];
        dx = x1[XX] - x2[XX] + shift[XX]; 
        dy = x1[YY] - x2[YY] + shift[YY];
        dz = x1[ZZ] - x2[ZZ] + shift[ZZ];
        rsq = dx*dx + dy*dy + dz*dz;
        
        realv4 cc = rcut2 - rsq;
        skipmask = simd_vselle(cc, zero_v4, skipmask);
        
        rsq += (1 - interact)*NBNXN_AVOID_SING_R2_INC;
        
        #ifdef GMX_DOUBLE
        rinv = 1.0 / simd_vsqrtd(rsq);
        #else
        rinv = (1.0 /simd_vsqrts(rsq));
        #endif
        
        rinv = rinv * skipmask;
        
        
        rinvsq = rinv*rinv;


        rinvsix = interact*rinvsq*rinvsq*rinvsq;
        FrLJ6 = c6*rinvsix;
        FrLJ12 = c12*rinvsix*rinvsix;
        frLJ = FrLJ12 - FrLJ6;
       
        
        fscal = frLJ*rinvsq;

        fx = fscal*dx;
        fy = fscal*dy;
        fz = fscal*dz;
        
        {
            realv4 zero_fv4 = 0;
            realv4 f0, f1, f2, f3;
            swapABCD(fx, fy, fz, zero_fv4, f0, f1, f2, f3);

            realv4 f0t, f1t, f2t;
            realv4 fi0, fi1, fi2;
            simd_load(fi0, fi + 0 + XX);
            simd_load(fi1, fi + 4 + XX);
            simd_load(fi2, fi + 8 + XX);
            f0t = simd_vshff(f1, f0, 0x4e);
            f0t = simd_vshff(f0t, f0, 0x84);
            f1t = simd_vshff(f2, f1, 0x49);
            f2t = simd_vshff(f3, f2, 0x4e);
            f2t = simd_vshff(f3, f2t, 0x98);

            fi0 += f0t;
            fi1 += f1t;
            fi2 += f2t;
            f0 += f1;
            f2 += f3;
            fj[j] = f0;
            fj[j] += f2;

            simd_store(fi0, fi + 0 + XX);
            simd_store(fi1, fi + 4 + XX);
            simd_store(fi2, fi + 8 + XX);
        }

    }

    {
        realv4 fj0, fj1, fj2;
        realv4 f0t, f1t, f2t;

        simd_load(fj0, f + 0 + XX);
        simd_load(fj1, f + 4 + XX);
        simd_load(fj2, f + 8 + XX);
        f0t = simd_vshff(fj[1], fj[0], 0x4e);
        f0t = simd_vshff(f0t, fj[0], 0x84);
        f1t = simd_vshff(fj[2], fj[1], 0x49);
        f2t = simd_vshff(fj[3], fj[2], 0x4e);
        f2t = simd_vshff(fj[3], f2t, 0x98);
        fj0 -= f0t;
        fj1 -= f1t;
        fj2 -= f2t;
        simd_store(fj0, f + 0 + XX);
        simd_store(fj1, f + 4 + XX);
        simd_store(fj2, f + 8 + XX);
    }
    lwpf_stop(COMPUTE3)
}

static void compute1_2_test_vec( atom *a, atom *b, real *fi, real *f,const real *nbfp,
                                  real rcut2, real k_rf2,
                                  int ci_sh, int ci, int cj, int ntype2,
                                  unsigned int excl)
{
    lwpf_start(COMPUTE2)
    int i;
    int j;

    real test[16 * 3], c6_array[4], c12_array[4];
    realv4 x1[4],x2[4], shift[4], xq;
    realv4 sk, inter, k_rf2_v4;
    real skipmask_array[4], interact_array[4];
    realv4 fj[4];
    
    k_rf2_v4 = k_rf2;

    realv4 zero_v4 = 0.0;
    simd_load(xq, a->q);
    shift[XX] = a->shiftvec[XX]; 
    shift[YY] = a->shiftvec[YY];
    shift[ZZ] = a->shiftvec[ZZ];
    simd_load(x2[0], &(a->x[0][0]));
    simd_load(x2[1], &(a->x[1][1]));
    simd_load(x2[2], &(a->x[2][2]));
    x1[0] = simd_vshff(x2[2], x2[1], (2 << 6) + (1 << 4) + (3 << 2) + 2);
    x1[2] = simd_vshff(x2[1], x2[0], (1 << 6) + (0 << 4) + (2 << 2) + 1);
    x1[1] = simd_vshff(x1[0], x1[2], (3 << 6) + (1 << 4) + (2 << 2) + 0);
    x1[0] = simd_vshff(x1[0], x2[0], (2 << 6) + (0 << 4) + (3 << 2) + 0);
    x1[2] = simd_vshff(x2[2], x1[2], (3 << 6) + (0 << 4) + (3 << 2) + 1);



    for(i = 0;i < 16 * 3;i++)
        test[i] = 0;
    
    
    for (j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
    {
        
        realv4 dx, dy, dz;
        realv4 rsq, rinv;
        realv4 rinvsq, rinvsix;
        realv4 c6, c12;
        realv4 FrLJ6 = 0, FrLJ12 = 0, frLJ = 0, VLJ = 0;
    
        realv4 qq;
        realv4 fcoul;
        realv4 fscal;
        realv4 fx, fy, fz;
        realv4 skipmask, interact;
        //int interact;
        //interact = ((excl>>(i*NBNXN_CPU_CLUSTER_I_SIZE + j)) & 1);
        //skipmask = !(cj == ci_sh && j <= i);

        for(i = 0;i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
        {
            int type_i_off;
            type_i_off = a->type[i] * ntype2;
    
            c6_array[i]       = nbfp[type_i_off+b->type[j]*2 ];
            c12_array[i]      = nbfp[type_i_off+b->type[j]*2+1];
        }
        simd_load(c6, c6_array);
        simd_load(c12, c12_array);

        if (cj == ci_sh) {
            skipmask = skipmask1[j];
        }
        else {
            skipmask = one_fv4;
        }
        if (excl == -1) {
            #pragma freq always
            interact = one_fv4;
        }
        else {
            for(i = 0;i < NBNXN_CPU_CLUSTER_I_SIZE; i++) {
                int inter = ((excl>>(i*NBNXN_CPU_CLUSTER_I_SIZE + j)) & 1);
                interact_array[i] = inter;
            }
            simd_load(interact, interact_array);
        }
        
        x2[XX] = b->x[j][XX]; 
        x2[YY] = b->x[j][YY];
        x2[ZZ] = b->x[j][ZZ];
        dx = x1[XX] - x2[XX] + shift[XX]; 
        dy = x1[YY] - x2[YY] + shift[YY];
        dz = x1[ZZ] - x2[ZZ] + shift[ZZ];
        rsq = dx*dx + dy*dy + dz*dz;
        
        realv4 cc = rcut2 - rsq;
        skipmask = simd_vselle(cc, zero_v4, skipmask);
        
        rsq += (1 - interact)*NBNXN_AVOID_SING_R2_INC;
        
        #ifdef GMX_DOUBLE
        rinv = 1.0 / simd_vsqrtd(rsq);
        #else
        rinv = (1.0 /simd_vsqrts(rsq));
        #endif
        
        rinv = rinv * skipmask;
        
        
        rinvsq = rinv*rinv;


        rinvsix = interact*rinvsq*rinvsq*rinvsq;
        FrLJ6 = c6*rinvsix;
        FrLJ12 = c12*rinvsix*rinvsix;
        frLJ = FrLJ12 - FrLJ6;
       
        realv4 yq = b->q[j];
        qq = skipmask * xq * yq;
        fcoul = qq*(interact*rinv*rinvsq - k_rf2_v4);
        
        
        fscal = frLJ*rinvsq + fcoul;

        fx = fscal*dx;
        fy = fscal*dy;
        fz = fscal*dz;
        
        {
            realv4 zero_fv4 = 0;
            realv4 f0, f1, f2, f3;
            swapABCD(fx, fy, fz, zero_fv4, f0, f1, f2, f3);

            realv4 f0t, f1t, f2t;
            realv4 fi0, fi1, fi2;
            simd_load(fi0, fi + 0 + XX);
            simd_load(fi1, fi + 4 + XX);
            simd_load(fi2, fi + 8 + XX);
            f0t = simd_vshff(f1, f0, 0x4e);
            f0t = simd_vshff(f0t, f0, 0x84);
            f1t = simd_vshff(f2, f1, 0x49);
            f2t = simd_vshff(f3, f2, 0x4e);
            f2t = simd_vshff(f3, f2t, 0x98);

            fi0 += f0t;
            fi1 += f1t;
            fi2 += f2t;
            f0 += f1;
            f2 += f3;
            fj[j] = f0;
            fj[j] += f2;

            simd_store(fi0, fi + 0 + XX);
            simd_store(fi1, fi + 4 + XX);
            simd_store(fi2, fi + 8 + XX);
        }

    }

    {
        realv4 fj0, fj1, fj2;
        realv4 f0t, f1t, f2t;

        simd_load(fj0, f + 0 + XX);
        simd_load(fj1, f + 4 + XX);
        simd_load(fj2, f + 8 + XX);
        f0t = simd_vshff(fj[1], fj[0], 0x4e);
        f0t = simd_vshff(f0t, fj[0], 0x84);
        f1t = simd_vshff(fj[2], fj[1], 0x49);
        f2t = simd_vshff(fj[3], fj[2], 0x4e);
        f2t = simd_vshff(fj[3], f2t, 0x98);
        fj0 -= f0t;
        fj1 -= f1t;
        fj2 -= f2t;
        simd_store(fj0, f + 0 + XX);
        simd_store(fj1, f + 4 + XX);
        simd_store(fj2, f + 8 + XX);
    }
    lwpf_stop(COMPUTE2)
}


static void compute1_1_test_vec( atom *a, atom *b, real *fi, real *f,const real *nbfp,
                                  real rcut2, real k_rf2,
                                  int ci_sh, int ci, int cj, int ntype2,
                                  unsigned int excl)
{
    lwpf_start(COMPUTE1)
    int i;
    int j;

    real test[16 * 3], c6_array[4], c12_array[4];
    realv4 x1[4],x2[4], shift[4], xq;
    realv4 sk, inter, k_rf2_v4;
    real skipmask_array[4], interact_array[4];
    realv4 fj[4];
    
    
    k_rf2_v4 = k_rf2;

    realv4 zero_v4 = 0.0;
    simd_load(xq, a->q);
    shift[XX] = a->shiftvec[XX]; 
    shift[YY] = a->shiftvec[YY];
    shift[ZZ] = a->shiftvec[ZZ];
    simd_load(x2[0], &(a->x[0][0]));
    simd_load(x2[1], &(a->x[1][1]));
    simd_load(x2[2], &(a->x[2][2]));
    x1[0] = simd_vshff(x2[2], x2[1], (2 << 6) + (1 << 4) + (3 << 2) + 2);
    x1[2] = simd_vshff(x2[1], x2[0], (1 << 6) + (0 << 4) + (2 << 2) + 1);
    x1[1] = simd_vshff(x1[0], x1[2], (3 << 6) + (1 << 4) + (2 << 2) + 0);
    x1[0] = simd_vshff(x1[0], x2[0], (2 << 6) + (0 << 4) + (3 << 2) + 0);
    x1[2] = simd_vshff(x2[2], x1[2], (3 << 6) + (0 << 4) + (3 << 2) + 1);

    for(i = 0;i < 16 * 3;i++)
        test[i] = 0;
    
    
    for (j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
    {
        
        realv4 dx, dy, dz;
        realv4 rsq, rinv;
        realv4 rinvsq, rinvsix;
        realv4 c6, c12;
        realv4 FrLJ6 = 0, FrLJ12 = 0, frLJ = 0, VLJ = 0;
    
        realv4 qq;
        realv4 fcoul;
        realv4 fscal;
        realv4 fx, fy, fz;
        realv4 skipmask, interact;


        for(i = 0;i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
        {
            int type_i_off;
            type_i_off = a->type[i] * ntype2;
    
            c6_array[i]       = nbfp[type_i_off+b->type[j]*2 ];
            c12_array[i]      = nbfp[type_i_off+b->type[j]*2+1];
        }
        simd_load(c6, c6_array);
        simd_load(c12, c12_array);

        if (cj == ci_sh) {
            skipmask = skipmask1[j];
        }
        else {
            skipmask = one_fv4;
        }
        if (excl == -1) {
            #pragma freq always
            interact = one_fv4;
        }
        else {
            for(i = 0;i < NBNXN_CPU_CLUSTER_I_SIZE; i++) {
                int inter = ((excl>>(i*NBNXN_CPU_CLUSTER_I_SIZE + j)) & 1);
                interact_array[i] = inter;
            }
            simd_load(interact, interact_array);
        }
        
        x2[XX] = b->x[j][XX]; 
        x2[YY] = b->x[j][YY];
        x2[ZZ] = b->x[j][ZZ];
        dx = x1[XX] - x2[XX] + shift[XX]; 
        dy = x1[YY] - x2[YY] + shift[YY];
        dz = x1[ZZ] - x2[ZZ] + shift[ZZ];
        rsq = dx*dx + dy*dy + dz*dz;
        
        realv4 cc = rcut2 - rsq;
        skipmask = simd_vselle(cc, zero_v4, skipmask);
        
        rsq += (1 - interact)*NBNXN_AVOID_SING_R2_INC;
        
        #ifdef GMX_DOUBLE
        rinv = 1.0 / simd_vsqrtd(rsq);
        #else
        rinv = (1.0 /simd_vsqrts(rsq));
        #endif
        
        rinv = rinv * skipmask;
        
        
        rinvsq = rinv*rinv;


        rinvsix = interact*rinvsq*rinvsq*rinvsq;
        FrLJ6 = c6*rinvsix;
        FrLJ12 = c12*rinvsix*rinvsix;
        frLJ = FrLJ12 - FrLJ6;
       
        realv4 yq = b->q[j];
        qq = skipmask * xq * yq;
        fcoul = qq*(interact*rinv*rinvsq - k_rf2_v4);
        
        
        fscal = frLJ*rinvsq + fcoul;
        fscal = simd_vshff( fcoul, fscal, (3 << 6) + (2 << 4) + (1 << 2) + 0);

        fx = fscal*dx;
        fy = fscal*dy;
        fz = fscal*dz;

        {
            realv4 zero_fv4 = 0;
            realv4 f0, f1, f2, f3;
            swapABCD(fx, fy, fz, zero_fv4, f0, f1, f2, f3);

            realv4 f0t, f1t, f2t;
            realv4 fi0, fi1, fi2;
            simd_load(fi0, fi + 0 + XX);
            simd_load(fi1, fi + 4 + XX);
            simd_load(fi2, fi + 8 + XX);
            f0t = simd_vshff(f1, f0, 0x4e);
            f0t = simd_vshff(f0t, f0, 0x84);
            f1t = simd_vshff(f2, f1, 0x49);
            f2t = simd_vshff(f3, f2, 0x4e);
            f2t = simd_vshff(f3, f2t, 0x98);

            fi0 += f0t;
            fi1 += f1t;
            fi2 += f2t;
            f0 += f1;
            f2 += f3;
            fj[j] = f0;
            fj[j] += f2;

            simd_store(fi0, fi + 0 + XX);
            simd_store(fi1, fi + 4 + XX);
            simd_store(fi2, fi + 8 + XX);
        }

    }

    {
        realv4 fj0, fj1, fj2;
        realv4 f0t, f1t, f2t;

        simd_load(fj0, f + 0 + XX);
        simd_load(fj1, f + 4 + XX);
        simd_load(fj2, f + 8 + XX);
        f0t = simd_vshff(fj[1], fj[0], 0x4e);
        f0t = simd_vshff(f0t, fj[0], 0x84);
        f1t = simd_vshff(fj[2], fj[1], 0x49);
        f2t = simd_vshff(fj[3], fj[2], 0x4e);
        f2t = simd_vshff(fj[3], f2t, 0x98);
        fj0 -= f0t;
        fj1 -= f1t;
        fj2 -= f2t;
        simd_store(fj0, f + 0 + XX);
        simd_store(fj1, f + 4 + XX);
        simd_store(fj2, f + 8 + XX);
    }
    lwpf_stop(COMPUTE1)
}


#if 0
static void compute1_1_test( atom *a, atom *b, real *fi, real *f,const real *nbfp,
                                  real rcut2, real k_rf2,
                                  int ci_sh, int ci, int cj, int ntype2,
                                  unsigned int excl)
{
    int i; 
    for (i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
    {
        int ai;
        int type_i_off;
        int j;
    
        ai = ci*NBNXN_CPU_CLUSTER_I_SIZE + i;
    
        type_i_off = a->type[i] * ntype2;
    
        for (j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
        {
            int aj;
            real dx, dy, dz;
            real rsq, rinv;
            real rinvsq, rinvsix;
            real c6, c12;
            real FrLJ6 = 0, FrLJ12 = 0, frLJ = 0, VLJ = 0;
    
            real qq;
            real fcoul;
            real fscal;
            real fx, fy, fz;
            real skipmask;
            int interact;
            interact = ((excl>>(i*NBNXN_CPU_CLUSTER_I_SIZE + j)) & 1);
            skipmask = !(cj == ci_sh && j <= i);
            aj = cj*NBNXN_CPU_CLUSTER_I_SIZE + j;
            dx = a->x[i][XX] - b->x[j][XX] + a->shiftvec[XX];
            dy = a->x[i][YY] - b->x[j][YY] + a->shiftvec[YY];
            dz = a->x[i][ZZ] - b->x[j][ZZ] + a->shiftvec[ZZ];
            rsq = dx*dx + dy*dy + dz*dz;
            skipmask = (rsq >= rcut2) ? 0 : skipmask;
            rsq += (1 - interact)*NBNXN_AVOID_SING_R2_INC;
            rinv = 1.0 / sqrt(rsq);
            rinv = rinv * skipmask;
            rinvsq = rinv*rinv;
            if (i < NBNXN_CPU_CLUSTER_I_SIZE/2)
            {
                c6  = nbfp[type_i_off+b->type[j]*2 ];
                c12 = nbfp[type_i_off+b->type[j]*2+1];
                rinvsix = interact*rinvsq*rinvsq*rinvsq;
                FrLJ6 = c6*rinvsix;
                FrLJ12 = c12*rinvsix*rinvsix;
                frLJ = FrLJ12 - FrLJ6;
            }
            qq = skipmask * a->q[i] * b->q[j];
            fcoul = qq*(interact*rinv*rinvsq - k_rf2);
            if (i < NBNXN_CPU_CLUSTER_I_SIZE/2)
            {
                fscal = frLJ*rinvsq + fcoul;
            }
            else
            {
                fscal = fcoul;
            }
            fx = fscal*dx;
            fy = fscal*dy;
            fz = fscal*dz;
            fi[i*3 +XX] += fx;
            fi[i*3 +YY] += fy;
            fi[i*3 +ZZ] += fz;
    
            f[j*3 +XX] -= fx;
            f[j*3 +YY] -= fy;
            f[j*3 +ZZ] -= fz;
        }
    }
}
#endif



#define ELEC_TYPE_SIZE 12
#define E_T_S ELEC_TYPE_SIZE
//static inline long rpcc(){
//  long rpc;
//  asm volatile("rcsr %0, 4" : "=r"(rpc));
//  return rpc;
//}

void nbnxn_kernel_ElecRF_VdwLJ_F_ref_slave(ElecRF_VdwLJ_F_ref *params)
{
  lwpf_enter(CASE2);

    lwpf_start(ALL)
    lwpf_start(PRE)
    
    long transport_bytes = 0;
    
    int i,j,k,l;
    ElecRF_VdwLJ_F_ref lp;
    volatile unsigned long get_reply, put_reply;
    int ntype, natoms, nci, ntype2;
    gmx_bool do_LJ, half_LJ, do_coul, do_self;
    atom *atoms;
    nbnxn_ci_t *nbln;
    real facel, rcut2, k_rf2, tabq_scale;
    atom a, b;
    int C_P_ = (1 << C_P) - 1;
    int C_T_ = (1 << C_T) - 1;
    
    athread_syn(ARRAY_SCOPE, 0xffff);
    transport_bytes += sizeof(ElecRF_VdwLJ_F_ref);
    if (_MYID == 0){
        get_reply = 0;
        athread_get(BCAST_MODE, params, &lp, sizeof(ElecRF_VdwLJ_F_ref), &get_reply,0xff,0,0);
        while(get_reply != 1);
    }
    athread_syn(ARRAY_SCOPE, 0xffff);
    real nbfp[E_T_S * E_T_S * 2];
    int get_update_f[64];
    //athread_syn(ARRAY_SCOPE, 0xffff);
    transport_bytes += sizeof(real) * 2 * lp.ntype * lp.ntype + sizeof(int) * 64;
    if (_MYID == 0){
        get_reply = 0;
        athread_get(BCAST_MODE, lp.nbfp, nbfp, sizeof(real) * 2 * lp.ntype * lp.ntype, &get_reply,0xff,0,0);
        athread_get(BCAST_MODE, lp.get_update_f, get_update_f, sizeof(int) * 64, &get_reply, 0xff, 0, 0);
        while(get_reply != 2);
    }
    athread_syn(ARRAY_SCOPE, 0xffff);
    //pe_syn();
    

    natoms        = lp.natoms; 
    nci           = lp.nci   ; 
    ntype         = lp.ntype ; 
    rcut2         = lp.rcut2 ; 
    facel         = lp.facel ;
    atoms         = lp.atoms ; 
    ntype2        = 2 * ntype;
    k_rf2         = lp.k_rf2 ;
    rank          = lp.rank;

    one_fv4       = 1.f;
    skipmask1[0]  = 0.f;
    skipmask1[1]  = simd_set_floatv4(1.f, 0.f, 0.f, 0.f);
    skipmask1[2]  = simd_set_floatv4(1.f, 1.f, 0.f, 0.f);
    skipmask1[3]  = simd_set_floatv4(1.f, 1.f, 1.f, 0.f);

    
    //int get_update_f[64];
    //volatile int *get_update_f = lp.get_update_f;
    //athread_syn(
    
    nbnxn_ci_t l_nbln[I_CLUSTER_SIZE];
    nbnxn_cj_t l_l_cj[J_CLUSTER_SIZE];
    real f_shift[I_CLUSTER_SIZE][3];

    for(i = 0; i < I_C_S; i++)
        for(j = 0; j < 3; j++)
            f_shift[i][j] = 0;
     
    
    int cach_f, cach_t, cach_p;
    int pn, cjind, cjind1, cjind0;
    int ish, ishf, ci_sh;
    int ci, n, np, d;
    int ip, jp, f_cache_update_ps;
    real shiftvec[(NBNXN_CI_SHIFT + 1) *3 * 3 + 4];

    get_reply = 0;
    athread_get(PE_MODE, lp.shiftvec, shiftvec, sizeof(real) * (NBNXN_CI_SHIFT + 1) * 3 * 3, &get_reply,0,0,0);
    transport_bytes += sizeof(real) * (NBNXN_CI_SHIFT + 1) * 3 * 3;
    while(get_reply != 1);
    pe_syn();   
      
    atom atoms_cache[1 << C_P][1 << C_T];
    update_f f_cache[1 << C_P][1 << C_T];
    int f_cache_ps[1 << C_P], cache_ps[1 << C_P];

    int cach_size = 1 << C_P;
    for(i = 0;i < cach_size;i++)
      cache_ps[i] = f_cache_ps[i] = -1;
    
    real fi[NBNXN_CPU_CLUSTER_I_SIZE * DIM], f[NBNXN_CPU_CLUSTER_I_SIZE * DIM];
    int get = 0, miss = 0;
    int np_st = get_update_f[ _MYID];
    int np_ed;
    if(_MYID == 63)
        np_ed = nci;
    else
        np_ed = get_update_f[ _MYID + 1];

    int cluster_number = ((((natoms + 3) >> 2) + 8)&(~7));

    update_f *update_data = &lp.update_date[_MYID * cluster_number];

    long *f_cache_bmap_h = (long*) f_cache_bmap_host;
    int ncache_line = (cluster_number + (1 << C_T) - 1) / (1 << C_T);
    int ncache_bitmap = (ncache_line + 63) >> 6;
    long f_cache_bmap[F_CACHE_BMAP_SIZE];
    // long cache_bitmap[((ncache_bitmap + 3) >> 2) << 2];
    {
        doublev4 zero_dv4;
        zero_dv4 = 0;
        #pragma unroll[4]
        for (i = 0; i < ncache_bitmap; i += 4) {
            simd_store(zero_dv4, f_cache_bmap + i + 0);
        }
    }

    lwpf_stop(PRE)
    
    
    for(np = np_st; np < np_ed; np += I_C_S)
    {
        int i_size;
        if(np_ed - np > I_C_S)
            i_size = I_C_S;
        else
            i_size = np_ed - np;
        get_reply = 0;
        athread_get(PE_MODE, lp.nbln + np, l_nbln, sizeof(nbnxn_ci_t) * i_size, &get_reply,0,0,0);
        transport_bytes += sizeof(nbnxn_ci_t) * i_size;
        while(get_reply != 1);
        pe_syn();

        for(n = 0; n < i_size;n++)
        {
            nbln    = &l_nbln[n];
            ish     = (nbln->shift & NBNXN_CI_SHIFT);
            ishf    = ish*DIM;
            cjind0  = nbln->cj_ind_start;
            cjind1  = nbln->cj_ind_end;
            ci      = nbln->ci;
            ci_sh   = (ish == CENTRAL ? ci : -1);
            do_LJ   = (nbln->shift & NBNXN_CI_DO_LJ(0));
            do_coul = (nbln->shift & NBNXN_CI_DO_COUL(0));
            half_LJ = ((nbln->shift & NBNXN_CI_HALF_LJ(0)) || !do_LJ) && do_coul;
            do_self = do_coul;
            
            get++;
            reset_realv4_f(fi);

            int ci_t = ci & ((1 << C_T) - 1);
            int ci_p = (ci >> C_T) & ((1 << C_P) - 1);
            if(cache_ps[ci_p] != ci - ci_t)
            {
                lwpf_start(GET_CACHE)
                miss++;
                get_reply = 0;
                athread_get(PE_MODE, atoms + ci - ci_t, atoms_cache[ci_p], (1 << C_T) * sizeof(atom), &get_reply,0,0,0);
                transport_bytes += (1 << C_T) * sizeof(atom);
                while(get_reply != 1); 
                pe_syn();
                cache_ps[ci_p] = ci - ci_t;
                lwpf_stop(GET_CACHE)
            }

            atom a = atoms_cache[ci_p][ci_t];
            for (d = 0; d < DIM; d++)
            {
                a.shiftvec[d] = shiftvec[ishf+d]; 
            }


            for(cjind = cjind0; cjind < cjind1; cjind+= J_C_S)
            {
                int j_size = J_C_S;
                if(J_C_S > cjind1 - cjind)
                    j_size = cjind1 - cjind;
                
                
                get_reply = 0;
                athread_get(PE_MODE, lp.l_cj + cjind, l_l_cj, sizeof(nbnxn_cj_t) * j_size, &get_reply,0,0,0);
                transport_bytes += sizeof(nbnxn_cj_t) * j_size;
                while(get_reply != 1);
                pe_syn();

                int cj_ps;
                for(cj_ps = 0; cj_ps < j_size; cj_ps++)
                {
                    
                    int cj = l_l_cj[cj_ps].cj;
                    int cj_t = cj & ((1 << C_T) - 1);
                    int cj_p = (cj >> C_T) & ((1 << C_P) - 1);
                    
                    get++;
                    if(cache_ps[cj_p] != cj - cj_t)
                    {
                        lwpf_start(GET_CACHE)
                        miss++;
                        get_reply = 0;
                        athread_get(PE_MODE, atoms + cj - cj_t, atoms_cache[cj_p], 
                                             (1 << C_T) * sizeof(atom), &get_reply,0,0,0);
                        transport_bytes += (1 << C_T) * sizeof(atom);
                        cache_ps[cj_p] = cj - cj_t;          
                        while(get_reply != 1); 
                        pe_syn();
                        lwpf_stop(GET_CACHE)
                    }
                    
                    atom b = atoms_cache[cj_p][cj_t];

                    reset_realv4_f(f);

                    static int flag = 0;
                    lwpf_start(CACU)
                    if(half_LJ)
                    {
                        compute1_1_test_vec( &a, &b, fi, f, nbfp,
                                          rcut2, k_rf2,
                                          ci_sh, ci, cj, ntype2,
                                          l_l_cj[cj_ps].excl);
            
                    } else if(do_coul)
                    {

                        compute1_2_test_vec( &a, &b, fi, f, nbfp,
                                          rcut2, k_rf2,
                                          ci_sh, ci, cj, ntype2,
                                          l_l_cj[cj_ps].excl);
            
            
                    } else
                    {
                        compute1_3_test_vec( &a, &b, fi, f, nbfp,rcut2, 
                                          ci_sh, ci, cj, ntype2,
                                          l_l_cj[cj_ps].excl);
            
                    }
                    lwpf_stop(CACU)


                    if(f_cache_ps[cj_p] != cj - cj_t)
                    {
                        if(f_cache_ps[cj_p] != -1)
                        {
                            lwpf_start(PUT_CACHE)
                            get_reply = 0;
                            athread_put(PE_MODE, f_cache[cj_p], &update_data[f_cache_ps[cj_p]], 
                                            (1 << C_T) * sizeof(update_f), &get_reply, 0, 0);
                            transport_bytes += (1 << C_T) * sizeof(update_f);
                            while(get_reply != 1);
                            pe_syn();
                            set_f_cache_bmap(f_cache_bmap, f_cache_ps[cj_p] >> C_T);
                            lwpf_stop(PUT_CACHE)

                            lwpf_start(GET_CACHE)
                            if (get_f_cache_bmap(f_cache_bmap, (cj - cj_t) >> C_T)) {
                                get_reply = 0;
                                athread_get(PE_MODE, &update_data[cj - cj_t], f_cache[cj_p], 
                                        (1 << C_T) * sizeof(update_f), &get_reply, 0, 0, 0);
                                transport_bytes += (1 << C_T) * sizeof(update_f);
                                while(get_reply != 1);
                                pe_syn();
                            }
                            else {
                                reset_f_cache(cj_p);
                            }
                            lwpf_stop(GET_CACHE)

                        } else
                        {
                            reset_f_cache(cj_p);
                        }
                        f_cache_ps[cj_p] = cj - cj_t;
                    }

                    lwpf_start(CACU2)
                    add_realv4_f(&(f_cache[cj_p][cj_t].f), f);
                    // for(i = 0;i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
                    // {
                    //     f_cache[cj_p][cj_t].f[i][XX] += f[i * 3 + XX];
                    //     f_cache[cj_p][cj_t].f[i][YY] += f[i * 3 + YY];
                    //     f_cache[cj_p][cj_t].f[i][ZZ] += f[i * 3 + ZZ];
                    // }
                    lwpf_stop(CACU2)
                }
            
            }
    
            if(f_cache_ps[ci_p] != ci - ci_t)
            {
               
                if(f_cache_ps[ci_p] != -1)
                {
                    lwpf_start(PUT_CACHE)
                    get_reply = 0;
                    athread_put(PE_MODE, f_cache[ci_p], &update_data[f_cache_ps[ci_p]], 
                                            (1 << C_T) * sizeof(update_f), &get_reply, 0, 0);
                    transport_bytes += (1 << C_T) * sizeof(update_f);
                    while(get_reply != 1);
                    pe_syn();
                    set_f_cache_bmap(f_cache_bmap, f_cache_ps[ci_p] >> C_T);
                    lwpf_stop(PUT_CACHE)

                    lwpf_start(GET_CACHE)
                    if (get_f_cache_bmap(f_cache_bmap, (ci - ci_t) >> C_T)) {
                        get_reply = 0;
                        athread_get(PE_MODE, &update_data[ci - ci_t], f_cache[ci_p], 
                                (1 << C_T) * sizeof(update_f), &get_reply, 0, 0, 0);
                        transport_bytes += (1 << C_T) * sizeof(update_f);
                        while(get_reply != 1);
                        pe_syn();
                    }
                    else {
                        reset_f_cache(ci_p);
                    }
                    lwpf_stop(GET_CACHE)

                } else
                {
                    reset_f_cache(ci_p);
                }
                f_cache_ps[ci_p] = ci - ci_t;
            }

            lwpf_start(CACU2)
            add_realv4_f(&(f_cache[ci_p][ci_t].f), fi);
            // for(i = 0;i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
            // {
            //     f_cache[ci_p][ci_t].f[i][XX] += fi[i * 3 + XX];     
            //     f_cache[ci_p][ci_t].f[i][YY] += fi[i * 3 + YY];     
            //     f_cache[ci_p][ci_t].f[i][ZZ] += fi[i * 3 + ZZ];     
            // }

            // add_realv4_f(f_shift[n], fi);
            for(i = 0;i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
            {
                f_shift[n][XX] += fi[i * 3 + XX];
                f_shift[n][YY] += fi[i * 3 + YY];
                f_shift[n][ZZ] += fi[i * 3 + ZZ];
            }
            lwpf_stop(CACU2)
        }

        lwpf_start(PUT_CACHE)
        get_reply = 0;
        athread_put(PE_MODE, f_shift, &lp.fshift_CPE[np * 3], 3 * i_size * sizeof(real), 
                             &get_reply,0,0);
        transport_bytes += 3 * i_size * sizeof(real);
        while(get_reply != 1);
        pe_syn();
        lwpf_stop(PUT_CACHE)

        reset_realv4_f(f_shift);
    }

    for(i = 0;i < (1 << C_P);i++)
    {
        if(f_cache_ps[i] != -1)
        {
            lwpf_start(PUT_CACHE)
            get_reply = 0;
            athread_put(PE_MODE, f_cache[i], &update_data[f_cache_ps[i]], 
                                    (1 << C_T) * sizeof(update_f), &get_reply, 0, 0);
            transport_bytes += (1 << C_T) * sizeof(update_f);
            set_f_cache_bmap(f_cache_bmap, f_cache_ps[i] >> C_T);
            while(get_reply != 1);
            pe_syn();
            lwpf_stop(PUT_CACHE)
        }
    }
    get_reply = 0;
    athread_put(PE_MODE, f_cache_bmap, f_cache_bmap_h + _MYID * F_CACHE_BMAP_SIZE, 
            F_CACHE_BMAP_SIZE * sizeof(long), &get_reply, 0, 0);
    transport_bytes += F_CACHE_BMAP_SIZE * sizeof(long);
    while(get_reply != 1);
    pe_syn();
    lwpf_stop(ALL)

  lwpf_exit(CASE2);
  

  //slave_count++;
  //if (slave_count == 100)
  //{
  //  if (lp.rank == 0)
  //    cal_locked_printf("_MYID: %d, transport_bytes: %ld\n", _MYID, transport_bytes);
  //  slave_count = 0;
  //}

}

#define update_f_size 32
#define U_F_S update_f_size
void update_data_f_CT3_slave(update_data_ *params)
{
    update_data_ lp;
    volatile unsigned long get_reply, put_reply;
    get_reply = 0;
    athread_get(PE_MODE, params, &lp, sizeof(update_data_), &get_reply,0,0,0);
    while(get_reply != 1);
    pe_syn();

    long f_cache_bmap_all[CPE_SIZE][F_CACHE_BMAP_SIZE];
    long *f_cache_bmap_h = (long*) f_cache_bmap_host;
    athread_syn(ARRAY_SCOPE, 0xffff);
    if (_MYID == 0) {
        get_reply = 0;
        athread_get(BCAST_MODE, f_cache_bmap_h, f_cache_bmap_all, CPE_SIZE * F_CACHE_BMAP_SIZE * sizeof(long), &get_reply, 0xff, 0, 0);
        while(get_reply != 1);
    }
    athread_syn(ARRAY_SCOPE, 0xffff);
    
    int natoms = lp.natoms;
    int cluster_number = ((((natoms + 3) >> 2) + 8)&(~7));
    int i;
    update_f data_f[CPE_SIZE][update_f_size >> 2];
    real f[3 * update_f_size];
    update_f *host_update_f = &lp.update_data[0];

    long mask_zero[64 / update_f_size];
    for (i = 0; i < 64 / update_f_size; ++i) {
        mask_zero[i] = ((1ul << update_f_size) - 1) << (i * update_f_size);
    }

    for(i = _MYID * update_f_size; i < lp.natoms; i += CPE_SIZE * update_f_size)
    {
        int i_size;
        if(natoms - i > update_f_size)
            i_size = U_F_S;
        else
            i_size = natoms - i;
        
        get_reply = 0;
        athread_get(PE_MODE, &lp.f[i * 3], &f[0], sizeof(real) * i_size * 3, &get_reply,0,0,0);
        while(get_reply != 1);
        pe_syn();

        int k;
        for(k = 0; k < CPE_SIZE; k++)
        if (f_cache_bmap_all[k][((i / NBNXN_CPU_CLUSTER_I_SIZE) >> C_T) >> 6] & mask_zero[(((i / NBNXN_CPU_CLUSTER_I_SIZE) >> C_T) % 64) / update_f_size])
        {
            get_reply = 0;
            athread_get(PE_MODE, &lp.update_data[(k * cluster_number) + (i >> 2)], &data_f[k][0], sizeof(update_f) * ((i_size + 3) >> 2), &get_reply,0,0,0);
            while(get_reply != 1);
            pe_syn();
        }
        int j = 0;
        for(k = 0; k < CPE_SIZE; k++)
        {
            for(j = 0; j < i_size; j += 4)
            {
                if (get_f_cache_bmap(f_cache_bmap_all[k], ((i + j) / NBNXN_CPU_CLUSTER_I_SIZE) >> C_T))
                {
                int id = j >> 2;
                int jd = j * 3;
                f[jd + 0 + XX] += data_f[k][id].f[0][XX];
                f[jd + 0 + YY] += data_f[k][id].f[0][YY];
                f[jd + 0 + ZZ] += data_f[k][id].f[0][ZZ];

                f[jd + 3 + XX] += data_f[k][id].f[1][XX];
                f[jd + 3 + YY] += data_f[k][id].f[1][YY];
                f[jd + 3 + ZZ] += data_f[k][id].f[1][ZZ];
                
                f[jd + 6 + XX] += data_f[k][id].f[2][XX];
                f[jd + 6 + YY] += data_f[k][id].f[2][YY];
                f[jd + 6 + ZZ] += data_f[k][id].f[2][ZZ];
                
                f[jd + 9 + XX] += data_f[k][id].f[3][XX];
                f[jd + 9 + YY] += data_f[k][id].f[3][YY];
                f[jd + 9 + ZZ] += data_f[k][id].f[3][ZZ];
                }
            }
        }

        get_reply = 0;
        athread_put(PE_MODE,  &f[0], &lp.f[i * 3], sizeof(real) * i_size * 3, &get_reply,0,0);
        while(get_reply != 1);
        pe_syn();
    }
}

#define P_A_S (64) // pack_atoms_size 

void pack_atoms_slave(update_data_ *params)
{
    update_data_ lp;
    volatile unsigned long get_reply, put_reply;
    athread_syn(ARRAY_SCOPE, 0xffff);
    if (_MYID == 0) {
        get_reply = 0;
        athread_get(BCAST_MODE, params, &lp, sizeof(update_data_), &get_reply,0xff,0,0);
        while(get_reply != 1);
    }
    athread_syn(ARRAY_SCOPE, 0xffff);

    int natoms = lp.natoms;
    int cluster_size = ((natoms + 3) >> 2);
    int i;
    real lx[P_A_S * NBNXN_CPU_CLUSTER_I_SIZE * DIM];
    real lq[P_A_S * NBNXN_CPU_CLUSTER_I_SIZE];
    int ltype[P_A_S * NBNXN_CPU_CLUSTER_I_SIZE];
    atom atoms[P_A_S];

    for(i = _MYID * P_A_S; i < cluster_size; i += CPE_SIZE * P_A_S)
    {
        int i_size;
        if(cluster_size - i > P_A_S)
            i_size = P_A_S;
        else
            i_size = cluster_size - i;
        
        get_reply = 0;
        athread_get(PE_MODE, lp.x + i * NBNXN_CPU_CLUSTER_I_SIZE * DIM, lx, sizeof(real) * i_size * NBNXN_CPU_CLUSTER_I_SIZE * DIM, &get_reply,0,0,0);
        while(get_reply != 1);
        pe_syn();
        get_reply = 0;
        athread_get(PE_MODE, lp.q + i * NBNXN_CPU_CLUSTER_I_SIZE, lq, sizeof(real) * i_size * NBNXN_CPU_CLUSTER_I_SIZE, &get_reply,0,0,0);
        while(get_reply != 1);
        pe_syn();
        get_reply = 0;
        athread_get(PE_MODE, lp.type + i * NBNXN_CPU_CLUSTER_I_SIZE, ltype, sizeof(int) * i_size * NBNXN_CPU_CLUSTER_I_SIZE, &get_reply,0,0,0);
        while(get_reply != 1);
        pe_syn();

        {
            int i, j;
            for(i = 0; i < i_size; ++i)
            {
                for(j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; ++j)
                {
                    int id = (i << 2) + j;
                    atoms[i].x[j][XX] = lx[id * 3 + XX];
                    atoms[i].x[j][YY] = lx[id * 3 + YY];
                    atoms[i].x[j][ZZ] = lx[id * 3 + ZZ];
                    atoms[i].q[j]     = lq[id] * lp.mul;
                    atoms[i].type[j]  = ltype[id];
                }
            }
        }

        get_reply = 0;
        athread_put(PE_MODE, atoms, lp.atoms + i, sizeof(atom) * i_size, &get_reply,0,0);
        while(get_reply != 1);
        pe_syn();

    }
}


static void compute2_2_3_test( atom *a, atom *b, real *fi, real *f,const real *nbfp,
                                  real rcut2, real cpot,
                                  int ci_sh, int ci, int cj, int ntype2,
                                  unsigned int excl, real *one_Vvdw_ci)
{
    int i;
    for (i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
    {
        int ai;
        int type_i_off;
        int j;
        ai = ci*NBNXN_CPU_CLUSTER_I_SIZE + i;
        type_i_off = a->type[i]*ntype2;
        for (j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
        {
            int aj;
            real dx, dy, dz;
            real rsq, rinv;
            real rinvsq, rinvsix;
            real c6, c12;
            real FrLJ6 = 0, FrLJ12 = 0, frLJ = 0, VLJ = 0;
            real qq;
            real fcoul;
            real vcoul;
            real fscal;
            real fx, fy, fz;
            real skipmask;
            int interact;
            interact = ((excl>>(i*NBNXN_CPU_CLUSTER_I_SIZE + j)) & 1);
            skipmask = !(cj == ci_sh && j <= i);
            aj = cj*NBNXN_CPU_CLUSTER_I_SIZE + j;
            dx = a->x[i][XX] - b->x[j][XX] + a->shiftvec[XX];
            dy = a->x[i][YY] - b->x[j][YY] + a->shiftvec[YY];
            dz = a->x[i][ZZ] - b->x[j][ZZ] + a->shiftvec[ZZ];
            rsq = dx*dx + dy*dy + dz*dz;
            skipmask = (rsq >= rcut2) ? 0 : skipmask;
            rsq += (1 - interact)*NBNXN_AVOID_SING_R2_INC;
            rinv = 1.0 / sqrt(rsq);
            rinv = rinv * skipmask;
            rinvsq = rinv*rinv;
            {
                c6  = nbfp[type_i_off+b->type[j]*2 ];
                c12 = nbfp[type_i_off+b->type[j]*2+1];
                rinvsix = interact*rinvsq*rinvsq*rinvsq;
                FrLJ6 = c6*rinvsix;
                FrLJ12 = c12*rinvsix*rinvsix;
                frLJ = FrLJ12 - FrLJ6;
                VLJ = (FrLJ12 + c12 * cpot)/12 -
                    (FrLJ6 + c6 * cpot)/6;
                VLJ = VLJ * interact;
                VLJ = VLJ * skipmask;
                *one_Vvdw_ci += VLJ;
            }
            {
                fscal = frLJ*rinvsq;
            }
            fx = fscal*dx;
            fy = fscal*dy;
            fz = fscal*dz;
            fi[i*3 +XX] += fx;
            fi[i*3 +YY] += fy;
            fi[i*3 +ZZ] += fz;
            f[j*3 +XX] -= fx;
            f[j*3 +YY] -= fy;
            f[j*3 +ZZ] -= fz;
        }
    }
}






static void compute2_2_2_test( atom *a, atom *b, real *fi, real *f,const real *nbfp,
                                  real rcut2, real c_rf, real k_rf2, real k_rf, real cpot,
                                  int ci_sh, int ci, int cj, int ntype2,
                                  unsigned int excl, real *one_Vc_ci, real *one_Vvdw_ci)
{
    int i;
    for (i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
    {
        int ai;
        int type_i_off;
        int j;
        ai = ci*NBNXN_CPU_CLUSTER_I_SIZE + i;
        type_i_off = a->type[i]*ntype2;
        for (j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
        {
            int aj;
            real dx, dy, dz;
            real rsq, rinv;
            real rinvsq, rinvsix;
            real c6, c12;
            real FrLJ6 = 0, FrLJ12 = 0, frLJ = 0, VLJ = 0;
            real qq;
            real fcoul;
            real vcoul;
            real fscal;
            real fx, fy, fz;
            real skipmask;
            int interact;
            interact = ((excl>>(i*NBNXN_CPU_CLUSTER_I_SIZE + j)) & 1);
            skipmask = !(cj == ci_sh && j <= i);
            aj = cj*NBNXN_CPU_CLUSTER_I_SIZE + j;
            dx = a->x[i][XX] - b->x[j][XX] + a->shiftvec[XX];
            dy = a->x[i][YY] - b->x[j][YY] + a->shiftvec[YY];
            dz = a->x[i][ZZ] - b->x[j][ZZ] + a->shiftvec[ZZ];
            rsq = dx*dx + dy*dy + dz*dz;
            skipmask = (rsq >= rcut2) ? 0 : skipmask;
            rsq += (1 - interact)*NBNXN_AVOID_SING_R2_INC;
            rinv = 1.0 / sqrt(rsq);
            rinv = rinv * skipmask;
            rinvsq = rinv*rinv;
            {
                c6  = nbfp[type_i_off+b->type[j]*2 ];
                c12 = nbfp[type_i_off+b->type[j]*2+1];
                rinvsix = interact*rinvsq*rinvsq*rinvsq;
                FrLJ6 = c6*rinvsix;
                FrLJ12 = c12*rinvsix*rinvsix;
                frLJ = FrLJ12 - FrLJ6;
                VLJ = (FrLJ12 + c12 * cpot)/12 -
                    (FrLJ6 + c6 * cpot)/6;
                VLJ = VLJ * interact;
                VLJ = VLJ * skipmask;
                *one_Vvdw_ci += VLJ;
            }
            qq = skipmask * a->q[i] * b->q[j];
            fcoul = qq *(interact*rinv*rinvsq - k_rf2);
            vcoul = qq *(interact*rinv + k_rf*rsq - c_rf);
            *one_Vc_ci += vcoul;
            {
                fscal = frLJ*rinvsq + fcoul;
            }
            fx = fscal*dx;
            fy = fscal*dy;
            fz = fscal*dz;
            fi[i*3 +XX] += fx;
            fi[i*3 +YY] += fy;
            fi[i*3 +ZZ] += fz;
            f[j*3 +XX] -= fx;
            f[j*3 +YY] -= fy;
            f[j*3 +ZZ] -= fz;
        }
    }
}






static void compute2_2_1_test( atom *a, atom *b, real *fi, real *f,const real *nbfp,
                                  real rcut2, real c_rf, real k_rf2, real k_rf, real cpot,
                                  int ci_sh, int ci, int cj, int ntype2,
                                  unsigned int excl, real *one_Vc_ci, real *one_Vvdw_ci)
{
    int i;
    for (i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
    {
        int ai;
        int type_i_off;
        int j;
        ai = ci*NBNXN_CPU_CLUSTER_I_SIZE + i;
        type_i_off = a->type[i]*ntype2;
        for (j = 0; j < NBNXN_CPU_CLUSTER_I_SIZE; j++)
        {
            int aj;
            real dx, dy, dz;
            real rsq, rinv;
            real rinvsq, rinvsix;
            real c6, c12;
            real FrLJ6 = 0, FrLJ12 = 0, frLJ = 0, VLJ = 0;
            real qq;
            real fcoul;
            real vcoul;
            real fscal;
            real fx, fy, fz;
            real skipmask;
            int interact;
            interact = ((excl>>(i*NBNXN_CPU_CLUSTER_I_SIZE + j)) & 1);
            skipmask = !(cj == ci_sh && j <= i);
            aj = cj*NBNXN_CPU_CLUSTER_I_SIZE + j;
            dx = a->x[i][XX] - b->x[j][XX] + a->shiftvec[XX];
            dy = a->x[i][YY] - b->x[j][YY] + a->shiftvec[YY];
            dz = a->x[i][ZZ] - b->x[j][ZZ] + a->shiftvec[ZZ];
            rsq = dx*dx + dy*dy + dz*dz;
            skipmask = (rsq >= rcut2) ? 0 : skipmask;
            rsq += (1 - interact)*NBNXN_AVOID_SING_R2_INC;
            rinv = 1.0 / sqrt(rsq);
            rinv = rinv * skipmask;
            rinvsq = rinv*rinv;
            if (i < NBNXN_CPU_CLUSTER_I_SIZE/2)
            {
                c6  = nbfp[type_i_off+b->type[j]*2 ];
                c12 = nbfp[type_i_off+b->type[j]*2+1];
                rinvsix = interact*rinvsq*rinvsq*rinvsq;
                FrLJ6 = c6*rinvsix;
                FrLJ12 = c12*rinvsix*rinvsix;
                frLJ = FrLJ12 - FrLJ6;
                VLJ = (FrLJ12 + c12 * cpot)/12 -
                    (FrLJ6 + c6 * cpot)/6;
                VLJ = VLJ * interact;
                VLJ = VLJ * skipmask;
                *one_Vvdw_ci += VLJ;
            }
            qq = skipmask * a->q[i] * b->q[j];
            fcoul = qq *(interact*rinv*rinvsq - k_rf2);
            vcoul = qq *(interact*rinv + k_rf*rsq - c_rf);
            *one_Vc_ci += vcoul;
            if (i < NBNXN_CPU_CLUSTER_I_SIZE/2)
            {
                fscal = frLJ*rinvsq + fcoul;
            }
            else
            {
                fscal = fcoul;
            }
            fx = fscal*dx;
            fy = fscal*dy;
            fz = fscal*dz;
            fi[i*3 +XX] += fx;
            fi[i*3 +YY] += fy;
            fi[i*3 +ZZ] += fz;
            f[j*3 +XX] -= fx;
            f[j*3 +YY] -= fy;
            f[j*3 +ZZ] -= fz;
        }
    }
}





void nbnxn_kernel_ElecRF_VdwLJ_VF_ref_slave(ElecRF_VdwLJ_VF_ref *params)
{

    int i,j,k,l;
    ElecRF_VdwLJ_VF_ref lp;
    volatile unsigned long get_reply, put_reply;
    int ntype, natoms, nci, ntype2;
    gmx_bool do_LJ, half_LJ, do_coul, do_self;
    atom *atoms;
    nbnxn_ci_t *nbln;
    real facel, rcut2, k_rf2, k_rf, c_rf;
    atom a, b;
    int C_P_ = (1 << C_P) - 1;
    int C_T_ = (1 << C_T) - 1;
    real one_Vc_ci   = 0.0;
    real one_Vvdw_ci = 0.0;

    get_reply = 0;
    athread_get(PE_MODE, params, &lp, sizeof(ElecRF_VdwLJ_VF_ref), &get_reply,0,0,0);
    while(get_reply != 1);
    pe_syn();
    
    real nbfp[E_T_S * E_T_S * 2];
 
    get_reply = 0;
    athread_get(PE_MODE, lp.nbfp, nbfp, sizeof(real) * 2 * lp.ntype * lp.ntype, &get_reply,0,0,0);
    while(get_reply != 1);
    pe_syn();   
 

    natoms        = lp.natoms; 
    nci           = lp.nci   ; 
    ntype         = lp.ntype ; 
    rcut2         = lp.rcut2 ; 
    facel         = lp.facel ;
    atoms         = lp.atoms ; 
    ntype2        = 2 * ntype;
    k_rf          = lp.k_rf;
    k_rf2         = lp.k_rf2 ;
    c_rf          = lp.c_rf;
    rank          = lp.rank;
    
    
    volatile int *get_update_f = lp.get_update_f;
    
    
    nbnxn_ci_t l_nbln[I_CLUSTER_SIZE];
    nbnxn_cj_t l_l_cj[J_CLUSTER_SIZE];
    real f_shift[I_CLUSTER_SIZE][5];

    for(i = 0; i < I_C_S; i++)
        for(j = 0; j < 5; j++)
            f_shift[i][j] = 0;
     
    
    int cach_f, cach_t, cach_p;
    int pn, cjind, cjind1, cjind0;
    int ish, ishf, ci_sh;
    int ci, n, np, d;
    int ip, jp, f_cache_update_ps;

    atom atoms_cache[1 << C_P][1 << C_T];
    update_f f_cache[1 << C_P][1 << C_T];
    int f_cache_ps[1 << C_P], cache_ps[1 << C_P];
    real shiftvec[(NBNXN_CI_SHIFT + 1) *3 * 3 + 4];

    get_reply = 0;
    athread_get(PE_MODE, lp.shiftvec, shiftvec, sizeof(real) * (NBNXN_CI_SHIFT + 1) * 3 * 3, 
            &get_reply,0,0,0);
    while(get_reply != 1);
    pe_syn();   
      
    int cach_size = 1 << C_P;
    for(i = 0;i < cach_size;i++)
      cache_ps[i] = f_cache_ps[i] = -1;
   




    real fi[NBNXN_CPU_CLUSTER_I_SIZE * DIM], f[NBNXN_CPU_CLUSTER_I_SIZE * DIM];
    int get = 0, miss = 0;
    int np_st = get_update_f[ _MYID];
    int np_ed;
    if(_MYID == 63)
        np_ed = nci;
    else
        np_ed = get_update_f[ _MYID + 1];

    int cluster_number = ((((natoms + 3) >> 2) + 8)&(~7));

    update_f *update_data = &lp.update_date[_MYID * cluster_number];

    long *f_cache_bmap_h = (long*) f_cache_bmap_host;
    int ncache_line = (cluster_number + (1 << C_T) - 1) / (1 << C_T);
    int ncache_bitmap = (ncache_line + 63) >> 6;
    long f_cache_bmap[F_CACHE_BMAP_SIZE];
    // long cache_bitmap[((ncache_bitmap + 3) >> 2) << 2];
    {
        doublev4 zero_dv4;
        zero_dv4 = 0;
        #pragma unroll[4]
        for (i = 0; i < ncache_bitmap; i += 4) {
            simd_store(zero_dv4, f_cache_bmap + i + 0);
        }
    }



    
    for(np = np_st; np < np_ed; np += I_C_S)
    {
        int i_size;
        if(np_ed - np > I_C_S)
            i_size = I_C_S;
        else
            i_size = np_ed - np;
        get_reply = 0;
        athread_get(PE_MODE, lp.nbln + np, l_nbln, sizeof(nbnxn_ci_t) * i_size, &get_reply,0,0,0);
        while(get_reply != 1);
        pe_syn();

        for(n = 0; n < i_size;n++)
        {
            nbln    = &l_nbln[n];
            ish     = (nbln->shift & NBNXN_CI_SHIFT);
            ishf    = ish*DIM;
            cjind0  = nbln->cj_ind_start;
            cjind1  = nbln->cj_ind_end;
            ci      = nbln->ci;
            ci_sh   = (ish == CENTRAL ? ci : -1);
            do_LJ   = (nbln->shift & NBNXN_CI_DO_LJ(0));
            do_coul = (nbln->shift & NBNXN_CI_DO_COUL(0));
            half_LJ = ((nbln->shift & NBNXN_CI_HALF_LJ(0)) || !do_LJ) && do_coul;
            do_self = do_coul;
            
            one_Vc_ci   = 0.0;   
            one_Vvdw_ci = 0.0;   
            get++;
            
            int ci_t = ci & ((1 << C_T) - 1);
            int ci_p = (ci >> C_T) & ((1 << C_P) - 1);
            if(cache_ps[ci_p] != ci - ci_t)
            {
                miss++;
                get_reply = 0;
                athread_get(PE_MODE, atoms + ci - ci_t, atoms_cache[ci_p], 
                                     (1 << C_T) * sizeof(atom), &get_reply,0,0,0);
                while(get_reply != 1); 
                pe_syn();
                cache_ps[ci_p] = ci - ci_t;
            }

            for (i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
            {
                for (d = 0; d < DIM; d++)
                {
                    fi[i*3 +d] = 0; 
                }
            }

            atom a = atoms_cache[ci_p][ci_t];
            for (d = 0; d < DIM; d++)
            {
                a.shiftvec[d] = shiftvec[ishf+d]; 
            }

            if (do_self)
            {
                real Vc_sub_self;
                Vc_sub_self = 0.5 * c_rf;
                if (lp.l_cj[nbln->cj_ind_start].cj == ci_sh)
                {
                    for (i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
                    {
                        one_Vc_ci -= a.q[i] * a.q[i] * Vc_sub_self;
                    }
                }
            }
            for(cjind = cjind0; cjind < cjind1; cjind+= J_C_S)
            {
                int j_size = J_C_S;
                if(J_C_S > cjind1 - cjind)
                    j_size = cjind1 - cjind;
                
                
                get_reply = 0;
                athread_get(PE_MODE, lp.l_cj + cjind, l_l_cj, sizeof(nbnxn_cj_t) * j_size, 
                        &get_reply,0,0,0);
                while(get_reply != 1);
                pe_syn();

                int cj_ps;
                for(cj_ps = 0; cj_ps < j_size; cj_ps++)
                {
                    
                    int cj = l_l_cj[cj_ps].cj;
                    int cj_t = cj & ((1 << C_T) - 1);
                    int cj_p = (cj >> C_T) & ((1 << C_P) - 1);
                    
                    get++;
                    if(cache_ps[cj_p] != cj - cj_t)
                    {
                        miss++;
                        get_reply = 0;
                        athread_get(PE_MODE, atoms + cj - cj_t, atoms_cache[cj_p], 
                                             (1 << C_T) * sizeof(atom), &get_reply,0,0,0);
                        cache_ps[cj_p] = cj - cj_t;          
                        while(get_reply != 1); 
                        pe_syn();
                    }
                    
                    atom b = atoms_cache[cj_p][cj_t];

                    for(i = 0; i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
                    {
                        f[i * 3 + XX] = 0;
                        f[i * 3 + YY] = 0;
                        f[i * 3 + ZZ] = 0;
                    }

                    if(half_LJ)
                    {
                        compute2_2_1_test( &a, &b, fi, f, nbfp, 
                                          rcut2, c_rf, k_rf2, k_rf, lp.cpot, 
                                          ci_sh, ci, cj, ntype2,
                                          l_l_cj[cj_ps].excl, &one_Vc_ci, &one_Vvdw_ci);
            
                    } else if(do_coul)
                    {

                        compute2_2_2_test( &a, &b, fi, f, nbfp, 
                                          rcut2, c_rf, k_rf2, k_rf, lp.cpot, 
                                          ci_sh, ci, cj, ntype2,
                                          l_l_cj[cj_ps].excl, &one_Vc_ci, &one_Vvdw_ci);
 
                    } else
                    {
                        compute2_2_3_test( &a, &b, fi, f, nbfp, 
                                          rcut2, lp.cpot, 
                                          ci_sh, ci, cj, ntype2,
                                          l_l_cj[cj_ps].excl, &one_Vvdw_ci);
             
                    }


                    if(f_cache_ps[cj_p] != cj - cj_t)
                    {
                        if(f_cache_ps[cj_p] != -1)
                        {
                            get_reply = 0;
                            athread_put(PE_MODE, f_cache[cj_p], &update_data[f_cache_ps[cj_p]], 
                                            (1 << C_T) * sizeof(update_f), &get_reply, 0, 0);
                            while(get_reply != 1);
                            pe_syn();
                            set_f_cache_bmap(f_cache_bmap, f_cache_ps[cj_p] >> C_T);

                            if (get_f_cache_bmap(f_cache_bmap, (cj - cj_t) >> C_T)) {
                                get_reply = 0;
                                athread_get(PE_MODE, &update_data[cj - cj_t], f_cache[cj_p], 
                                        (1 << C_T) * sizeof(update_f), &get_reply, 0, 0, 0);
                                while(get_reply != 1);
                                pe_syn();
                            }
                            else {
                                reset_f_cache(cj_p);
                            }

                        } else
                        {
                            reset_f_cache(cj_p);
                        }
                        f_cache_ps[cj_p] = cj - cj_t;
                    }

                    add_realv4_f(&(f_cache[cj_p][cj_t].f), f);
                    // for(i = 0;i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
                    // {
                    //     f_cache[cj_p][cj_t].f[i][XX] += f[i * 3 + XX];
                    //     f_cache[cj_p][cj_t].f[i][YY] += f[i * 3 + YY];
                    //     f_cache[cj_p][cj_t].f[i][ZZ] += f[i * 3 + ZZ];
                    // }
                }
            
            }
    
            if(f_cache_ps[ci_p] != ci - ci_t)
            {
               
                if(f_cache_ps[ci_p] != -1)
                {
                    get_reply = 0;
                    athread_put(PE_MODE, f_cache[ci_p], &update_data[f_cache_ps[ci_p]], 
                                            (1 << C_T) * sizeof(update_f), &get_reply, 0, 0);
                    while(get_reply != 1);
                    pe_syn();
                    set_f_cache_bmap(f_cache_bmap, f_cache_ps[ci_p] >> C_T);

                    if (get_f_cache_bmap(f_cache_bmap, (ci - ci_t) >> C_T)) {
                        get_reply = 0;
                        athread_get(PE_MODE, &update_data[ci - ci_t], f_cache[ci_p], 
                                (1 << C_T) * sizeof(update_f), &get_reply, 0, 0, 0);
                        while(get_reply != 1);
                        pe_syn();
                    }
                    else {
                        reset_f_cache(ci_p);
                    }

                } else
                {
                    reset_f_cache(ci_p);
                }
                f_cache_ps[ci_p] = ci - ci_t;
            }

            add_realv4_f(&(f_cache[ci_p][ci_t].f), fi);
            // for(i = 0;i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
            // {
            //     f_cache[ci_p][ci_t].f[i][XX] += fi[i * 3 + XX];     
            //     f_cache[ci_p][ci_t].f[i][YY] += fi[i * 3 + YY];     
            //     f_cache[ci_p][ci_t].f[i][ZZ] += fi[i * 3 + ZZ];     
            // }


            for(i = 0;i < NBNXN_CPU_CLUSTER_I_SIZE; i++)
            {
                f_shift[n][XX] += fi[i * 3 + XX];
                f_shift[n][YY] += fi[i * 3 + YY];
                f_shift[n][ZZ] += fi[i * 3 + ZZ];
            }
            f_shift[n][3] = one_Vc_ci;
            f_shift[n][4] = one_Vvdw_ci;
        }

        get_reply = 0;
        athread_put(PE_MODE, &f_shift[0][0], &lp.fshift_CPE[np * 5], 5 * i_size * sizeof(real), 
                             &get_reply,0,0);
        while(get_reply != 1);
        pe_syn();

        for(i  = 0;i < i_size;i++)
        {
            f_shift[i][XX] = 0;
            f_shift[i][YY] = 0;
            f_shift[i][ZZ] = 0;
        }
    }

    for(i = 0;i < (1 << C_P);i++)
    {
        if(f_cache_ps[i] != -1)
        {
            get_reply = 0;
            athread_put(PE_MODE, f_cache[i], &update_data[f_cache_ps[i]], 
                                    (1 << C_T) * sizeof(update_f), &get_reply, 0, 0);
            set_f_cache_bmap(f_cache_bmap, f_cache_ps[i] >> C_T);
            while(get_reply != 1);
            pe_syn();
        }
    }
    get_reply = 0;
    athread_put(PE_MODE, f_cache_bmap, f_cache_bmap_h + _MYID * F_CACHE_BMAP_SIZE, 
            F_CACHE_BMAP_SIZE * sizeof(long), &get_reply, 0, 0);
    while(get_reply != 1);
    pe_syn();

}


