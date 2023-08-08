#include <slave.h>
#include <simd.h>
#define NCPE (64)
typedef float real;

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
typedef struct memcpy_ref{
    float *out;
    uint64_t *recv_ptr;
    int message_number,numprocs;
} memcpy_ref;
void memcpy_float_slave(memcpy_ref *hparam){
#define DLEN (1024)
    float buffer[DLEN];
    float *xout;
    float *xin;
    memcpy_ref s_param;

    uint64_t recv_ptr[512];
    int numprocs;
    int i, j, k;
    int myid = _MYID;

    bcast_get(hparam, &s_param, sizeof(memcpy_ref));

    xout     = s_param.out;
    numprocs = s_param.numprocs;

    int cnt = s_param.message_number;

    bcast_get(s_param.recv_ptr, recv_ptr, sizeof(uint64_t)*numprocs);

    // return;

    for(k = myid / 4; k < numprocs; k += 16) {
        xin = recv_ptr[k];
        for(i = (myid % 4) * DLEN; i < cnt; i += 4 * DLEN) {
            int i_size = i + DLEN < cnt ? DLEN : cnt - i;

            pe_get(xin + i, buffer, i_size * sizeof(float));
            pe_put(buffer, xout + k*cnt + i, i_size * sizeof(float));

        }
    }
#undef DLEN
}
