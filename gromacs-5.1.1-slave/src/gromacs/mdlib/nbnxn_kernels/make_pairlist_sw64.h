#define NCPE (64)
#define MAX_NCI (20000)
#define MAX_NCJ (2500000)
#define MAX_1NCJ (512)
#define MAX_LNCI (256)

#define MAX_NCI_1SLAVE (128 * 200)
#define MAX_NCJ_1SLAVE (128 * 2000)
#define MAX_NCI_64SLAVE (64 * MAX_NCI_1SLAVE)
#define MAX_NCJ_64SLAVE (64 * MAX_NCJ_1SLAVE)

#define NBNXN_INTERACTION_MASK_ALL        0xffffffffU
#define NBNXN_INTERACTION_MASK_DIAG       0x08ceU
#define NBNXN_CPU_CLUSTER_I_SIZE 4
#define STRIDE_XYZ   3
#define NNBSBB_C         4
#define NNBSBB_D         2
#define BB_X  0
#define BB_Y  1
#define BB_Z  2
#define DIM   3

struct nbnxn_ci{
    int ci;
    int shift;
    int cj_ind_start;
    int cj_ind_end;
};
struct nbnxn_cj{
    int          cj;
    unsigned int excl;
};
 struct nbnxn_bb {
    float lower[NNBSBB_C];
    float upper[NNBSBB_C];
};

extern int nbl_nci_host[], nbl_ncj_host[];
extern int nbl_old_nci, nbl_old_ncj;
extern struct nbnxn_ci nbl_ci_host[];
extern struct nbnxn_cj nbl_cj_host[];

struct mklist_ind_t{
    int ci, ci_x, ci_y;
    int n;
};
typedef struct make_pairlist_ref{
    /* func params */
    const nbnxn_grid_t *gridi;
    const nbnxn_grid_t *gridj;
    const int *flags_i;
    const nbnxn_bb_t *bb_i;
    const nbnxn_atomdata_t *nbat;
    nbnxn_pairlist_t *nbl;
    const t_blocka *excl;
    int conv_i, nth, th, ci_block, cell0_i, nb_kernel_type;
    int gridj_flag_shift, gridi_flag_shift;
    int  na_cj_2log, nsubpair_max, progBal;
    float rbb2;
    const float *bbcz_i;
    const float *bbcz_j;
    float nsubpair_tot_est;
    real rl_fep2;
    real rl2;
    const struct nbnxn_search* nbs;
    ivec shp;

    /* parallel params */
    struct mklist_ind_t *ind_start;
    
    /* debug params */
    int mpi_rank;
}make_pairlist_ref;

typedef struct mklist_gather_ref{
    struct nbnxn_ci *ci, *ci64;
    struct nbnxn_cj *cj, *cj64;
    int nci[NCPE], ncj[NCPE];
    int old_nci, old_ncj;

    /* debug params */
    int mpi_rank;
}mklist_gather_ref;

