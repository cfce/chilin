#define MAXMOTIFS 100000
#define MINSCORE  0     // minimum PSSM matrix score to be considered as a possible motif hit  
#define MEAN_OPTION 0 // return the mean probability 
#define MAX_OPTION 1 // return the max probability 
#define CUTOFF_OPTION 2 // return the probability for all sites above CUTOFF
#define CUTOFF 100   // threshold to return motifs

struct motif_T{
    int iseq;
    int istart;
    int iend;
    int lenseq;
    int orient;
    float score;
};

void motifscan_subseq( int **seq, int **locations, int nlocations, float **pssm_p, int motiflen, int bgmkv, struct motif_T *motif_record, int nscore, int *bgprob, int *cond2g1, int *cond3g2, int score_option );

