/* 
Jun S. Song 2005
Cliff A. Meyer 2010
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include "genomescan_Markov0_3.h"
//#include "mt64.h"

#define MINCOUNT   20        // minimum number of hits required to calculate centrality zscore 
#define TABLESIZE  1048576   // Size of each lookup table
#define SEQSIZE    10        // Size of subsequence in the lookup table
#define BITSHIFT   20
#define FWD    0     // + DNA strand
#define RVS    1     // - DNA strand
#define MIN(x,y) (x) < (y) ? (x) : (y)
#define max(x,y) (x) > (y) ? (x) : (y)
#define MONTECARLOITERS 1000
#define MINFREQ 0.01 // minimum frequency in PSSM matrix

/*  FUNCTIONS TO LOAD FILES */
//extern int getNOfFiles( char *dirname);
//extern char **getfiles( char *dirname, int nOfFiles);
//extern void usage(void);
//extern void initialize8910(void);

/* HEADER INFORMATION */
//struct head {
//  unsigned int start ;  // start of bed sequence; will be set to 0 if not given
//  unsigned int end ;    // end of bed sequence; will be set to 0 if not given 
//  char chrname[20];     // 
//} header;

/*
struct motif_T{
    int iseq;
    int istart;
    int iend;
    float score;
    int lenseq;
    int orient;
};
*/

/* GLOBAL VARIABLES */
int (*pssm)[4];             // PSSM Matrix
int (*lookup)[TABLESIZE];   // PSSM score look up table in 10 basepairs
int motiflength = 0; 
int hashlength;             // Number of hash tabled needed
int *motifbg;               // position specific background prob along the motif

/* PROTOTYPE */
void buildHash(int nOfHash);
int  loadPSSM(char *filename, float *entropy, int *bgprob);
//void parseHeader(char *headerline);

/* LOAD PSSM MATRICES */
int passPSSM( float **pssm_p, int length,  float *entropy, int *bgprob){
    register unsigned int i,j;  // counting variables
    float entrpy = 0.0, ftot = 0;

    /* Allocate Memory for the PSSM matrix. */
    pssm = malloc(length*4*sizeof(int));

    /* Verify allocation */
    if(!pssm){
        printf("Memory request failed.\n");
        exit(1);
    }

    // adjust minimum frequency
    for(i=0; i<length; i++){
        ftot = 0.0;
        for( j = 0; j < 4; j++ ){
            if( pssm_p[i][j] <= 0.0 ){
                pssm_p[i][j] = MINFREQ;
            }
            ftot += pssm_p[i][j];
        }

        for(j=0; j<4; j++){
            /*pssm_p[i][j] = 100000*log( pssm_p[i][j]/ftot );*/
            pssm[i][j] = (int)(100000*log(pssm_p[i][j]/ftot));
        }
        
    }

  /* Load the matrix elements */
  for(i=0; i < length; i ++){
    entrpy += exp(*(pssm[i]  )/100000.0)* ( *(pssm[i])   - bgprob[0]);
    entrpy += exp(*(pssm[i]+1)/100000.0)* ( *(pssm[i]+1) - bgprob[1]);
    entrpy += exp(*(pssm[i]+2)/100000.0)* ( *(pssm[i]+2) - bgprob[2]);
    entrpy += exp(*(pssm[i]+3)/100000.0)* ( *(pssm[i]+3) - bgprob[3]);
  }
  
  *entropy = entrpy/100000.0;

  return length;
}

/* LOAD PSSM MATRICES */
int loadPSSM(char *filename, float *entropy, int *bgprob){
  FILE *pssmfile;
  int length;  // length of the PSSM matrix
  register unsigned int i;  // couting variable
  float entrpy = 0.0;
  pssmfile = fopen(filename,"r");
  if (pssmfile == NULL) perror ("Error opening PSSM file");
  
  fscanf(pssmfile, "%d", &length);
  
  /* Allocate Memory for the PSSM matrix. */
  pssm = malloc(length*4*sizeof(int));
  /* Verify allocation */
  if(!pssm){
    printf("Memory request failed.\n");
    exit(1);
  }

  /* Load the matrix elements */
  for(i=0; i < length; i ++){
    fscanf(pssmfile, "%d%d%d%d", pssm[i], pssm[i]+1, pssm[i]+2, pssm[i]+3);
    entrpy += exp(*(pssm[i])/100000.0)* ( *(pssm[i]) - bgprob[0]);
    entrpy += exp(*(pssm[i]+1)/100000.0)* ( *(pssm[i]+1) - bgprob[1]);
    entrpy += exp(*(pssm[i]+2)/100000.0)* ( *(pssm[i]+2) - bgprob[2]);
    entrpy += exp(*(pssm[i]+3)/100000.0)* ( *(pssm[i]+3) - bgprob[3]);
  }
  fclose(pssmfile);

  *entropy = entrpy/100000.0;
  return length;
}

/* BUILD HASH TABLES */
void buildHash(int nOfHash){
  register unsigned int i, j, k,m, n,p,q, bp, score;
  /* allocate memory */
  lookup = malloc(nOfHash*TABLESIZE*sizeof(int));
  if(!lookup){
    printf("Memory request failed: lookup.\n");
    exit(1);
  }

  for (i = 0; i < nOfHash; i++){
    k = (i==(nOfHash-1))? (motiflength % SEQSIZE) : SEQSIZE;
    if (k  ==0) k = SEQSIZE;
    q = SEQSIZE * i;
    m = pow(4,k);
    
    for(j = 0; j < m; j ++) {
      p = j;
      score = 0;
      for( n =0; n < k; n ++){
     bp = (p & 0x3);
     score += pssm[n+q][bp];
     p = p >> 2;
      }
      
      lookup[i][j] = score;
      //printf("%d %d %d\n", i,j, lookup[i][j]);
    } 
  }

  motifbg = malloc(motiflength*sizeof(int));
  if(!motifbg){
    printf("Memory request failed: motifbg.\n");
    exit(1);
  }
}


float calc_score( unsigned int bp2int, int bgmkv, int reset, int *bgprob, int *cond2g1, int *cond3g2 ){

    static unsigned int i, st1, st2, st3, st4, st5, st6, st7, st8, st9, st10;
    static unsigned int  MAX, count, pos, shift, offset;
    static unsigned long long int motif, tmpmotif;
    static int bgscore, pssmscore;
    static int UPPERLIMIT;
    float combinedscore = 0.0;

    MAX = motiflength -1;
    UPPERLIMIT = max(bgmkv,MAX);

    if ( reset == 1 ){
        i=0, st1=0, st2=0, st3=0, st4=0, st5=0, st6=0, st7=0, st8=0, st9=0, st10=0;
        count=0, pos=0, shift=0, offset=0;
        motif=0LL, tmpmotif = 0LL;
        bgscore = 0, pssmscore = 0;
        combinedscore = 0.0;
        MAX = motiflength -1;
        UPPERLIMIT = max(bgmkv,MAX);
        motifbg[0] = 0;
        shift = MAX*2;  // amount by which to shift the motif
        return 0;
    }
 
    motif = (motif >> 2) + ((long long int) bp2int << shift);

    st2  = (st2>>2)  + (bp2int<<2);
    st3  = (st3>>2)  + (bp2int<<4);
    st4  = (st4>>2)  + (bp2int<<6);
    st5  = (st5>>2)  + (bp2int<<8);
    st6  = (st6>>2)  + (bp2int<<10);
    st7  = (st7>>2)  + (bp2int<<12);
    st8  = (st8>>2)  + (bp2int<<14);
    st9  = (st9>>2)  + (bp2int<<16);
    st10 = (st10>>2) + (bp2int<<18);

    //fprintf(stdout,"%d\n", bgmkv);
    if ( count < UPPERLIMIT){
        count++;
        if(count <= MAX){
            switch ( count){
            case 1: offset += (bp2int << (2*(count-1))); motifbg[1] =  bgprob[bp2int]; bgscore += motifbg[1]; break;
            case 2: offset += (bp2int << (2*(count-1))); 
                if (bgmkv == 0) motifbg[2] = bgprob[bp2int];
                else            motifbg[2] = cond2g1[offset]; 
                bgscore += motifbg[2]; break;
            case 3: offset += (bp2int << (2*(count-1))); 
                if (bgmkv == 0)    motifbg[3] = bgprob[bp2int];
                else if(bgmkv ==1) motifbg[3] = cond2g1[st2];
                else               motifbg[3] = cond3g2[offset];
                bgscore += motifbg[3]; break;
            default:
                switch(bgmkv){
                    //case 3: motifbg[count] = cond4g3[st4]; break;
                    case 2: motifbg[count] = cond3g2[st3]; break;
                    case 1: motifbg[count] = cond2g1[st2]; break;
                    case 0: motifbg[count] = bgprob[bp2int]; break;  
                }
                bgscore += motifbg[count]; break;

            }
            /*continue;*/
        }
        else {  /* if motif length is less than 10 */
            bgscore -= motifbg[0];
            for (i = 0; i < MAX; i ++){
                motifbg[i] = motifbg[i+1];
            }
            offset += (bp2int << (2*(count-1)));
            switch(count){
            case 2: motifbg[MAX] = cond2g1[offset]; break;    
            case 3: motifbg[MAX] = cond3g2[offset]; break;
            }
            bgscore += motifbg[MAX];
        }
    }
    else{
        bgscore -= motifbg[0];
        for (i = 0; i < MAX; i ++){
            motifbg[i] = motifbg[i+1];
        }
        switch(bgmkv){
        //case 3: motifbg[MAX] = cond4g3[st4]; break;
        case 2: motifbg[MAX] = cond3g2[st3]; break;
        case 1: motifbg[MAX] = cond2g1[st2]; break;
        case 0: motifbg[MAX] = bgprob[bp2int]; break;
        }
        bgscore += motifbg[MAX];
    }

    pssmscore = 0; 
    tmpmotif = motif;
    for (i = 0; i < hashlength; i ++){
        int subseq = (int)(tmpmotif & 0xfffff);
        pssmscore += lookup[i][subseq];
        tmpmotif = (tmpmotif >> BITSHIFT);
    }

    /*fprintf( stdout, "%d\t%d\n", pssmscore, bgscore );*/
    combinedscore = (pssmscore - bgscore)*1.0/100000;
    return combinedscore;
}

/*
  GET PSSM SCORES FROM SEQUENCE FILES
  array of encoded sequences 
  indices of subsequences to be searched: seqindex start end
*/

void getPSSMScores_pp( int **totdata, int **locations, int nlocations, int bgmkv, struct motif_T *motif_record, int *bgprob, int *cond2g1, int *cond3g2, int score_option ){

    //FILE *seqfile; FILE *outfile;
    //register unsigned int i;
    //st1=0, st2=0, st3=0, st4=0, st5=0, st6=0, st7=0, st8=0, st9=0, st10=0;
    register unsigned int bp2int, MAX, count=0, shift, offset=0;
    //register unsigned long long int motif=0LL, tmpmotif = 0LL;
    register int bgscore = 0;
    int UPPERLIMIT;
    int i,j,k;
    int ir;
    int idx;
    int iseq, istart, iend, lenseq;
    unsigned int loc;
    float combinedscore = 0.0;
    float norm = 0.0; //norm const
    float eps  = 1e-5;

    MAX = motiflength-1;
    UPPERLIMIT = max(bgmkv,MAX);

    loc = 0;
    shift = MAX*2;  // amount by which to shift the motif
    motifbg[0] = 0;

    for ( idx=0, k=0; k < nlocations; k++ ){ 

        iseq   = locations[k][0];
        istart = locations[k][1];
        iend   = locations[k][2];
        lenseq = locations[k][3];

        if ( score_option == MEAN_OPTION || score_option == MAX_OPTION ){
            (motif_record[k]).iseq = (int)iseq;
            (motif_record[k]).lenseq = lenseq;
            (motif_record[k]).score = 0.0;
        }

        //(motif_record[k]).pscore = 0.0;
        norm = 0.0;
        //pnorm = 0.0;

        //fprintf( stdout, "---- idx, k, iseq, istart, iend %d\t%d\t%d\t%d\t%d\n", idx, k, iseq, istart, iend );
        /* Forward */
        calc_score( 0, 0, 1, bgprob, cond2g1, cond3g2 );
        for( j = istart; j < iend; j++){

            bp2int = (unsigned int)( totdata[iseq][j] ); 
            if ( bp2int >= 0 && bp2int < 4 ){
                combinedscore = calc_score( bp2int, bgmkv, 0, bgprob, cond2g1, cond3g2 );
                //fprintf( stdout, "FWD%d\t%d\t%f\n",k ,j, exp(combinedscore) );
                //remove 200: if ( j - 200) >= (motiflength - 1)){      
                if ( j >= (motiflength - 1)){      
                // do not use CUTOFF 
                //if( combinedscore > CUTOFF && idx < MAXMOTIFS - 1 ){                
                    if ( score_option == MEAN_OPTION ){

                        (motif_record[k]).score  += exp( combinedscore );
 
                    }else if ( score_option == MAX_OPTION ){

                        if ( motif_record[(((idx-1)>(k))?(idx-1):(k))].score < exp(combinedscore) + eps && idx < MAXMOTIFS ) {
                            (motif_record[idx]).score  = exp(combinedscore);
                            (motif_record[idx]).istart = j - motiflength + 1; 
                            (motif_record[idx]).iend   = j+1; 
                            (motif_record[idx]).orient = FWD;
                            (motif_record[idx]).iseq = (int)iseq;
                            (motif_record[idx]).lenseq = lenseq;
                            idx++;
                        }

                    }else if( score_option == CUTOFF_OPTION ){

                        if ( exp(combinedscore) > CUTOFF  && idx < MAXMOTIFS){
                            (motif_record[idx]).iseq = (int)iseq;
                            (motif_record[idx]).score = exp(combinedscore);
                            (motif_record[idx]).istart = j - motiflength + 1; 
                            (motif_record[idx]).iend   = j+1; 
                            (motif_record[idx]).orient = FWD;
                            idx++;
                        } 
                    }
 
                    norm += 1.0;
                }
                               
            }else{

                count = 0; motifbg[0] = 0; bgscore = 0; offset = 0;

            }



        }

        /* reverse */
        /*count=0; bgscore =0; pssmscore = 0;motif=0LL; tmpmotif = 0LL; combinedscore = 0.0; motifbg[0] = 0;
        offset = 0;*/

        calc_score( 0, 0, 1, bgprob, cond2g1, cond3g2 );

        for( j = ( iend - 1 ); j > ( istart - 1 ); j-- ){

            bp2int = (unsigned int)( 3 - totdata[iseq][j] );

            if ( bp2int >= 0 && bp2int < 4 ){
                combinedscore = calc_score( bp2int, bgmkv, 0, bgprob, cond2g1, cond3g2 );

                if( j <= (lenseq - motiflength)){ 
                // do not use CUTOFF
                //if( combinedscore > CUTOFF && idx < MAXMOTIFS - 1 ){      
                    if ( score_option == MEAN_OPTION ){

                        (motif_record[k]).score  += exp( combinedscore );
                    
                    }else if ( score_option == MAX_OPTION ){

                        if ( motif_record[(((idx-1)>(k))?(idx-1):(k))].score < exp(combinedscore + eps ) && idx < MAXMOTIFS ) {
                            (motif_record[idx]).score  = exp( combinedscore );
                            (motif_record[idx]).istart = j; 
                            (motif_record[idx]).iend   = j + motiflength; 
                            (motif_record[idx]).orient = RVS;
                            (motif_record[idx]).iseq   = (int)iseq;
                            (motif_record[idx]).lenseq = lenseq;
                            idx++;
                        }

                    }else if( score_option == CUTOFF_OPTION ){
 
                        if ( exp(combinedscore) > CUTOFF  && idx < MAXMOTIFS){
                            (motif_record[idx]).iseq   = (int)iseq;
                            (motif_record[idx]).lenseq = lenseq;
                            (motif_record[idx]).score  = exp(combinedscore);
                            (motif_record[idx]).istart = j; 
                            (motif_record[idx]).iend   = j + motiflength; 
                            (motif_record[idx]).orient = RVS;
                            idx++;
                        }

                    }
                    norm += 1.0;
                } 
                
            }else{
                count = 0; motifbg[0] = 0; bgscore = 0; offset = 0;
            }
        }

        if (score_option == MEAN_OPTION ){

            (motif_record[k]).score /= norm;
        
        }
        else if ( score_option == MAX_OPTION ){

            // If max option is being used a random best motif needs to be selected
            // Scores for idx > k need to be reset 
            for ( i=idx-1, j=idx-1; j>k-1; j-- ){ 
                // find lowest index in set of max scores
                if( (motif_record[j]).score > ( (motif_record[idx-1]).score - eps ) ){
                    i = j;
                }
                // fprintf( stdout, "score %d %d %d %d %4.2f %4.2f\n", k, idx, j, i, (motif_record[idx-1]).score, (motif_record[j]).score );
            }

            // randomly pick an index between the lowest and highest indices 
            ir = (int)((idx-i)*((double)(rand())/((double)(RAND_MAX)+(double)(1.0))));

            // fprintf( stdout, "random int %d %d %d %d %4.2f %4.2f\n", k, idx-i, ir, i, (motif_record[idx-1]).score, (motif_record[i+ir]).score );

            (motif_record[k]).iseq   = (motif_record[i+ir]).iseq;
            (motif_record[k]).score  = (motif_record[i+ir]).score;
            (motif_record[k]).istart = (motif_record[i+ir]).istart; 
            (motif_record[k]).iend   = (motif_record[i+ir]).iend; 
            (motif_record[k]).orient = (motif_record[i+ir]).orient;

            //reset records beyond current sequence
            for ( j=k+1; j<idx; j++ ){ 
                (motif_record[j]).iseq   = -1; 
                (motif_record[j]).score  =  0;
                (motif_record[j]).istart = -1;
                (motif_record[j]).iend   = -1;
                (motif_record[j]).orient = -1;
                (motif_record[j]).lenseq = -1;
            }
    
            // reset index for new records
            idx = k+1;

        }
        //(motif_record[k]).pscore /= pnorm;
    }

    //( motif_record[idx] ).istart = -1; ( motif_record[idx] ).score = -1;   
    //( motif_record[idx] ).istart = -1; ( motif_record[idx] ).score = 0;   
 
    //qsort( motif_record, MAXMOTIFS, sizeof(struct motif_T), (void *)comp_position );

    //for ( i = 0; i < MAXMOTIFS-1; i++ ){
    //   if ( motif_record[i+1].iseq == motif_record[i].iseq ){
    //      if ( motif_record[i+1].istart < motif_record[i].istart + motiflength ){
    //         if ( motif_record[i].score > motif_record[i+1].score )
    //            motif_record[i].score = 0.0; 
    //         else
    //            motif_record[i+1].score = 0.0;
    //      }
    //   }
    //}

    //qsort( motif_record, MAXMOTIFS, sizeof(struct motif_T), (void *)comp_motif );
    return;

}

//void motifscan_subseq( int **seq, int **locations, int nlocations, float **pssm_p, int motiflen, int bgmkv, float **score, int nscore, int *bgprob, int *cond2g1, int *cond3g2, int score_option ){
void motifscan_subseq( int **seq, int **locations, int nlocations, float **pssm_p, int motiflen, int bgmkv, struct motif_T *motif_record, int nscore, int *bgprob, int *cond2g1, int *cond3g2, int score_option ){
    // float CUTOFF = 0.0;
    float entropy; 

    motiflength = motiflen;
    passPSSM( pssm_p, motiflength, &entropy, bgprob );

    // CUTOFF = 0.1*entropy;
    hashlength = (int) ceil( motiflength*1.0/10 );
    buildHash(hashlength);
    getPSSMScores_pp( seq, locations, nlocations, bgmkv, motif_record, bgprob, cond2g1, cond3g2, score_option );

    free(lookup);
    free(pssm);
    free(motifbg);
    return;

} 
