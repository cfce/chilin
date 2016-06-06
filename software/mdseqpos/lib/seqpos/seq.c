#include <Python.h>
/*#include <Numeric/arrayobject.h>*/
#include <numpy/arrayobject.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include "genomescan_Markov0_3.h"
#define SEQSCAN_DOC "Scan sequence for occurrence of motif. seqscan( list_of_sequences, PSSM_matrix, backgroup probability matrices, markov order, option ).\
PSSM_matrix is a motifLength x 4 Float Numeric array where columns are A, C, G, T fractions. No base at any position may have a frequency less that 0.01, if this occurs the frequencies are adjusted to meet this requirement.   \
background prob matrix computed using count.py \
Markov order 0,1,2\
probability option 0 mean 1 max\
returns numpy arrays: sequence index, start, end, orientation, score" 

//#define CUTOFF 0
//#define bgmkv  2

static PyObject *pyError;

static PyObject *seqscan(PyObject *self, PyObject *args)
{

    int i,j;
    double tmp = 0;
    int numLines;      /* how many lines we passed for parsing */
    int numHits;       /* how many motif hits found */
    int nbg;           /* number of background parameters passed */
    int numChars;      /* how many characters we passed for parsing */
    int **seqarray;    /* integer array of encoded nucleotide strings */
    int **locations;   /* integer array of regions of seqarray to be scanned for motifs */    
    int bp2int;
    char *line;        /* pointer to the line as a string */
    char tmpstr[100];
    PyObject *listObj; /* the list of strings */
    PyObject *strObj;  /* one string in the list */
    PyObject *listBG; /* the list of background frequencies */
    PyObject *orderBG; /* order of the Markov model */
    /*PyObject *pssmObj;  python pssm  */
    PyObject *Py_ret;  /* motif location array and motif score array */
    PyArrayObject *pssmArrayObj; /* Numeric pssm matrix */
    /*PyArrayObject *Py_MotifScore;  Numeric array of motif scores */
    PyArrayObject *motifIdx;    /* index of sequence*/
    PyArrayObject *motifStart;  /* start of hit within seq        */
    PyArrayObject *motifEnd;    /* end of hit within seq +1       */
    PyArrayObject *motifOrient; /* strand on which motif is found */
    PyArrayObject *motifScore;  /* score of hit                   */
    int motifLength;   /* length of motif */
    int motifBases;    /* number of bases motif */
    float ** pssm_p;   
    struct motif_T * motif_record;
    int n_motif_record; /* space allocated for motif_record */
    int dimensions[2];
    npy_intp dims[2];
    npy_intp dim[1];
    int bgmkv;
    int prob_option;
    int bgprob[4] = { -121897 ,  -158735 ,  -158735 ,  -121897 };
    int cond2g1[16] = { -110466 , -103555 , -123712 , -150393 , 
                        -176887 , -136752 , -156749 , -160550 , 
                        -144109 , -303335 , -136752 , -140393 , 
                        -134139 , -107272 , -140049 , -110466 };
    int cond3g2[64] = { -93492 , -134318 , -110188 , -114941 , 
                        -91580 , -103850 , -108575 , -110842 , 
                        -115123 , -149673 , -121660 , -130999 , 
                        -132085 , -169047 , -149085 , -154868 , 
                        -190443 , -157522 , -183796 , -175579 , 
                        -146660 , -137741 , -127711 , -134700 , 
                        -161014 , -142165 , -147708 , -161769 , 
                        -175577 , -142454 , -167459 , -160272 , 
                        -159051 , -127639 , -126014 , -162764 , 
                        -299958 , -293920 , -288751 , -329295 , 
                        -137105 , -127338 , -137741 , -137048 , 
                        -143697 , -123923 , -121029 , -164246 , 
                        -136786 , -137444 , -149166 , -115832 , 
                        -114167 , -107624 , -111537 , -98683 , 
                        -146944 , -136673 , -149957 , -128074 , 
                        -113112 , -125407 , -124022 , -93492 };

    if (! PyArg_ParseTuple( args, "O!O!O!O!i", &PyList_Type, &listObj, &PyArray_Type, &pssmArrayObj, &PyList_Type, &listBG, &PyInt_Type, &orderBG, &prob_option ) )
        return NULL;

    if ( pssmArrayObj == NULL )
        return NULL;        
 
    /* get the background frequencies from listBG */
    bgmkv = PyInt_AsLong(orderBG);
    //fprintf(stdout, "%d\n", bgmkv);    

    /* give new background probabilities */
    nbg = PyList_Size(listBG);
    if( nbg > 0 ){
        
        for( i = 0; i < 4; i++){ 
            strObj = PyList_GetItem(listBG, i);
            bgprob[i] = PyInt_AsLong( strObj );
        }
        
        for( i = 4; i < 20; i++){ 
            strObj = PyList_GetItem(listBG, i);
            cond2g1[i - 4] = PyInt_AsLong( strObj );
        }
        
        for( i = 20; i < 84; i++){ 
            strObj = PyList_GetItem(listBG, i);
            cond3g2[i - 20] = PyInt_AsLong( strObj );
        }

    }

    /* get the number of lines passed to us */
    numLines = PyList_Size(listObj);

    /* should raise an error here. */
    if ( numLines < 1 || numLines > MAXMOTIFS ){
        printf( tmpstr, "argument 1 must be a list of length > 1 and < %d", MAXMOTIFS );
        PyErr_SetString( PyExc_ValueError, tmpstr );
        return 0;    // trigger exception
    }

    /*pssmArrayObj = (PyArrayObject *) PyArray_ContiguousFromObject( pssmObj, PyArray_FLOAT, 0, 3 );*/
    if ( pssmArrayObj == NULL )
        return NULL;

    motifLength = pssmArrayObj->dimensions[0];
    motifBases  = pssmArrayObj->dimensions[1];

    if ( motifLength == 1 )
        return PyErr_Format( PyExc_RuntimeError, "Incorrect motif length (%d)", motifLength );
    if ( pssmArrayObj->dimensions[1] != 4 )
        //return PyErr_Format( PyExc_RuntimeError, "Incorrect PSSM dimension not 4 (%ld)", PyInt_AsLong(pssmArrayObj->dimensions[1]) );
        return PyErr_Format( PyExc_RuntimeError, "Incorrect PSSM dimension not 4 (%d)", motifBases);

    pssm_p = ( float ** )malloc( motifLength*sizeof( float * ) );
    for( i = 0; i < motifLength; i++ ){
        pssm_p[i] = ( float * )malloc( 4*sizeof( float ) );
        for( j = 0; j < 4; j++ ){
            pssm_p[i][j] = (float)( *( double * )( pssmArrayObj->data + i*pssmArrayObj->strides[0] + j*pssmArrayObj->strides[1] ));
            tmp += pssm_p[i][j]; 
            //fprintf( stdout, "%e\t", *( double * )( pssmArrayObj->data + i*pssmArrayObj->strides[0] + j*pssmArrayObj->strides[1] ));
        }
    }

    seqarray  = (int ** )malloc((numLines)*sizeof(int *));
    locations = (int ** )malloc((numLines)*sizeof(int *));
    //locations = (int ** )malloc((MAXMOTIFS)*sizeof(int *));
    //motif_record = ( struct motif_T * )malloc( MAXMOTIFS*sizeof( struct motif_T ));

    /* iterate over items of the list, grabbing strings, and parsing for numbers */
    for (i=0; i<numLines; i++){

        /* grab the string object from the next element of the list */
        strObj = PyList_GetItem(listObj, i); /* Can't fail */
 
        /* make it a string */
        line = PyString_AsString( strObj );

        for( numChars = 0; line[numChars]; numChars++ ){}
        seqarray[i]  = (int *)malloc((numChars)*sizeof(int));
        locations[i] = (int *)malloc((4)*sizeof(int));

        locations[i][0] = i; 
        locations[i][1] = 0;
        locations[i][2] = numChars;
        locations[i][3] = numChars;

        for( j=0; j<numChars; j++ ){

            switch ( line[j] ){
            case 'A': case 'a': bp2int = 0;
                break;
            case 'C': case 'c': bp2int = 1;
                break;
            case 'G': case 'g': bp2int = 2;
                break;
            case 'T': case 't': bp2int = 3;
                break;
            default: bp2int = 4;
                break;
            }
 
            seqarray[i][j] = bp2int;
            tmp += (double)bp2int;

        } 
    }

    // allocate and initialize motif_record 
    //motif_record = ( struct motif_T * )malloc( numLines*sizeof( struct motif_T ));
    n_motif_record = MAXMOTIFS+numChars;
    motif_record = ( struct motif_T * )malloc( (n_motif_record)*sizeof( struct motif_T ));
    if(!motif_record) {    
        printf("Memory request fails for motif_ score!\n");
        exit(1);
    }

    for( i = 0; i < n_motif_record; i++ ){
        motif_record[i].iseq   = -1;
        motif_record[i].istart = -1; 
        motif_record[i].iend   = -1;
        motif_record[i].orient = -1;        
        motif_record[i].lenseq = -1;
        motif_record[i].score  = 0.0;
    }
 
    motifscan_subseq( seqarray, locations, numLines, pssm_p, motifLength, bgmkv, motif_record, 0, bgprob, cond2g1, cond3g2, prob_option );
    dimensions[0] = numLines;

    for( numHits = 0, i = 0; i < MAXMOTIFS && (float)( motif_record[i].score ) > 0; numHits += 1, i++ ){}
    //dims[0] = numLines;
    dims[0] = numHits;
    dims[1] = 1;
 
    motifIdx    = (PyArrayObject *) PyArray_SimpleNew( 1, dims, PyArray_INT );
    motifStart  = (PyArrayObject *) PyArray_SimpleNew( 1, dims, PyArray_INT );
    motifEnd    = (PyArrayObject *) PyArray_SimpleNew( 1, dims, PyArray_INT );
    motifOrient = (PyArrayObject *) PyArray_SimpleNew( 1, dims, PyArray_INT );

    //dim[0] = numLines;
    dim[0] = numHits;
 
    //motifScore  = ( PyArrayObject * )PyArray_SimpleNew( 1, dim , PyArray_DOUBLE );
    motifScore  = ( PyArrayObject * )PyArray_SimpleNew( 1, dim , PyArray_FLOAT );

    //for( i = 0; i < numLines; i++ ){
    for( i = 0; i < numHits; i++ ){
        *( int * )( motifIdx->data     + i*motifIdx->strides[0] )    = (int)( motif_record[i].iseq );
        *( int * )( motifStart->data   + i*motifStart->strides[0] )  = (int)( motif_record[i].istart );
        *( int * )( motifEnd->data     + i*motifEnd->strides[0] )    = (int)( motif_record[i].iend );
        *( int * )( motifOrient->data  + i*motifOrient->strides[0] ) = (int)( motif_record[i].orient );
        *( float * )( motifScore->data + i*motifScore->strides[0] )  = (float)( motif_record[i].score );
    }

    for (i = 0; i < numLines; i++){
        free( seqarray[i] );
        free( locations[i] );
        //free( score[i] );
    }

    //free( score );
    free( seqarray  );
    free( locations );
    free( motif_record );

    for( i=0; i<motifLength; i++ ){
        free( pssm_p[i] );
    } 
    
    free( pssm_p );

    /*Py_DECREF( pssmArrayObj );*/
    Py_ret = ( PyObject * )Py_BuildValue( "OOOOO", motifIdx, motifStart, motifEnd, motifOrient, motifScore );
    Py_DECREF( motifIdx );
    Py_DECREF( motifStart );
    Py_DECREF( motifEnd );
    Py_DECREF( motifOrient );
    Py_DECREF( motifScore );

    return Py_ret; 

}

PyMethodDef myMethods[] = {
    { "seqscan", seqscan, METH_VARARGS, SEQSCAN_DOC },
    { NULL, NULL, 0, NULL }
};

void init_seq(void)
{       
    PyObject *m,*d;
    PyObject *tmp;
    m=Py_InitModule("_seq", myMethods);
    import_array();
    pyError = PyErr_NewException("seq.error", NULL, NULL);
    Py_INCREF(pyError);
    PyModule_AddObject(m, "error", pyError);

    d = PyModule_GetDict(m);

    tmp = PyInt_FromLong(MAX_OPTION);
    PyDict_SetItemString(d,"MAX_OPTION",tmp);
    Py_DECREF(tmp);

    tmp = PyInt_FromLong(MEAN_OPTION);
    PyDict_SetItemString(d,"MEAN_OPTION",tmp);
    Py_DECREF(tmp);

    tmp = PyInt_FromLong(CUTOFF_OPTION);
    PyDict_SetItemString(d,"CUTOFF_OPTION",tmp);
    Py_DECREF(tmp);

    tmp = PyInt_FromLong(MINSCORE);
    PyDict_SetItemString(d,"MINSCORE",tmp);
    Py_DECREF(tmp);

    tmp = PyInt_FromLong(MAXMOTIFS);
    PyDict_SetItemString(d,"MAXMOTIFS",tmp);
    Py_DECREF(tmp);

}
