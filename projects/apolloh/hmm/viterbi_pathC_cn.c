#include "mex.h"
#include <string.h>
#include <math.h>
/* This function is a C implementation by Guillaume Alain
 * of viterbi_path_SS.m from Sohrab Shah, whose prototype is
 *  [path loglik seg] = viterbi_path_SS(prior, transmat, obslik).
 *
 * It supports non-stationary transition matrices by allowing the second
 * argument to be of size either (K,K) or (K,K,T-1).
 *
 * The size of the other arguments have to be
 *    prior : (K,1) or (1,K)
 *    obslik : (K,T)
 *
 * This function WORKS WITH LOGARITHMS. You need to take logs of all the
 * arguments that you'd usually put into Sohrab's viterbi_path_SS function.
 * It can usually be done usually easily by taking logs as the function is
 * called, but you have to be careful about takings logs of zero.
 * If you don't like having to call
 *      viterbi_path(log(prior+eps), log(transmat+eps), log(obslik))
 * you can always write your own wrapper function that does that. You might
 * want to be careful also about using 'eps' as the lowest possible number,
 * because it might not always be small enough compared to the values of
 * the log-likelihoods for the observations that can be big when they come
 * from a gaussian, for example.
 *
 * The decision to work with logarithms was because hmmmix-soft was using
 * them a lot and it would be stable without having to do renormalizing
 * trick all the time. The downside is that you can never had true zero
 * transition probabilities.
 *
 * Output arguments that are not specified in Matlab will not be computed.
 * This function does not compute the "Bayes segment factor" found in the
 * fourth column of the "seg" output argument. It just puts zeros in there.
 */

void copyVector(double *, double *, unsigned int);
void addVectors(double *, double *, double *, unsigned int);
void setVectorToValue(double *, double, unsigned int);
void setVectorToValue_int(int *, int, unsigned int);
void maxVectorInPlace(double *, int *, double *, unsigned int);
void preparePositionSpecificMatrix(double * transSlice, double * C, unsigned int, double * ZS, double, double, unsigned int);
double createCopyNumberMatrix(double * C, unsigned int, double * CNS, double, double, unsigned int);
void deepCopy(double * transSlice, double* transmat, unsigned int K);
void logMatrixInPlace(double * A, unsigned int K);
void outputMatrix(double * A, unsigned int K);
double distanceTransitionFunction(double, double, double);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double * prior, * transmat, * obslik, * transSlice, * C;
    int K, T, tmp, tmp1, transitionMatrixStrides;
    double * copyNumber, * CNS, * ZS, * posn, * txnExpLen;
    double * delta, *changes;
    int changesCounter;
    int * psi, * path;
    int t, k, j;
    double loglik, curC, rho;
    
    double * d; /* buffer */
    
    double *outputToolPtr;
    
    if ( !( (nrhs==8) && (nlhs==3) ) )
        mexErrMsgTxt("viterbi_path_MatlabC: requires 8 inputs and 3 outputs");
    
    prior=mxGetPr(prhs[0]);
    transmat=mxGetPr(prhs[1]);  /* base matrix - stationary values, KxK */
    obslik=mxGetPr(prhs[2]);
    copyNumber=mxGetPr(prhs[3]);/* copy number vector, Nx1 */
    CNS=mxGetPr(prhs[4]);/* vector of indices, 1x6 */
    ZS=mxGetPr(prhs[5]);/* 1x20 */
    posn=mxGetPr(prhs[6]); /* vector of posns, Nx1 */
    txnExpLen=mxGetPr(prhs[7]); /* integer */
    
    /* Check size of initial state distribution */
    K = mxGetN(prhs[0]); /* Number of columns; number of states */
    tmp = mxGetM(prhs[0]); /* Number of rows */
    /*mexPrintf("fwd_backC: Initial state vector is %d-by-%d\n",tmp,K);*/
    
    /* Check size of obslik */
    obslik = mxGetPr(prhs[2]);
    T = mxGetN(prhs[2]); /* Number of columns in obslik; # data points */
    tmp = mxGetM(prhs[2]); /* Number of rows in obslik; should be K */
    if (tmp != K)
        mexErrMsgTxt("fwd_backC: The obslik must have K rows.");
    /*mexPrintf("fwd_backC: Obslik is %d-by-%d\n",tmp,T);*/
    
    tmp = mxGetM(prhs[3]); /* Number of rows in copynumber vector */  
    tmp1 = mxGetN(prhs[3]);
    if (tmp != T)
        mexErrMsgTxt("fwd_backC: The copy number vector must be length T.");
    /*mexPrintf("fwd_backC: Copy number vector is %d-by-%d\n",tmp,tmp1);*/
    
    tmp = mxGetM(prhs[1]);
    tmp1 = mxGetN(prhs[1]);
    if (tmp != K || tmp1 != K)
        mexErrMsgTxt("fwd_backC: The transition matrix must be of size KxK.");
    /*mexPrintf("Base txn matrix is %d-by-%d\n",tmp,tmp1);*/
    
    tmp = mxGetM(prhs[6]);
    tmp1 = mxGetN(prhs[6]);
    if (tmp != T)
        mexErrMsgTxt("fwd_backC: The positions vector must be of size Tx1.");
    /*mexPrintf("Positions vector is %d-by-%d\n",tmp,tmp1);*/
    
    
    delta = mxMalloc(K*T*sizeof(double));
    psi = mxMalloc(K*T*sizeof(int));
    path = mxMalloc(T*sizeof(int));
    
    t = 0;
    addVectors(delta + t*K, prior, obslik + t*K, K);
    setVectorToValue_int(psi + t*K, 0, K);
    
    d = mxMalloc(K*sizeof(double));
    transSlice = (double *)mxMalloc(K*K*sizeof(double));
    C = mxMalloc(K*K*sizeof(double));
    
    /* forward */
    for(t=1;t<T;++t) { /* position */
        /* Each iteration, we overwrite transSlice with transmat
         * to start over when at a new probe */
        deepCopy(transSlice, transmat, K);
        deepCopy(C, transmat, K);
        /* modify transSlice inplace by adding position-specific probs */
        rho = 1.0 - distanceTransitionFunction(posn[t-1],posn[t],txnExpLen[0]);    
        curC = createCopyNumberMatrix(C,K,CNS,copyNumber[t-1],copyNumber[t],0);        
        preparePositionSpecificMatrix(transSlice, C, K, ZS, rho, curC, 0);   
        logMatrixInPlace(transSlice, K);
        /* mexPrintf("viterbi_path_MatlabC: Normalized Log2 Matrix: \n"); */
        /* outputMatrix(transSlice, K); */
        for(j=0;j<K;++j) { /* column */
            /* addVectors(d, delta + (t-1)*K, transmat + j*K + transitionMatrixStrides*(t-1), K); */
            addVectors(d, delta + (t-1)*K, transSlice + j*K , K);
            maxVectorInPlace(delta + j + t*K, psi + j + t*K, d, K);
            delta[j+t*K] += obslik[j+t*K];
        }
   
    }
    
    
    /* backward */
    t = T-1;
    maxVectorInPlace(d, path + t, delta + t*K, K); /* using the first value of d to store junk */
    loglik = d[0];
    
    for(t=T-2;t>=0;--t) {
        path[t] = psi[path[t+1] + (t+1)*K];
        /*mexPrintf("Setting path[%d]=%d\n", t, path[t]); */
    }
    
    changes = mxMalloc(4*T*sizeof(double));
    changesCounter = 0;
    changes[changesCounter + 0*T] = 0;
    changes[changesCounter + 1*T] = 0; /* overwritten */
    changes[changesCounter + 2*T] = path[0];
    changes[changesCounter + 3*T] = 0;
    changesCounter = 1;
    
    for(t=1;t<T;++t) {
        if (path[t] != path[t-1]) {
            changes[changesCounter + 0*T] = t;
            changes[(changesCounter-1) + 1*T] = t-1;
            changes[changesCounter + 2*T] = path[t];
            changes[changesCounter + 3*T] = 0; /* that computeSegmentBayesFactor */
            changesCounter++;
        }
    }
    changes[(changesCounter-1) + 1*T] = T-1;
    
    plhs[0] = mxCreateDoubleMatrix(1, T, mxREAL);
    outputToolPtr = mxGetPr(plhs[0]);
    /* Be careful to add +1 to path values. This is because C starts from 0
     * and Matlab starts from 1. I figured it would be easier to think in C
     * anb convert only at the end. Hence, the path[t] + 1.
     */
    for(t=0;t<T;++t)
        outputToolPtr[t] = (double)(path[t]+1);
    
    /* A junk value for the loglik for backward compatibility.
     * We're not scaling so we don't get a loglik value.
     */
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    outputToolPtr = mxGetPr(plhs[1]);
    outputToolPtr[0] = loglik;
    
    plhs[2] = mxCreateDoubleMatrix(changesCounter, 4, mxREAL);
    outputToolPtr = mxGetPr(plhs[2]);
    for(t=0;t<changesCounter; ++t) {
        /* +1 because of the Matlab offset from C */
        outputToolPtr[t + 0*changesCounter] = changes[t + 0*T] + 1;
        outputToolPtr[t + 1*changesCounter] = changes[t + 1*T] + 1;
        outputToolPtr[t + 2*changesCounter] = changes[t + 2*T] + 1;
        outputToolPtr[t + 3*changesCounter] = changes[t + 3*T];
    }
    
    
    mxFree(delta); mxFree(psi); mxFree(path);
    mxFree(d);
    mxFree(changes);
    
    return;
}



/* Method to assign transSlice as the position specific matrix that includes cnv freq prior
 preparePositionSpecificMatrix(double * transSlice, double * C, unsigned int, double * ZS, double, unsigned int) */
void preparePositionSpecificMatrix(double * transSlice, double * C, unsigned int K, double * ZS, double rho, double curC, unsigned int boolTest) {
    unsigned int i, j;
    double sum;
  
    /* Add the distance to our output matrix, in place 
     * Also multiple by copy number values */       
    for (i = 0; i < K; i++){ /* rows */
        for (j = 0; j < K; j++){ /* columns */           
            if (ZS[i]==ZS[j]){             
                transSlice[i + j*K] = rho * C[i+j*K];                   
            }else{
                transSlice[i + j*K] = ((1.0-rho)/curC) * C[i+j*K];
            }            
        }
    }
       
  
    if (boolTest) {
        mexPrintf("fwd_backC: Rho:\t%f\n",rho);
        mexPrintf("fwd_backC: Raw Matrix: \n");        
        outputMatrix(transSlice, K);
    }
    
    /* Normalize matrix by rows */        
    for (i=0;i<K;i++) { /* rows */
        sum = 0;
        for (j=0;j<K;j++) { /* columns */
            sum += transSlice[i+j*K];
        }        
        if (sum > 0){
            for (j=0;j<K;j++) { /* columns */            
                transSlice[i+j*K] /= sum;
            }
        }
        /*mexPrintf("Sum: %f\n",sum);*/   
    }
    
    if (boolTest) {
        mexPrintf("fwd_backC: Normalized Matrix: \n");
        outputMatrix(transSlice, K);
    }
    
}

/* Create a copy number matrix by modifying C in place
 *
 */
double createCopyNumberMatrix(double * C, unsigned int K, double * CNS, double prevCN, double curCN, unsigned int boolTest) {
    unsigned int i, j;
    double curC;
    /*mexPrintf("prevCN=%d\tcurCN=%d\n",(int)prevCN,(int)curCN);*/
    for (i = CNS[(int)prevCN-1]-1; i < CNS[(int)prevCN]-1; i++){ /* columns */
        for (j = CNS[(int)curCN-1]-1; j < CNS[(int)curCN]-1; j++){  /* rows */              
            C[i+j*K] = 1.0;
        }
    }      
    if (boolTest) {
        mexPrintf("fwd_backC: Copy number Matrix: \n");
        outputMatrix(C, K);
    }
    curC = CNS[(int)curCN] - CNS[(int)curCN-1] + 1;   
    return curC;
}

/* Position specific distance used in transition matrix
 * returns double
 */
double distanceTransitionFunction(double prevPosn, double curPosn, double L) {
    double distance = 0;
    double rho = 0;
    distance = curPosn - prevPosn + 1.0; /* won't encounter next chr */
    rho = (1.0/2.0)*(1.0-exp(-distance/(2.0*L)));
    return rho;
}


/* Method to make a deep copy of a K-by-K matrix */
void deepCopy(double * transSlice, double* transmat, unsigned int K) {
    unsigned int i, j;
    
    /* deep copy transmat to transSlice */
    for (i=0;i<K;i++) /* rows */
        for (j=0;j<K;j++) /* columns */
            transSlice[i + j*K] = transmat[i + j*K];
}


/* logs each element in a K-by-K matrix, A */
void logMatrixInPlace(double * A, unsigned int K) {
    unsigned int i, j;
    for (i=0;i<K;i++) /* rows */
        for (j=0;j<K;j++) /* columns */
            A[i + j*K] = log(A[i + j*K]);
}

/* Output matrix for debugging purposes */
void outputMatrix(double * A, unsigned int K) {
    unsigned int i, j;
    
    for (i=0;i<K;i++) { /*rows*/
        for (j=0;j<K;j++) { /*cols*/
            mexPrintf("%f\t", A[K*j+i]);
        }
        mexPrintf("\n");
    }
}

void multiplyMatrixInPlace(double * result, double * trans, double * v, unsigned int K) {
    
    unsigned int i, d;
    
    for(d=0;d<K;++d) {
        result[d] = 0;
        for (i=0;i<K;++i){
            result[d] += trans[d + i*K] * v[i];
        }
    }
    return;
}

void transposeSquareInPlace(double * out, double * in, unsigned int K) {
    
    unsigned int i, j;
    
    for(i=0;i<K;++i){
        for(j=0;j<K;++j){
            out[j+i*K] = in[i+j*K];
        }
    }
    return;
}

void outerProductUVInPlace(double * Out, double * u, double * v, unsigned int K) {
    unsigned int i, j;
    
    for(i=0;i<K;++i){
        for(j=0;j<K;++j){
            Out[i + j*K] = u[i] * v[j];
        }
    }
    return;
}



void copyVector(double * Out, double * In, unsigned int L) {
    unsigned int i;
    
    for(i=0;i<L;++i)
        Out[i] = In[i];
    
    return;
}

void addVectors(double * Out, double * u, double * v, unsigned int L) {
    unsigned int i;
    
    for(i=0;i<L;++i)
        Out[i] = u[i] + v[i];
    
    return;
    
}

void setVectorToValue(double * A, double value, unsigned int L) {
    unsigned int i;
    
    for(i=0;i<L;++i)
        A[i] = value;
    
    return;
}

void setVectorToValue_int(int * A, int value, unsigned int L) {
    unsigned int i;
    
    for(i=0;i<L;++i)
        A[i] = value;
    
    return;
}


void maxVectorInPlace(double * Out_value, int * Out_index, double * A, unsigned int L) {
    unsigned int i;
    double maxvalue;
    int index;
    
    maxvalue = A[0];
    index = 0;
    
    for(i=1;i<L;++i) {
        if (maxvalue < A[i]) {
            index = i;
            maxvalue = A[i];
        }
    }
    
    *Out_value = maxvalue;
    *Out_index = index;
    
    return;
}
