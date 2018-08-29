#include <math.h>
#include "mex.h"
#include <string.h>

void deepCopy(double * transSlice, double* transmat, unsigned int, unsigned int);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double * cn, * CNS, * C; /* inputs */
    double * cy, * outputToolPtr; /* output */
    unsigned int tmp, tmp1, T, i, j, k, t, K;

    if (nrhs!=3 || !(nlhs==1))
        mexErrMsgTxt("fwd_backC: requires 3 inputs and 1 outputs");

    C = mxGetPr(prhs[0]);
    K = mxGetM(prhs[0]);
    T = mxGetN(prhs[0]);
    /*mexPrintf("C matrix is %d-by-%d\n",K,T);*/
    
    cn = mxGetPr(prhs[1]);
    tmp1 = mxGetN(prhs[1]); /* Number of columns; number of positions */
    tmp = mxGetM(prhs[1]); /* Number of rows */
    /*mexPrintf("CN vector is %d-by-%d\n",tmp,tmp1);*/
    
    CNS = mxGetPr(prhs[2]);
    tmp1 = mxGetN(prhs[2]); /* Number of columns; number of states */
    tmp = mxGetM(prhs[2]); /* Number of rows */
    /*mexPrintf("CNS vector is %d-by-%d\n",tmp,tmp1);*/
    
    cy = mxMalloc(K*T*sizeof(double));
    
    deepCopy(cy,C,K,T);
    
    for(t=0; t<T; t++){
        for(i = CNS[(int)cn[t]-1]-1; i < CNS[(int)cn[t]]-1; i++){  
            /*mexPrintf("T=%f\tt=%d\ti=%d\t%f:%f\n",T,t,i,CNS[(int)cn[t]-1]-1,CNS[(int)cn[t]]-1);*/
            cy[i+t*K] = 1.0;       
        }
    }
  
    plhs[0] = mxCreateDoubleMatrix(K, T, mxREAL);
    outputToolPtr = mxGetPr(plhs[0]);
    memcpy(outputToolPtr, cy, K*T*sizeof(double));
    
    mxFree(cy);
    return;
}


/*
 * Create a deep copy of a matrix
 */
void deepCopy(double * transSlice, double* transmat, unsigned int K, unsigned int T) {
    unsigned int i, j;
    
    /* deep copy transmat to transSlice */
    for (i=0;i<K;i++) /* rows */
        for (j=0;j<T;j++) /* columns */
            transSlice[i + j*K] = transmat[i + j*K];
}
