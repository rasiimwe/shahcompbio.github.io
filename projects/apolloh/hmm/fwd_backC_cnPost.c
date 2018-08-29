#include <math.h>
#include "mex.h"
#include <string.h>

double normalizeInPlace(double *, unsigned int, unsigned int);
void multiplyInPlace(double *, double *, double *, unsigned int);
void multiplyMatrixInPlace(double *, double *, double *, unsigned int);
void transposeSquareInPlace(double *, double *, unsigned int);
void outerProductUVInPlace(double *, double *, double *, unsigned int);
void componentVectorMultiplyInPlace(double *, double *, double *, unsigned int);
void preparePositionSpecificMatrix(double * transSlice, double * C, unsigned int, double * ZS, double, unsigned int, unsigned int);
void createCopyNumberMatrix(double * C, unsigned int, double * CNS, double * CN, unsigned int);
void deepCopy(double * transSlice, double* transmat, unsigned int K);
void outputMatrix(double * A, unsigned int K);
double distanceTransitionFunction(double, double, double);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double * init_state_distrib, * transmat, * obslik;
    double * copyNumber, * CNS, * ZS, * posn, * txnExpLen;
    int K, T, tmp, tmp1, cnStates; 
    
    /* the tranposed version of transmat*/
    double * transmatT, * transSlice, * C;
    
    double * scale, * alpha, * beta, * gamma;
    int t, d, i, j;
    double loglik = 0.0;
    double rho = 0.0;
  
    double *m, *b; 
    double *outputToolPtr;
    
    if (nrhs!=8 || !(nlhs==4))
        mexErrMsgTxt("fwd_backC: requires 8 inputs and 4 outputs");
    
    init_state_distrib=mxGetPr(prhs[0]);
    transmat=mxGetPr(prhs[1]); /* base matrix - stationary values, KxK */
    obslik=mxGetPr(prhs[2]); /* obslik, KxT*/
    copyNumber=mxGetPr(prhs[3]);/* copy number vector, 6xT */
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
    
    cnStates = mxGetM(prhs[3]); /* Number of rows in copynumber vector */  
    tmp1 = mxGetN(prhs[3]); /* Number of columns in copynumber vector */
    if (tmp1 != T)
        mexErrMsgTxt("fwd_backC: The copy number posteriors must have T rows.");
    /*if (tmp != 6)
        mexErrMsgTxt("fwd_backC: The copy number posteriors must have 6 columns.");
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
    
    scale = mxMalloc(T*sizeof(double));    
    alpha = mxMalloc(K*T*sizeof(double));
    beta = mxMalloc(K*T*sizeof(double));
    gamma = mxMalloc(K*T*sizeof(double));
    
    /* Use transSlice as the temporary txn matrix to modify
     * with position-specific; each iteration, we overwrite it with transmat
     * to start over at a new probe */
    transSlice = mxMalloc(K*K*sizeof(double));    
    transmatT = mxMalloc(K*K*sizeof(double));
    C = mxMalloc(K*K*sizeof(double));
    
    /********* Forward. ********/
   t = 0;
   multiplyInPlace(alpha + t*K, init_state_distrib, obslik + t*K, K);
   /*mexPrintf("Normalize alpha at t=%d\n",t);*/
   scale[t] = normalizeInPlace(alpha + t*K, K, t);
    
   m = mxMalloc(K*sizeof(double));
    
    for(t=1;t<T;++t){
        /* Each iteration, we overwrite transSlice with transmat
         * to start over when at a new probe */
        deepCopy(transSlice, transmat, K);
        deepCopy(C, transmat, K);
        /* modify transSlice inplace by adding position-specific probs */
        rho = 1.0 - distanceTransitionFunction(posn[t-1],posn[t],txnExpLen[0]);       
        createCopyNumberMatrix(C,K,CNS,copyNumber+t*cnStates,0); 
        preparePositionSpecificMatrix(transSlice, C, K, ZS, rho, t, 0);   
        transposeSquareInPlace(transmatT, transSlice, K);
        multiplyMatrixInPlace(m, transmatT, alpha + (t-1)*K, K);
        multiplyInPlace(alpha + t*K, m, obslik + t*K, K); 
        scale[t] = normalizeInPlace(alpha + t*K, K, t); 
       
    }
    
    loglik = 0;
    for(t=0;t<T;++t)
        loglik += log(scale[t]);
    
    /********* Backward. ********/
    
    
    t = T-1;
    /* I don't think we need to initialize beta to all zeros. */
    for(d=0;d<K;++d) {
        beta[d + t*K] = 1;
        gamma[d + t*K] = alpha[d + t*K];
    }
    
    b = mxMalloc(K*sizeof(double));/*mxCreateDoubleMatrix(K,1,mxREAL);*/
    /*    eta = mxMalloc(K*K*T*sizeof(double));
     * squareSpace = mxMalloc(K*K*sizeof(double));
     */
    /* Put the last slice of eta as zeros, to be compatible with Sohrab and Gavin's code.
     * There are no values to put there anyways. This means that you can't normalise the
     * last matrix in eta, but it shouldn't be used. Note the d<K*K range.
     */
    /*
     * for(d=0;d<(K*K);++d) {
     * mexPrintf("setting *(eta + %d) = 0 \n", d+t*K*K);
     *(eta + d + t*K*K) = 0;(double)7.0f;
     * }
     */
    /* We have to remember that the 1:T range in Matlab is 0:(T-1) in C. */
    for(t=(T-2);t>=0;--t) {        
        /* setting beta */
        multiplyInPlace(b, beta + (t+1)*K, obslik + (t+1)*K, K);
        /* Using "m" again instead of defining a new temporary variable.
         * We using a lot of lines to say
         * beta(:,t) = normalize(transmat * b);
         */
        deepCopy(transSlice, transmat, K);
        deepCopy(C, transmat, K);
        /* modify transSlice inplace by adding position-specific probs */
        rho = 1.0 - distanceTransitionFunction(posn[t],posn[t+1],txnExpLen[0]);
        createCopyNumberMatrix(C,K,CNS,copyNumber+(t+1)*cnStates,0);
        preparePositionSpecificMatrix(transSlice, C, K, ZS, rho,t,0);        
        multiplyMatrixInPlace(m, transSlice, b, K);
        normalizeInPlace(m, K, t);
        for(d=0;d<K;++d) { beta[d + t*K] = m[d]; }
        /* using "m" again as valueholder */
        
        /* setting eta, whether we want it or not in the output */
        /*	    outerProductUVInPlace(squareSpace, alpha + t*K, b, K);
         * componentVectorMultiplyInPlace(eta + t*K*K, transmat + K*K*t, squareSpace, K*K);
         * normalizeInPlace(eta + t*K*K, K*K);
         */
        
        /* setting gamma */        
        multiplyInPlace(m, alpha + t*K, beta + t*K, K);
        normalizeInPlace(m, K, t);
        for(d=0;d<K;++d) { gamma[d + t*K] = m[d]; } 
       
    }
    
    
    
    plhs[0] = mxCreateDoubleMatrix(K, T, mxREAL);
    outputToolPtr = mxGetPr(plhs[0]);
    memcpy(outputToolPtr, gamma, K*T*sizeof(double));
    
    plhs[1] = mxCreateDoubleMatrix(K, T, mxREAL);
    outputToolPtr = mxGetPr(plhs[1]);
    memcpy(outputToolPtr, alpha, K*T*sizeof(double));
    
    plhs[2] = mxCreateDoubleMatrix(K, T, mxREAL);
    outputToolPtr = mxGetPr(plhs[2]);
    memcpy(outputToolPtr, beta, K*T*sizeof(double));
    
    /* This handles the two possible cases for outputs, based on the number of outputs.
     * It's either
     * [gamma,alpha,beta,eta,loglik]
     * or
     * [gamma,alpha,beta,loglik].
     */
    if (nlhs == 4) {
        plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
        outputToolPtr = mxGetPr(plhs[3]);
        outputToolPtr[0] = loglik;
        /*  } else if (nlhs == 5) {
         * eta_dims[0] = K;
         * eta_dims[1] = K;
         * eta_dims[2] = T;
         * plhs[3] = mxCreateNumericArray(eta_ndim, eta_dims, mxDOUBLE_CLASS, mxREAL);
         * outputToolPtr = mxGetPr(plhs[3]);
         * memcpy(outputToolPtr, eta, K*K*T*sizeof(double));
         *
         *
         * plhs[4] = mxCreateDoubleMatrix(1,1,mxREAL);
         * outputToolPtr = mxGetPr(plhs[4]);
         * outputToolPtr[0] = loglik;
         */
    }
    
    
    mxFree(b); mxFree(m); /*mxFree(squareSpace);*/
    mxFree(scale); mxFree(transmatT);
    mxFree(alpha); mxFree(beta); mxFree(gamma); /*mxFree(eta); mxFree(eta_dims);*/
    mxFree(C);
    return;
}

/* Method to assign transSlice as the position specific matrix that includes cnv freq prior
 preparePositionSpecificMatrix(double * transSlice, double * C, unsigned int, double * ZS, double, unsigned int) */
void preparePositionSpecificMatrix(double * transSlice, double * C, unsigned int K, double * ZS, double rho, unsigned int col, unsigned int boolTest) {
    unsigned int i, j;
    double sum;
  
    /* Add the distance to our output matrix, in place 
     * Also multiple by copy number values */       
    for (i = 0; i < K; i++){ /* rows */
        for (j = 0; j < K; j++){ /* columns */           
            if (ZS[i]==ZS[j]){             
                transSlice[i + j*K] = rho * C[i+j*K];                   
            }else{
                transSlice[i + j*K] = ((1.0-rho)/(double)K)*2.0 * C[i+j*K];
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
        if (sum == 0){
            mexPrintf("fwd_backC: Row:%d\tPosn:%d\n",i,col);
            mexErrMsgTxt("Normalizing genotype/copy number matrix - column containing only zeros.");
        }else{
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
 * CN is a vector of 5 values!
 */
void createCopyNumberMatrix(double * C, unsigned int K, double * CNS, double * CN, unsigned int boolTest) {
    unsigned int i, j;
    
    for (i=0; i<K; i++){ /* columns */
        for (j=0; j<K; j++){ /* rows */
            C[i+j*K] = CN[(int)CNS[j]-1];
        }
    }
          
    if (boolTest) {
        mexPrintf("CN[0]=%f\tCN[1]=%f\tCN[2]=%f\tCN[3]=%f\tCN[4]=%f\tCN[5]=%f\n",CN[0],CN[1],CN[2],CN[3],CN[4],CN[5]);
        mexPrintf("fwd_backC: Copy number Matrix: \n");
        outputMatrix(C, K);
    }

}

/*
 * Create a deep copy of a matrix
 */
void deepCopy(double * transSlice, double* transmat, unsigned int K) {
    unsigned int i, j;
    
    /* deep copy transmat to transSlice */
    for (i=0;i<K;i++) /* rows */
        for (j=0;j<K;j++) /* columns */
            transSlice[i + j*K] = transmat[i + j*K];
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

/* And returns the normalization constant used.
 * I'm assuming that all I'll want to do is to normalize columns
 * so I don't need to include a stride variable.
 */
double normalizeInPlace(double * A, unsigned int N, unsigned int col) {
    unsigned int n;
    double sum = 0;
    
    for(n=0;n<N;++n) {
        sum += A[n];
        if (A[n] < 0) {
            mexErrMsgTxt("We don't want to normalize if A contains a negative value. This is a logical error.");
        }
    }
    
    if (sum == 0){
        mexPrintf("Failed at position: %d\n",col);
        mexErrMsgTxt("We are asked to normalize a section of a vector containing only zeros.");
    }else {
        for(n=0;n<N;++n)
            A[n] /= sum;
    }
    return sum;
}

void multiplyInPlace(double * result, double * u, double * v, unsigned int K) {
    unsigned int n;
    
    for(n=0;n<K;++n)
        result[n] = u[n] * v[n];
    
    return;
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

/* this works for matrices also if you just set the length "L" to be the right value,
 * often K*K, instead of just K in the case of vectors
 */
void componentVectorMultiplyInPlace(double * Out, double * u, double * v, unsigned int L) {
    unsigned int i;
    
    for(i=0;i<L;++i)
        Out[i] = u[i] * v[i];
    
    return;
}
