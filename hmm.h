/*
 *	Author: Kan Ouivirach
 *	File: hmm.h
 *	Purpose: Define HMM class
 *	Credit: Tapas Kanungo, kanungo@cfar.umd.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <iostream>
#include <string.h>
#include <vector>
#include <string>
#include <sstream>

#include "matrix.h"

#ifndef _NO_NAMESPACE
using namespace std;
using namespace math;
#define STD std
#else
#define STD
#endif

#ifndef _NO_TEMPLATE
typedef matrix<double> Matrix;
#else
typedef matrix Matrix;
#endif

#ifndef _NO_EXCEPTION
#  define TRYBEGIN()	try {
#  define CATCHERROR()	} catch (const STD::exception& e) { \
						cerr << "Error: " << e.what() << endl; }
#else
#  define TRYBEGIN()
#  define CATCHERROR()
#endif

class HMM {
private:
	// properties
	int N;				// number of states; Q = { 1, 2, ..., N }
	int M;				// number of observation symbols; V = { 1, 2, ..., M }
	double **A;		// transition probability 
	double **B;		// emission probability 
	double *pi;		// initial state distribution 

	// sufficient statistics for the HMM
	double **numeratorA_t;
	double *denominatorA_t;
	double **numeratorB_t;
	double *denominatorB_t;
	double *pi_t;
	double logprobfinal;	// final log-likelihood after batch training
	int NS;								// the number of sequences used in batch training

	int dim;					// number of dimensions (continuous HMMs)
	double **meanMat;	// mean matrix
	double ***covMat;	// covariance matrix

	double loglik;		// log likelihood of the HMM
	
	// methods
	void computeGamma( int T, double **alpha, double **beta, double **gamma );
	void computeXi( int T, int *O, double **alpha, double **beta, double ***xi );
	
	void computeGamma( int T, double **alpha, double **gamma );			// only used in the IBW algorithm
	void computeXi( int T, int *O, double **alpha, double ***xi );	// only used in the IBW algorithm
	
	void computeContXi( int T, double **O, double **alpha, double **beta, double ***xi );

	double computeGaussian( int state, double *X );

	int genInitialState();
	int genNextState( int q_t );
	int genSymbol( int q_t );

public:
	// methods
	HMM() {
		this->N      = 0;
		this->M      = 0;
		this->dim    = 0;
		this->loglik = 0;
	};
	~HMM() {};

	// methods to return HMM's properties
	int getNumStates()  { return this->N; }
	int getNumSymbols() { return this->M; }
	int getNumDim()     { return this->dim; }
	double getLoglik()  { return this->loglik; }
	double *getPi()     { return this->pi; }
	double **getA()     { return this->A; }
	double **getB()     { return this->B; }

	double *getPi_t()           { return this->pi_t; }
	double **getNumeratorA_t()  { return this->numeratorA_t; }
	double *getDenominatorA_t() { return this->denominatorA_t; }
	double **getNumeratorB_t()  { return this->numeratorB_t; }
	double *getDenominatorB_t() { return this->denominatorB_t; }
	double getLogProbFinal()    { return this->logprobfinal; }
	int getNS()                 { return this->NS; }

	// methods to set HMM's properties
	void setPi( double *pi );
	void setA( double **A );
	void setB( double **B );

	void readHMM( FILE *fp );
	void readContHMM( FILE *fp );
	void readHMMWithSuffStats( FILE *fp );
	
	void printHMM();
	void printHMM( char *filename );
	void printHMM( FILE *fp );
	void printHMMWithSuffStats( FILE *fp );
	void printContHMM();
	void initHMM( int seed, int N, int M, int type );	// randomly initialize a HMM ( type 0: normal, type 1: left-to-right )
	void printSuffStats();

	void readSequence( FILE *fp, int *pT, int **pO );
	void readContSequence( FILE *fp, int *pT, double ***pO );
	void printSequence( FILE *fp, int T, int *O );
	void printContSequence( FILE *fp, int T, double **O );
	void genSequenceArray( int seed, int T, int *O, int *p );	// randomly generate a sequence given a HMM

	// allocate Xi
	double ***allocXi( int T, int N );

	// forward algorithm
	void forward( int T, int *O, double **alpha );
	void forwardWithScale( int T, int *O, double **alpha, double *scale, int type );	// run the forward algorithm given a HMM ( type 0: normal, type 1: left-to-right )
	void forwardWithScale( int T, vector<int> O, double **alpha, double *scale );
	void forwardContWithScale( int T, double **O, double **alpha, double *scale );

	// backward algorithm
	void backward( int T, int *O, double **beta );
	void backwardWithScale( int T, int *O, double **beta, double *scale, int type );	// run the backward algorithm given a HMM ( type 0: normal, type 1: left-to-right )
	void backwardContWithScale( int T, double **O, double **beta, double *scale );

	// Viterbi algorithm
	void viterbi( int T, int *O );
	
	// conventional Baum-Welch algorithm
	void baumwelch( int T, int *O, double **alpha, double **beta, double **gamma, int *pniter, double *plogprobinit, double *plogprobfinal );
	void train_bw( char *dir_name, int T, int *O, double **alpha, double **beta, double **gamma, int *pniter, double *plogprobinit, double *plogprobfinal );
	void train_bw_cont( int NS, int T, double **O, double **alpha, double **beta, double **gamma, int *pniter, double *plogprobinit, double *plogprobfinal );

	// Baum-Welch algorithm with the sufficient statistics (incremental maximum-likelihood algorithm)
	
	void train_iml( int T, int *O, double **alpha, double **beta, double **gamma, int niter, double *plogprobinit, double *plogprobfinal, 
				double *pi_t, double **numeratorA_t, double *denominatorA_t, double **numeratorB_t, double *denominatorB_t, 
				int NS, ostringstream &os_seq_name, double forgetting_factor, int actual_obs_length, int num_of_seqs_trained, double logprob_trained );

	// Baum-Welch algorithm with the sufficient statistics (incremental Baum-Welch algorithm)
	void train_ibw( int T, int *O, double **alpha, double **beta, double **gamma, int *pniter, double *plogprobinit, double *plogprobfinal );
};


