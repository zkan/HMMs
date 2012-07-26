#include "hmm.h"

int main( int argc, char **argv ) {
	FILE *fp;
	int T;
	HMM *hmm;
	int *O;		/* observations */
	double **alpha, **beta, **gamma;
	double logprobinit, logprobfinal;
	int niter;	/* number of iterations */

	if( argc != 3 ) {
		cout << "Usage error" << endl;
		cout << "Usage: xxx <model.hmm> <sequence.seq>" << endl;
		exit( 1 );
	}
	
	fp = fopen( argv[1], "r" );

	if( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found\n", argv[1] );
		exit( 1 );
	}

	hmm = new HMM();
	hmm->readHMM( fp );
	hmm->printHMM();
	fclose( fp );

	fp = fopen( argv[2], "r" );

	if( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found\n", argv[1] );
		exit( 1 );
	}
	
	hmm->readSequence( fp, &T, &O );

	alpha = (double **) calloc( (unsigned) ( T + 1 ), sizeof(double*) );
	for( int i = 1; i <= T; i++ ) {
		alpha[i] = (double *) calloc( (unsigned) ( hmm->getNumStates() + 1 ), sizeof(double) );
	}
	beta = (double **) calloc( (unsigned) ( T + 1 ), sizeof(double*) );
	for( int i = 1; i <= T; i++ ) {
		beta[i] = (double *) calloc( (unsigned) ( hmm->getNumStates() + 1 ), sizeof(double) );
	}
	gamma = (double **) calloc( (unsigned) ( T + 1 ), sizeof(double*) );
	for( int i = 1; i <= T; i++ ) {
		gamma[i] = (double *) calloc( (unsigned) ( hmm->getNumStates() + 1 ), sizeof(double) );
	}
	
	hmm->baumwelch( T, O, alpha, beta, gamma, &niter, &logprobinit, &logprobfinal );
	fprintf( stdout, "Number of iterations: %d\n", niter );
	fprintf( stdout, "Log prob(O|init model): %E\n", logprobinit );
	fprintf( stdout, "Log prob(O|estimated model): %E\n", logprobfinal );

	hmm->printHMM();

	return 1;
}



