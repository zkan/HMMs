#include "hmm.h"
#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdio.h>
#include <ctime>

int main( int argc, char **argv ) {
	FILE *fp;
	int T;
	HMM *hmm;
	int *O;					// observations
	double **alpha, **beta, **gamma;
	double logprobinit, logprobfinal;
	int niter;			// number of iterations

	int c;
	extern char *optarg;
	extern int optind, opterr, optopt;

//	cout << argc << endl;
//	cout << optind << endl;

	if( ( argc - optind ) != 2 ) {
		cout << "Usage error" << endl;
		cout << "Usage: " << argv[0] << " <model.hmm> <the sequence directory>" << endl;
		cout << "Note: This program needs at least one sequence to train." << endl;
		cout << "      The sequence file extension should be '.seq'" << endl;
		return EXIT_FAILURE;
	}
	
	fp = fopen( argv[optind], "r" );
	if( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found\n", argv[optind] );
		return EXIT_FAILURE;
	}
	
	hmm = new HMM();
	hmm->readHMM( fp );
	cout << "### Initial HMM ###" << endl;
	hmm->printHMM();
	if( fp ) {
		fclose( fp );
		fp = NULL;
	}
	
//	ostringstream os_seq_name;
//	os_seq_name << argv[optind + 1];
//	cout << os_seq_name.str() << endl;
	
	char *dir_name = argv[optind + 1];

	T = 50000;
	alpha = ( double** ) calloc( ( unsigned ) ( T + 1 ), sizeof( double* ) );
	for( int i = 1; i <= T; i++ ) {
		alpha[i] = ( double* ) calloc( ( unsigned ) ( hmm->getNumStates() + 1 ), sizeof( double ) );
	}
	beta = ( double** ) calloc( ( unsigned ) ( T + 1 ), sizeof( double* ) );
	for( int i = 1; i <= T; i++ ) {
		beta[i] = ( double* ) calloc( ( unsigned ) ( hmm->getNumStates() + 1 ), sizeof( double ) );
	}
	gamma = ( double** ) calloc( ( unsigned ) ( T + 1 ), sizeof( double* ) );
	for( int i = 1; i <= T; i++ ) {
		gamma[i] = ( double* ) calloc( ( unsigned ) ( hmm->getNumStates() + 1 ), sizeof( double ) );
	}
	
	timeval start, end;
  long mtime, seconds, useconds;
  gettimeofday(&start, NULL);

	hmm->train_bw( dir_name, T, O, alpha, beta, gamma, &niter, &logprobinit, &logprobfinal );
	
	gettimeofday(&end, NULL);
	seconds = end.tv_sec - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;
	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
	printf("Elapsed time: %ld milliseconds\n", mtime);
	
	fprintf( stdout, "Number of iterations: %d\n", niter );
	fprintf( stdout, "Log prob(O|init model): %E\n", logprobinit );
	fprintf( stdout, "Log prob(O|estimated model): %E\n", logprobfinal );

	double **tmp_A;		// temporary current transition probability
	double **tmp_B;		// temporary current emission probability
	double *tmp_pi;		// temporary current initial state distribution
	tmp_A = ( double** ) calloc( ( unsigned ) ( hmm->getNumStates() + 1 ), sizeof( double* ) );
	for( int i = 1; i <= hmm->getNumStates(); i++ ) {
		tmp_A[i] = ( double* ) calloc( ( unsigned ) ( hmm->getNumStates() + 1 ), sizeof( double ) );
	}
	tmp_B = ( double** ) calloc( ( unsigned ) ( hmm->getNumStates() + 1 ), sizeof( double* ) );
	for( int i = 1; i <= hmm->getNumStates(); i++ ) {
		tmp_B[i] = ( double* ) calloc( ( unsigned ) ( hmm->getNumSymbols() + 1 ), sizeof( double ) );
	}
	tmp_pi = ( double* ) calloc( ( unsigned ) ( hmm->getNumStates() + 1 ), sizeof( double ) );

	tmp_A = hmm->getA();
	tmp_B = hmm->getB();
	tmp_pi = hmm->getPi();

	cout << "### Checksum ###" << endl;
	float checkSum = 0;
	cout << "A: " << endl;
	for( int i = 1; i <= hmm->getNumStates(); i++ ) {
		cout << "sum of row " << i << ": ";
		for( int j = 1; j <= hmm->getNumStates(); j++ ) {
			checkSum += tmp_A[i][j];
		}
		cout << checkSum << endl;
		checkSum = 0;
	}
	cout << "B: " << endl;
	for( int i = 1; i <= hmm->getNumStates(); i++ ) {
		cout << "sum of row " << i << ": ";
		for( int j = 1; j <= hmm->getNumSymbols(); j++ ) {
			checkSum += tmp_B[i][j];
		}
		cout << checkSum << endl;
		checkSum = 0;
	}
	cout << "Pi: " << endl;
	cout << "sum of Pi: ";
	for( int i = 1; i <= hmm->getNumStates(); i++ ) {
		checkSum += tmp_pi[i];
	}
	cout << checkSum << endl;
	checkSum = 0;

	cout << "### Final HMM ###" << endl;
	hmm->printHMM();
	cout << "## Sufficient statistics ##" << endl;
	hmm->printSuffStats();

	// print an HMM to file
//*
	fp = fopen( argv[optind], "w" );
	if( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found\n", argv[optind] );
		exit( 1 );
	}
	hmm->printHMMWithSuffStats( fp );
	fclose( fp );
//*/
	return EXIT_SUCCESS;
}



