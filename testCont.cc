#include "hmm.h"

int main( int argc, char **argv ) {
	FILE *fp;
	int T;
	int NS = 1000;
	HMM *hmm;
	double **O;		/* observations */
	double **alpha, **beta, **gamma, *scale;
	double logprobinit, logprobfinal;
	int niter;	/* number of iterations */

	int c;
	extern char *optarg;
	extern int optind, opterr, optopt;
//*
	while( ( c = getopt( argc, argv, "N:" ) ) != EOF ) {
		switch( c ) {
			case 'N':
				sscanf( optarg, "%d", &NS );
				break;
		}
	}

	if( ( argc - optind ) != 2 ) {
		cout << "Usage error" << endl;
		cout << "Usage: xxx -N <number of sequences> <model.hmm> <the first sequence: name.#.seq>" << endl;
		exit( 1 );
	}
//*/
//	fp = fopen( argv[1], "r" );	
//	fp = fopen( "model/x.hmm", "r" );
//	fp = fopen( "model/cont.hmm", "r" );
//	fp = fopen( "model/walkin_3s5d.hmm", "r" );
//	fp = fopen( "model/walkin_trained_cont.hmm", "r" );
//	fp = fopen( "model/3dim.hmm", "r" );
	//fp = fopen( "../data/model/bicyclein/cont/bicyclein_3s6d.hmm", "r" );
	//fp = fopen( "../data/tmp.hmm", "r" );
	//fp = fopen( "model/model20s_cont.hmm", "r" );

//*	
	if( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found\n", "" );
		exit( 1 );
	}

	fp = fopen( argv[optind], "r" );

	hmm = new HMM();
	hmm->readContHMM( fp );
	hmm->printContHMM();
	fclose( fp );

//	fp = fopen( argv[2], "r" );
//	fp = fopen( "sequence/x1.seq", "r" );
//	fp = fopen( "sequence/test/test1.seq", "r" );
//	fp = fopen( "sequence/data_3s3d/data1.seq", "r" );
	//fp = fopen( "../data/bicyclein/cont/1.seq", "r" );
	//fp = fopen( "sequence/all_cont/1.seq", "r" );

	fp = fopen( argv[optind + 1], "r" );
	if( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found\n", "" );
		exit( 1 );
	}

	hmm->readContSequence( fp, &T, &O );
//	hmm->printContSequence( stdout, T, O );

	T = 300;

	alpha = (double **) calloc( (unsigned) ( T + 100 ), sizeof(double*) );
	for( int i = 1; i <= T + 100; i++ ) {
		alpha[i] = (double *) calloc( (unsigned) ( hmm->getNumStates() + 100 ), sizeof(double) );
	}
	beta = (double **) calloc( (unsigned) ( T + 100 ), sizeof(double*) );
	for( int i = 1; i <= T + 100; i++ ) {
		beta[i] = (double *) calloc( (unsigned) ( hmm->getNumStates() + 100 ), sizeof(double) );
	}
	gamma = (double **) calloc( (unsigned) ( T + 100 ), sizeof(double*) );
	for( int i = 1; i <= T + 100; i++ ) {
		gamma[i] = (double *) calloc( (unsigned) ( hmm->getNumStates() + 100 ), sizeof(double) );
	}
	scale = (double *) calloc ( (unsigned) ( T + 100 ), sizeof( double ) );
//*/
//	double X[] = { 0.1, 0.2, 0.15, 0.6, 0.22 };
//	cout << "Gaussian probability at X: " << hmm->computeGaussian( 1, X ) << endl;

//	hmm->forwardCont( T, O, alpha );
//	hmm->forwardContWithScale( T, O, alpha, scale );
//	fprintf( stdout, "Log prob(O|model) using Forward with scale = %E\n", hmm->getLoglik() );

//	hmm->readContSequence( fp, &T, &O );

//*
	hmm->trainContModel( NS, T, O, alpha, beta, gamma, &niter, &logprobinit, &logprobfinal );
	fprintf( stdout, "Number of iterations: %d\n", niter );
	fprintf( stdout, "Log prob(O|init model): %E\n", logprobinit );
	fprintf( stdout, "Log prob(O|estimated model): %E\n", logprobfinal );

	hmm->printContHMM();
//*/

/*
	T = 300;
	alpha = (double **) calloc( (unsigned) ( T + 100 ), sizeof(double*) );
	for( int i = 1; i <= T + 100; i++ ) {
		alpha[i] = (double *) calloc( (unsigned) ( 50 + 100 ), sizeof(double) );
	}
	scale = (double *) calloc ( (unsigned) ( T + 100 ), sizeof( double ) );

//	cout << optind << endl;

	fp = fopen( "../data/model/trained/cont/walkin_5s6d.hmm", "r" );
	hmm = new HMM();
	hmm->readContHMM( fp );
	fp = fopen( argv[optind], "r" );
	hmm->readContSequence( fp, &T, &O );
	hmm->forwardContWithScale( T, O, alpha, scale );
	fprintf( stdout, "%E\t", hmm->getLoglik() );

	fp = fopen( "../data/model/trained/cont/walkout_5s6d.hmm", "r" );
	hmm = new HMM();
	hmm->readContHMM( fp );
	fp = fopen( argv[optind], "r" );
	hmm->readContSequence( fp, &T, &O );
	hmm->forwardContWithScale( T, O, alpha, scale );
	fprintf( stdout, "%E\t", hmm->getLoglik() );

	fp = fopen( "../data/model/trained/cont/bicyclein_5s6d.hmm", "r" );
	hmm = new HMM();
	hmm->readContHMM( fp );
	fp = fopen( argv[optind], "r" );
	hmm->readContSequence( fp, &T, &O );
	hmm->forwardContWithScale( T, O, alpha, scale );
	fprintf( stdout, "%E\t", hmm->getLoglik() );

	fp = fopen( "../data/model/trained/cont/bicycleout_5s6d.hmm", "r" );
	hmm = new HMM();
	hmm->readContHMM( fp );
	fp = fopen( argv[optind], "r" );
	hmm->readContSequence( fp, &T, &O );
	hmm->forwardContWithScale( T, O, alpha, scale );
	fprintf( stdout, "%E\n", hmm->getLoglik() );
*/

	return 1;
}



