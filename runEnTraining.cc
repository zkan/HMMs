/*
 *	Author: Kan Ouivirach
 *	File: enTraining.cc
 *	Purpose: Implement the Ensemble Training
 *	Credit: Cavalin: Evaluation of Incremental Learning Algorithms for An HMM-Based Handwritten Isolated Digits Recognizer in ICFHR'08
 */

#include "hmm.h"

string func_itos( int i ) {
	stringstream s;
	s << i;
	return s.str();
}

int main( int argc, char **argv ) {
	FILE *fp;
	int T;					// length of an observation sequence
	HMM *curr_hmm, *new_hmm;		// a current model and a new model
	int *O;					// observation sequence
	double **alpha, **beta, **gamma;
	double logprobinit, logprobfinal;
	int niter;				// number of iterations
	float curr_weight = 0;			// the current weight for the current model
	float new_weight = 1;			// a new weight for a new model

	double **tmp_curr_A;		/* Temporal Current Transition probability */
	double **tmp_curr_B;		/* Temporal Current Emission probability */
	double *tmp_curr_pi;		/* Temporal Current Initial state distribution */

	double **tmp_new_A;		/* Temporal New Transition probability */
	double **tmp_new_B;		/* Temporal New Emission probability */
	double *tmp_new_pi;		/* Temporal New Initial state distribution */

	int c;
	extern char *optarg;
	extern int optind, opterr, optopt;

	int num_of_models = 1;

	while( ( c = getopt( argc, argv, "N:" ) ) != EOF ) {
		switch( c ) {
			case 'N':
				sscanf( optarg, "%d", &num_of_models );
				break;
		}
	}

//	cout << argc << endl;
//	cout << optind << endl;

	if(  ( argc - optind ) != 2 ) {
		cout << "Usage error" << endl;
		cout << "Usage: xxx <curr_model.hmm> <sequence.seq>" << endl;
		exit( 1 );
	}
	
	fp = fopen( argv[argc - 2], "r" );
	if( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found\n", argv[1] );
		exit( 1 );
	}
	curr_hmm = new HMM();
	curr_hmm->readHMM( fp );
//	curr_hmm->printHMM();
	fclose( fp );

	for( int n = 1; n <= num_of_models; n++ ) {
		// use the same initial model for every sequence
		// 
		fp = fopen( argv[argc - 2], "r" );
		if( fp == NULL ) {
			fprintf( stderr, "Error: File %s not found\n", argv[1] );
			exit( 1 );
		}
	//	int seed = (int)getpid();
		new_hmm = new HMM();
		new_hmm->readHMM( fp );
		fclose( fp );
	//	new_hmm->initHMM( seed, curr_hmm->getNumStates(), curr_hmm->getNumSymbols() );
	//	new_hmm->printHMM();

		// train the new model with a new sequence
		// 
		string ss;
		char tmp[1000] = "";
		char filename[100] = "/home/zkan/Desktop/accv-video/final/activity-group/5/";
		ss = func_itos( n );
		strcpy( tmp, ss.c_str() );
		strcat( filename, strcat( tmp, ".seq" ) );

//		fp = fopen( filename, "r" );
		fp = fopen( filename, "r" );
		if( fp == NULL ) {
			fprintf( stderr, "Error: File %s not found\n", filename );
			exit( 1 );
		}
		new_hmm->readSequence( fp, &T, &O );
		fclose( fp );

		// set weight = the observation length
		//
		new_weight = T;
//		T = 100000;

		alpha = (double **) calloc( (unsigned) ( T + 1 ), sizeof(double*) );
		for( int i = 1; i <= T; i++ ) {
			alpha[i] = (double *) calloc( (unsigned) ( new_hmm->getNumStates() + 1 ), sizeof(double) );
		}
		beta = (double **) calloc( (unsigned) ( T + 1 ), sizeof(double*) );
		for( int i = 1; i <= T; i++ ) {
			beta[i] = (double *) calloc( (unsigned) ( new_hmm->getNumStates() + 1 ), sizeof(double) );
		}
		gamma = (double **) calloc( (unsigned) ( T + 1 ), sizeof(double*) );
		for( int i = 1; i <= T; i++ ) {
			gamma[i] = (double *) calloc( (unsigned) ( new_hmm->getNumStates() + 1 ), sizeof(double) );
		}

		new_hmm->baumWelch( T, O, alpha, beta, gamma, &niter, &logprobinit, &logprobfinal );
		fprintf( stdout, "Number of iterations: %d\n", niter );
		fprintf( stdout, "Log prob(O|init model): %E\n", logprobinit );
		fprintf( stdout, "Log prob(O|estimated model): %E\n", logprobfinal );
	//	new_hmm->printHMM();

		int N = new_hmm->getNumStates();	// the number of states
		int M = new_hmm->getNumSymbols();	// the number of symbols

		// initialize the temp variables
		//
		tmp_curr_A = (double **) calloc( (unsigned) ( N + 1 ), sizeof(double*) );
		for( int i = 1; i <= N; i++ ) {
			tmp_curr_A[i] = (double *) calloc( (unsigned) ( N + 1 ), sizeof(double) );
		}
		tmp_curr_B = (double **) calloc( (unsigned) ( N + 1 ), sizeof(double*) );
		for( int i = 1; i <= N; i++ ) {
			tmp_curr_B[i] = (double *) calloc( (unsigned) ( M + 1 ), sizeof(double) );
		}
		tmp_curr_pi = (double *) calloc( (unsigned) ( N + 1 ), sizeof( double ) );

		tmp_curr_A = curr_hmm->getA();
		tmp_curr_B = curr_hmm->getB();
		tmp_curr_pi = curr_hmm->getPi();
		
		tmp_new_A = (double **) calloc( (unsigned) ( N + 1 ), sizeof(double*) );
		for( int i = 1; i <= N; i++ ) {
			tmp_new_A[i] = (double *) calloc( (unsigned) ( N + 1 ), sizeof(double) );
		}
		tmp_new_B = (double **) calloc( (unsigned) ( N + 1 ), sizeof(double*) );
		for( int i = 1; i <= N; i++ ) {
			tmp_new_B[i] = (double *) calloc( (unsigned) ( M + 1 ), sizeof(double) );
		}
		tmp_new_pi = (double *) calloc( (unsigned) ( N + 1 ), sizeof( double ) );

		tmp_new_A = new_hmm->getA();
		tmp_new_B = new_hmm->getB();
		tmp_new_pi = new_hmm->getPi();

		// merge the current model and the new trained model
		// 
		float tmp_weight = curr_weight + new_weight;
		for( int i = 1; i <= N; i++ ) {
			for( int j = 1; j <= N; j++ ) {
				tmp_new_A[i][j] = ( ( curr_weight * tmp_curr_A[i][j] ) + ( new_weight * tmp_new_A[i][j] ) ) / tmp_weight;
			}
		}
		for( int i = 1; i <= N; i++ ) {
			for( int j = 1; j <= M; j++ ) {
				tmp_new_B[i][j] = ( ( curr_weight * tmp_curr_B[i][j] ) + ( new_weight * tmp_new_B[i][j] ) ) / tmp_weight;
			}
		}
		for( int i = 1; i <= N; i++ ) {
			tmp_new_pi[i] = ( ( curr_weight * tmp_curr_pi[i] ) + ( new_weight * tmp_new_pi[i] ) ) / tmp_weight;
		}
		curr_hmm->setPi( tmp_new_pi );
		curr_hmm->setA( tmp_new_A );
		curr_hmm->setB( tmp_new_B );
//*
		float checkSum = 0;
		cout << "A: " << endl;
		for( int i = 1; i <= curr_hmm->getNumStates(); i++ ) {
			cout << "sum of row " << i << ": ";
			for( int j = 1; j <= curr_hmm->getNumStates(); j++ ) {
				checkSum += tmp_new_A[i][j];
			}
			cout << checkSum << endl;
			checkSum = 0;
		}
		cout << "B: " << endl;
		for( int i = 1; i <= curr_hmm->getNumStates(); i++ ) {
			cout << "sum of row " << i << ": ";
			for( int j = 1; j <= curr_hmm->getNumSymbols(); j++ ) {
				checkSum += tmp_new_B[i][j];
			}
			cout << checkSum << endl;
			checkSum = 0;
		}
		cout << "Pi: " << endl;
		cout << "sum of Pi: ";
		for( int i = 1; i <= curr_hmm->getNumStates(); i++ ) {
			checkSum += tmp_new_pi[i];
		}
		cout << checkSum << endl;
		checkSum = 0;
//*/
		curr_hmm->printHMM();
	//	new_hmm->printHMM();
		curr_weight = tmp_weight;
		cout << "new total weight: " << curr_weight << endl;

//		for( int i = 1; i <= T; i++ ) {
//			free( alpha[i] );
//		}
		free( alpha );
//		for( int i = 1; i <= T; i++ ) {
//			free( beta[i] );
//		}
		free( beta );
//		for( int i = 1; i <= T; i++ ) {
//			free( gamma[i] );
//		}
		free( gamma );

	}

	return 1;
}
