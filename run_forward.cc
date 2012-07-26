#include "hmm.h"

#define HMM_METHOD 0
#define SVM_OR_PCA_METHOD 1

int main( int argc, char **argv ) {
	FILE *fp;
	int T;
	HMM *hmm;
	int *O;		/* observations */
	double **alpha, *scale;
	double proba;

	int c;
	extern char *optarg;
	extern int optind, opterr, optopt;
//*
	int tmp_T = 0;
	int window_size = 5;
	int method = 0;	// 0 for HMM, 1 for PCA or SVM
	int type = 0;
	while( ( c = getopt( argc, argv, "T:m:w:" ) ) != EOF ) {
		switch( c ) {
			case 'T':
				sscanf( optarg, "%d", &tmp_T );
				break;
			case 'm':
				sscanf( optarg, "%d", &method );
				break;
		   case 't':
				sscanf( optarg, "%d", &type );
				break;
			case 'w':
				sscanf( optarg, "%d", &window_size );
				break;
		}
	}
//	cout << "argc: " << argc << endl;
//	cout << "optind: " << optind << endl;

	if( ( argc - optind ) != 2 ) {
		cout << "Usage error" << endl;
		cout << "Usage: xxx -T <num of obs> <model.hmm> <sequence.seq>" << endl;
		exit( 1 );
	}
	
	fp = fopen( argv[optind], "r" );
	if( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found\n", argv[optind] );
		exit( 1 );
	}

	hmm = new HMM();
	hmm->readHMM( fp );
	//hmm->printHMM();
	fclose( fp );

	fp = fopen( argv[optind + 1], "r" );
	if( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found\n", argv[optind + 1] );
		exit( 1 );
	}

	hmm->readSequence( fp, &T, &O );
	if( tmp_T ) {
		if( tmp_T > T ) {
			cout << "Usage error" << endl;
			cout << "Cannot set the number of observations bigger than the actual one" << endl;
			exit( 1 );
		}
		T = tmp_T;
	}
//*/
	alpha = ( double ** ) calloc( ( unsigned ) ( T + 1 ), sizeof( double* ) );
	for( int i = 1; i <= T; i++ ) {
		alpha[i] = ( double * ) calloc( ( unsigned ) ( hmm->getNumStates() + 1 ), sizeof( double ) );
	}
	scale = ( double * ) calloc ( ( unsigned ) ( T + 1 ), sizeof( double ) );
	
	for( int i = 1; i <= T; i++ ) {
		for( int j = 1; j <= hmm->getNumStates(); j++ ) {
			alpha[i][j] = 0;
		}
		scale[i] = 0;
	}
		
	if( method == HMM_METHOD ) {
	//	hmm->forward( T, O, alpha );
	//	fprintf( stdout, "Log prob(O|model) using Forward = %E\n", hmm->getLoglik() );

		hmm->forwardWithScale( T, O, alpha, scale, type );
	//	fprintf( stdout, "Log prob(O|model) using Forward with Scale = %E\n", hmm->getLoglik() );
	//	fprintf( stdout, "Log prob per observation (O|model) using Forward with Scale = %E\n", hmm->getLoglik()/T );

		if( isnan( hmm->getLoglik() ) ) {
			fprintf( stdout, "-50\n" );
		}
		else {
			fprintf( stdout, "%lf\n", hmm->getLoglik()/T );
		}

		// Test the forward algorithm for L-to-R model
		/*
		for( int i = 1; i <= T; i++ ) {
			for( int j = 1; j <= hmm->getNumStates(); j++ ) {
				alpha[i][j] = 0;
			}
			scale[i] = 0;
		}		
		hmm->forwardWithScale( T, O, alpha, scale, 1 );
		if( isnan( hmm->getLoglik() ) ) {
			fprintf( stdout, "-1000000" );
		}
		else {
			fprintf( stdout, "%lf\n", hmm->getLoglik()/T );
		}
		*/
	}
	else if( method == SVM_OR_PCA_METHOD ) {
		vector<double> loglik;
		int t = 1;
		
		do {
			hmm->forwardWithScale( t, O, alpha, scale, 0 );
			
			if( isnan( hmm->getLoglik() ) ) {
				fprintf( stdout, "-50\n" );
				
				loglik.push_back( -50 );
			}
			else {				
			//	fprintf( stdout, "%lf\n", hmm->getLoglik()/t );
				loglik.push_back( hmm->getLoglik() );
			}
			
			if( loglik.size() == window_size ) {
				double mean_loglik_in_window = 0;
				for( int i = 0; i < loglik.size(); i++ ) {
					mean_loglik_in_window += loglik[i];
				}
				mean_loglik_in_window /= loglik.size();
				cout << mean_loglik_in_window << endl;
				
			//	loglik.erase( loglik.begin() );
				loglik.clear();
			}
			
			t++;
		}
		while( t <= T );
	}

//	if( hmm->getLoglik()/T < -1.32718 ) cout << "ALARM!!!" << endl;

//	if( hmm->getLoglik()/T < -0.82524416352941 ) cout << "ALARM!!!" << endl;

//	if( isnan( hmm->getLoglik()/T ) ) cout << "ALARM!!!" << endl;

/*
	T = 300;
	alpha = (double **) calloc( (unsigned) ( T + 1 ), sizeof(double*) );
	for( int i = 1; i <= T; i++ ) {
		alpha[i] = (double *) calloc( (unsigned) ( 50 + 1 ), sizeof(double) );
	}
	scale = (double *) calloc ( (unsigned) (T + 1), sizeof( double ) );

	fp = fopen( "../data/model/trained/walkin_5s5s.hmm", "r" );
	hmm = new HMM();
	hmm->readHMM( fp );
	fp = fopen( argv[optind], "r" );
	hmm->readSequence( fp, &T, &O );
	hmm->forwardWithScale( T, O, alpha, scale );
	fprintf( stdout, "%E\t", hmm->getLoglik() );

	fp = fopen( "../data/model/trained/walkout_5s5s.hmm", "r" );
	hmm = new HMM();
	hmm->readHMM( fp );
	fp = fopen( argv[optind], "r" );
	hmm->readSequence( fp, &T, &O );
	hmm->forwardWithScale( T, O, alpha, scale );
	fprintf( stdout, "%E\t", hmm->getLoglik() );

	fp = fopen( "../data/model/trained/bicyclein_5s5s.hmm", "r" );
	hmm = new HMM();
	hmm->readHMM( fp );
	fp = fopen( argv[optind], "r" );
	hmm->readSequence( fp, &T, &O );
	hmm->forwardWithScale( T, O, alpha, scale );
	fprintf( stdout, "%E\t", hmm->getLoglik() );

	fp = fopen( "../data/model/trained/bicycleout_5s5s.hmm", "r" );
	hmm = new HMM();
	hmm->readHMM( fp );
	fp = fopen( argv[optind], "r" );
	hmm->readSequence( fp, &T, &O );
	hmm->forwardWithScale( T, O, alpha, scale );
	fprintf( stdout, "%E\n", hmm->getLoglik() );
*/

//	cout << optind << endl;

/*
	T = 10000;
	alpha = (double **) calloc( (unsigned) ( T + 1 ), sizeof(double*) );
	for( int i = 1; i <= T; i++ ) {
		alpha[i] = (double *) calloc( (unsigned) ( T + 1 ), sizeof(double) );
	}
	scale = (double *) calloc ( (unsigned) (T + 1), sizeof( double ) );

	fp = fopen( "/home/zkan/Desktop/Pat_work/data/dental_expertModel.hmm", "r" );
	hmm = new HMM();

	hmm->readHMM( fp );
	fp = fopen( argv[optind], "r" );
	hmm->readSequence( fp, &T, &O );
	hmm->forwardWithScale( T, O, alpha, scale );
	fprintf( stdout, "Log likelood from the expert model: %E\n", hmm->getLoglik() );

	fp = fopen( "/home/zkan/Desktop/Pat_work/data/dental_noviceModel.hmm", "r" );
	hmm = new HMM();
	hmm->readHMM( fp );
	fp = fopen( argv[optind], "r" );
	hmm->readSequence( fp, &T, &O );
	hmm->forwardWithScale( T, O, alpha, scale );
	fprintf( stdout, "Log likelood from the novice model: %E\n", hmm->getLoglik() );
//*/

	return 1;
}







