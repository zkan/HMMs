#include "hmm.h"
#include "utility.h"

#include <unistd.h>
#include <sys/time.h>
#include <stdio.h>
#include <ctime>

string func_itos( int i ) {
	stringstream s;
	s << i;
	return s.str();
}

int main( int argc, char **argv ) {
	FILE *fp;
	int T;
	HMM *hmm;
	int *O;		// observations
	double **alpha, **beta, **gamma;
	double logprobinit;
	int niter = 100;	// number of iterations

	// forgetting factor for sufficient statistics
	double forgetting_factor = 1;

	int c;
	extern char *optarg;
	extern int optind, opterr, optopt;

	while( ( c = getopt( argc, argv, "f:l:" ) ) != EOF ) {
		switch( c ) {
			case 'f':
				sscanf( optarg, "%lf", &forgetting_factor );
				break;
			case 'l':
				sscanf( optarg, "%d", &niter );
				break;
		}
	}
	if( forgetting_factor > 1 ) {
		cerr << "Usage error: the forgetting factor is less than or equal to 1!" << endl;
		return EXIT_FAILURE;
	}

	if( ( ( argc - optind ) != 2 ) || ( argc == 2 ) ) {
		cout << "Usage error" << endl;
		cout << "Usage: " << argv[0] << " <model.hmm> <the sequence directory>" << endl;
		cout << "Option: " << endl;
		cout << "-f <forgetting factor>" << endl;
		cout << "-l <number of iterations>" << endl;
		return EXIT_FAILURE;
	}
	
	fp = fopen( argv[optind], "r" );
	if( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found\n", argv[1] );
		return EXIT_FAILURE;
	}

	hmm = new HMM();
	hmm->readHMMWithSuffStats( fp );
//	cout << "### Initial HMM ###" << endl;
//	hmm->printHMM();
//	hmm->printSuffStats();
	if( fp ) {
		fclose( fp );
		fp = NULL;
	}
	
	char *dir_name = argv[optind + 1];
//	cout << dir_name << endl;
//	exit( 1 );

	// load the sufficient statistics
	double *pi_t = hmm->getPi_t();
	double **numeratorA_t = hmm->getNumeratorA_t();
	double *denominatorA_t = hmm->getDenominatorA_t();
	double **numeratorB_t = hmm->getNumeratorB_t();
	double *denominatorB_t = hmm->getDenominatorB_t();
	double logprobfinal = hmm->getLogProbFinal();
	double logprob_trained = hmm->getLogProbFinal();
	int num_of_seqs_trained = hmm->getNS();
//	hmm->printSuffStats();
//	exit( 1 );
/*
	// the sufficient statistics
	double ***xi_t;
	double **gamma_t;
	xi_t = hmm->allocXi( T + 1, hmm->getNumStates() + 1 );
	gamma_t = ( double** ) calloc( ( unsigned ) ( T + 1 ), sizeof( double* ) );
	for( int i = 1; i <= T; i++ ) {
		gamma_t[i] = ( double* ) calloc( ( unsigned ) ( hmm->getNumStates() + 1 ), sizeof( double ) );
	}
	// set the sufficient statistics
	for( int t = 1; t <= T; t++ ) {
		for( int i = 1; i <= hmm->getNumStates(); i++ ) {
			gamma_t[t][i] = 0;
		}
	}
	for( int t = 1; t <= T - 1; t++ ) {
		for( int i = 1; i <= hmm->getNumStates(); i++ ) {
			for( int j = 1; j <= hmm->getNumStates(); j++ ) {
				xi_t[t][i][j] = 0;
			}
		}
	}
*/

/*
	for( int i = 0; files[i]; ++i ) {
		// ignore . and ..
		if( files[i][0] != '.' ) {
//			printf( "%s\n", files[i] );
			os_seq_name_tmp << dir_name << files[i];

			FILE *fp;
			fp = fopen( os_seq_name_tmp.str().c_str(), "r" );
			if( fp == NULL ) {
				fprintf( stderr, "Error: File %s not found\n", os_seq_name_tmp.str().c_str() );
				exit( 1 );
			}
			readSequence( fp, &T, &O );
			if( fp ) {
				fclose( fp );
				fp = NULL;
			}
			os_seq_name_tmp.str( "" );
		}
	}
*/
	char **files = get_all_files( dir_name );
	if( !files ) {
		cerr << "Error: cannot get all files under the directory!" << endl;
		return EXIT_FAILURE;
	}
	ostringstream os_seq_name;
//	cout << "### Sequences ###" << endl;
	for( int i = 0; files[i]; ++i ) {
		// ignore . and ..
		if( files[i][0] != '.' ) {
			os_seq_name << dir_name << files[i];
//			cout << os_seq_name.str() << endl;
			os_seq_name.str( "" );
		}
	}
	
	// avoid the situation when the next sequence is shorter or longer
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

	int actual_obs_length = 0;
//	cout << "### Training process ###" << endl;
	for( int i = 0; files[i]; ++i ) {
		// ignore . and ..
		if( files[i][0] != '.' ) {
			os_seq_name << dir_name << files[i];
			
			FILE *fp;
			fp = fopen( os_seq_name.str().c_str(), "r" );
			if( fp == NULL ) {
				fprintf( stderr, "Error: File %s not found\n", os_seq_name.str().c_str() );
				exit( 1 );
			}
			hmm->readSequence( fp, &actual_obs_length, &O );
//			ostringstream os_seq_name;
//			os_seq_name << argv[optind + 1];
			if( fp ) {
				fclose( fp );
				fp = NULL;
			}
//			cout << os_seq_name.str() << endl; exit(1);
			
//			cout << "number of seqs trained: " << num_of_seqs_trained << endl;

    timeval start, end;
    long mtime, seconds, useconds;
    gettimeofday(&start, NULL);

//      cout << "seq: " << os_seq_name.str() << endl;
    hmm->train_iml( T, O, alpha, beta, gamma, niter, &logprobinit, &logprobfinal, pi_t, numeratorA_t, denominatorA_t, numeratorB_t, denominatorB_t, 1, os_seq_name, forgetting_factor, actual_obs_length, num_of_seqs_trained, logprob_trained );
			
    gettimeofday(&end, NULL);
    seconds = end.tv_sec - start.tv_sec;
    useconds = end.tv_usec - start.tv_usec;
    mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
//    printf("Elapsed time: %ld milliseconds\n", mtime);
						
			logprob_trained = logprobfinal;
			num_of_seqs_trained++;
						
			os_seq_name.str( "" );
//			fprintf( stdout, "Number of iterations: %d\n", niter );
//			fprintf( stdout, "Log prob(O|init model): %E\n", logprobinit );
//			fprintf( stdout, "Log prob(O|estimated model): %E\n", logprobfinal );
			
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
	/*
			cout << endl << "## Checksum ##" << endl;
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
			cout << "number of seqs trained: " << num_of_seqs_trained << endl;
			for( int i = 1; i <= hmm->getNumStates(); i++ ) {
				cout <<  pi_t[i] << " ";
			}
			cout << endl;
			
	//		pi_t = hmm->getPi_t();
	//		numeratorA_t = hmm->getNumeratorA_t();
	//		denominatorA_t = hmm->getDenominatorA_t();
	//		numeratorB_t = hmm->getNumeratorB_t();
	//		denominatorB_t = hmm->getDenominatorB_t();
	
		//*/
	//		cout << "### Final HMM ###" << endl;
//			hmm->printHMM();
	//		cout << "## Sufficient statistics ##" << endl;
//			hmm->printSuffStats();
//			cout << endl;
		}
	//	if( i == 2 ) exit( 1 );
	}
	
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
	
/*
	exit( 1 );
	for( int n = 1; n <= 1; n++ ) {
//		scanf( "%s", filename );
//		if( !strcmp( filename, "q" ) ) {
//			break;
//		}

//		fp = fopen( argv[argc - 1], "r" );

		// Read the sequence name
		// 
	//	cout << func_itos( n ) << endl;

		string ss;
		char tmp[1000] = "";
		char filename[100] = "/home/zkan/AIT-PhD-Stuffs/src/HMMs/sequence/train-1-seq/";
		ss = func_itos( n );
		strcpy( tmp, ss.c_str() );
		strcat( filename, strcat( tmp, ".seq" ) );

		fp = fopen( filename, "r" );
		if( fp == NULL ) {
			fprintf( stderr, "Error: File %s not found\n", filename );
			exit( 1 );
		}		
		hmm->readSequence( fp, &actual_obs_length, &O );
	//	cout << T << endl;
	
		ostringstream os_seq_name;
		os_seq_name << argv[optind + 1];
		
//		dir_name = argv[optind + 1];

		// avoid the situation when the next sequence is shorter or longer
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

		hmm->train_iml( dir_name, T, O, alpha, beta, gamma, niter, &logprobinit, &logprobfinal, 
					pi_t, numeratorA_t, denominatorA_t, numeratorB_t, denominatorB_t, n, os_seq_name, 
					forgetting_factor, actual_obs_length, 0, 0 );

	//	cout << "Gamma: " << endl;
	//	for( int t = 1; t <= T; t++ ) {
	//		for( int s = 1; s <= hmm->getNumStates(); s++ ) {
	//			cout << gamma[t][s] << " ";
	//		}
	//		cout << endl;
	//	}
				
		fprintf( stdout, "Number of iterations: %d\n", niter );
		fprintf( stdout, "Log prob(O|init model): %E\n", logprobinit );
		fprintf( stdout, "Log prob(O|estimated model): %E\n", logprobfinal );

		hmm->printHMM();

		// -------------------------------------------------------------------------------------------------
		
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

		cout << endl << "## Checksum ##" << endl;
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

		// -------------------------------------------------------------------------------------------------

		checkSum = 0;
		cout << "Gamma:" << endl;
		for( int t = 1; t <= actual_obs_length; t++ ) {
			cout << "T= " << t << endl;
			for( int i = 1; i <= hmm->getNumStates(); i++ ) {
				checkSum += gamma_t[t][i];
			}
			cout << checkSum << endl;
			checkSum = 0;
		}

		cout << "Xi:" << endl;
		for( int t = 1; t <= actual_obs_length - 1; t++ ) {
			for( int i = 1; i <= hmm->getNumStates(); i++ ) {
				for( int j = 1; j <= hmm->getNumSymbols(); j++ ) {
					cout << xi_t[t][i][j] << " ";
				}
				cout << endl;
			}
			cout << endl;
		}
		cout << endl;

		cout << "### Final HMM ###" << endl;
		hmm->printHMM();
		cout << "## Sufficient statistics ##" << endl;
		hmm->printSuffStats();

//		free( alpha );
//		free( beta );
//		free( gamma );
//		free( gamma_t );
//		free( xi_t );

//		fprintf( stdout, "Log prob(O|init model): %E\n", logprobinit );
//		fprintf( stdout, "Log prob(O|estimated model): %E\n", logprobfinal );
	}
*/

//	fclose( fp );

//	return EXIT_SUCCESS;
}



