/*
 *	Author: Kan Ouivirach
 *	File: genseq.cc
 *	Purpose: driver for generating a sequence of observation symbols.
 *	Credit: Tapas Kanungo, kanungo@cfar.umd.edu
 */

#include "hmm.h"
#include <sys/types.h>
#include <unistd.h> 

using namespace std;

void Usage( char *name );

int main( int argc, char **argv ) {
	int T; 		// length of observation sequence (default)
	int	*O;					// the observation sequence O[1..T]
	int	*q; 				// the state sequence q[1..T] 
	FILE *fp;				// HMM parameters in this file 
	int	sflg = 0, tflg = 0, errflg = 0;
	int	seed;				// random number seed 
	int	c;	
	extern char *optarg;
	extern int optind, opterr, optopt;

	char outputfile[10000] = "";

	while ( ( c = getopt( argc, argv, "S:T:o:" ) ) != EOF ) {
		switch( c ) {
		case 'S':
			// set random number generator seed 
			if( sflg )
				errflg++;
			else {	
				sflg++;
				sscanf( optarg, "%d", &seed );
			}
			break;
		case 'T':
			// set sequence length 
			if( tflg )
				errflg++;
			else {	
				tflg++;
				sscanf( optarg, "%d", &T );
			}
			break;
		case 'o':
			sscanf( optarg, "%s", outputfile );
			tflg++;
			break;
		case '?':
			errflg++;
		}
	}
	
	if ( ( argc - optind ) != 1 ) errflg++; // number or arguments not OK
/*
	cout << "argc: " << argc << endl;
	cout << "optind: " << optind << endl;
	cout << "tflg: " << tflg << endl;
	cout << "errflg: " << errflg << endl;
*/
	if ( errflg || !tflg ) {
		Usage( argv[0] );
		exit( 1 );
	}

	// read HMM file 
	fp = fopen( argv[optind], "r" );		
	if ( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found \n", argv[optind] );
		exit ( 1 );
	}
	HMM *hmm = new HMM();
	hmm->readHMM( fp );
	// close the file, if any 
	if( fp ) {
		fclose( fp );
		fp = NULL;
	}

	// length of observation sequence, T 
	O = (int *) calloc ( ( unsigned ) ( T + 1 ), sizeof( int ) );
	q = (int *) calloc ( ( unsigned ) ( T + 1 ), sizeof( int ) );

	if( !sflg ) seed = ( int ) getpid();

//	fprintf( stderr, "RandomSeed: %d\n", seed );
	hmm->genSequenceArray( seed, T, O, q );
	
//	hmm->printSequence( stdout, T, O );
//	cout << "generating a sequence.." << endl;

	FILE *fw;	// file to write
//	cout << "name: " << name << endl;
	fw = fopen( outputfile, "w" );
	if ( fw != NULL ) {
		hmm->printSequence( fw, T, O );
	}
	else {
		cerr << "Cannot open the file to print." << endl;
	}
	// close the file, if any 
	if( fw ) {
		fclose( fw );
		fw = NULL;
	}

	delete hmm;
}

void Usage( char *name ) {
	printf( "Usage error \n" );
	printf( "Usage: %s -T <sequence length> <model.hmm> \n", name );
	printf( "Usage: %s -S <seed> -T <sequence length> <model.hmm> \n", name );
	printf( "Usage: %s -S <seed> -T <sequence length> -o <output sequence file name> <model.hmm> \n", name );
	printf( "  T = length of sequence\n" );
	printf( "  S = random number seed \n" );
	printf( "  o = output sequence file name\n" );
	printf( "  model.hmm is a file with HMM parameters\n" );
}




