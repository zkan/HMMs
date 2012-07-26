#include "hmm.h"

int main( int argc, char **argv ) {
	FILE *fp;
	int T;
	HMM *hmm;
	int *O;		/* observations */

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
	hmm->viterbi( T, O );
	fprintf( stdout, "Log prob(O|model) using Viterbi = %E\n", hmm->getLoglik() );
	fclose( fp );

	return 1;
}
