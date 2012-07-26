/*
 *	Author: Kan Ouivirach
 *	File: genseq.cc
 *	Purpose: test the function of reading a HMM
 */

#include "hmm.h"
#include <sys/types.h>
#include <unistd.h> 

using namespace std;

void Usage( char *name );

int main( int argc, char **argv ) {
	FILE *fp;				// HMM parameters in this file 
	extern char *optarg;
	extern int optind, opterr, optopt;
	
	HMM *hmm = new HMM();

/*
	// read HMM file 
	fp = fopen( argv[optind], "r" );		
	if ( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found \n", argv[optind] );
		Usage( argv[0] );
		exit ( 1 );
	}
	hmm->readHMM( fp );
	hmm->printHMM();
	// close the file, if any 
	if( fp ) {
		fclose( fp );
		fp = NULL;
	}
//*/
	// read HMM file
	fp = fopen( argv[optind], "r" );		
	if ( fp == NULL ) {
		fprintf( stderr, "Error: File %s not found \n", argv[optind] );
		Usage( argv[0] );
		exit ( 1 );
	}
	hmm->readHMMWithSuffStats( fp );
	hmm->printHMM();
	hmm->printSuffStats();
	// close the file, if any 
	if( fp ) {
		fclose( fp );
		fp = NULL;
	}
	
	delete hmm;
}

void Usage( char *name ) {
	printf( "Usage error \n" );
	printf( "Usage: %s <model.hmm> \n", name );
}




