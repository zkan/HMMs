#include "hmm.h"
#include <sys/types.h>
#include <unistd.h>

int main( int argc, char **argv ) {
	FILE *fp;
	int T;
	HMM *hmm;
	int *O;		/* observations */
	double **alpha, **beta, **gamma;
	double logprobinit, logprobfinal;
	int niter;	/* number of iterations */

	int N = 1;	/* number of states */
	int M = 1;	/* number of symbols */
	int type = 0;	// 0 - full-connected model, 1 - left-to-right model
	int seed;

	int sflg = 0;

	int c;
	extern char *optarg;
	extern int optind, opterr, optopt;

	while( ( c = getopt( argc, argv, "S:N:M:t:" ) ) != EOF ) {
		switch( c ) {
			case 'S':
				sscanf( optarg, "%d", &seed );
				sflg = 1;
				break;
			case 'N':
				sscanf( optarg, "%d", &N );
				break;
			case 'M':
				sscanf( optarg, "%d", &M );
				break;
			case 't':
				sscanf( optarg, "%d", &type );
				break;	
		}
	}

//	cout << argc << endl;
//	cout << optind << endl;

	if( ( argc - optind ) != 0 ) {
		cout << "Usage error" << endl;
		cout << "Usage: xxx -N <number of states> -M <number of symbols>" << endl;
		cout << "Usage: xxx -S <seed> -N <number of states> -M <number of symbols>" << endl;
		exit( 1 );
	}

	if( !sflg ) seed = (int)getpid();
	
	hmm = new HMM();
	hmm->initHMM( seed, N, M, type );
	hmm->printHMM();
//	hmm->printSuffStats();

//	cout << "######################################################################" << endl;
	
//	hmm->initHMM( seed, N, M, 1 );
//	hmm->printHMM();

	return 1;
}



