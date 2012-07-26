/*
 *	Author: Kan Ouivirach
 *	File: hmm.cc
 *	Purpose: Implement HMM classes
 *	Credit: Tapas Kanungo, kanungo@cfar.umd.edu
 */

#include "hmm.h"
#define PI 3.141592653589793

void HMM::setPi( double *pi ) {
	for( int i = 1; i <= this->N; i++ ) {
		this->pi[i] = pi[i];
	}
}

void HMM::setA( double **A ) {
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			this->A[i][j] = A[i][j];
		}
	}
}

void HMM::setB( double **B ) {
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->M; j++ ) {
			this->B[i][j] = B[i][j];
		}
	}
}

void HMM::readHMM( FILE *fp ) {
	fscanf( fp, "M= %d\n", &( this->M ) );
	fscanf( fp, "N= %d\n", &( this->N ) );

	fscanf( fp, "A:\n");
	this->A = ( double** )  calloc( ( unsigned ) ( this->N + 1 ), sizeof( double* ) );
	for( int i = 1; i <= this->N; i++ ) {
		this->A[i] = ( double* ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	}
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			fscanf( fp, "%lf", &( this->A[i][j] ) );
		}
		fscanf( fp, "\n" );
	}

	fscanf( fp, "B:\n" );
	this->B = ( double** ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double* ) );
	for( int i = 1; i <= this->N; i++ ) {
		this->B[i] = ( double * ) calloc( ( unsigned ) ( this->M + 1 ), sizeof( double ) );
	}
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->M; j++ ) {
			fscanf( fp, "%lf", &( this->B[i][j] ) );
		}
		fscanf( fp, "\n" );
	}

	fscanf( fp, "pi:\n" );
	this->pi = ( double* ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	for( int i = 1; i <= this->N; i++ ) {
		fscanf( fp, "%lf", &( this->pi[i] ) );
	}

	// initialize the sufficient statistics
	this->pi_t = ( double* ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	for( int i = 1; i <= this->N; i++ ) {
		this->pi_t[i] = 0;
	}

	this->numeratorA_t = ( double** ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double* ) );
	for( int i = 1; i <= this->N; i++ ) {
		this->numeratorA_t[i] = ( double* ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	}
	this->denominatorA_t = ( double* ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	for( int i = 1; i <= this->N; i++ ) {
		this->denominatorA_t[i] = 0;
		for( int j = 1; j <= this->N; j++ ) {
			this->numeratorA_t[i][j] = 0;
		}
	}

	this->numeratorB_t = ( double** ) calloc( (unsigned) ( this->N + 1 ), sizeof( double* ) );
	for( int i = 1; i <= this->N; i++ ) {
		this->numeratorB_t[i] = ( double* ) calloc( (unsigned) ( this->M + 1 ), sizeof( double ) );
	}
	this->denominatorB_t = ( double* ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	for( int i = 1; i <= this->N; i++ ) {
		this->denominatorB_t[i] = 0;
		for( int j = 1; j <= this->M; j++ ) {
			this->numeratorB_t[i][j] = 0;
		}
	}
	this->logprobfinal = 0;
	this->NS = 0;
}

void HMM::readHMMWithSuffStats( FILE *fp ) {
	fscanf( fp, "M= %d\n", &(this->M) );
	fscanf( fp, "N= %d\n", &(this->N) );

	fscanf( fp, "A:\n");
	this->A = ( double ** ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double* ) );
	for( int i = 1; i <= this->N; i++ ) {
		this->A[i] = ( double * ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	}
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			fscanf( fp, "%lf", &( this->A[i][j] ) );
		}
		fscanf( fp, "\n" );
	}

	fscanf( fp, "B:\n" );
	this->B = ( double** ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double* ) );
	for( int i = 1; i <= this->N; i++ ) {
		this->B[i] = ( double* ) calloc( ( unsigned ) ( this->M + 1 ), sizeof( double ) );
	}
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->M; j++ ) {
			fscanf( fp, "%lf", &( this->B[i][j] ) );
		}
		fscanf( fp, "\n" );
	}

	fscanf( fp, "pi:\n" );
	this->pi = ( double* ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	for( int i = 1; i <= this->N; i++ ) {
		fscanf( fp, "%lf", &( this->pi[i] ) );
	}
	fscanf( fp, "\n" );

	// to store the sufficient statistics
	fscanf( fp, "numeratorA_t:\n" );
	this->numeratorA_t = ( double** ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double* ) );
	for( int i = 1; i <= this->N; i++ ) {
		this->numeratorA_t[i] = ( double* ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	}
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			fscanf( fp, "%lf", &( this->numeratorA_t[i][j] ) );
		}
		fscanf( fp, "\n" );
	}
	fscanf( fp, "denominatorA_t:\n" );
	this->denominatorA_t = ( double* ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	for( int i = 1; i <= this->N; i++ ) {
		fscanf( fp, "%lf", &( this->denominatorA_t[i] ) );
	}
	fscanf( fp, "\n" );
	
	fscanf( fp, "numeratorB_t:\n" );
	this->numeratorB_t = ( double** ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double* ) );
	for( int i = 1; i <= this->N; i++ ) {
		this->numeratorB_t[i] = ( double* ) calloc( ( unsigned ) ( this->M + 1 ), sizeof( double ) );
	}
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->M; j++ ) {
			fscanf( fp, "%lf", &( this->numeratorB_t[i][j] ) );
		}
		fscanf( fp, "\n" );
	}
	fscanf( fp, "denominatorB_t:\n" );
	this->denominatorB_t = ( double* ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	for( int i = 1; i <= this->N; i++ ) {
		fscanf( fp, "%lf", &( this->denominatorB_t[i] ) );
	}
	fscanf( fp, "\n" );
	
	fscanf( fp, "pi_t:\n" );
	this->pi_t = ( double* ) calloc( ( unsigned ) ( this->N + 1 ), sizeof( double ) );
	for( int i = 1; i <= this->N; i++ ) {
		 fscanf( fp, "%lf", &( this->pi_t[i] ) );
	}
	fscanf( fp, "\n" );
	
	fscanf( fp, "logprobfinal:\n" );
	fscanf( fp, "%lf", &( this->logprobfinal ) );
	fscanf( fp, "\n" );
	
	fscanf( fp, "NS:\n" );
	fscanf( fp, "%d", &( this->NS ) );
}

void HMM::readContHMM( FILE *fp ) {
	fscanf( fp, "N= %d\n", &(this->N) );
	fscanf( fp, "dim= %d\n", &(this->dim) );

	fscanf( fp, "A:\n");
	this->A = (double **) calloc( (unsigned) ( this->N + 1 ), sizeof(double*) );
	for( int i = 1; i <= this->N; i++ ) {
		this->A[i] = (double *) calloc( (unsigned) ( this->N + 1 ), sizeof(double) );
	}
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			fscanf( fp, "%lf", &(this->A[i][j]) );
		}
		fscanf( fp, "\n" );
	}

	fscanf( fp, "B:\n" );
	/* initialize mean matrix */
	this->meanMat = (double **) calloc( (unsigned) ( this->N + 1 ), sizeof(double*) );
	for( int i = 1; i <= this->N; i++ ) {
		this->meanMat[i] = (double *) calloc( (unsigned) ( this->dim + 1 ), sizeof(double) );
	}
	/* initialize cov matrix */
	this->covMat = (double ***) calloc( (unsigned) ( this->N + 1 ), sizeof(double**) );
	for( int i = 1; i <= this->N; i++ ) {
		this->covMat[i] = (double **)calloc( (unsigned) (this->dim + 1), sizeof(double*) );
		for( int j = 1; j <= this->dim; j++ ) {
			this->covMat[i][j] = (double *)calloc( (unsigned) (this->dim + 1), sizeof(double) );
		}
	}
	for( int i = 1; i <= this->N; i++ ) {
		fscanf( fp, "mean:\n" );
		for( int j = 1; j <= this->dim; j++ ) {
			fscanf( fp, "%lf", &(this->meanMat[i][j]) );
		}
		fscanf( fp, "\n" );
		fscanf( fp, "cov:\n" );
		for( int j = 1; j <= this->dim; j++ ) {
			for( int k = 1; k <= this->dim; k++ ) {
				fscanf( fp, "%lf", &(this->covMat[i][j][k]) );
			}
			fscanf( fp, "\n" );
		}
	}

	fscanf( fp, "pi:\n" );
	this->pi = (double *) calloc( (unsigned) ( this->N + 1 ), sizeof( double ) );
	for( int i = 1; i <= this->N; i++ ) {
		fscanf( fp, "%lf", &( this->pi[i] ) );
	}
}

void HMM::printHMM() {
	cout << "M= " << this->M << endl;
	cout << "N= " << this->N << endl;

	cout << "A:" << endl;
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			printf( "%.12g ", this->A[i][j] );
		}
		cout << endl;
	}
	cout << "B:" << endl;
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->M; j++ ) {
			printf( "%.12g ", this->B[i][j] );
		}
		cout << endl;
	}
	cout << "pi:" << endl;
	for( int i = 1; i <= this->N; i++ ) {
		printf( "%.12g ", this->pi[i] );
	}
	cout << endl;
}

void HMM::printHMM( char *filename ) {
	FILE *fp = fopen( filename, "w" );
	if( fp == NULL ) {
		cerr << "Error: Cannot create file " << filename << endl;
		exit( 1 );
	}

	fprintf( fp, "M= %d\n", this->M );
	fprintf( fp, "N= %d\n", this->N );

	fprintf( fp, "A:\n" );
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			fprintf( fp, "%.12g ", this->A[i][j] );
		}
		fprintf( fp, "\n" );
	}

	fprintf( fp, "B:\n" );
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->M; j++ ) {
			fprintf( fp, "%.12g ", this->B[i][j] );
		}
		fprintf( fp, "\n" );
	}

	fprintf( fp, "pi:\n" );
	for( int i = 1; i <= this->N; i++ ) {
		fprintf( fp, "%.12g ", this->pi[i] );
	}
	fprintf( fp, "\n" );
	
	if( fp ) {
		fclose( fp );
		fp = NULL;
	}
}

void HMM::printHMM( FILE *fp ) {
	fprintf( fp, "M= %d\n", this->M );
	fprintf( fp, "N= %d\n", this->N );

	fprintf( fp, "A:\n" );
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			fprintf( fp, "%.12g ", this->A[i][j] );
		}
		fprintf( fp, "\n" );
	}

	fprintf( fp, "B:\n" );
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->M; j++ ) {
			fprintf( fp, "%.12g ", this->B[i][j] );
		}
		fprintf( fp, "\n" );
	}

	fprintf( fp, "pi:\n" );
	for( int i = 1; i <= this->N; i++ ) {
		fprintf( fp, "%.12g ", this->pi[i] );
	}
	fprintf( fp, "\n" );
}

void HMM::printSuffStats() {
	cout << "numeratorA_t:" << endl;
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			cout << this->numeratorA_t[i][j] << " ";
		}
		cout << endl;
	}
	cout << "denominatorA_t:" << endl;
	for( int i = 1; i <= this->N; i++ ) {
		cout << this->denominatorA_t[i] << " ";
	}
	cout << endl;
	cout << "numeratorB_t:" << endl;
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->M; j++ ) {
			cout << this->numeratorB_t[i][j] << " ";
		}
		cout << endl;
	}
	cout << "denominatorB_t:" << endl;
	for( int i = 1; i <= this->N; i++ ) {
		cout << this->denominatorB_t[i] << " ";
	}
	cout << endl;
	cout << "pi_t:" << endl;
	for( int i = 1; i <= this->N; i++ ) {
		cout << this->pi_t[i] << " ";
	}
	cout << endl;
	cout << "logprobfinal:" << endl;
	cout << this->logprobfinal << endl;
	cout << "NS:" << endl;
	cout << this->NS << endl;
}

void HMM::printHMMWithSuffStats( FILE *fp ) {
	fprintf( fp, "M= %d\n", this->M );
	fprintf( fp, "N= %d\n", this->N );

	fprintf( fp, "A:\n" );
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			fprintf( fp, "%.12g ", this->A[i][j] );
		}
		fprintf( fp, "\n" );
	}
	fprintf( fp, "B:\n" );
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->M; j++ ) {
			fprintf( fp, "%.12g ", this->B[i][j] );
		}
		fprintf( fp, "\n" );
	}
	fprintf( fp, "pi:\n" );
	for( int i = 1; i <= this->N; i++ ) {
		fprintf( fp, "%.12g ", this->pi[i] );
	}
	fprintf( fp, "\n" );
	
	fprintf( fp, "numeratorA_t:\n" );
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			fprintf( fp, "%.12g ", this->numeratorA_t[i][j] );
		}
		fprintf( fp, "\n" );
	}
	fprintf( fp, "denominatorA_t:\n" );
	for( int i = 1; i <= this->N; i++ ) {
		fprintf( fp, "%.12g ", this->denominatorA_t[i] );
	}
	fprintf( fp, "\n" );
	fprintf( fp, "numeratorB_t:\n" );
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->M; j++ ) {
			fprintf( fp, "%.12g ", this->numeratorB_t[i][j] );
		}
		fprintf( fp, "\n" );
	}
	fprintf( fp, "denominatorB_t:\n" );
	for( int i = 1; i <= this->N; i++ ) {
		fprintf( fp, "%.12g ", this->denominatorB_t[i] );
	}
	fprintf( fp, "\n" );
	fprintf( fp, "pi_t:\n" );
	for( int i = 1; i <= this->N; i++ ) {
		fprintf( fp, "%.12g ", this->pi_t[i] );
	}
	fprintf( fp, "\n" );
	fprintf( fp, "logprobfinal:\n" );
	fprintf( fp, "%.12g\n", this->logprobfinal );
	fprintf( fp, "NS:\n" );
	fprintf( fp, "%d\n", this->NS );
}

// for continuous HMMs
void HMM::printContHMM() {
	cout << "N= " << this->N << endl;
	cout << "dim= " << this->dim << endl;

	cout << "A:" << endl;
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			printf( "%.12g ", this->A[i][j] );
		}
		cout << endl;
	}

	cout << "B:" << endl;
	for( int i = 1; i <= this->N; i++ ) {
		cout << "mean:" << endl;
		for( int j = 1; j <= this->dim; j++ ) {
			printf( "%.12g ", this->meanMat[i][j] );
		}
		cout << endl;
		cout << "cov:" << endl;
		for( int j = 1; j <= this->dim; j++ ) {
			for( int k = 1; k <= this->dim; k++ ) {
				printf( "%.12g ", this->covMat[i][j][k] );
			}
			cout << endl;
		}
	}
	
	cout << "pi:" << endl;
	for( int i = 1; i <= this->N; i++ ) {
		printf( "%.12g ", this->pi[i] );
	}
	cout << endl;
}

// initialize a HMM at random ( type 0: normal, type 1: left-to-right )
void HMM::initHMM( int seed, int N, int M, int type ) {
	double sum;

	srand( seed );

	this->M = M;
	this->N = N;

	this->A = (double **) calloc( (unsigned) ( this->N + 1 ), sizeof(double*) );
	for( int i = 1; i <= this->N; i++ ) {
		this->A[i] = (double *) calloc( (unsigned) ( this->N + 1 ), sizeof(double) );
	}
	this->B = (double **) calloc( (unsigned) ( this->N + 1 ), sizeof(double*) );
	for( int i = 1; i <= this->N; i++ ) {
		this->B[i] = (double *) calloc( (unsigned) ( this->M + 1 ), sizeof(double) );
	}
	this->pi = (double *) calloc( (unsigned) ( this->N + 1 ), sizeof( double ) );
	
	// set all of the values to 0
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			this->A[i][j] = 0;
		}
		for( int k = 1; k <= this->M; k++ ) {
			this->B[i][k] = 0;
		}
		this->pi[i] = 0;
	}
	
	if( type == 0 ) {
		for( int i = 1; i <= this->N; i++ ) {
			sum = 0.0;
			for( int j = 1; j <= this->N; j++ ) {
				this->A[i][j] = (double)rand()/RAND_MAX;
				sum += this->A[i][j];
			}
			for( int j = 1; j <= this->N; j++ ) {
				this->A[i][j] /= sum;
			}
		}
		for( int j = 1; j <= this->N; j++ ) {
			sum = 0.0;
			for( int k = 1; k <= this->M; k++ ) {
				this->B[j][k] = (double)rand()/RAND_MAX;
				sum += this->B[j][k];
			}
			for( int k = 1; k <= this->M; k++ ) {
				this->B[j][k] /= sum;
			}
		}
		sum = 0.0;
		for( int i = 1; i <= this->N; i++ ) {
			this->pi[i] = (double)rand()/RAND_MAX;
			sum += this->pi[i];
		}
		for( int i = 1; i <= this->N; i++ ) {
			this->pi[i] /= sum;
		}
	}
	else if( type == 1 ) {
		for( int i = 1; i <= this->N; i++ ) {
			sum = 0.0;
			for( int j = 1; j <= this->N; j++ ) {
				if( ( i > j ) || ( j - i >= 2 ) ) {
					continue;
				}
				this->A[i][j] = (double)rand()/RAND_MAX;
				sum += this->A[i][j];
			}
			for( int j = 1; j <= this->N; j++ ) {
				this->A[i][j] /= sum;
			}
		}
		for( int j = 1; j <= this->N; j++ ) {
			sum = 0.0;
			for( int k = 1; k <= this->M; k++ ) {
				this->B[j][k] = (double)rand()/RAND_MAX;
				sum += this->B[j][k];
			}
			for( int k = 1; k <= this->M; k++ ) {
				this->B[j][k] /= sum;
			}
		}
		this->pi[1] = 1;
		for( int i = 2; i <= this->N; i++ ) {
			this->pi[i] = 0;
		}
/*
		sum = 0.0;
		for( int i = 1; i <= this->N; i++ ) {
			this->pi[i] = (double)rand()/RAND_MAX;
			sum += this->pi[i];
		}
		for( int i = 1; i <= this->N; i++ ) {
			this->pi[i] /= sum;
		}
//*/
	}
}

void HMM::readSequence( FILE *fp, int *pT, int **pO ) {
	int *O;
	fscanf( fp, "T= %d\n", pT );
	O = ( int* ) calloc ( ( unsigned ) ( *pT + 1 ), sizeof( int ) );
	for( int i = 1; i <= *pT; i++ ) {
		fscanf( fp, "%d", &O[i] );
	}
	*pO = O;
}

void HMM::readContSequence( FILE *fp, int *pT, double ***pO ) {
	double **O;
	fscanf( fp, "T= %d\n", pT );
	O = ( double** ) calloc( ( unsigned ) ( *pT + 1 ), sizeof( double* ) );
	for( int i = 1; i <= *pT; i++ ) {
		O[i] = ( double* ) calloc ( ( unsigned ) ( this->dim + 1 ), sizeof( double ) );
	}
	for( int i = 1; i <= *pT; i++ ) {
		for( int j = 1; j <= this->dim; j++ ) {
			fscanf( fp, "%lf", &O[i][j] );
		}
		fscanf( fp, "\n" );
	}
	*pO = O;
}

void HMM::printSequence( FILE *fp, int T, int *O ) {
	fprintf( fp, "T= %d\n", T );
	for( int i = 1; i <= T; i++ ) {
		fprintf( fp,"%d ", O[i] );
	}
	fprintf( fp, "\n" );
}

void HMM::printContSequence( FILE *fp, int T, double **O ) {
	fprintf( fp, "T= %d\n", T );
	for( int i = 1; i <= T; i++ ) {
		for( int j = 1; j <= this->dim; j++ ) {
			fprintf( fp, "%lf ", O[i][j] );
		}
		fprintf( fp, "\n" );
	}
}

int HMM::genInitialState() {
	double val, accum;
	int i, q_t;
 
	val = (double)rand()/RAND_MAX;

//	cout << "val: " << val << endl;
//	cout << "RAND_MAX: " << RAND_MAX << endl;

	accum = 0.0;
	q_t = this->N;	/* number of states */
	for( i = 1; i <= this->N; i++ ) {
		if( val < this->pi[i] + accum ) {
    	q_t = i;
    	break;
    }
  	else {
	    accum += this->pi[i];
  	}
	}
 
	return q_t;
}

int HMM::genNextState(int q_t) {
	double val, accum;
	int j, q_next;

	val = ( double )rand()/RAND_MAX;
	accum = 0.0;
	q_next = this->N;
	for( j = 1; j <= this->N; j++ ) {
		if( val < this->A[q_t][j] + accum ) {
			q_next = j;
			break;
		}
		else {
			accum += this->A[q_t][j];
		}
	}

	return q_next;
}

int HMM::genSymbol(int q_t) {
	double val, accum;
	int j, o_t;

	val = ( double )rand()/RAND_MAX;
	accum = 0.0;
	o_t = this->M;
	for (j = 1; j <= this->M; j++) {
		if ( val < this->B[q_t][j] + accum ) {
			o_t = j;
			break;
		}
		else {
			accum += this->B[q_t][j];
		}
	}
	return o_t;
}

void HMM::genSequenceArray( int seed, int T, int *O, int *q ) {
	int t = 1;
//	int q_t, o_t;

	srand( seed );
 
	q[1] = genInitialState();
	O[1] = genSymbol( q[1] );
 
	for( t = 2; t <= T; t++ ) {
		q[t] = genNextState( q[t-1] );
		O[t] = genSymbol( q[t] );
	}
}

double HMM::computeGaussian( int state, double *X ) {
	double prob;
/*
	for( int i = 1; i <= this->dim; i++ ) {
		cout << X[i] << " ";
	}
	cout << endl;
//*/
	Matrix m_mean_T( 1, this->dim );
	for( int i = 0; i < this->dim; i++ ) {
		m_mean_T( 0, i ) = X[i + 1] - this->meanMat[state][i + 1];
	}
	Matrix m_mean = ~m_mean_T;

	Matrix m_cov( this->dim, this->dim );
	for( int i = 0; i < this->dim; i++ ) {
		for( int j = 0; j < this->dim; j++ ) {
			m_cov( i, j ) = this->covMat[state][i + 1][j + 1];
		}
	}

	// check if the cov matrix is invertable. if not, do SVD of cov matrix
/*
	double var = 0;
	double alpha = 0.01;
	double cond = m_cov( 0, 0 )/m_cov( this->dim - 1, this->dim - 1 );
	if( cond > alpha ) {
		var = ( ( alpha * m_cov( this->dim - 1, this->dim - 1 ) ) - m_cov( 0, 0 ) )/( 1 - alpha );
		for( int i = 0; i < this->dim; i++ ) {
			for( int j = 0; j < this->dim; j++ ) {
				if( i == j ) {
					m_cov( i, j ) += var;
				}
			}
		}
		//exit(1);
	}
//*/
	Matrix m_cov_inv = !m_cov;
	double det = m_cov.Det();

//	cout << m_cov << endl;
//	cout << m_mean_T << endl;
//	cout << m_mean << endl;
//	cout << m_cov << endl;
//	cout << m_cov_inv << endl;
//	cout << m_mean_T * m_cov_inv << endl;
//	cout << ( m_mean_T * m_cov_inv ) * m_mean << endl;

	Matrix tmp( 1, 1 );
	tmp = ( m_mean_T * m_cov_inv ) * m_mean;

	double e = exp( -0.5 * tmp( 0, 0 ) );
//*
	if( e == 0 ) {
		e = 1e-10;
	}
	else if( e > 1e+10 ) {
		e = 1e+10;
	}
//*/	
	prob = 1/sqrt( pow( 2 * PI, this->dim ) * abs( det ) ) * e;
//	cout << "det: " << abs( det ) << endl;
//	cout << "prob: " << prob << endl;
//	cout << pow( 2 * PI, this->dim ) * abs( det ) << endl;
//	cout << sqrt( pow( 2 * PI, this->dim ) * abs( det ) ) << endl;
	//cout << "exp: " << exp( -0.5 * tmp( 0, 0 ) ) << endl;
	//cout << "tmp: " << tmp( 0, 0 ) << endl;

//	exit( 1 );

	return prob;
}

// ---------------- Used in training ----------------
void HMM::computeGamma( int T, double **alpha, double **beta, double **gamma ) {
	double denominator;
	for( int t = 1; t <= T; t++ ) {
		denominator = 0.0;
		for( int j = 1; j <= this->N; j++ ) {
			gamma[t][j] = alpha[t][j] * beta[t][j];
		//	cout << "beta: " << beta[t][j] << " ";
		//	cout << "alpha: " << alpha[t][j] << " ";
			denominator += gamma[t][j];
		//	cout << j << "-denom: " << denominator << endl;
		}
//		cout << endl;
		for( int i = 1; i <= this->N; i++ ) {
		//	cout << i << "-denom: " << denominator << endl;
		//	if( denominator == 0 ) exit( 1 );
			gamma[t][i] = gamma[t][i]/denominator;
		}
	}
/*
	cout << "Gamma while computing" << endl;
	for( int t = 1; t <= T; t++ ) {
		for( int i = 1; i <= this->N; i++ ) {
			cout << gamma[t][i] << " ";
		}
		cout << endl;
	}
//*/
}

void HMM::computeXi( int T, int *O, double **alpha, double **beta, double ***xi ) {
	double sum;
	for( int t = 1; t <= T - 1; t++ ) {
		sum = 0.0;
		for( int i = 1; i <= this->N; i++ ) {
			for( int j = 1; j <= this->N; j++ ) {
				xi[t][i][j] = alpha[t][i] * beta[t+1][j] * (this->A[i][j]) * (this->B[j][O[t+1]]);
				sum += xi[t][i][j];
			}
		}
		for( int i = 1; i <= this->N; i++ ) {
			for( int j = 1; j <= this->N; j++ ) {
				if( sum == 0 ) {
					sum = 1e-80;
				}
				xi[t][i][j] /= sum;
			}
		}
	}
}

double ***HMM::allocXi( int T, int N ) {
	double ***xi;
	xi = (double ***)malloc( ( T + 5 ) * sizeof(double **) );
	xi--;
	for( int t = 1; t <= T; t++ ) {
		xi[t] = (double **) calloc( (unsigned) ( N + 5 ), sizeof(double*) );
		for( int i = 1; i <= N; i++ ) {
			xi[t][i] = (double *) calloc( (unsigned) ( N + 5 ), sizeof(double) );
		}
	}
	for( int t = 1; t <= T; t++ ) {
		for( int i = 1; i <= N; i++ ) {
			for( int j = 1; j <= N; j++ ) {
				xi[t][i][j] = 0;
			}
		}
	}
	return xi;
}


















