/*
 *	Author: Kan Ouivirach
 *	File: viterbi.cc
 *	Purpose: Implement HMM class: Viterbi algorithm
 *	Credit: Tapas Kanungo, kanungo@cfar.umd.edu
 */

#include "hmm.h"

#define VITHUGE 100000000000.0

void HMM::viterbi( int T, int *O ) {
	int maxvalind;
	double maxval, val;

	double **delta;
	int **psi;
	int *q;

	double **biot;

	cout << "Observation: ";
	for( int i = 1; i <= T; i++ ) {
		cout << O[i] << " ";
	}
	cout << endl;

	q = (int *) calloc( (unsigned) ( T + 1 ), sizeof(int*) );

	delta = (double **) calloc( (unsigned) ( T + 1 ), sizeof(double*) );
	for( int i = 1; i <= T; i++ ) {
		delta[i] = (double *) calloc( (unsigned) ( this->N + 1 ), sizeof(double) );
	}

	psi = (int **) calloc( (unsigned) ( T + 1 ), sizeof(int*) );
	for( int i = 1; i <= T; i++ ) {
		psi[i] = (int *) calloc( (unsigned) ( this->N + 1 ), sizeof(int) );
	}
	
	biot = (double **) calloc( (unsigned) ( this->N + 1 ), sizeof(double*) );
	for( int i = 1; i <= this->N; i++ ) {
		biot[i] = (double *) calloc( (unsigned) ( T + 1 ), sizeof(double) );
	}

	/* 0. Preprocessing */
	for( int i = 1; i <= this->N; i++ ) {
		this->pi[i] = log( this->pi[i] );
	}
	for( int i = 1; i <= this->N; i++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			this->A[i][j] = log( this->A[i][j] );
		}
	}
	for( int i = 1; i <= this->N; i++ ) {
		for( int t = 1; t <= T; t++ ) {
			biot[i][t] = log( this->B[i][O[t]] );
		}
	}

	/* 1. Initialization */
	for( int i = 1; i <= this->N; i++ ) {
		delta[1][i] = this->pi[i] + biot[i][1];
		psi[1][i] = 0;
	}
	/* 2. Recursion */
	for( int t = 2; t <= T; t++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			maxval = -VITHUGE;
			maxvalind = 1;
			for( int i = 1; i <= this->N; i++ ) {
				val = delta[t-1][i] + this->A[i][j];
				if( val > maxval ) {
					maxval = val;
					maxvalind = i;
				}
			}
//			alpha[t+1][j] = sum * this->B[j][O[t+1]];
			delta[t][j] = maxval + biot[j][t];
			psi[t][j] = maxvalind;
		}
	}
	/* 3. Termination */
	this->loglik = -VITHUGE;
	q[T] = 1;
	for( int i = 1; i <= this->N; i++ ) {
//		this->prob += alpha[T][i];
		if( delta[T][i] > this->loglik ) {
			this->loglik = delta[T][i];
			q[T] = i;
		}
	}

	/* 4. Path (state sequence) backtracking */
	for( int t = T - 1; t >= 1; t-- ) {
		q[t] = psi[t+1][q[t+1]];
	}

	cout << "Best path: ";
	for( int t = 1; t <= T; t++ ) {
		cout << q[t] << " ";
	}
	cout << endl;

}
