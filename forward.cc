/*
 *	Author: Kan Ouivirach
 *	File: forward.cc
 *	Purpose: Implement HMM class: Forward algorithm
 *	Credit: Tapas Kanungo, kanungo@cfar.umd.edu
 */

#include "hmm.h"

void HMM::forward( int T, int *O, double **alpha ) {
	double sum;
/*
	cout << "Observation: ";
	for( int i = 1; i <= T; i++ ) {
		cout << O[i] << " ";
	}
	cout << endl;
*/
	// 1. Initialization 
	for( int i = 1; i <= this->N; i++ ) {
		alpha[1][i] = this->pi[i] * this->B[i][O[1]];
//		cout << alpha[1][i] << " ";
	}
	// 2. Induction 
	for( int t = 1; t < T; t++ ) {
		for( int j = 1; j <= this->N; j++ ) {
			sum = 0.0;
			for( int i = 1; i <= this->N; i++ ) {
				sum += alpha[t][i] * this->A[i][j];
			}
			alpha[t+1][j] = sum * this->B[j][O[t+1]];
//			cout << alpha[t+1][j] << " ";
		}
	}
	// 3. Termination 
	this->loglik = 0.0;
	for( int i = 1; i <= this->N; i++ ) {
		this->loglik += alpha[T][i];
	}
	this->loglik = log( this->loglik );
/*	
	cout << "Alpha (no scaling):" << endl;
	for( int t = 1; t <= T; t++ ) {
		for( int i = 1; i <= this->N; i++ ) {
			cout << alpha[t][i] << " ";
		}
		cout << endl;
	}
	cout << endl;
//*/
}

void HMM::forwardWithScale( int T, int *O, double **alpha, double *scale, int type = 0 ) {
	if( type == 0 ) {
		double sum;
	/*
		cout << "Observation: ";
		for( int i = 1; i <= T; i++ ) {
			cout << O[i] << " ";
		}
		cout << endl;
	//*/
		// 1. Initialization 
		scale[1] = 0.0;
		for( int i = 1; i <= this->N; i++ ) {
	//		cout << "T: " << T << endl;
	//		cout << "pi: " << this->pi[i] << ", B: " << this->B[i][O[1]] << endl;
			alpha[1][i] = this->pi[i] * this->B[i][O[1]];
	//		if( alpha[1][i] == 0 ) {
	//			alpha[1][i] = 1e-80;
	//		}
			scale[1] += alpha[1][i];
	//		cout << alpha[1][i] << " ";
	//		cout << this->pi[i] << " x " << this->B[i][O[1]] << ", ";
	//		cout << this->pi[i] << "-";
	//		cout << scale[1] << " ";
	//		cout << "Obs: " << O[1] << endl;
		}
	//	if( scale[1] == 0 ) {
	//		scale[1] = 1e-80;
	//	}
		// Normalize the alpha
		//
		for( int i = 1; i <= this->N; i++ ) {
			alpha[1][i] /= scale[1];
		//	cout << alpha[1][i] << " ";
		}
		// 2. Induction 
		for( int t = 1; t < T; t++ ) {
			scale[t+1] = 0.0;
			for( int j = 1; j <= this->N; j++ ) {
				sum = 0.0;
				for( int i = 1; i <= this->N; i++ ) {
					sum += alpha[t][i] * this->A[i][j];
				}
				alpha[t+1][j] = sum * this->B[j][O[t+1]];
			//	cout << alpha[t+1][j] << " ";
				scale[t+1] += alpha[t+1][j];
			//	if( alpha[t+1][j] == 0 ) {
			//		alpha[t+1][j] = 1e-80;
			//	}
			}
		//	cout << endl;
		//	if( scale[t+1] == 0 ) {
		//		scale[t+1] = 1e-80;
		//	}
			for( int j = 1; j <= this->N; j++ ) {
				alpha[t+1][j] /= scale[t+1];
		//		cout << "scale: " << scale[t+1] << " ";
			}
	//		cout << endl;
		}
		// 3. Termination 
		this->loglik = 0.0;
		for( int t = 1; t <= T; t++ ) {
			if( scale[t] != 0 ) {
				this->loglik += log( scale[t] );
			//	cout << "log scale: " << log( scale[t] ) << " ";
			}
		//	cout << endl;
		//	cout << "LOG LIK: " << this->loglik << " ";
		}
	/*
		cout << "Alpha (with scaling):" << endl;
		for( int t = 1; t <= T; t++ ) {
			for( int i = 1; i <= this->N; i++ ) {
				cout << alpha[t][i] << " ";
			}
			cout << endl;
		}
		cout << endl;
	//*/
	}
	else {
		double sum;
		// 1. Initialization 
		scale[1] = 0.0;
		for( int i = 1; i <= this->N; i++ ) {
			alpha[1][i] = this->pi[i] * this->B[i][O[1]];
			scale[1] += alpha[1][i];
		}
		// Normalize the alpha
		for( int i = 1; i <= this->N; i++ ) {
			alpha[1][i] /= scale[1];
		}
		// 2. Induction 
		for( int t = 1; t < T; t++ ) {
			scale[t+1] = 0.0;
			for( int j = 1; j <= this->N; j++ ) {
				sum = 0.0;
				for( int i = 1; i <= this->N; i++ ) {
					if( j < i ) {
						continue;
					}
					sum += alpha[t][i] * this->A[i][j];
				}
				alpha[t+1][j] = sum * this->B[j][O[t+1]];
				scale[t+1] += alpha[t+1][j];
			}
			for( int j = 1; j <= this->N; j++ ) {
				alpha[t+1][j] /= scale[t+1];
			}
		}
		// 3. Termination 
		this->loglik = 0.0;
		for( int t = 1; t <= T; t++ ) {
			if( scale[t] != 0 ) {
				this->loglik += log( scale[t] );
			}
		}
	}
}

void HMM::forwardWithScale( int T, vector<int> O, double **alpha, double *scale ) {
	double sum;
/*
	cout << "Observation: ";
	for( int i = 1; i <= T; i++ ) {
		cout << O[i] << " ";
	}
	cout << endl;
*/
	// 1. Initialization 
	scale[1] = 0.0;
	for( int i = 1; i <= this->N; i++ ) {
//		cout << "pi: " << this->pi[i] << ", B: " << this->B[i][O[1]] << endl;
//cout << "before" << endl;
		alpha[1][i] = this->pi[i] * this->B[i][O[1]];
//cout << "after" << endl;
		scale[1] += alpha[1][i];
//		cout << alpha[1][i] << " ";
	}
//	cout << "???" << endl;
	if( scale[1] == 0 ) {
		scale[1] = 1e-80;
	}
	for( int i = 1; i <= this->N; i++ ) {
		alpha[1][i] /= scale[1];
	}
//	cout << "??" << endl;
	// 2. Induction 
	for( int t = 1; t < T; t++ ) {
		scale[t+1] = 0.0;
		for( int j = 1; j <= this->N; j++ ) {
			sum = 0.0;
			for( int i = 1; i <= this->N; i++ ) {
				sum += alpha[t][i] * this->A[i][j];
		//		cout << "N: " << this->N << endl;
		//		cout << "alpha: " << alpha[t][i] << endl;
			}
//			cout << "sum: " << sum;
			alpha[t+1][j] = sum * this->B[j][O[t+1]];
			scale[t+1] += alpha[t+1][j];
		}
//		cout << endl;
		if( scale[t+1] == 0 ) {
			scale[t+1] = 1e-80;
		}
		for( int j = 1; j <= this->N; j++ ) {
			alpha[t+1][j] /= scale[t+1];
//			cout << "scale: " << scale[t+1] << " ";
		}
//		cout << endl;
	}
	// 3. Termination 
	this->loglik = 0.0;
	for( int t = 1; t <= T; t++ ) {
		if( scale[t] != 0 ) {
			this->loglik += log( scale[t] );
		}
//		cout << this->prob << " ";
	}
}

void HMM::forwardContWithScale( int T, double **O, double **alpha, double *scale ) {
	double sum;
	// 1. Initialization 
	scale[1] = 0.0;
	for( int i = 1; i <= this->N; i++ ) {
//		alpha[1][i] = this->pi[i] * this->B[i][O[1]];
		alpha[1][i] = this->pi[i] * computeGaussian( i, O[1] );
		scale[1] += alpha[1][i];
//		cout << "alpha: " << alpha[1][i] << endl;
//		cout << "gaussian: " << computeGaussian( i, O[1] ) << endl;
	}
	if( scale[1] == 0 ) {
		scale[1] = 1e-80;
	}
	for( int i = 1; i <= this->N; i++ ) {
		alpha[1][i] /= scale[1];
	}
	// 2. Induction 
	for( int t = 1; t < T; t++ ) {
		scale[t+1] = 0.0;
		for( int j = 1; j <= this->N; j++ ) {
			sum = 0.0;
			for( int i = 1; i <= this->N; i++ ) {
				sum += alpha[t][i] * this->A[i][j];
			}
//			alpha[t+1][j] = sum * this->B[j][O[t+1]];
			alpha[t+1][j] = sum * computeGaussian( j, O[t+1] );
			scale[t+1] += alpha[t+1][j];
		//	cout << computeGaussian( j, O[t+1] ) << " ";
		}
		if( scale[t+1] == 0 ) {
			scale[t+1] = 1e-80;
		}
		for( int j = 1; j <= this->N; j++ ) {
			alpha[t+1][j] /= scale[t+1];
		}
	}
	// 3. Termination 
	this->loglik = 0.0;
	for( int t = 1; t <= T; t++ ) {
		if( scale[t] != 0 ) {
			this->loglik += log( scale[t] );
		}
	}

//	for( int t = 1; t <= T; t++ ) {
//		cout << scale[t] << " ";
//	}
//	cout << endl;
/*
	for( int t = 1; t <= T; t++ ) {
	//	cout << scale[t] << " : ";
		for( int i = 1; i <= this->N; i++ ) {
			cout << alpha[t][i] << " ";
//			if( alpha[t][i] == 0 ) {
//				alpha[t][i] = 1e-30;
//			}
		}
		cout << endl;
	}
//*/

}




