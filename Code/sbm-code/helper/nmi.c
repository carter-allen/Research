#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Random.h>
#include <R_ext/Rdynload.h>

#define LOGEPS (-37) /* roughly log(2^(-53)), log of machine precision */

double entropy (double *p, int K) {
  double ent = 0.0;
  int i;
	for (i=0; i < K; i++){
    if (p[i] != 0.0) ent += (p[i] * log(p[i]) );
	}
  return(-ent);
}

static double log1pe (double x) {
  double d = (x > 0) ? -x : x;
  d = (d < LOGEPS) ? 0 : log1p(exp(d));
  if (x > 0) d += x;
  return d;
}

int gindex (int a, int b, int K){
  int l, t, k2 = K * (K - 1) / 2;
  if (b < a) { 
    t = b; b = a; a = t;
  }
  if(a==b) {l = k2 + a;} else {l = a * (K - 1) - a * (a + 1) / 2 + b - 1;}
  return(l) ;
}

SEXP LHOOD (SEXP Am, SEXP sigmav, SEXP gammav, SEXP etav, SEXP Zv, SEXP ns){
  SEXP lhoods;

  double xi;
  int  i, j, n;


  double * gamma;
  double * eta;
  double * lhood;
  int * A;
  int * sigma;
  int * Z;

  A = INTEGER(Am);
  sigma = INTEGER(sigmav);
  gamma = REAL(gammav);
  eta = REAL(etav);
  Z = INTEGER(Zv);
  n = INTEGER(ns)[0];

  PROTECT(lhoods = allocVector(REALSXP, 1));
  lhood = REAL(lhoods);

  *lhood = 0.0;

  for (i = 0; i < (n-1); i++){
		for (j = (i + 1); j < n; j++){
      xi = (sigma[i]==sigma[j]) * gamma[(sigma[i] - 1)] + eta[(Z[i] - 1)] + eta[(Z[j] - 1)];
      *lhood += ( A[j*n + i] * xi - log1pe(xi));
    }
	}

	UNPROTECT(1);
	return(lhoods);
}

SEXP LHOOD_E (SEXP Em, SEXP sigmav, SEXP gammav, SEXP etav, SEXP Zv, SEXP ns, SEXP es){
  SEXP lhoods;

  double xi;
  int  i, j, k, n, e;

  double * gamma;
  double * eta;
  double * lhood;
  int * E;
  int * sigma;
  int * Z;

  E = INTEGER(Em);
  sigma = INTEGER(sigmav);
  gamma = REAL(gammav);
  eta = REAL(etav);
  Z = INTEGER(Zv);
  n = INTEGER(ns)[0];
  e = INTEGER(es)[0];

  PROTECT(lhoods = allocVector(REALSXP, 1));
  lhood = REAL(lhoods);

  *lhood = 0.0;

  for (i = 0; i < (n-1); i++){
    for (j = (i + 1); j < n; j++){
      xi = (sigma[i]==sigma[j]) * gamma[(sigma[i] - 1)] + eta[(Z[i] - 1)] + eta[(Z[j] - 1)];
      *lhood -= log1pe(xi);
    }
	}

  for (k = 0; k < e; k++){
      i = E[k] - 1;
      j = E[e + k] - 1;
      xi = (sigma[i]==sigma[j]) * gamma[(sigma[i] - 1)] + eta[(Z[i] - 1)] + eta[(Z[j] - 1)];
      *lhood += xi ;
  }

  UNPROTECT(1);

  return(lhoods);
}

SEXP LHOOD_G (SEXP Em, SEXP sigmav, SEXP gammav, SEXP etav, SEXP Zv, SEXP ns, SEXP es, SEXP Ks, SEXP Ps){
  SEXP lhoods;

  double xi;
  int  i, j, k, n, e, m, K, P, ind1, ind2;

  double * gamma;
  double * eta;
  double * lhood;
  int * E;
  int * sigma;
  int * Z;

  E = INTEGER(Em);
  sigma = INTEGER(sigmav);
  gamma = REAL(gammav);
  eta = REAL(etav);
  Z = INTEGER(Zv);
  n = INTEGER(ns)[0];
  e = INTEGER(es)[0];
  K = INTEGER(Ks)[0];
  P = INTEGER(Ps)[0];

  int p1 = P * (P - 1) / 2;
  int p2 = p1 + P;

  double * lwz = malloc(K * p2 * sizeof(double));
  double * lz = malloc(p2 * sizeof(double));

  int * nw = malloc(K * p2 * sizeof(int));
  int * nb = malloc(p2 * sizeof(int));

  PROTECT(lhoods = allocVector(REALSXP, 1));
  lhood = REAL(lhoods);

  *lhood = 0.0;

	for (i =0; i < (K * p2); i++)	nw[i] = 0;
	for (i =0; i < p2; i++)	nb[i] = 0;

/*
  for(j = 0; j < (P-1), j++){
    for(k = j; k < (P-1), k++){
      for(i = 0; i < K, i++){
        lwz[i * p2 + gindex(j,k,P)] = log1pe(gamma[i] + eta[j] + eta[k]);
      }
      lz[gindex(j,k,P)] = log1pe(eta[j] + eta[k]);
    }
  }
*/

//int count = 0;

  for (i = 0; i < (n-1); i++){
		for (j = (i + 1); j < n; j++){
      if(sigma[i]==sigma[j]) {nw[(sigma[i]) * p2 + gindex(Z[i],Z[j],P)] ++;} else {nb[gindex(Z[i],Z[j],P)] ++;}
		}
	}

  for(j = 0; j < P; j++){
    for(k = j; k < P; k++){
      ind2 = gindex(j,k,P);
      for(i = 0; i < K; i++){
        ind1 = i * p2 + ind2;
        lwz[ind1] = log1pe(gamma[i] + eta[j] + eta[k]);
        *lhood -= (lwz[ind1] * nw[ind1]); 
//count += nw[ind1];
      }
      lz[ind2] = log1pe(eta[j] + eta[k]);
      *lhood -= (lz[ind2] * nb[ind2]); 
//count += nb[ind2];
    }
  }

//printf("wtf: %d", count);

  for (k = 0; k < e; k++){
      i = E[k] - 1;
      j = E[e + k] - 1;
      xi = (sigma[i]==sigma[j]) * gamma[sigma[i]] + eta[Z[i]] + eta[Z[j]];
      *lhood += xi ;
  }

  UNPROTECT(1);
  free(lwz);
  free(lz);
  free(nw);
  free(nb);

  return(lhoods);
}



SEXP NMI0 (SEXP refv, SEXP estv, SEXP ns, SEXP Ks) {

	SEXP nmis;
	double * nmi, ent_est, ent_ref, mi = 0.0;

	int i, j, k, l, n, K, count=0;
	int * ref;
	int * est;

	ref = INTEGER(refv);
	est = INTEGER(estv);
	n = INTEGER(ns)[0];
	K = INTEGER(Ks)[0];

	double * pest = malloc(K * sizeof(double));
	double * pref = malloc(K * sizeof(double));
  double * pp = malloc(K * K * sizeof(double));

	PROTECT(nmis = allocVector(REALSXP, 1));
	nmi = REAL(nmis);

	*nmi = 0.0;

// compute p
  for (j = 0; j < K; j++) pest[j] = 0;

  for (i = 0; i < n; i++) { if (est[i] <= K) (* (pest + est[i] - 1))++;}
 
  for (j = 0; j < K; j++) {pest[j] = pest[j] / n;}// printf("%lf ", pest[j]);}

  for (j = 0; j < K; j++) pref[j] = 0;

  for (i = 0; i < n; i++)   (* (pref + ref[i] - 1))++;
 
  for (j = 0; j < K; j++) {pref[j] = pref[j] / n;}// printf("%lf ", pest[j]);}
/*
  for (j = 0; j < K; j++){
    pref[j] = 0;
    for (i = 0; i < n; i++){
			pref[j] += *(ref + j*n + i);
		}
		pref[j] = pref[j] / n;
		// printf("%lf ", pref[j]);
  }
*/
	ent_est = entropy(pest, K);
	ent_ref = entropy(pref, K);

  //printf("\n est: %lf, ref: %lf \n", ent_est, ent_ref);

// compute pp
  for (j=0; j<K; j++){
		for (k=0; k<K; k++){
			count = 0;
			for (i=0; i<n; i++){
				if ((est[i] == j + 1) & (ref[i] == k + 1)) count ++; 
			}			
			pp[j * K + k] = (double)(count) / n;
			// printf("p[ %d, %d ]= %lf, ", j, k, pp[j * K + k]);
    } 
	}

// compute MI
  for (j=0; j<K; j++){
		for (k=0; k<K; k++){
			if ( (pp[j*K + k]  * pest[j] * pref[k]) != 0.0) mi += pp[j*K + k] * (log(pp[j*K + k]) - log(pest[j]) - log(pref[k]));
		}
	}

// compute NMI
  *nmi = 2 * mi / (ent_est + ent_ref);


	UNPROTECT(1);

  free(pref);
  free(pest);
  free(pp);
	return(nmis);
}




SEXP NMI (SEXP refm, SEXP estv, SEXP ns, SEXP Ks) {

	SEXP nmis;
	double * nmi, ent_est, ent_ref, mi = 0.0;

	int i, j, k, l, n, K, count=0;
	int * ref;
	int * est;

	ref = INTEGER(refm);
	est = INTEGER(estv);
	n = INTEGER(ns)[0];
	K = INTEGER(Ks)[0];

	double * pest = malloc(K * sizeof(double));
	double * pref = malloc(K * sizeof(double));
  double * pp = malloc(K * K * sizeof(double));

	PROTECT(nmis = allocVector(REALSXP, 1));
	nmi = REAL(nmis);

	*nmi = 0.0;

// compute p
  for (j = 0; j < K; j++) pest[j] = 0;

	for (i = 0; i < n; i++)  { if (est[i] <= K) (* (pest + est[i] - 1))++;}
 
  for (j = 0; j < K; j++) {pest[j] = pest[j] / n;}// printf("%lf ", pest[j]);}

  for (j = 0; j < K; j++){
    pref[j] = 0;
    for (i = 0; i < n; i++){
			pref[j] += *(ref + j*n + i);
		}
		pref[j] = pref[j] / n;
		// printf("%lf ", pref[j]);
  }

	ent_est = entropy(pest, K);
	ent_ref = entropy(pref, K);

  //printf("\n est: %lf, ref: %lf \n", ent_est, ent_ref);

// compute pp
  for (j=0; j<K; j++){
		for (k=0; k<K; k++){
			count = 0;
			for (i=0; i<n; i++){
				if (est[i] == j + 1) count += ref[k*n + i]; 
			}			
			pp[j * K + k] = (double)(count) / n;
			// printf("p[ %d, %d ]= %lf, ", j, k, pp[j * K + k]);
    } 
	}

// compute MI
  for (j=0; j<K; j++){
		for (k=0; k<K; k++){
			if ( (pp[j*K + k]  * pest[j] * pref[k]) != 0.0) mi += pp[j*K + k] * (log(pp[j*K + k]) - log(pest[j]) - log(pref[k]));
		}
	}

// compute NMI
  *nmi = 2 * mi / (ent_est + ent_ref);


	UNPROTECT(1);

  free(pref);
  free(pest);
  free(pp);
	return(nmis);
}






SEXP ADJRAND (SEXP refm, SEXP estv, SEXP ns, SEXP Ks) {

	SEXP RandIndexs;
	double * RandIndex;

	int i, j, k, l, n, K, count=0;
	int * ref;
	int * est;
  int x = 0, y = 0, z = 0, n2;

	ref = INTEGER(refm);
	est = INTEGER(estv);
	n = INTEGER(ns)[0];
	K = INTEGER(Ks)[0];
	n2 = n * (n - 1) / 2;

  int * nn = malloc(K * K * sizeof(int));
  int * a = malloc(K * sizeof(int));
  int * b = malloc(K * sizeof(int));

	PROTECT(RandIndexs = allocVector(REALSXP, 1));
	RandIndex = REAL(RandIndexs);
  * RandIndex = 0.0;
// compute nn
  for (j=0; j<K; j++){
		for (k=0; k<K; k++){
			count = 0;
			for (i=0; i<n; i++){
				if (est[i] == j + 1) count += ref[k*n + i]; 
			}			
			nn[j * K + k] = count;
			x += (nn[j * K + k] * (nn[j * K + k] - 1) / 2); 
			//printf("nn[ %d, %d ]= %d, ", j, k, nn[j * K + k]);
    } 
	}

  for (j=0; j<K; j++){
    a[j] = 0; b[j] = 0;
		for (k=0; k<K; k++){
			a[j] += nn[j * K + k];		
		  b[j] += nn[k * K + j];	
    } 
		y += (a[j] * (a[j] - 1) / 2); 
		z += (b[j] * (b[j] - 1) / 2); 
		//printf("a [%d] = %d, b [%d] = %d ], ", j, a[j], j, b[j]);
	}


// compute Rand Index
  *RandIndex = (x * n2 - y * z) / (n2 * (y + z) / 2.0 - y * z);


	UNPROTECT(1);

	free(a);
	free(b);
  free(nn);
	return(RandIndexs);
}

