#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Random.h>
#include <R_ext/Rdynload.h>

#include <math.h>
#include <stdlib.h>
#include <igraph.h>

#define LOGEPS (-37) /* roughly log(2^(-53)), log of machine precision */

#define COMM_FIELD "weight"

double rpg (int n, double z); /* from pg.c */
SEXP sweep (SEXP sn, SEXP sm, SEXP sa, SEXP sind); /* from sweep.c */

/* from rinterface.c [igraph] */

static int R_SEXP_to_vector_copy (SEXP sv, igraph_vector_t *v) {
  return igraph_vector_init_copy(v, REAL(sv), length(sv));
}

static int R_SEXP_to_vector (SEXP sv, igraph_vector_t *v) {
  v->stor_begin=REAL(sv);
  v->stor_end=v->stor_begin+length(sv);
  v->end=v->stor_end;
  return 0;
}

static int R_SEXP_to_igraph (SEXP graph, igraph_t *res) {
  res->n=REAL(VECTOR_ELT(graph, 0))[0];
  res->directed=LOGICAL(VECTOR_ELT(graph, 1))[0];
  R_SEXP_to_vector(VECTOR_ELT(graph, 2), &res->from);
  R_SEXP_to_vector(VECTOR_ELT(graph, 3), &res->to);
  R_SEXP_to_vector(VECTOR_ELT(graph, 4), &res->oi);
  R_SEXP_to_vector(VECTOR_ELT(graph, 5), &res->ii);
  R_SEXP_to_vector(VECTOR_ELT(graph, 6), &res->os);
  R_SEXP_to_vector(VECTOR_ELT(graph, 7), &res->is);
  
  /* attributes */
  REAL(VECTOR_ELT(VECTOR_ELT(graph, 8), 0))[0] = 1; /* R objects refcount */
  REAL(VECTOR_ELT(VECTOR_ELT(graph, 8), 0))[1] = 0; /* igraph_t objects */
  res->attr=VECTOR_ELT(graph, 8);
  
  return 0;
}

static SEXP R_igraph_getListElement(SEXP list, const char *str) {
  SEXP names = getAttrib(list, R_NamesSymbol);
  int i;
  for (i = 0; i < length(list); i++)
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0)
      return VECTOR_ELT(list, i);
  return R_NilValue;
}

/* numeric edge attributes for the whole graph */
static int R_igraph_attribute_get_edge_attr (const igraph_t *graph,
    const char *name, igraph_vector_t *value) {
  SEXP eal = VECTOR_ELT(graph->attr, 3);
  SEXP ea = R_igraph_getListElement(eal, name);
  igraph_vector_t newvalue;

  if (ea == R_NilValue)
    IGRAPH_ERROR("No such attribute", IGRAPH_EINVAL);
  PROTECT(ea = coerceVector(ea, REALSXP));
  R_SEXP_to_vector_copy(coerceVector(ea, REALSXP), &newvalue);
  igraph_vector_destroy(value);
  *value = newvalue;
  UNPROTECT(1);
  return 0;
}


/* auxiliary */

static double log1pe (double x) {
  double d = (x > 0) ? -x : x;
  d = (d < LOGEPS) ? 0 : log1p(exp(d));
  if (x > 0) d += x;
  return d;
}

static double lse (double w1, double w2) {
  double d, w;
  if (!isfinite(w1)) return w2;
  if (!isfinite(w2)) return w1;
  if (w1 > w2) {
    d = w2 - w1; w = w1;
  }
  else {
    d = w1 - w2; w = w2;
  }
  if (d < LOGEPS) return w;
  return w + log1p(exp(d));
}

static int lsample (int n, double *t) {
  int i;
  double u, z, c = 0;
  /* normalize */
  z = t[0]; for (i = 1; i < n; i++) z = lse(z, t[i]);
  /* sample */
  GetRNGstate(); u = unif_rand(); PutRNGstate();
  for (i = 0; i < n && u >= c; i++)
    c += exp(t[i] - z);
  return i - 1;
}


/* group index in beta: returns negative index if g1 and g2 are invalid */
static int gindex (int g1, int g2, int k) {
  if (g1 == g2) return -1; /* within */
  if (g2 < g1) { /* swap? */
    int t = g2; g2 = g1; g1 = t;
  }
  return g1 * (k - 1) - g1 * (g1 + 1) / 2 + g2 - 1; /* zero based */
}

/* TODO: optimize by avoiding allocations at each call */
static void remap (int n, int k, int *sigma, double *pi, double *beta,
    int *card, double *sigmap) {
  int c = 0; /* #labels in sigma */
  int i, j, l;
  int *revind, *pcard;
  double *ppi, *pbeta;
  revind = Calloc(k, int);
  for (l = 0; l < k; l++) {
    card[l] = 0;
    revind[l] = -1;
  }

  /* [ find order ] */
  for (i = 0; i < n; i++) {
    l = sigma[i];
    if (card[l] == 0) /* first visit? */
      revind[l] = c++; /* reverse order of appearance of `l` */
    card[l]++;
  }
  for (l = 0; l < k; l++) {
    if (revind[l] == -1) /* no appearance? */
      revind[l] = c++;
  }

  /* [ permute ] */
  for (i = 0; i < n; i++)
    sigma[i] = revind[sigma[i]];

  pcard = Calloc(k, int);
  ppi = Calloc(k, double);
  pbeta = Calloc(k * (k - 1) / 2, double);
  for (i = 0; i < k; i++) { /* copy */
    ppi[i] = pi[i];
    pcard[i] = pcard[i];
    for (j = i + 1; j < k; j++) {
      int gk = gindex(i, j, k);
      pbeta[gk] = beta[gk];
    }
  }
  for (i = 0; i < k; i++) { /* permute: l -> revind[l] */
    int l = revind[i];
    pi[l] = ppi[i];
    card[l] = pcard[i];
    for (j = i + 1; j < k; j++) {
      int gk = gindex(i, j, k);
      int pk = gindex(l, revind[j], k);
      if (gk < 0 || pk < 0) error("invalid index");
      beta[pk] = pbeta[gk];
    }
  }

  if (sigmap) { /* posterior means? */
    double *psigma = Calloc(k, double);
    for (i = 0; i < n; i++) {
      for (j = 0; j < k; j++) psigma[j] = sigmap[i + j * n];
      for (j = 0; j < k; j++) {
        int l = revind[j];
        sigmap[i + l * n] = psigma[j];
      }
    }
    Free(psigma);
  }

  Free(revind); Free(pcard); Free(ppi); Free(pbeta);
}



/* [ R API ] */

SEXP rtruncnorm (SEXP smu, SEXP ssigma) {
  double *mu = REAL(smu);
  double *sigma = REAL(ssigma);
  int i, n = length(smu);
  double u, *x;
  SEXP sx;
  PROTECT(sx = allocVector(REALSXP, n));
  x = REAL(sx);
  GetRNGstate();
  for (i = 0; i < n; i++) {
    if (mu[i] >= 0) { /* rejection by support? */
      do
        u = mu[i] + sigma[i] * norm_rand();
      while (u <= 0);
    }
    else {
      double ms = -mu[i] / sigma[i];
      do
        u = ms * exp_rand(); /* ~ Exp(ms) */
      while (unif_rand() >= exp(-.5 * u * u));
      u *= sigma[i];
    }
    x[i] = u;
  }
  PutRNGstate();
  UNPROTECT(1);
  return sx;
}


/* if `weight` == 0, treat as binary case */
SEXP comm_beta_mcmc_bin (SEXP sg, SEXP ssigma, SEXP sk, SEXP stau2,
    SEXP sbeta, SEXP sweight, SEXP sV) {
  SEXP sm;
  igraph_t g;
  igraph_vector_t count;
  double tau2, *beta, *V, weight, *m;
  int i, j, n, e, ne, k, gk, ck, w, *sigma;

  R_SEXP_to_igraph(sg, &g);
  sigma = INTEGER(ssigma);
  k = INTEGER(sk)[0];
  tau2 = REAL(stau2)[0];
  beta = REAL(sbeta);
  weight = REAL(sweight)[0];
  w = (weight > 0) ? (int) weight : 1;
  V = REAL(sV);
  n = igraph_vcount(&g);

  /* initialize */
  ck = k * (k - 1) / 2;
  PROTECT(sm = allocVector(REALSXP, n + ck));
  m = REAL(sm);

  /* compute mean */
  ne = igraph_ecount(&g);
  if (weight > 0) {
    igraph_vector_init(&count, ne);
    R_igraph_attribute_get_edge_attr(&g, COMM_FIELD, &count);
  }
  for (i = 0; i < n + ck; i++) m[i] = 0;
  for (e = 0; e < ne; e++) { /* for all edges */
    double y = (weight > 0) ? VECTOR(count)[e] / weight : 1.;
    igraph_edge(&g, e, &i, &j);
    gk = gindex(sigma[i], sigma[j], k);
    if (gk >= 0) m[gk] += y;
    m[i + ck] += y;
    m[j + ck] += y;
  }
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      int gk = gindex(sigma[i], sigma[j], k);
      if (gk >= 0) m[gk] -= .5;
    }
    m[i + ck] -= (n - 1) / 2.0;
  }
  if (weight > 0)
    for (i = 0; i < n + ck; i++) m[i] *= weight;
  if (weight > 0) igraph_vector_destroy(&count);
  
  /* V = Sigma^(-1) = I_n / tau2 */
  for (i = 0; i < n + ck; i++) {
    V[i * (n + ck + 1)] = 1 / tau2; /* V[i][i] */
    for (j = i + 1; j < n + ck; j++)
      V[i + (n + ck) * j] = 0.; /* V[i][j] */
  }

  /* [ compute statistics ] */
  for (i = 0; i < n; i++) { /* for each pair i < j */
    int s = sigma[i];
    for (j = i + 1; j < n; j++) {
      double wij = beta[i + ck] + beta[j + ck];
      int gk = gindex(sigma[j], s, k);
      if (gk >= 0) wij += beta[gk];
      wij = rpg(w, wij);
      /* update entries */
      if (gk >= 0) {
        V[gk * (n + ck + 1)] += wij; /* V[gk][gk] */
        V[gk + (i + ck) * (n + ck)] += wij; /* V[gk][i + ck] */
        V[gk + (j + ck) * (n + ck)] += wij; /* V[gk][j + ck] */
      }
      V[(n + ck + 1) * (i + ck)] += wij; /* V[i + ck][i + ck] */
      V[(n + ck + 1) * (j + ck)] += wij; /* V[j + ck][j + ck] */
      V[(i + ck) + (n + ck) * (j + ck)] = wij; /* V[i + ck][j + ck] */
    }
  }
  UNPROTECT(1);
  return sm;
}


/* if `weight` == 0, treat as binary case; if `k` == 1, no community effect */
SEXP comm_beta_map_bin (SEXP sg, SEXP ssigma, SEXP sk, SEXP stau2,
    SEXP sbeta, SEXP sweight, SEXP sV) {
  SEXP sm;
  igraph_t g;
  igraph_vector_t count;
  double tau2, *beta, *V, weight, *m;
  int i, j, n, e, ne, k, gk, ck, *sigma;

  R_SEXP_to_igraph(sg, &g);
  sigma = isNull(ssigma) ? NULL : INTEGER(ssigma);
  k = INTEGER(sk)[0];
  tau2 = REAL(stau2)[0];
  beta = REAL(sbeta);
  weight = REAL(sweight)[0];
  V = REAL(sV);
  n = igraph_vcount(&g);

  /* [ initialize ] */
  ck = k * (k - 1) / 2;
  PROTECT(sm = allocVector(REALSXP, n + ck));
  m = REAL(sm);

  /* compute mean: m = X^T * y - Sigma^(-1) * beta */
  ne = igraph_ecount(&g);
  if (weight > 0) {
    igraph_vector_init(&count, ne);
    R_igraph_attribute_get_edge_attr(&g, COMM_FIELD, &count);
  }
  if (weight > 0)
    for (i = 0; i < n + ck; i++) m[i] = -beta[i] / tau2 / weight;
  else
    for (i = 0; i < n + ck; i++) m[i] = -beta[i] / tau2;
  for (e = 0; e < ne; e++) { /* for all edges */
    double y = (weight > 0) ? VECTOR(count)[e] / weight : 1.;
    igraph_edge(&g, e, &i, &j);
    if (sigma != NULL && k > 1) {
      gk = gindex(sigma[i], sigma[j], k);
      if (gk >= 0) m[gk] += y;
    }
    m[i + ck] += y;
    m[j + ck] += y;
  }
  if (weight > 0) igraph_vector_destroy(&count);

  /* [ compute statistics ] */
  /* V = Sigma^(-1) = I_n / tau2 */
  for (i = 0; i < n + ck; i++) {
    if (weight > 0)
      V[i * (n + ck + 1)] = 1. / tau2 / weight; /* V[i][i] */
    else
      V[i * (n + ck + 1)] = 1. / tau2; /* V[i][i] */
    for (j = i + 1; j < n + ck; j++)
      V[i + (n + ck) * j] = 0.; /* V[i][j] */
  }

  for (i = 0; i < n; i++) { /* for all pairs i < j */
    for (j = i + 1; j < n; j++) {
      double muij, muetaij;
      double etaij = beta[i + ck] + beta[j + ck];
      int gk = -1;
      if (sigma != NULL && k > 1) {
        gk = gindex(sigma[j], sigma[i], k);
        if (gk >= 0) etaij += beta[gk];
      }
      /* MAIN */
      muij = 1. / (1 + exp(-etaij)); /* inverse logit (link) */
      muetaij = muij * (1 - muij); /* variance */
      /* MAIN */
      /* update entries */
      if (gk >= 0) {
        m[gk] -= muij; /* from X^T * mu */
        V[gk * (n + ck + 1)] += muetaij; /* V[gk][gk] */
        V[gk + (i + ck) * (n + ck)] += muetaij; /* V[gk][i + ck] */
        V[gk + (j + ck) * (n + ck)] += muetaij; /* V[gk][j + ck] */
      }
      m[i + ck] -= muij;
      m[j + ck] -= muij;
      V[(n + ck + 1) * (i + ck)] += muetaij; /* V[i + ck][i + ck] */
      V[(n + ck + 1) * (j + ck)] += muetaij; /* V[j + ck][j + ck] */
      V[(i + ck) + (n + ck) * (j + ck)] = muetaij; /* V[i + ck][j + ck] */
    }
  }
  UNPROTECT(1);
  return sm;
}



/* samples/fits sigma in-place */
SEXP comm_sigma_bin (SEXP sg, SEXP ssigma, SEXP slogpi, SEXP sbeta, SEXP sweight,
    SEXP ssample) {
  igraph_t g;
  double *logpi, *beta, weight, w, lhood;
  int *sigma, *card, *isneighbor;
  int i, j, e, s, k, ck, n, sample;
  igraph_vector_t adj;
  igraph_vector_t lp;
  igraph_vector_t count;

  R_SEXP_to_igraph(sg, &g);
  sigma = INTEGER(ssigma);
  logpi = REAL(slogpi);
  beta = REAL(sbeta);
  weight = REAL(sweight)[0];
  w = (weight > 0) ? weight : 1.;
  sample = INTEGER(ssample)[0];
  k = length(slogpi); /* number of groups */
  n = igraph_vcount(&g);
  igraph_vector_init(&adj, n);
  igraph_vector_init(&lp, k);
  if (weight > 0) {
    igraph_vector_init(&count, igraph_ecount(&g));
    R_igraph_attribute_get_edge_attr(&g, COMM_FIELD, &count);
  }

  card = Calloc(k, int);
  isneighbor = Calloc(n, int);
  for (i = 0; i < k; i++) card[i] = 0;
  for (i = 0; i < n; i++) card[sigma[i]]++;

  ck = k * (k - 1) / 2;
  lhood = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) isneighbor[j] = 0;
    igraph_neighbors(&g, &adj, i, IGRAPH_ALL);
    for (j = 0; j < igraph_vector_size(&adj); j++)
      isneighbor[(int) VECTOR(adj)[j]] = 1;

    for (s = 0; s < k; s++) {
      double lps = logpi[s];
      for (j = 0; j < n; j++) {
        if (i != j) {
          double xij = beta[i + ck] + beta[j + ck];
          int gk = gindex(sigma[j], s, k);
          double yij = isneighbor[j];
          if (gk >= 0) xij += beta[gk];
          if (weight > 0 && isneighbor[j]) {
            igraph_get_eid(&g, &e, i, j, 0, 0);
            yij = VECTOR(count)[e];
          }
          lps += yij * xij - w * log1pe(xij); /* MAIN */
        }
      }
      VECTOR(lp)[s] = lps;
    }
    lhood += VECTOR(lp)[sigma[i]]; /* FIXME: move to after sigma is assigned? */

    /* check if sigma[i] needs to be sampled: */
    if (card[sigma[i]] == 1) continue; /* single node in group */
    s = sample ? lsample(k, VECTOR(lp)) : igraph_vector_which_max(&lp);
    if (s != sigma[i]) { /* update card? */
      card[sigma[i]]--;
      card[s]++;
    }
    sigma[i] = s;
  }
  remap(n, k, sigma, logpi, beta, card, NULL);

  igraph_vector_destroy(&adj);
  igraph_vector_destroy(&lp);
  if (weight > 0) igraph_vector_destroy(&count);
  Free(card); Free(isneighbor);
  return ScalarReal(lhood);
}


/* Interface */

static const R_CallMethodDef callMethods[] = {
  /* auxiliar */
  {"rtruncnorm", (DL_FUNC) &rtruncnorm, 2},
  {"sweep", (DL_FUNC) &sweep, 4},
  /* logistic */
  {"comm_beta_mcmc_bin", (DL_FUNC) &comm_beta_mcmc_bin, 7},
  {"comm_beta_map_bin", (DL_FUNC) &comm_beta_map_bin, 7},
  {"comm_sigma_bin", (DL_FUNC) &comm_sigma_bin, 6},
  {NULL, NULL, 0}
};

void R_init_sbmlogit (DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

