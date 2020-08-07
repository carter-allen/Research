#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>

#define M_4_PI (2.0 * M_2_PI)
#define M_2_PI2 (M_2_PI * M_1_PI)
/* terms in `a' series */
#define TRUNC (0.64)
#define A1N(x) (K * exp(-0.5 * K * K * (x))) /* x > TRUNC, tail */
#define A2N(x) (K * exp(-1.5 * log(M_PI_2 * (x)) - M_2_PI2 * K * K / (x)))
#define AN(x) ((x) > TRUNC) ? A1N(x) : A2N(x)
#define LPNORM(x) pnorm((x),0,1,1,1)


/* Polya-Gamma deviates */

/* computes p/(p+q) where 
    p  = 0.5 * __PI * exp(-1.0 * fz * TRUNC) / fz
    q  = 2 * exp(-1.0 * Z) * pigauss(TRUNC, Z)
  and pigauss is the CDF of the inverse normal */
static double probexp (double z, double fz) {
  double b = sqrt(1.0 / TRUNC) * (TRUNC * z - 1);
  double a = -sqrt(1.0 / TRUNC) * (TRUNC * z + 1);
  double x0 = log(fz) + fz * TRUNC;
  double xb = x0 - z + LPNORM(b);
  double xa = x0 + z + LPNORM(a);
  return 1.0 / (1.0 + M_4_PI * (exp(xb) + exp(xa)));
}

/* assumes z >= 0 */
static double rtigauss (double z) {
  double x = TRUNC + 1.0;
  if (z * TRUNC < 1.0) { /* mu > t? */
    double alpha = 0.0;
    while (runif(0, 1) > alpha) {
      double E1 = rexp(1);
      double E2 = rexp(1);
      while (E1 * E1 > 2 * E2 / TRUNC) {
        E1 = rexp(1); E2 = rexp(1);
      }
      x = 1 + E1 * TRUNC;
      x = TRUNC / (x * x);
      alpha = exp(-0.5 * z * z * x);
    }
  }
  else {
    double mu = 1.0 / z;
    while (x > TRUNC) {
      double y = rnorm(0, 1);
      double muy2 = 0.5 * mu * y;
      muy2 *= muy2;
      x = mu + 2 * muy2 * (1 - sqrt(1 + mu / muy2));
      if (runif(0, 1) > mu / (mu + x)) x = mu * mu / x;
    }
  }
  return x;
}

/* assumes z >= 0 */
static double rpgaux (double z, double fz, double pq) {
  for (;;) {
    double x, s, y;
    int n = 0; /* iteration counter */
    int within = 1; /* flag for y <= s */
    double K = M_PI_2; /* n = 0, K = pi / 2 */
    if (runif(0, 1) < pq) /* truncated Exp(1)? */
      x = TRUNC + rexp(1) / fz;
    else /* truncated inverse Gaussian */
      x = rtigauss(z);
    /* FIXME: check if x is close to zero */

    s = AN(x);
    y = runif(0, 1) * s;
    while (within) {
      n++; K += M_PI;
      if (n % 2 == 1) {
        s -= AN(x);
        if (y <= s) return 0.25 * x;
      }
      else {
        s += AN(x);
        within = y <= s;
      }
    }
  }
}


double rpg (int n, double z) {
  double fz, pq, r;
  int i;
  z = 0.5 * fabs(z);
  fz = 0.5 * (M_PI_2 * M_PI_2 + z * z);
  pq = probexp(z, fz);
  GetRNGstate();
  if (n == 1)
    r = rpgaux(z, fz, pq);
  else {
    r = 0;
    for (i = 0; i < n; i++) r += rpgaux(z, fz, pq);
  }
  PutRNGstate();
  return r;
}

double epg (double n, double z) {
  if (z == 0) return 1;
  return n / (2 * z) * tanh(z / 2);
}

