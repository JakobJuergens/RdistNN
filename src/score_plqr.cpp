// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

struct PlqrCtx {
  NumericMatrix X;      // n x p
  NumericMatrix D;      // n x p_d
  NumericMatrix mD;     // n x p_d
  NumericVector Y;      // n
  NumericVector qx;     // n
  NumericVector x0;     // p
  NumericVector beta;   // p_d

  int p;                // dim X
  int p_d;              // dim D
  int j;                // component index (1-based from R)
  std::string localizer;
  int k;
  double h;             // <=0 => auto
  double tau;
};

static inline void weights_dnn1(const NumericVector& dx, NumericVector& w) {
  int m = which_min(dx);
  std::fill(w.begin(), w.end(), 0.0);
  w[m] = 1.0;
}

static inline void weights_dnnk(const NumericVector& dx, NumericVector& w, int k) {
  std::vector<int> ord(dx.size());
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int a, int b){
    if (dx[a] < dx[b]) return true;
    if (dx[a] > dx[b]) return false;
    return a < b; // tie-break by index to mimic R's order()
  });
  std::fill(w.begin(), w.end(), 0.0);
  int kk = std::min(k, (int)dx.size());
  for (int t=0; t<kk; ++t) w[ ord[t] ] = 1.0 / kk;
}

static inline void weights_kernel(const NumericVector& dx, NumericVector& w, double h) {
  std::fill(w.begin(), w.end(), 0.0);
  double hh = h;
  if (!(hh > 0)) {
    std::vector<double> pos; pos.reserve(dx.size());
    for (int i=0; i<dx.size(); ++i) if (dx[i] > 0 && R_finite(dx[i])) pos.push_back(dx[i]);
    if (!pos.empty()) {
      std::sort(pos.begin(), pos.end());
      size_t m = pos.size();
      hh = (m % 2u) ? pos[m/2] : 0.5 * (pos[m/2 - 1] + pos[m/2]);  // true median
    } else {
      hh = 1.0;
    }
  }
  double ssum = 0.0;
  for (int i=0; i<dx.size(); ++i) { double g = std::exp(-0.5 * (dx[i]/hh)*(dx[i]/hh)); w[i] = g; ssum += g; }
  if (ssum <= 0) for (int i=0; i<dx.size(); ++i) w[i] = 1.0/dx.size();
  else           for (int i=0; i<dx.size(); ++i) w[i] /= ssum;
}

// [[Rcpp::export]]
SEXP rustat_plqr_ctx_create(NumericMatrix X, NumericVector Y,
                            NumericMatrix D, NumericVector qx, NumericMatrix mD,
                            NumericVector x0,
                            std::string localizer, int k, double h,
                            double tau, int j, NumericVector beta) {
  if (X.nrow() != D.nrow() || X.nrow() != mD.nrow() || X.nrow() != Y.size() || X.nrow() != qx.size())
    stop("X, D, mD, Y, qx must share n.");
  if (X.ncol() != x0.size()) stop("x0 length must equal ncol(X).");
  if (D.ncol() != beta.size()) stop("beta length must equal ncol(D).");
  if (!(j >= 1 && j <= D.ncol())) stop("j out of range.");

  XPtr<PlqrCtx> ctx(new PlqrCtx(), true);
  ctx->X = X; ctx->D = D; ctx->mD = mD;
  ctx->Y = Y; ctx->qx = qx; ctx->x0 = x0; ctx->beta = beta;
  ctx->p = X.ncol(); ctx->p_d = D.ncol();
  ctx->j = j; ctx->localizer = localizer; ctx->k = k; ctx->h = h; ctx->tau = tau;
  return ctx;
}

// [[Rcpp::export]]
NumericVector rustat_plqr_eval_batch(IntegerMatrix subs, SEXP ctx_) {
  XPtr<PlqrCtx> ctx(ctx_);
  int K = subs.nrow(), s = subs.ncol();
  NumericVector out(K);
  NumericVector dx(s), w(s);

  const int j0 = ctx->j - 1; // 0-based column

  for (int kk=0; kk<K; ++kk) {
    // distances to x0 for subset rows
    for (int t=0; t<s; ++t) {
      int i = subs(kk, t) - 1;
      double acc = 0.0;
      for (int p=0; p<ctx->p; ++p) {
        double dif = ctx->X(i,p) - ctx->x0[p];
        acc += dif*dif;
      }
      dx[t] = std::sqrt(acc);
    }
    // weights
    if (ctx->localizer == "dnn1")      weights_dnn1(dx, w);
    else if (ctx->localizer == "dnnk") weights_dnnk(dx, w, std::max(1, std::min(ctx->k, s)));
    else                                weights_kernel(dx, w, ctx->h);

    // compute value
    double val = 0.0;
    for (int t=0; t<s; ++t) {
      int i = subs(kk, t) - 1;
      // threshold D_i^T beta + qx_i
      double db = 0.0;
      for (int pd=0; pd<ctx->p_d; ++pd) db += ctx->D(i,pd) * ctx->beta[pd];
      double thr = db + ctx->qx[i];
      int ind = (ctx->Y[i] <= thr) ? 1 : 0;
      double resid = ctx->tau - ind;
      double rDj = ctx->D(i, j0) - ctx->mD(i, j0);
      val += w[t] * resid * rDj;
    }
    out[kk] = val;
  }
  return out;
}
