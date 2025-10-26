// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
using namespace Rcpp;

// Minimal context just for CATE; kept inside this .cpp
struct CateCtx {
  NumericMatrix X;      // n x p
  NumericVector rho;    // length n (pseudo-outcome)
  NumericVector x0;     // length p (target point)
  int p;
  std::string localizer; // "dnn1" or "dnnk" or "kernel"
  int k;                // for dnnk
  double h;             // bandwidth for kernel (<=0 => auto)
};

// helper: assign weights for dnn1
static inline void weights_dnn1(const NumericVector& dx, NumericVector& w) {
  int m = which_min(dx);
  std::fill(w.begin(), w.end(), 0.0);
  w[m] = 1.0;
}

// helper: assign weights for dnnk
static inline void weights_dnnk(const NumericVector& dx, NumericVector& w, int k) {
  IntegerVector ord = seq(0, dx.size() - 1);
  std::sort(ord.begin(), ord.end(), [&](int i, int j){ return dx[i] < dx[j]; });
  std::fill(w.begin(), w.end(), 0.0);
  int kk = std::min(k, (int)dx.size());
  for (int t = 0; t < kk; ++t) w[ ord[t] ] = 1.0 / kk;
}

// helper: Gaussian kernel weights, auto-bandwidth if needed
static inline void weights_kernel(const NumericVector& dx, NumericVector& w, double h) {
  std::fill(w.begin(), w.end(), 0.0);
  double hh = h;
  if (!(hh > 0)) {
    // crude auto h: median of positive distances, else 1.0
    std::vector<double> pos; pos.reserve(dx.size());
    for (int i=0;i<dx.size();++i) if (dx[i] > 0 && R_finite(dx[i])) pos.push_back(dx[i]);
    if (!pos.empty()) {
      std::nth_element(pos.begin(), pos.begin() + pos.size()/2, pos.end());
      hh = pos[pos.size()/2];
    } else {
      hh = 1.0;
    }
  }
  double s = 0.0;
  for (int i=0;i<dx.size();++i) {
    double g = std::exp(-0.5 * (dx[i]/hh)*(dx[i]/hh));
    w[i] = g; s += g;
  }
  if (s <= 0) for (int i=0;i<dx.size();++i) w[i] = 1.0/dx.size();
  else for (int i=0;i<dx.size();++i) w[i] /= s;
}

// Create a context external pointer
// [[Rcpp::export]]
SEXP rustat_cate_ctx_create(NumericMatrix X, NumericVector rho, NumericVector x0,
                            std::string localizer, int k, double h) {
  if (X.nrow() != rho.size()) stop("X and rho dimensions mismatch");
  if (x0.size() != X.ncol())  stop("x0 length must equal ncol(X)");
  XPtr<CateCtx> ctx(new CateCtx(), true);
  ctx->X = X;
  ctx->rho = rho;
  ctx->x0 = x0;
  ctx->p = X.ncol();
  ctx->localizer = localizer;
  ctx->k = k;
  ctx->h = h;
  return ctx;
}

// Evaluate a batch of subsets: subs is K x s (1-based indices from R)
// [[Rcpp::export]]
NumericVector rustat_cate_eval_batch(IntegerMatrix subs, SEXP ctx_) {
  XPtr<CateCtx> ctx(ctx_);
  int K = subs.nrow(), s = subs.ncol();
  NumericVector out(K);
  NumericVector dx(s), w(s);

  for (int kk=0; kk<K; ++kk) {
    // distances for each element of the subset to x0
    for (int t=0; t<s; ++t) {
      int i = subs(kk, t) - 1; // R is 1-based
      double acc = 0.0;
      for (int j=0; j<ctx->p; ++j) {
        double dif = ctx->X(i,j) - ctx->x0[j];
        acc += dif * dif;
      }
      dx[t] = std::sqrt(acc);
    }
    // weights
    if (ctx->localizer == "dnn1") {
      weights_dnn1(dx, w);
    } else if (ctx->localizer == "dnnk") {
      weights_dnnk(dx, w, std::max(1, ctx->k));
    } else { // "kernel"
      weights_kernel(dx, w, ctx->h);
    }
    // weighted average of rho in the subset
    double val = 0.0;
    for (int t=0; t<s; ++t) {
      int i = subs(kk, t) - 1;
      val += w[t] * ctx->rho[i];
    }
    out[kk] = val;
  }
  return out;
}
