#ifndef RDISTNN_COVARIATE_TIES_HPP
#define RDISTNN_COVARIATE_TIES_HPP

#include <algorithm>
#include <cstddef>

// Used by coordinate-aware R native metric rankers. X is the original finite,
// column-major n x p design, before scaling or whitening. Refine only exactly
// equal represented distances; primary distance arithmetic stays unchanged.
namespace rdistnn_ordering {

inline bool covariate_rows_less(const double *X, int n, int p, int i, int j) {
  for (int column = 0; column < p; ++column) {
    const std::size_t offset = static_cast<std::size_t>(column) * n;
    const double a = X[offset + i], b = X[offset + j];
    if (a < b)
      return true;
    if (b < a)
      return false;
  }
  return i < j;
}

inline void refine_distance_ties(const double *distance, const double *X, int n, int p,
                                 int *order) {
  for (int begin = 0; begin < n;) {
    int end = begin + 1;
    while (end < n && distance[order[end]] == distance[order[begin]])
      ++end;
    if (end - begin > 1) {
      std::sort(order + begin, order + end,
                [=](int i, int j) { return covariate_rows_less(X, n, p, i, j); });
    }
    begin = end;
  }
}

} // namespace rdistnn_ordering
#endif
