// src/fast_order.hpp
#pragma once
#include <cstddef>

namespace rustat {

// Returns true iff all x[i] are finite and >= 0
bool all_finite_nonneg(const double *x, int n) noexcept;

// Stable argsort into out_idx (0-based indices, ascending).
// Finite/>=0 -> radix path; otherwise stable comparison sort fallback.
void argsort_indices(const double *x, int n, int *out_idx);

// Stable 1-based ranks from values (ties by original index).
// ranks_out[i] = rank of x[i] in {1..n}.
void ranks_from_values(const double *x, int n, int *ranks_out);

// Returns k (clamped to [0,n]). Writes sorted (ascending) 0-based indices
// of the k smallest values into out_idx_sorted[0..k-1].
int topk_indices(const double *x, int n, int k, int *out_idx_sorted);

} // namespace rustat
