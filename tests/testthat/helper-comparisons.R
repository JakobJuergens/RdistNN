# Exact comparisons for retained state and bit-sensitive numeric fixtures.
# On failure, avoid recursively diffing closure environments or formatting
# entire raw/large numeric vectors. This changes diagnostics, not tolerance.
.expect_identical_fast <- function(object, expected, info = NULL) {
  same <- identical(object, expected)
  message <- if (same) "Objects are identical." else sprintf(
    "Expected `%s` to be identical to `%s`.",
    paste(deparse(substitute(object)), collapse = " "),
    paste(deparse(substitute(expected)), collapse = " ")
  )
  testthat::expect(same, message, info = info)
  invisible(object)
}
