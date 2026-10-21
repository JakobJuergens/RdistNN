test_that("signed table and evaluator extensions share compensated positive mass", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  table_name <- "test_compensated_signed_table"
  eval_name <- "test_compensated_signed_eval"
  on.exit(rm(list = c(table_name, eval_name), envir = .rank_registry), add = TRUE)
  tab <- c(2^100, 1, -2^100)
  register_rank_weight(table_name, list(table = function(n) tab))
  register_rank_weight(eval_name, as_rank_weight(function(ranks) tab[ranks]))
  table_rule <- get_rank_weight(table_name)
  eval_rule <- get_rank_weight(eval_name)
  ranks <- cbind(1:3, 3:1)
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    if (native) expect_false(is.null(.rdistnn_compiled_fun("rustat_compensated_sum")))
    expect_identical(rank_weight_eval(table_rule, 1:3), tab)
    expect_identical(rank_weight_eval_batch(table_rule, ranks),
      rank_weight_eval_batch(eval_rule, ranks))
    table_weights <- weights_matrix(matrix(1:3), X0 = matrix(c(0, 4)),
      ranking = "l2", weight = table_rule)$W
    eval_weights <- weights_matrix(matrix(1:3), X0 = matrix(c(0, 4)),
      ranking = "l2", weight = eval_rule)$W
    expect_identical(table_weights, eval_weights)
  }
})

test_that("signed table mass validation recovers intermediate overflow without losing tiny terms", {
  old <- options(RdistNN.use_compiled = TRUE)
  on.exit(options(old), add = TRUE)
  name <- "test_range_signed_table"
  on.exit(rm(list = name, envir = .rank_registry), add = TRUE)
  M <- .Machine$double.xmax
  tiny <- .Machine$double.xmin * .Machine$double.eps
  tab <- numeric()
  register_rank_weight(name, list(table = function(n) tab))
  rule <- get_rank_weight(name)
  accepted <- list(c(M, M, -M), c(M, M, -M, -M, 1),
    c(M, M, -M, -M, 1e-300), c(M, M, -M, -M, tiny))
  rejected <- list(c(M, M), c(M, M, -M, -M),
    c(M, M, -M, -M, -tiny), c(-M, -M, M))
  for (native in c(TRUE, FALSE)) {
    options(RdistNN.use_compiled = native)
    if (native) expect_false(is.null(.rdistnn_compiled_fun("rustat_compensated_sum")))
    for (value in accepted) {
      tab <- value
      expect_identical(rank_weight_eval(rule, seq_along(tab)), tab)
    }
    for (value in rejected) {
      tab <- value
      expect_error(rank_weight_eval(rule, seq_along(tab)), "positive finite sum")
    }
  }
})
