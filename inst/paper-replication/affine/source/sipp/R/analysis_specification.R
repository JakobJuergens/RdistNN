# Economic/reporting choices, independent of nuisance fitting.
make_sipp_401k_analysis_spec <- function(adjustment = "baseline",
  reporting = "age_income", profile_incomes = NULL, metric_scales = NULL,
  asset_cap = 100000) {
  if (!adjustment %in% c("baseline", "demographic", "ira"))
    stop("Unknown adjustment specification.", call. = FALSE)
  if (!reporting %in% c("age_income", "income", "age40"))
    stop("Unknown reporting specification.", call. = FALSE)
  if (!is.null(profile_incomes) && (length(profile_incomes) < 1L ||
      any(!is.finite(profile_incomes)) || anyDuplicated(profile_incomes)))
    stop("Profile incomes must be distinct finite values.", call. = FALSE)
  dim <- if (reporting == "age_income") 2L else 1L
  if (!is.null(metric_scales) && (length(metric_scales) != dim ||
      any(!is.finite(metric_scales)) || any(metric_scales <= 0)))
    stop("Metric scales must be positive, finite and match reporting dimension.", call. = FALSE)
  if (length(asset_cap) != 1L || !is.finite(asset_cap) || asset_cap <= 0)
    stop("Asset cap must be a positive finite scalar.", call. = FALSE)
  list(adjustment = adjustment, reporting = reporting,
       profile_incomes = profile_incomes, metric_scales = metric_scales,
       asset_cap = asset_cap)
}

sipp_401k_adjustment_columns <- function(adjustment) {
  common <- c("age", "age2", "asinh_inc", "asinh_inc2", "educ", "fsize",
              "marr", "male", "twoearn")
  education <- c("nohs", "hs", "smcol", "col")
  switch(adjustment, demographic = c(common, education),
         baseline = c(common, "db", "hown", education),
         ira = c(common, "db", "hown", education, "pira"),
         stop("Unknown adjustment specification.", call. = FALSE))
}

sipp_401k_outcome_definitions <- function(cap) {
  list(raw_dollars = list(column = "Y", units = "1991 dollars", transformation = "identity"),
       asinh_dollars = list(column = "Y_asinh", units = "asinh of 1991 dollars",
                            transformation = "asinh(Y)"),
       positive_assets = list(column = "Y_positive", units = "probability",
                              transformation = "1{Y>0}", bounds = c(0, 1)),
       capped_dollars = list(column = "Y_capped", units = "1991 dollars",
                            transformation = "max(-cap,min(cap,Y))",
                            cap = cap, bounds = c(-cap, cap)))
}
