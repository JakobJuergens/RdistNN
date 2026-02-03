#!/usr/bin/env bash
set -euo pipefail

if grep -R -n "\[\[Rcpp::plugins(cpp11)\]\]" src; then
  echo "Found cpp11 plugin annotation in src/. Please use cpp17."
  exit 1
fi

files=$(find src inst/include -type f \( -name "*.cpp" -o -name "*.hpp" \))
if [[ -n "${files}" ]]; then
  clang-format --dry-run --Werror ${files}
fi
