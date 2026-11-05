#!/bin/sh
# Expected failures for the Gaussian-independent quadrature API.
set -eu
for pair in 10:4 11:6 12:8; do
  test_number=${pair%:*}
  last_case=${pair#*:}
  case_number=1
  while [ "$case_number" -le "$last_case" ]; do
    if "./unitTest${test_number}" "$case_number" >/dev/null 2>&1; then
      echo "unitTest${test_number}: invalid case ${case_number} was accepted"
      exit 1
    fi
    case_number=$((case_number + 1))
  done
done
echo 'All 18 invalid quadrature cases rejected.'
