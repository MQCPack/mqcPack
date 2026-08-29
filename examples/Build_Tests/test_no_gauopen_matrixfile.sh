#!/bin/sh

output_file="${TMPDIR:-/tmp}/mqc-no-gauopen-$$.log"
trap 'rm -f "$output_file"' EXIT HUP INT TERM

if ./unitTest06 >"$output_file" 2>&1; then
  echo "unitTest06 unexpectedly completed successfully"
  cat "$output_file"
  exit 1
fi

expected_message="MQCPack was configured without GauOpen support; MatrixFile operations are unavailable."
if ! grep -F "$expected_message" "$output_file" >/dev/null; then
  echo "unitTest06 did not report the expected no-GauOpen diagnostic"
  cat "$output_file"
  exit 1
fi

echo "PASS: FormChk-only MatrixFile operation reports unsupported feature"
exit 0
