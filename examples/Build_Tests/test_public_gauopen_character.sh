#!/bin/sh

output_file=unitTest03-public-gauopen.out
faf_file=unitTest03-character.faf
expected_message="Character FAF I/O requires the GauOpen frontier TypeA API; the configured public GauOpen API does not provide character record support."

if ./unitTest03 >"$output_file" 2>&1; then
  echo "unitTest03 unexpectedly completed with public GauOpen"
  rm -f "$output_file" "$faf_file"
  exit 1
fi

if ! grep -F "$expected_message" "$output_file" >/dev/null 2>&1; then
  echo "unitTest03 did not report the expected public-GauOpen character diagnostic"
  cat "$output_file"
  rm -f "$output_file" "$faf_file"
  exit 1
fi

rm -f "$output_file" "$faf_file"
echo "PASS: public GauOpen character FAF operation reports unsupported feature"
