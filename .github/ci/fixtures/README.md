# Native output fixtures

These small RDS fixtures contain synthetic sequences and coordinates, not patient
data. They were generated on 2026-09-17 using the preserved upstream Longcellsrc
1.1.0 implementation with R 4.4.3. The reviewed 1.1.1 backend and integrated
LongcellPre 1.3.0 both matched them exactly in the local validation.

- `native-baseline.rds`: 35 UMI graph comparisons (seed 7301), plus barcode outputs
  for 300 tags and 1,500 whitelist entries at cosine thresholds 0, 0.25 and 0.75.
- `iso-baseline.rds`: 50 isoform-match tables (seed 429).

The accompanying scripts fix the RNG kind and seeds, reconstruct the inputs, and
compare actual values with `identical()`. They also test boundary handling, empty
inputs, weight overflow and the native graph edge limit. These CI scripts only
read the reference fixtures; they cannot regenerate them from the candidate under
test. Intentional changes to scientific semantics need explicit review and a
separate justification for any fixture replacement.
