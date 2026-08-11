# General utilities

Source authority: [`src/mqc_general.F03`](../../src/mqc_general.F03) and
[`src/mqc_general_lapack.F03`](../../src/mqc_general_lapack.F03).

`MQC_General` is the lowest reusable MQCPack layer. It owns constants, version
and error helpers, intrinsic-value printing and conversion, small array
kernels, packed/full storage conversions, sorting and contractions, and the
project's BLAS/LAPACK wrappers.

## Preferred entry points

- Versioning: `mqc_version`, `mqc_version_check`, `mqc_version_print`.
- Errors: `mqc_error` and the typed `mqc_error_*` helpers. Hard invalid input
  or invalid state should use these rather than an unstructured `stop`.
- Output: `mqc_print`, `mqc_print_scalar`, `mqc_print_vector`,
  `mqc_print_matrix`, and `mqc_print_r4tensor` for intrinsic character,
  integer, real, complex, and supported logical arrays.
- Numeric conversion: `mqc_float`, `num2char`, `integer2character`,
  `real2character`, and `complex2character`.
- Character kernels: `mqc_adjustl`, `mqc_adjustr`, `mqc_len_trim`, and scalar
  `mqc_trim`. These operate on intrinsic values; `MQC_Algebra2` supplies the
  corresponding `MQC_Variable` layer.
- Array/storage helpers: `mqc_packedDiagonalMatrix2FullMatrix`,
  `mqc_matrixSymm2Full`, `matrixOrderedColumns`, `flatten`, and `contraction`.
- Linear solves: the `mqc_dgesv` generic dispatches matrix and vector right-hand
  sides. Other BLAS/LAPACK-backed kernels are declared in this module and
  implemented partly through `mqc_general_lapack.F03`.

Consult the [procedure index](PROCEDURE_INDEX.md) for all specific overloads.

## Numerical rules

- Use `mqc_small` for numerical zero and threshold comparisons.
- Use `real64`, `int64`, and `int32` explicitly in new interfaces even though
  supported build configurations promote default kinds.
- Use `mqc_float(...)` when converting integers or literals in mixed numerical
  formulas.
- Treat packed/full conversion routines as storage transformations. Confirm the
  expected packed ordering before interchanging their output with external
  software.

## Error and output behavior

The `iOut` convention is established throughout MQCPack; unit 6 is the normal
screen output. Preserve caller-supplied output units. Public output labels may
be consumed by tests or downstream scripts and should not be casually changed.

Character printing supports scalars and rank-1 fixed-width vectors through the
intrinsic kernels. `MQC_Variable%print` delegates to these routines, so fixes at
this layer should benefit both intrinsic and object use.

Logical rank-1 vectors and rank-2 matrices are available through `mqc_print`.
Logical vectors are also available through the rank-specific
`mqc_print_vector` generic.

## Use-association caution

`MQC_General` exports a broad set of names and several later modules re-export
them through use association. When a downstream unit uses multiple MQCPack
modules, prefer `only:` lists if a generic becomes ambiguous. Do not bypass an
object-layer generic merely because its intrinsic implementation routine is
visible.
