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
- Geometry kernels: `mqc_distance(point1,point2)` returns the Euclidean
  distance between equal-length intrinsic real vectors, and
  `mqc_distance_matrix(points)` returns all pairwise distances for points
  stored by column in `points(nDimensions,nPoints)`.
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
- Distance results retain the coordinate units supplied by the caller. The
  distance routines do not assume bohr, angstrom, or three dimensions.
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

## Element-radius lookups

The element-radius API is defined directly in `src/mqc_general.F03` and
requires no molecule or algebra objects:

- `mqc_element_has_bragg_slater_radius(atomicNumber)` reports whether a
  tabulated value exists. Atomic number zero identifies a ghost center and
  returns false.
- `mqc_element_bragg_slater_radius(atomicNumber)` returns the Slater empirical
  atomic radius in bohr, including Slater's 0.25-angstrom hydrogen value.
- `mqc_element_becke_1988_radius(atomicNumber)` returns the radius used for the
  original Becke molecular partition in bohr. It uses the same table except
  for Becke's 0.35-angstrom hydrogen value.

The current table covers atomic numbers 1 through
`MQC_BRAGG_SLATER_MAX_ATOMIC_NUMBER` (currently 86). Missing values terminate
through `mqc_error`; callers must not silently assign a physical radius to a
ghost or unsupported element. Select radii using true atomic numbers, not
effective nuclear charges. The radius procedures are implemented directly in
`src/mqc_general.F03`; they are not generated or included from another file.

These routines intentionally do not define a generic, context-free “atomic
radius.” Covalent, van der Waals, ionic, and integration radii have different
scientific meanings and must use separate, source-identified APIs if they are
added later. The present API is the minimum atomic-size boundary required for
an original-Becke reference implementation.

The table and procedures moved from `MQC_Molecule` without numerical or
unit changes. Callers should use `MQC_General` directly. There are no new
wrappers in `MQC_Molecule`; its existing unrestricted use association may
still expose the same entities transitively.
