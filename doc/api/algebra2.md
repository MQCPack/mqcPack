# `MQC_Variable`

Source authority: [`src/mqc_algebra2.F03`](../../src/mqc_algebra2.F03).

`MQC_Algebra2` provides `MQC_Variable`, a flexible container whose numerical
rank can change at runtime. It supports integer, real, and complex numerical
data and character scalars/fixed-width vectors.

## State and access

Preferred type-bound entry points include:

- `init`, `clear`, `put`, `putMQC`, and `getVal`;
- `getRank`, `getType`, and `getCharacterLength`;
- `print` and `isConformable`;
- numerical methods `trimZero`, `trace`, `iPower`, `rPower`, `column`,
  `eigen`, `diag`, `det`, and `subMatrix`;
- character methods `change_case`, `adjustl`, `adjustr`, `len_trim`, and
  scalar `trim`.

Generic interfaces include `MQC`, `RANK`, `SIZE`, `RESHAPE`, `LEN`,
`ADJUSTL`, `ADJUSTR`, `LEN_TRIM`, `TRIM`, `String_Change_Case`, `INT`,
`FLOAT`, `Contraction`, `Dot_Product`, `MatMul`, `Transpose`, assignment, and
the supported arithmetic/operator families.

The module also exposes public `dataType` and `storageFormat` metadata, but
consumers should prefer the query and initialization interfaces rather than
manually constructing an inconsistent object.

## Character contract

Character `MQC_Variable` values currently support only:

- rank zero: one character scalar;
- rank one: a vector whose elements all have the same declared character
  length.

The backing store is a deferred-length allocatable character array. Scalar
character data uses a one-element backing array while the object reports rank
zero. Preserve that distinction and the fixed-width vector invariant.

Character objects support assignment, `RANK`, `SIZE`, `LEN`, reshape within
the supported ranks, `MQC`, `put`, `getVal`, copying, printing, case changes,
`ADJUSTL`, `ADJUSTR`, `LEN_TRIM`, and scalar `TRIM`. Numerical operators must
not silently accept or promote character type.

FAF/MatrixFile word sizes and `TypeA` rules do not belong in this module.
Standalone character variables may have widths that cannot be encoded as an
exact MatrixFile record; encoding validation occurs in `MQC_Gaussian` and
`MQC_MatWrapper`.

## Intrinsic/object layering

Reusable character algorithms operate on intrinsic values in `MQC_General`:
`mqc_adjustl`, `mqc_adjustr`, `mqc_len_trim`, `mqc_trim`, and the character
printing kernels. `MQC_Algebra2` validates object state and rank, then delegates
to those kernels. Follow the same design for extensions so intrinsic and object
callers receive consistent behavior.

## Test reference

`unitTests/src/unitTest02.f03` is the focused behavioral reference for character
scalar/vector assignment, rank/size/length, reshape, put/get, clear, copy,
printing, case conversion, adjustment, trimmed length, and scalar trimming.
Use it when determining whether a proposed change preserves the current
contract.
