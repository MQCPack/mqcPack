# Algebra objects

Source authority: [`src/mqc_algebra.F03`](../../src/mqc_algebra.F03).

`MQC_Algebra` provides strongly shaped numerical objects and overloaded
intrinsic-like operations. Use it when rank and storage are conceptually stable:
`MQC_Scalar`, `MQC_Vector`, `MQC_Matrix`, and `MQC_R4Tensor`.

## Primary types

### `MQC_Scalar`

Holds one integer, real, or complex value. Preferred type-bound operations
include `print`, `rval`, `ival`, `cval`, `abs`, `phase`, `sqrt`, `exp`, `log`,
and `random`. Assignment and arithmetic are overloaded for intrinsic and
MQCPack scalar operands.

### `MQC_Vector`

Holds a rank-1 integer, real, or complex array and tracks row/column
orientation. Important methods include:

- construction/state: `init`, `set`, `size`, `print`;
- access/mutation: `at`, `vat`, `put`, `vput`, `getArrayI`, `getArrayR`,
  `getArrayC`;
- sequence operations: `push`, `unshift`, `pop`, `shift`;
- algebra: `norm`, `transpose`, `dagger`, `sum`, `product`, `power`;
- element operations: `maxval`, `minval`, `maxloc`, `minloc`, `argsort`,
  `sort`, `sqrt`, `exp`, `abs`;
- diagonal construction: `diag`/`diagf`.

### `MQC_Matrix`

Holds rank-2 integer, real, or complex data with explicit storage metadata.
Important methods include:

- construction/state: `init`, `identity`, `set`, `print`;
- access/mutation: `at`, `vat`, `mat`, `put`, `vput`, `mput`, `mpad`;
- transformations: `transpose`, `dagger`, `sqrt`, `exp`, `power`;
- decompositions/solves: `diag`, `svd`, `eigensys`, `inv`, `pinv`;
- invariants/inspection: `norm`, `det`, `minor`, `cofactor`, `trace`,
  `diagonal`, `rmsmax`, `sum`, `psum`.

### `MQC_R4Tensor`

Holds rank-4 integer, real, or complex data. Its public methods cover
initialization, printing, scalar/vector/matrix/tensor access (`at`, `vat`,
`mat`, `tat`) and corresponding mutation (`put`, `vput`, `mput`, `tput`).
Tensor contractions and outer products are also available through module
generics.

## Generics and operators

The module overloads assignment, arithmetic, comparisons, `MatMul`,
`Dot_Product`, `Transpose`, `Dagger`, `Contraction`, `PartialContraction`,
`Reshape`, `Size`, `Allocated`, `Sum`, `Sqrt`, `Exp`, `Cmplx`, `real`,
`aimag`, and several trigonometric functions. It also defines MQCPack-specific
operators for element-wise, dot, outer, cross-outer, and contraction-related
operations.

Always inspect the generic block for the exact supported operand pair. A
familiar intrinsic name does not imply that every intrinsic rank/type
combination is implemented for MQCPack objects.

## Storage conventions

Matrix and tensor storage metadata is part of the object contract. Established
matrix values include full, symmetric, diagonal, antisymmetric, Hermitian, and
anti-Hermitian forms (spelled with the `Stor...` names used in the source).
Do not directly infer a private backing-array shape or change storage metadata
without using the object's initialization/conversion routines.

The source notes that lower-triangular matrix data is stored row by row. It also
permits diagonal elements in some symmetry storage modes for practical storage
purposes, so do not assume that a storage flag alone proves every mathematical
property of the represented matrix.

## Choosing between algebra layers

Use `MQC_Algebra` when code benefits from a fixed scalar/vector/matrix/tensor
type and storage-aware operations. Use [`MQC_Variable`](algebra2.md) when rank
must change dynamically, a single container must hold several intrinsic types,
or character scalar/vector data is required.

Do not reach into private vector, matrix, or tensor arrays. Use `getArray*`,
`at`, and the documented mutation methods so allocation, type, orientation, and
storage invariants remain consistent.
