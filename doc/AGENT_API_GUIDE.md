# MQCPack agent API guide

This guide is the starting point for agents and developers working in a
repository that consumes MQCPack. It explains where the supported abstractions
live, which source files define them, and which scientific and binary-format
conventions must be preserved.

The source remains authoritative. This guide supplies navigation and semantics;
it does not replace checking the exact declaration used by the MQCPack version
or Git commit on which a downstream repository depends.

## Before reading or changing dependent code

1. Resolve the MQCPack source root. In many local environments it is available
   as `$mqcroot`; do not assume that variable exists everywhere.
2. Record the MQCPack version and, when reproducibility matters, the Git commit.
   `MQC_General` provides `mqc_version`, `mqc_version_check`, and
   `mqc_version_print`; the source commit distinguishes worktrees with the same
   released version.
3. Read the relevant curated page below, then inspect the exact type, generic
   interface, and implementation in `src/`.
4. Compile against `.mod`, object, GauOpen, and library files produced by one
   compiler family and compatible flags. Never mix GNU, Intel, and NVidia
   Fortran artifacts.
5. Preserve units, basis, rank, storage, spin-block ordering, record labels, and
   integer precision as API data—not incidental implementation details.

## Where to start

| Downstream task | Read first | Source authority |
| --- | --- | --- |
| Constants, errors, printing, conversions, packed arrays, BLAS/LAPACK wrappers | [General utilities](api/general.md) | `src/mqc_general.F03`, `src/mqc_general_lapack.F03` |
| Scalars, vectors, matrices, rank-4 tensors, operators, eigensystems, SVD | [Algebra objects](api/algebra.md) | `src/mqc_algebra.F03` |
| Rank-changing numerical or character values | [MQC_Variable](api/algebra2.md) | `src/mqc_algebra2.F03` |
| SCF/post-SCF intermediates, spin blocks, determinants, two-electron integrals | [Electronic-structure objects](api/est.md) | `src/mqc_est.F03` |
| FChk or FAF/MatrixFile reading and writing | [Gaussian interfaces](api/gaussian.md) | `src/mqc_gaussian.F03`, `src/mqc_matwrapper_direct.F03` |
| Text files, linked lists, basis functions, molecule data, C interoperability | [Supporting modules](api/files-and-integrals.md) | corresponding files under `src/` |
| Find an exact declaration or type-bound implementation | [Generated procedure index](api/PROCEDURE_INDEX.md) | all principal `src/*.F03` files and `src/mqc_util.c` |

## Layering and dependency direction

The intended direction is:

```text
MQC_General
  ├─ MQC_Binary / MQC_DataStructures / MQC_Files / MQC_Integrals
  ├─ MQC_Algebra
  │    ├─ MQC_Algebra2
  │    ├─ MQC_Molecule
  │    └─ MQC_EST
  └─ MQC_MatWrapper
       └─ MQC_Gaussian
            └─ downstream programs and higher-level wavefunction handling
```

Reusable operations on intrinsic Fortran values belong in `MQC_General`.
Algebra modules validate object state, dispatch on type/rank/storage, and then
delegate where appropriate. Gaussian word counts, `TypeA`, `Len4L`, record
labels, and GauOpen semantics stay at the MatrixFile boundary.

## What counts as the public API

Prefer, in order:

- documented derived types and their public type-bound procedures;
- generic interfaces such as `mqc_print`, `MatMul`, `Allocated`, `MQC`, and
  overloaded intrinsic or operator names;
- named module procedures that are explicitly used by examples or documented
  as entry points.

Many modules do not set a module-wide `private` default. Consequently, some
implementation routines are technically accessible through `use` association.
That does not make every such routine a stable downstream contract. Avoid a
specific implementation name when a generic or type-bound entry point exists.

The generated procedure index is deliberately called an inventory. It helps an
agent locate code, but does not upgrade internal routines into supported API.

## Scientific invariants

Before calling or changing MQCPack code, identify all applicable invariants:

- numerical kind: new interfaces should use `int32`, `int64`, and `real64`
  explicitly;
- units: especially Hartree versus eV, length units, and Gaussian Cartesian
  coordinates;
- basis: AO, MO, spin orbital, or generalized spin;
- spin blocks: alpha, beta, alpha-beta, and beta-alpha are distinct;
- storage: full, symmetric, diagonal, Hermitian, antisymmetric, packed, or
  specialized MatrixFile storage;
- shape and rank: do not infer them only from element count;
- file labels: FChk and FAF labels are compatibility keys;
- numerical thresholds: use `mqc_small` for zero/threshold decisions;
- resource ownership: an open file unit is connection state, not ordinary
  copyable metadata.

## Common downstream patterns

### Algebra objects

```fortran
program matrix_example
  use iso_fortran_env, only: int64
  use MQC_Algebra
  implicit none

  type(MQC_Matrix) :: overlapAO

  call overlapAO%identity(3_int64,3_int64)
  call overlapAO%print(6_int64,'AO overlap matrix')
end program matrix_example
```

Check the exact overload before supplying optional arguments; many established
interfaces also accept intrinsic integer, real, complex, or MQCPack scalar
values.

### MatrixFile character data

```fortran
program read_faf_title
  use MQC_Algebra2
  use MQC_Gaussian
  implicit none

  type(MQC_Gaussian_Unformatted_Matrix_File) :: faf
  type(MQC_Variable) :: title

  call faf%load('input.faf')
  call faf%getArray('TITLE',mqcVarOut=title)
  call mqc_print(title,header='TITLE')
  call faf%closeFile()
end program read_faf_title
```

MatrixFile support must be present in the configured MQCPack build. Preserve
the exact record label and verify binary-format changes with both the public API
and an independent Gaussian-aware reader such as `dumpbaf` when available.

## Build and ABI compatibility

MQCPack's Fortran `.mod` files and objects are compiler-specific. Default
integer/real promotion flags are also part of the effective ABI. A dependent
repository should record at least:

- MQCPack version and preferably source commit;
- compiler family and version;
- `--without-gauopen` or the exact `--with-gauopen=DIR` configuration;
- GauOpen integer ABI, currently `--with-gauopen-integer-bytes=8`;
- integer/real-size and OpenMP flags;
- BLAS/LAPACK link line;
- installation prefix providing the matching `libmqc.a` and `.mod` files.

Changing a public derived type changes its compiled module layout. Rebuild all
dependent objects after such a change; a stale object can compile cleanly and
then fail at runtime.

MQCPack has one tracked Autotools configuration. MatrixFile builds are selected
with `--with-gauopen=DIR`; FormChk-only builds use `--without-gauopen`. Do not
select a build variant by replacing or symlinking `configure.ac` or
`Makefile.am` files.

MatrixFile builds compile the external, unmodified GauOpen `qcmatrix.F` and
`qcmatrixio.F` sources directly. Dependency-only legacy flags are kept
separate from MQCPack's Fortran flags. The supported GauOpen integer ABI is
currently selected at configure time and limited to 8 bytes; one installed
library does not dispatch between 4-byte and 8-byte MatrixFiles at runtime.

In a FormChk-only build, `MQC_Gaussian` and its formatted-checkpoint API remain
linkable. The MatrixFile type remains declared for compile-time source
compatibility, but MatrixFile operations terminate with a clear unsupported-
feature diagnostic. Record the configuration explicitly instead of testing
MatrixFile availability by calling an operation.

## Guidance to place in a dependent repository

```text
MQCPack source is available at $mqcroot.

Before changing MQCPack-dependent code, read:
  $mqcroot/doc/AGENT_API_GUIDE.md

Use its module guide to locate the relevant interface, then inspect the exact
declaration under $mqcroot/src. Preserve units, basis, rank, storage format,
spin-block ordering, Gaussian labels, and compiler compatibility. Record the
MQCPack version or Git commit used by this repository.
```

If the dependent project pins a particular installation or commit, replace the
generic `$mqcroot` language with that exact path and revision.

## Maintaining this documentation

Update a curated page whenever a public type, generic, scientific convention,
file layout, or recommended use changes. Regenerate the inventory after source
procedure, interface, type, or type-bound binding changes:

```sh
python3 doc/api/generate_procedure_index.py
python3 doc/api/generate_procedure_index.py --check
```

The first command updates the file; the second is suitable for a consistency
check. The scanner has no third-party dependencies.
