# Supporting modules

This page covers file abstractions, linked data, contracted Gaussian basis
objects, molecule data, interoperability, and the status of the higher-level
wavefunction module.

## `MQC_Files`

Source: [`src/mqc_files.F03`](../../src/mqc_files.F03).

`MQC_FileInfo` is the abstract base carrying filename, open state, and unit
number; `IsOpen` queries connection state. `MQC_Text_FileInfo` adds a text
buffer, cursor and parsing controls, with bindings `OpenFile`, `CloseFile`,
`Rewind`, `LoadBuffer`, `GetBuffer`, `GetNextString`, `GetNextInteger`, and
`WriteLine`.

Use `MQC_getUnitNumber` and the module's unit registry rather than inventing a
fixed unit in reusable library code. Keep the registry synchronized when a
derived file type closes or replaces a connection.

Text parsing distinguishes loaded-buffer state, quote handling, repeated
delimiters, and whether a blank line means EOF. Configure those controls via
the existing parsing helpers rather than duplicating tokenization downstream.

## `MQC_DataStructures` and `MQC_Binary`

Sources: [`src/mqc_datastructures.F03`](../../src/mqc_datastructures.F03) and
[`src/mqc_binary.F03`](../../src/mqc_binary.F03).

`MQC_LinkedList` stores polymorphic values and provides typed return helpers for
integer, real, character, and `MQC_Bits` data. Because storage uses `class(*)`,
callers must preserve the expected dynamic type when retrieving values.

`MQC_Bits` represents bit strings used especially by determinant machinery.
It extends familiar bit operations through generics including `btest`,
`ibset`, `ibclr`, `iand`, `ieor`, `ior`, `popcnt`, `mvBits`, `iBits`, and
`iShft`, and provides `mqc_bits` construction and `mqc_print`/`%print`.

## `MQC_Integrals`

Source: [`src/mqc_integrals.F03`](../../src/mqc_integrals.F03).

`MQC_CGTF` represents one contracted Gaussian-type shell/function object. Its
public bindings are:

- `init`: set angular momentum, center, contraction coefficients, and primitive
  exponents;
- `print`: report the object;
- `shell2nBasis`: return the number of basis functions for the shell;
- `primitiveSelfOverlap`: compute primitive self-overlap information.

`MQC_basisSet` owns a collection of `MQC_CGTF` shells and tracks basis/shell
counts. Use `init`, `addShell`, and `print` to maintain its allocation and count
invariants.

Angular momentum, Cartesian center, coefficient/exponent pairing, primitive
normalization, and shell ordering are scientific data. Preserve them when
mapping to Gaussian basis records or EST basis matrices.

## `MQC_Molecule`

Source: [`src/mqc_molecule.F03`](../../src/mqc_molecule.F03).

`MQC_Molecule_Data` holds atom count, atomic numbers, masses, nuclear charges,
and Cartesian coordinates. Its public methods are `print`, `getNucRep`,
`getNumAtoms`, and `updateMolData`. Coordinate units and the orientation of the
3-by-`NAtoms` coordinate matrix must be explicit at call boundaries.

The nuclear-repulsion result is meaningful only when coordinates and charges
use the expected units/conventions. Do not infer units from a variable named
only `coordinates` in downstream code.

## C/Fortran interoperability

Sources: [`src/mqc_interface.F03`](../../src/mqc_interface.F03) and
[`src/mqc_util.c`](../../src/mqc_util.c).

These files bridge fixed-size Fortran character buffers and C strings for error
reporting, file-list creation, FormChk/MatrixFile discovery, file/executable
checks, scratch handling, and MatrixFile guidance. Preserve `bind(C)` names,
null termination, buffer lengths, integer conversions, and the established
Fortran-owned output path.

Prefer the Fortran wrapper routines from Fortran callers. Calling C helpers
directly creates additional responsibility for string termination and C kind
compatibility.

## `MQC_FullWavefunction`

Source: [`src/mqc_FullWavefunction.F03`](../../src/mqc_FullWavefunction.F03).

This higher-level module currently contains substantial historical and
work-in-progress code, much of it commented with active `hph` context. Do not
treat commented type sketches as public API and do not remove them as casual
cleanup. For stable downstream work, prefer the active `MQC_Gaussian` and
`MQC_EST` objects unless an active `MQC_FullWavefunction` entry point is
verified in the exact source revision.
