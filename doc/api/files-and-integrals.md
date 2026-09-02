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

- `init`: set angular momentum, center, contraction coefficients, primitive
  exponents, and an optional angular representation;
- `print`: report the object;
- `shell2nBasis`: return the exposed/requested number of basis functions;
- `shell2nCartesian`: return the complete Cartesian working dimension;
- `getAngularRepresentation`: return `MQC_CGTF_CARTESIAN` or
  `MQC_CGTF_REAL_PURE`;
- `getCartesianToBasis`: return the `nCartesian`-by-`nBasis` angular
  transformation;
- `primitiveValues`: evaluate one individually normalized stored primitive in
  the shell's requested representation, excluding both its contraction
  coefficient and the contracted-shell renormalization;
- `primitiveSelfOverlap`: compute primitive self-overlap information.

`MQC_gtoBasisSet` owns a collection of `MQC_CGTF` shells and tracks basis/shell
counts. Its public `nBasis` count is the exposed dimension and `nCartesian` is
the sum of complete Cartesian working dimensions. Use `init`, `addShell`, and
`print` to maintain its allocation and count invariants. The optional trailing
representation argument to `addShell` permits different shells in one basis
set to use different representations.

For a Gaussian MatrixFile, use
`MQC_Gaussian_Unformatted_Matrix_File%loadGTOBasisSet` to construct this
evaluable object from the file's signed shell types, primitive data, and shell
coordinates. Gaussian record labels and SP-shell expansion remain owned by
`MQC_Gaussian`; `MQC_Integrals` owns the resulting normalized shell objects.

Each shell is homogeneous: all of its exposed components are Cartesian or all
are Gaussian-convention real pure. Cartesian shells use an identity
transformation. For real-pure shells, values are formed as
`transpose(T) * cartesianValues`, where `T` is returned by
`getCartesianToBasis`. Complete Gaussian-ordered Cartesian `lArrays` remain
available as the working representation. `MQC_Value_CGFT` and
`basisSetValuesList` apply the transformation for contracted values at a point;
`MQC_Value_Primitive_Radial` and `MQC_Value_Primitive_Angular` retain their
low-level radial and scalar Cartesian-monomial meanings.

Pure-shell analytical overlap transformation is not part of this value-only
API increment. `MQC_Overlap_CGFT` therefore rejects real-pure D and higher
shells explicitly instead of interpreting pure indices as Cartesian indices.

Angular momentum, Cartesian center, coefficient/exponent pairing, primitive
normalization, and shell ordering are scientific data. Preserve them when
mapping to Gaussian basis records or EST basis matrices.

## `MQC_Molecule`

Source: [`src/mqc_molecule.F03`](../../src/mqc_molecule.F03).

`MQC_Molecule_Data` holds atom count, atomic numbers, masses, nuclear charges,
and Cartesian coordinates. Its public methods are `print`, `getNucRep`,
`getNumAtoms`, and `updateMolData`. Coordinate units and the orientation of the
3-by-`NAtoms` coordinate matrix must be explicit at call boundaries.

Element-property lookup procedures are module-level APIs rather than stored
components of every molecule object:

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
`src/mqc_molecule.F03`; they are not generated or included from another file.

These routines intentionally do not define a generic, context-free “atomic
radius.” Covalent, van der Waals, ionic, and integration radii have different
scientific meanings and must use separate, source-identified APIs if they are
added later. The present API is the minimum atomic-size boundary required for
an original-Becke reference implementation.

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
