# Electronic-structure objects

Source authority: [`src/mqc_est.F03`](../../src/mqc_est.F03).

`MQC_EST` builds quantum-chemistry intermediates on the algebra layer. Its
array labels, basis interpretation, spin layout, and storage conventions are
scientific API.

## Principal types

### `MQC_SCF_Integral`

Represents one-electron operators or related matrices using up to four spin
blocks. The source defines the generalized block arrangement as:

```text
| alpha       beta-alpha |
| alpha-beta  beta       |
```

The corresponding object components are named `Alpha`, `Beta`, `AlphaBeta`,
and `BetaAlpha`. Never interchange the off-diagonal blocks merely because their
dimensions happen to match.

Public methods cover printing, type/dimension/label queries, initialization,
block extraction, identity construction, diagonalization, SVD, generalized
eigensystems, inverse, trace, determinant, norm, powers, pseudoinverse,
element access/mutation, energy-list management, padding, sums, orbital
selection/reordering/combination/update, off-diagonal block swapping, outer
products, and diagonal extraction. See the generated index for exact binding
names and implementations.

`Array_Type` distinguishes established representations such as space, spin,
and general. `Array_Name` is a semantic label. Preserve both through
transformations.

### `MQC_SCF_Eigenvalues`

Stores alpha/beta eigenvalue vectors plus array name/type metadata. Its public
surface covers initialization, labels, block access, power, mutation,
append/access, outer products, extrema/locations, and absolute values.

### Wavefunction and basis containers

- `MQC_Wavefunction` groups MO coefficients, MO energies and symmetries, core
  Hamiltonian, Fock, density and SCF density matrices, overlap, electron and
  basis counts, charge, multiplicity, basis/symmetry labels, spin type, and a
  complex/real flag.
- `MQC_Basis_Set` stores the Gaussian MatrixFile-oriented shell-to-atom map,
  shell types, primitive counts/exponents, contraction coefficients, and shell
  coordinates as MQCPack matrices.
- `MQC_PSCF_Wavefunction` extends the single-reference wavefunction with
  core/valence/active/frozen counts, amplitudes, and energies.

### Determinants and two-electron integrals

- `MQC_Determinant` stores occupation strings through `MQC_Bits`, determinant
  counts/order, basis/core/virtual counts, and substitution metadata.
- `mqc_twoERIs` stores spin-resolved rank-4 tensors with integral-type and
  storage-type metadata. Use `getBlock`, `blockSize`, `type`, `at`, and
  `print` rather than reaching into private spin blocks.
- `mqc_twoERISet` owns multiple two-electron-integral objects and provides
  `erinum`, `eris`, `transform`, `addtoset`, and `print`.

## Generic operations

The module extends `mqc_print`, `Allocated`, `MatMul`, arithmetic and
comparison operators, contractions, outer products, and other algebraic
interfaces to EST types. Inspect the interface block before relying on an
operand pairing: supported spin structures and return types are more
constrained than an intrinsic matrix expression might suggest.

## Checklist for EST changes

For every input and result, state or determine:

- AO, MO, spin-orbital, or generalized-spin basis;
- restricted, unrestricted, or general wavefunction convention;
- real or complex data;
- active spin blocks and their row/column dimensions;
- full/symmetric or specialized integral storage;
- whether energies are Hartree and whether indices include core/frozen spaces;
- semantic label expected by Gaussian or downstream code.

Broad changes in this module can propagate through many overloads. Prefer small
changes and validate a representative restricted, unrestricted, and general
case whenever the changed code can reach all three.
