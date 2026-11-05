# Quadrature grids and rules

`MQC_Integrals` defines `quadrature_grid` and its extension
`molecular_quadrature_grid`. Both own allocatable storage. Every component is
private; implementation helpers are explicitly private in the host module.
There are no public pointer views, callbacks exposing storage, or mutable aliases. Intrinsic assignment
makes an independent deep copy; allocatable storage is automatically finalized.

Source authority: [`src/mqc_integrals.F03`](../../src/mqc_integrals.F03)
for declarations/types/bindings and
[`src/mqc_integrals_quadrature.F03`](../../src/mqc_integrals_quadrature.F03)
for contained procedures and Lebedev tables. The latter is a Fortran INCLUDE,
not a separately compiled module. The migration preserves the previously
validated gaussianIntegrals algorithms; it sets no convergence defaults.

This API entered the source release in version 26.9.3 on 2026-09-19.
Earlier validated migration installations reported 26.9.2; older builds with
that version do not necessarily include it. Use compiler-matched installations
and record source provenance.
The mathematical API works with or without GauOpen; Gaussian file loading and
PAD grid-selection policy belong to downstream callers.

## Construction

Constructors are module procedures whose first argument is the output object.
They replace its previous contents. Base constructors require the exact base
type, preventing accidental reconstruction of just the inherited portion of a
molecular object. Use `select type` when constructing a polymorphic variable.
All counts and dimensions use `integer(int64)` and coordinates/weights use
`real(real64)`.

| Procedure (arguments following `grid`) | Meaning |
| --- | --- |
| `gridTrapezoid1D(n,step,origin,units)` | Explicit uniform interval; endpoint half weights |
| `gridGaussLegendre1D(n,lower,upper,units)` | Finite interval; exactness degree `2*n-1` |
| `gridMappedRadial(n,scale,units)` | Semi-infinite mapped Gauss–Legendre; weights integrate `dr`, without `r**2` |
| `gridCartesian(n,step,origin,units)` | Existing equal-spacing, equal-count three-axis trapezoid product; origin is lower corner |
| `gridLebedev(n)` | Exact supported point count; unit directions and solid-angle weights summing to `4*pi` |
| `gridSphericalProduct(origin,radialPoints,radialWeights,angularPoints,angularWeights,units)` | Explicit supplied rules; applies `r**2` once |
| `gridFiniteSpherical(originBohr,rMaxBohr,nRadial,nAngular,radialRule)` | Existing finite Lebedev spherical convenience constructor |
| `gridBeckeMolecular(centersBohr,radiiBohr,radialPointsBohr,radialWeightsBohr,angularPoints,angularWeights)` | Molecular output type; unpruned concatenation and original-Becke partition |
| `gridFromArrays(points,weights,integrationDimension,coordinateUnits,weightUnits,rule,normalization)` | Copy external data and explicit metadata |
| `gridTakeArrays(points,weights,integrationDimension,coordinateUnits,weightUnits,rule,normalization)` | Transfer allocated arrays with `move_alloc`; source arrays become unallocated |

Points have shape `(coordinateDimension,nPoints)`, weights `(nPoints)`.
Imported rules require nonempty, finite, conformable arrays and explicit,
nonblank metadata strings. Imports preserve the supplied measure; they do not
renormalize weights or convert units. Negative and zero weights are valid.
Metadata supplied by an importer is its responsibility.

Element lookup is available in the same module. For example, use
`MQC_Integrals::becke1988_radii_from_atomic_numbers` and then call
`gridBeckeMolecular` with the returned bohr radii. True atomic numbers, not
ECP effective charges, identify elements. The convenience procedure obtains radii from `MQC_General`.
The quadrature procedures require only intrinsic facilities and MQC_General;
they do not use molecule, basis, or file objects.

## Operations and access

- `%isInitialized()` and `%pointCount()` are valid on an empty object; they
  return false and zero. `%clear()` releases owned storage and resets metadata,
  including molecular storage when called through a polymorphic base.
- `%validate()` checks an initialized object. Uninitialized scientific
  operations, invalid indices/shapes, and nonfinite data fail via `mqc_error`.
- `%getMetadata(coordinateDimension,integrationDimension,coordinateUnits,
  weightUnits,rule,normalization,exactnessDegree)` has optional output arguments.
  String outputs are allocatable characters and independent copies. A Lebedev
  surface has coordinate dimension 3 and integration dimension 2. Degree `-1`
  means unknown/not claimed; no Cartesian polynomial exactness is inferred for
  mapped, supplied, molecular, or finite spherical grids.
- `%copyBlock(first,points[,weights])` fills caller-owned buffers. The second
  extent of `points` determines the requested count. Indices are one-based and
  blocks must be nonempty and within the grid. Buffers cannot modify the object.
- `%integrate(values)` returns the real weighted sum. Values have one entry per
  point and include the integrand only, not another Jacobian or weight.
- `%overlapMatrix(values,overlap)` integrates products of real function columns,
  with values shaped `(nPoints,nFunctions)`. The allocated output is symmetric.

Molecular objects additionally expose `%atomCount()`, `%copyMolecule(centersBohr,
radiiBohr)`, `%copyAtomicRule(radialPointsBohr,radialWeightsBohr,angularPoints,
angularWeights)`, `%partitionRange(minimum,maximum)`, and
`%copyPartitionBlock(first,generatingAtoms,baseWeightsBohr3,partitionWeights)`.
The last operation fills caller buffers with equal lengths. Generating atoms
are one-based input-center indices, not atomic numbers or exclusive point
ownership. Final weights equal base volume weights times partition weights.
The exact supplied atomic rules and radii are retained for reproducibility;
no named radial mapping is inferred from arbitrary supplied arrays.

## Consumer example

```fortran
use iso_fortran_env, only: real64,int64
use MQC_Integrals
implicit none
type(quadrature_grid) :: grid
real(real64) :: points(1,4), integral
call gridGaussLegendre1D(grid,4_int64,-1.0_real64,1.0_real64,'bohr')
call grid%copyBlock(1_int64,points)
integral = grid%integrate(points(1,:)**2)
```

For large grids, allocate a bounded point buffer once and iterate over blocks.
The downstream gaussianIntegrals procedure
`gbs_mod::evaluateMOsOnGridObject` demonstrates this with at most 4096 points per
block, forwarding them to the existing intrinsic-array MO evaluator. It accepts
`class(quadrature_grid)` and requires three-dimensional coordinates in bohr.
The complete MO-value matrix is still retained for overlap integration, as it
was before this API. This change bounds coordinate scratch storage; it does
not introduce a fully streamed overlap algorithm.

The downstream `quadratureTest.exe` constructs the molecular object directly for Becke
and takes ownership of existing PAD-builder arrays for other choices. It uses
the block consumer and the object's overlap integration without altering CLI,
scientific output fields, or PAD production policy. For imported PAD grids,
detailed provenance remains in the existing setup output; the object does not
infer the rule selected by `auto`.

## Validation and current limits

`unitTest12` covers exact moments, procedural parity, negative Lebedev
weights, measure metadata, ownership transfer, independent copies, lifecycle,
polymorphic molecular integration, retained construction data, and rejection
of invalid input. It is registered with MQCPack `make check`, together with `unitTest10` (rules),
`unitTest11` (partition and molecular construction), and invalid-input/privacy
scripts. gaussianIntegrals retains its unit/PAD tests and FAF-backed
`make test-molecular-grid` convergence checks as downstream integration coverage.

This initial API supports real sampled integration and real overlap matrices.
Complex integration, pruning, general unequal-axis Cartesian products, and
production default selection remain future work. Finite/supplied spherical
base objects retain the resulting points and weights but not a separate full
construction recipe; molecular objects retain their supplied atomic rules.
The two provided types are supported by `%clear`; additional downstream type
extensions would require extending that lifecycle implementation.


## Procedural API

The owning types are optional; intrinsic-array constructors remain public.
All counts/indices are `int64` and real data are `real64`. For exact signatures
and intents, see the source/index. Fixed-size outputs must be allocated by the
caller; allocatable `intent(out)` results replace their previous contents.

| Entry point | Output and conventions |
| --- | --- |
| `setup_quadrature_trapezoid1d(nPoints,stepSize,origin,quadGrid,quadWeights)` | Caller-sized vectors, positive spacing, endpoints included |
| `setup_quadrature_trapezoid3d(nPoints,stepSize,origin,quadGrid,quadWeights)` | Caller-sized `(3,nPoints**3)` and weight vector; lower-corner origin |
| `setup_quadrature_gauss_legendre1d(nPoints,lowerBound,upperBound,quadGrid,quadWeights)` | Caller-sized finite-interval rule |
| `setup_quadrature_gauss_legendre_semi_infinite(nPoints,radialScale,quadGrid,quadWeights)` | Caller-sized vectors; `r=R*(1+x)/(1-x)`, weights include `dr/dx`, not `r**2` |
| `setup_lebedev_angular_grid(lebedevPoints,angularGrid,angularWeights)` | Allocated unit directions and solid-angle weights |
| `setup_quadrature_spherical_product(origin,radialGrid,radialWeights,angularGrid,angularWeights,quadGrid,quadWeights)` | Allocated points and volume weights; angular directions must be unit length |
| `setup_quadrature_lebedev_spherical(originBohr,rMaxBohr,nRadial,lebedevOrder,radialRule,quadGrid,quadWeights)` | Finite convenience grid in bohr; trapezoid or finite Gauss-Legendre radial rule |
| `becke1988_partition_weights` | Dimensionless `(nAtoms,nPoints)` partition weights for explicit centers/radii and points in bohr |
| `apply_becke1988_partition` | Partition and final weights for caller-supplied generating-atom indices and base weights |
| `setup_quadrature_becke_molecular` | Unpruned concatenated atomic grids, generating atoms, base/partition/final weights |
| `becke1988_radii_from_atomic_numbers` | Allocated radii in bohr from true atomic numbers through MQC_General |
| `setup_quadrature_becke_molecular_from_atomic_numbers` | Convenience construction with element-radius lookup |

`QUAD_RADIAL_TRAPEZOID=0`, `QUAD_RADIAL_GAUSS_LEGENDRE=1`, and
`QUAD_RADIAL_MAPPED_GAUSS_LEGENDRE=2` identify radial rules. The finite spherical
convenience constructor accepts only the first two; use an explicit mapped
radial rule with the supplied-rule product or molecular constructor.

`N_LEBEDEV_RULES`, `LEBEDEV_POINT_COUNTS`, and `LEBEDEV_EXACTNESS_DEGREES`
provide the 23-rule catalogue (6 through 2030 points, degrees 3 through 77).
`lebedevRuleIsSupported` queries exact availability; `lebedevRuleDegree` returns
the supported rule's degree. `lebedevSelectPointCount(requestedPoints,
roundUp,capAtMaximum)` selects a supported rule, rounding up by default;
requests above the maximum fail unless explicitly capped. Preserve signed
weights; some valid Lebedev rules contain negative weights.

Low-level procedural builders retain their established validation contracts.
The object constructors add finite-data and some count-overflow checks; this
migration does not claim identical defensive checks for every raw-array path.
