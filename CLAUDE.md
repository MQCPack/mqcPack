# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**MQCPack** (Merced Quantum Chemistry Package) is an object-oriented Fortran 2003+ library for quantum chemistry methodology development. It provides algebra objects, electronic structure theory (EST) data structures, and interfaces to quantum chemistry codes (primarily Gaussian).

## Build & Install

```bash
# Interactive installation (prompts for compiler, install path, BLAS/LAPACK, GauOpen)
./mqc_install

# After configuration, standard autotools build
make
make install

# Run tests
make check
```

Supported compilers: `gfortran` (v8+), `pgfortran`/`nvfortran` (v19+).  
Required dependencies: autotools, BLAS, LAPACK.  
Optional: GauOpen (Gaussian matrix file I/O), Gaussian 16+.

Two configure variants exist:
- `configure_FormCHK` — default, uses formatted checkpoint files
- `configure_MatrixFile` — adds GauOpen support for unformatted matrix files

## Source Layout

All Fortran sources live in `src/`. The library compiles to `src/libmqc.a`. Compilation order is controlled by `src/Makefile.am`. Working examples are in `examples/`, unit tests in `unitTests/`.

## Module Architecture

Modules are layered; lower layers must be compiled first:

**Foundation (Layer 0)**
- `mqc_general.F03` — constants, math utilities, error handling
- `mqc_binary.F03` — bit manipulation
- `mqc_datastructures.F03` — linked lists, READONLY lists

**Mathematics (Layer 1)**
- `mqc_algebra.F03` (~34K lines) — core algebra objects: `MQC_Scalar`, `MQC_Vector`, `MQC_Matrix`, `MQC_R4Tensor`; wraps LAPACK; handles packed/unpacked storage transparently
- `mqc_algebra2.F03` — dynamic-rank `MQC_Variable` arrays
- `mqc_integrals.F03` — basis function integrals (Obara-Saika recursion)

**Quantum Chemistry (Layer 2)**
- `mqc_est.F03` (~31K lines) — electronic structure objects: `MQC_Wavefunction`, `MQC_PSCF_Wavefunction`, `MQC_SCF_Integral`, `MQC_Determinant`, `MQC_TwoERIs`/`MQC_TwoERIset`; auto-detects restricted/unrestricted/general spin frameworks
- `mqc_gaussian.F03` — reads/writes Gaussian FChk files and matrix files
- `mqc_molecule.F03` — nuclear structure, molecular data
- `mqc_fullwavefunction.F03` — complete wavefunction representation

**Interface (Layer 3)**
- `mqc_interface.F03` — C-to-Fortran bindings
- `mqc_matwrapper.F03` — matrix file wrapper utilities
- `mqc_files.F03` — file I/O abstractions
- `mqc_util.c` — C utility functions

**GauOpen** (optional, `src/GauOpen/`): `qcmatrix_4/8.F` and `qcmatrixio_4/8.F` compiled as objects linked into `libmqc.a`.

## Key Design Patterns

- **Operator overloading**: arithmetic operators (`+`, `-`, `*`, `/`), comparison, and assignment are overloaded for all `MQC_*` types via Fortran's `INTERFACE` blocks in `mqc_algebra.F03`.
- **Spin framework transparency**: EST routines detect R/U/G frameworks at runtime; callers use the same API regardless of spin treatment.
- **Array packing**: `MQC_Matrix` transparently handles packed (triangular) vs. full storage; routines like `MQC_Matrix_Pack`/`MQC_Matrix_Unpack` are internal.
- **Two-electron integral storage**: `MQC_TwoERIs` stores ERIs in multiple formats (full, packed, sparse). `MQC_TwoERIset` is a container for multiple ERI objects (e.g., alpha/beta blocks). The `mqc_est.F03` module is currently being updated to ensure subroutines work across all storage types.
- **C interop**: selected routines use `ISO_C_BINDING` via `mqc_interface.F03`.

## CI/CD

GitHub Actions run on push/PR:
- `.github/workflows/linux.yml` — Ubuntu, gfortran + nvfortran, downloads GauOpen v2
- `.github/workflows/macos.yml` — macOS, gfortran via Homebrew

## Current Development Branch

Active branch `mmorato-dev` has changes in `src/mqc_est.F03` related to updating `MQC_TwoERIs` storage and iteration types. Subroutines are still being modified to handle different storage cases correctly.
