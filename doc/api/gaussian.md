# Gaussian file interfaces

Source authority: [`src/mqc_gaussian.F03`](../../src/mqc_gaussian.F03), the
GauOpen backend in
[`src/mqc_matwrapper_direct.F03`](../../src/mqc_matwrapper_direct.F03),
and the FormChk-only boundary in
[`src/mqc_matwrapper_stub.F03`](../../src/mqc_matwrapper_stub.F03).

`MQC_Gaussian` is the public boundary for Gaussian formatted-checkpoint and
unformatted MatrixFile/FAF data. `MQC_MatWrapper` adapts MQCPack-prefixed
operations to unmodified external GauOpen sources and centralizes raw record
layout.

`MQC_MatWrapper` is an internal backend boundary, not a supported downstream
API. Its prefixed adapter procedures are public only so `MQC_Gaussian` and
focused MQCPack tests can use the selected real or stub implementation. The
obsolete unprefixed wrapper procedures and copied GauOpen helper routines are
not compatibility interfaces. Downstream code should use the file types and
type-bound operations in `MQC_Gaussian`.

Configure MatrixFile support explicitly with
`--with-gauopen=/path/to/gauopen --with-gauopen-integer-bytes=8`. Configure
detects one of two supported source interfaces and defines the corresponding
adapter profile:

- public: version-2 files, logical `ASym`, `Close_MatF`, and the older
  `Rd_RInd` signature;
- frontier: version-3 files, integer `TypeA`, `Close_FAF`, the newer
  `Rd_RInd` signature, and `Rd_ChBuf`/`Wr_ChBuf`.

`--with-gauopen-api=auto` is the default. An explicitly verified source tree
can be selected with `--with-gauopen-api=public` or
`--with-gauopen-api=frontier`; a requested profile that conflicts with a
recognized source interface is rejected. Use `--without-gauopen` for the
FormChk-only configuration. The integer ABI and API profile are part of the
installed library's compatibility contract; do not mix the resulting library
or module files with artifacts from another compiler or configuration. The
supported build currently uses one configured 8-byte GauOpen ABI; it does not
dispatch between 4-byte and 8-byte MatrixFiles at runtime. GauOpen's legacy
source flags are isolated from the stricter flags used for MQCPack sources.

The FormChk-only library remains linkable through `MQC_Gaussian`, including
`MQC_Gaussian_FChk_File` and the FChk search routines. The MatrixFile type is
still declared so shared source can compile, but any MatrixFile read, write,
or layout operation terminates with an explicit diagnostic that MQCPack was
configured without GauOpen support. Callers must not use the MatrixFile type
as a runtime capability probe.

## Public file types

### `MQC_Gaussian_FChk_File`

Extends `MQC_Text_FileInfo` and records the title, job type, method, and basis.
Its `OpenFile` binding opens and initializes formatted-checkpoint parsing.
`Find_FChk_Entry` and related module procedures locate typed scalar or array
entries. Treat FChk entry tags as stable compatibility keys.

### `MQC_Gaussian_Unformatted_Matrix_File`

Extends `MQC_FileInfo`. Preferred type-bound operations are:

- lifecycle: `OpenFile`, `CloseFile`, `load`, `create`, `updateHeader`;
- wavefunction classification: `isRestricted`, `isUnrestricted`, `isGeneral`,
  `isComplex`;
- molecule/basis access: `getAtomicNumbers`, `getAtomCarts`, `getAtomWeights`,
  `getAtomInfo`, `getBasisInfo`, `getBasisArray`, `getBasisData`, `getMolData`;
- typed data: `getVal`, `getValReal`, `getArray`, `getESTObj`, `get2ERIs`;
- writing: `writeArray`, `writeArray2`, `writeBasisData`, `writeESTObj`,
  `write2ERIs`.

Use the exact signature in the source because several routines have multiple
optional destinations or dispatch according to the supplied MQCPack object.

## Connection ownership and assignment

Defined assignment copies persistent MatrixFile header, molecular, layout, and
cached-scalar metadata into independent destination storage. It deliberately
does not copy a live GauOpen connection. The destination is left disconnected:
blank filename, closed state, zero unit, blank read/write mode, and unset header
I/O flags.

Assigning over an open destination closes that destination first. It must not
alter the source connection. Code that needs a connected destination must open
or create it explicitly after assignment.

## Character FAF records

Character data is exposed through `getArray(...,mqcVarOut=...)` and
`writeArray2`. Logical values are character scalars or rank-1 fixed-width
vectors represented by `MQC_Variable`. These operations require the frontier
GauOpen profile. Public GauOpen does not export character-buffer operations;
attempted character FAF layout or buffer operations terminate with a clear
unsupported-feature diagnostic rather than writing a partial character
record. Numerical MatrixFile operations remain supported in the public
profile for its version-2 files. Because the legacy label reader discards raw
positive `TypeA`, the public profile rejects version-3 FAF files rather than
risk silently treating a frontier character record as an integer array.

Do not interpret raw positive `TypeA` as the element width without
normalization:

```text
payloadCharacters = Len4L * NTot
elementWidth       = payloadCharacters  when raw TypeA = 1
                     raw TypeA           when raw TypeA > 1
elementCount       = payloadCharacters / elementWidth
bufferCharacters   = Len4L * LenBuf
```

Raw `TypeA=1` is Gaussian's legacy marker for one scalar spanning the complete
payload. Values greater than one are explicit fixed element widths. The
translation belongs in
`MQC_MatrixFile_Get_Character_Record_Layout`; outgoing translation belongs in
`MQC_MatrixFile_Get_Character_Write_Layout`.

When scanning past an unrequested character record, the reader must consume it
with `Rd_ChBuf`, not a numeric skip routine, or later label reads become
desynchronized.

Outgoing exact-shape encoding currently requires:

- element width greater than one, because `TypeA=1` has legacy scalar meaning;
- positive element count and buffer word count;
- total payload characters divisible by the target file's `Len4L`.

These are file-format constraints, not general `MQC_Variable` restrictions.
When packing a deferred-length character buffer, initialize the full declared
length; assignment of a one-character blank can reallocate it to length one.
For real or complex array metadata, frontier `TypeA=-1` means
antisymmetric/anti-Hermitian while `TypeA=1` means symmetric/Hermitian; do not
collapse those signs when writing externally consumed records.

## MatrixFile layout types

`MQC_MatrixFile_Layout` stores the active integer width, file version, header
record count, `Len12L`, `Len4L`, and raw-format flag. Its `set` and
`setForWrite` methods establish per-file metadata.

`MQC_MatrixFile_Character_Record_Layout` holds logical element count/width and
payload/buffer sizes in both characters and file words. Keep word arithmetic
at this boundary rather than leaking it into algebra objects.

## Validation standard

For binary encoding changes, passing an MQCPack round trip alone is not enough.
Validate both:

1. public MQCPack write/read behavior, including skip and rewind/rescan paths;
2. independent Gaussian `dumpbaf` inspection when it is available.

Initialize the Gaussian/GDV/MQCPack shell environment with the repository's
documented `gdvcode2026` setup before build or Gaussian work. Never combine
GauOpen objects or `.mod` files from different compilers.
