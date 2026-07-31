#!/usr/bin/env python3
"""Generate MQCPack's source-level API inventory.

This intentionally uses only Python's standard library.  It is a lightweight
Fortran declaration scanner, not a compiler: the curated API pages remain the
authority for semantics and recommended downstream use.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass, field
from pathlib import Path


FORTRAN_SOURCES = (
    "src/mqc_general.F03",
    "src/mqc_general_lapack.F03",
    "src/mqc_binary.F03",
    "src/mqc_datastructures.F03",
    "src/mqc_files.F03",
    "src/mqc_integrals.F03",
    "src/mqc_algebra.F03",
    "src/mqc_algebra2.F03",
    "src/mqc_molecule.F03",
    "src/mqc_est.F03",
    "src/mqc_matwrapper_direct.F03",
    "src/mqc_matwrapper_stub.F03",
    "src/mqc_gaussian.F03",
    "src/mqc_FullWavefunction.F03",
    "src/mqc_interface.F03",
)


@dataclass
class Entry:
    name: str
    line: int
    detail: str = ""


@dataclass
class SourceInventory:
    path: str
    modules: list[Entry] = field(default_factory=list)
    types: list[Entry] = field(default_factory=list)
    bindings: list[Entry] = field(default_factory=list)
    interfaces: list[Entry] = field(default_factory=list)
    procedures: list[Entry] = field(default_factory=list)


def logical_statements(path: Path) -> list[tuple[int, str]]:
    """Join ordinary free-form continuation lines and discard comments."""
    statements: list[tuple[int, str]] = []
    buffer = ""
    start_line = 0
    for line_number, physical_line in enumerate(
        path.read_text(encoding="utf-8", errors="replace").splitlines(), 1
    ):
        code = physical_line.split("!", 1)[0].rstrip()
        if not code.strip():
            continue
        if not buffer:
            start_line = line_number
        piece = code.strip()
        if buffer and piece.startswith("&"):
            piece = piece[1:].lstrip()
        continued = piece.endswith("&")
        if continued:
            piece = piece[:-1].rstrip()
        buffer = f"{buffer} {piece}".strip()
        if not continued:
            statements.append((start_line, buffer))
            buffer = ""
    if buffer:
        statements.append((start_line, buffer))
    return statements


def parse_fortran(repo_root: Path, relative_path: str) -> SourceInventory:
    inventory = SourceInventory(relative_path)
    type_name: str | None = None
    interface_depth = 0

    module_re = re.compile(r"^\s*module\s+(?!procedure\b)([a-z]\w*)", re.I)
    type_re = re.compile(
        r"^\s*type\s*(?!\()(?!(?:is|default)\b)"
        r"(?:,[^:]*)?(?:::|\s+)\s*([a-z]\w*)\s*$",
        re.I,
    )
    interface_re = re.compile(r"^\s*interface(?:\s+(.+?))?\s*$", re.I)
    procedure_re = re.compile(
        r"^\s*(?:(?:pure|elemental|recursive|impure|module)\s+)*"
        r"(?:(?:integer|real|complex|logical|character|type\s*\([^)]*\)|"
        r"class\s*\([^)]*\)|double\s+precision)"
        r"(?:\s*\([^)]*\))?(?:\s*,[^:]*)?\s+)?"
        r"(subroutine|function)\s+([a-z]\w*)\b",
        re.I,
    )
    binding_re = re.compile(
        r"^\s*procedure(?:\s*\([^)]*\))?(?:\s*,[^:]*)?::\s*(.+)$", re.I
    )
    generic_re = re.compile(r"^\s*generic(?:\s*,[^:]*)?::\s*(.+)$", re.I)

    statements = logical_statements(repo_root / relative_path)
    has_module = any(module_re.match(statement) for _, statement in statements)
    module_contains = not has_module

    for line_number, statement in statements:
        lower = statement.lower().strip()

        module_match = module_re.match(statement)
        if module_match:
            inventory.modules.append(Entry(module_match.group(1), line_number))
            continue

        if re.match(r"^\s*end\s*type\b", statement, re.I):
            type_name = None
            continue

        if type_name is None:
            type_match = type_re.match(statement)
            if type_match:
                type_name = type_match.group(1)
                inventory.types.append(Entry(type_name, line_number))
                continue
        else:
            binding_match = binding_re.match(statement)
            generic_match = generic_re.match(statement)
            if generic_match:
                declaration = generic_match.group(1).strip()
                if "=>" in declaration:
                    binding, target = (
                        part.strip() for part in declaration.split("=>", 1)
                    )
                else:
                    binding, target = declaration, declaration
                inventory.bindings.append(
                    Entry(f"{type_name}%{binding}", line_number, target)
                )
            elif binding_match:
                declaration = binding_match.group(1)
                for item in declaration.split(","):
                    item = item.strip()
                    if not item:
                        continue
                    if "=>" in item:
                        binding, target = (part.strip() for part in item.split("=>", 1))
                    else:
                        binding, target = item, item
                    inventory.bindings.append(
                        Entry(f"{type_name}%{binding}", line_number, target)
                    )
            continue

        if re.match(r"^\s*end\s*interface\b", statement, re.I):
            interface_depth = max(0, interface_depth - 1)
            continue

        interface_match = interface_re.match(statement)
        if interface_match:
            interface_depth += 1
            if interface_depth == 1 and interface_match.group(1):
                inventory.interfaces.append(
                    Entry(interface_match.group(1).strip(), line_number)
                )
            continue

        if lower == "contains" and interface_depth == 0:
            module_contains = True
            continue

        if module_contains and interface_depth == 0 and not lower.startswith("end"):
            procedure_match = procedure_re.match(statement)
            if procedure_match:
                inventory.procedures.append(
                    Entry(procedure_match.group(2), line_number, procedure_match.group(1))
                )

    return inventory


def parse_c(repo_root: Path) -> SourceInventory:
    relative_path = "src/mqc_util.c"
    inventory = SourceInventory(relative_path)
    text = (repo_root / relative_path).read_text(encoding="utf-8", errors="replace")
    definition_re = re.compile(
        r"(?m)^[ \t]*(?!if\b|for\b|while\b|switch\b)"
        r"(?:static\s+)?(?:void|int|long|off_t|size_t|char\s*\*|FILE\s*\*)"
        r"\s+([a-zA-Z_]\w*)\s*\([^;{}]*?\)\s*\{"
    )
    for match in definition_re.finditer(text):
        line_number = text.count("\n", 0, match.start()) + 1
        inventory.procedures.append(Entry(match.group(1), line_number, "C function"))
    return inventory


def source_link(relative_path: str, line: int) -> str:
    return f"../../{relative_path}#L{line}"


def render_entries(entries: list[Entry], relative_path: str, include_detail: bool) -> list[str]:
    rows: list[str] = []
    seen: set[tuple[str, int, str]] = set()
    for entry in entries:
        key = (entry.name.lower(), entry.line, entry.detail.lower())
        if key in seen:
            continue
        seen.add(key)
        detail = f" — `{entry.detail}`" if include_detail and entry.detail else ""
        rows.append(
            f"- `{entry.name}`{detail} "
            f"([source]({source_link(relative_path, entry.line)}))"
        )
    return rows


def render(inventories: list[SourceInventory]) -> str:
    procedure_count = sum(len(item.procedures) for item in inventories)
    binding_count = sum(len(item.bindings) for item in inventories)
    lines = [
        "# MQCPack source procedure index",
        "",
        "> Generated by `python3 doc/api/generate_procedure_index.py`. Do not edit by hand.",
        "",
        "This is a source inventory, not a stability guarantee. Prefer the interfaces",
        "described in [AGENT_API_GUIDE.md](../AGENT_API_GUIDE.md) and the curated API",
        "pages. Many MQCPack modules do not declare a module-wide `private` default, so",
        "a technically accessible implementation routine is not automatically a supported",
        "downstream interface.",
        "",
        f"The current inventory contains {procedure_count} Fortran/C procedure definitions",
        f"and {binding_count} type-bound procedure or generic bindings.",
        "",
        "The scanner recognizes free-form Fortran declarations and ordinary C function",
        "definitions. Conditional compilation can change which entries are built, and",
        "generic overload membership must still be checked in the source interface block.",
        "",
    ]

    for inventory in inventories:
        lines.extend((f"## `{inventory.path}`", ""))
        sections = (
            ("Modules", inventory.modules, False),
            ("Derived types", inventory.types, False),
            ("Type-bound bindings", inventory.bindings, True),
            ("Generic and defined interfaces", inventory.interfaces, False),
            ("Procedure definitions", inventory.procedures, True),
        )
        for title, entries, include_detail in sections:
            if not entries:
                continue
            lines.extend((f"### {title}", ""))
            lines.extend(render_entries(entries, inventory.path, include_detail))
            lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def main() -> None:
    default_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parent / "PROCEDURE_INDEX.md",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail if the existing output differs instead of rewriting it",
    )
    args = parser.parse_args()

    repo_root = args.root.resolve()
    inventories = [parse_fortran(repo_root, path) for path in FORTRAN_SOURCES]
    inventories.append(parse_c(repo_root))
    generated = render(inventories)

    if args.check:
        current = args.output.read_text(encoding="utf-8") if args.output.exists() else ""
        if current != generated:
            raise SystemExit(
                f"{args.output} is stale; run doc/api/generate_procedure_index.py"
            )
        return

    args.output.write_text(generated, encoding="utf-8")


if __name__ == "__main__":
    main()
