#!/usr/bin/env python3
"""Generate a Markdown table documenting samplesheet columns from an
nf-core style assets/schema_input.json file."""
import json
import sys

def main(schema_path, out_path=None):
    with open(schema_path) as f:
        schema = json.load(f)

    props = schema.get("items", {}).get("properties", {})
    required = set(schema.get("items", {}).get("required", []))

    lines = []
    lines.append(f"# {schema.get('title', 'Samplesheet input')}\n")
    if schema.get("description"):
        lines.append(f"{schema['description']}\n")
    lines.append("| Column | Description | Type | Required | Pattern |")
    lines.append("|---|---|---|---|---|")

    for col, spec in props.items():
        desc = spec.get("description") or spec.get("errorMessage", "")
        ftype = spec.get("format", spec.get("type", ""))
        is_req = "Yes" if col in required else "No"
        pattern = f"`{spec['pattern']}`" if "pattern" in spec else ""
        desc = desc.replace("|", "\\|")
        lines.append(f"| `{col}` | {desc} | {ftype} | {is_req} | {pattern} |")

    md = "\n".join(lines) + "\n"

    if out_path:
        with open(out_path, "w") as f:
            f.write(md)
        print(f"Written to {out_path}")
    else:
        print(md)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: schema_input_to_md.py <path/to/schema_input.json> [output.md]")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
