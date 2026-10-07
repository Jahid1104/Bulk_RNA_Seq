##Bulk RNA Seq Data Analysis Workflow by Md Jahid Hasan Jone##

# Collect the percent of reads mapped from the Salmon output (run after 2_Salmon.sh).
#
# Expects the folder structure written by 2_Salmon.sh:
#
#     5_Results/2_Salmon/<sample_name>/aux_info/meta_info.json
#
# For each sample folder, this pulls "percent_mapped" (and the library type
# Salmon detected) out of meta_info.json and writes everything to one CSV:
#
#     5_Results/2_Salmon/salmon_mapping_rates.csv
#
# Needs only the Python standard library.

import csv
import json
import re
import sys
from pathlib import Path


# =====================================================================
# PATHS  <- EDIT THESE (replace /.../.../ with your path)
# =====================================================================
salmon_dir = Path("/.../.../Bulk_RNA_seq/5_Results/2_Salmon")
output_csv = Path("/.../.../Bulk_RNA_seq/5_Results/2_Salmon/salmon_mapping_rates.csv")


# Natural sort so mixed sample numbering (1F, 2F, ..., 10F, 1L, 2L, ...) orders
# sensibly instead of plain alphabetical order (which would put "10F" before "2F").
def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def collect_mapping_rates(base_dir: Path) -> list:
    rows = []
    sample_dirs = sorted(
        (d for d in base_dir.iterdir() if d.is_dir()),
        key=lambda d: natural_key(d.name),
    )

    for sample_dir in sample_dirs:
        meta_path = sample_dir / "aux_info" / "meta_info.json"
        if not meta_path.exists():
            print(f"  Missing meta_info.json for {sample_dir.name}, skipping")
            continue

        with open(meta_path) as f:
            meta = json.load(f)

        percent_mapped = meta.get("percent_mapped")
        if percent_mapped is None:
            print(f"  No percent_mapped in {meta_path}, skipping")
            continue

        # 2_Salmon.sh uses --libType A, so this is the type Salmon detected (e.g. IU = paired end)
        library_type = ",".join(meta.get("library_types", []))

        rows.append({
            "sample_name": sample_dir.name,
            "library_type": library_type,
            "percent_mapped": round(percent_mapped, 2),
        })

    return rows


def main():
    if "..." in str(salmon_dir):
        sys.exit("Edit the paths at the top of this script (replace /.../.../ with your path).")
    if not salmon_dir.is_dir():
        sys.exit(f"Salmon folder not found: {salmon_dir}")

    rows = collect_mapping_rates(salmon_dir)

    if not rows:
        print("No samples found. Check salmon_dir.")
        return

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample_name", "library_type", "percent_mapped"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} samples to {output_csv}")


if __name__ == "__main__":
    main()
