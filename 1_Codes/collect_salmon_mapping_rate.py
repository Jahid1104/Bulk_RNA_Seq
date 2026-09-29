"""
Collect percent_mapped from Salmon meta_info.json files.

Expects a folder structure like:

    <base_dir>/Salmon_PE/<sample_name>/aux_info/meta_info.json
    <base_dir>/Salmon_SE/<sample_name>/aux_info/meta_info.json

For each sample folder found, this pulls "percent_mapped" out of
meta_info.json and writes everything to one CSV.
"""

import json
from pathlib import Path
import csv

# ---- EDIT THIS ----
BASE_DIR = Path(r"R:\Md_Jahid_Hasan_Jone\1_Experiments_and_Data\5_RNA_seq\New_Name\Results\2_Salmon")
OUTPUT_CSV = Path(r"R:\Md_Jahid_Hasan_Jone\1_Experiments_and_Data\5_RNA_seq\New_Name\Results\2_Salmon\salmon_mapping_rates.csv")
# -------------------

SALMON_FOLDERS = ["Salmon_PE", "Salmon_SE"]


def collect_mapping_rates(base_dir: Path) -> list[dict]:
    rows = []
    for folder_name in SALMON_FOLDERS:
        folder = base_dir / folder_name
        if not folder.is_dir():
            print(f"Skipping {folder} (not found)")
            continue

        for sample_dir in sorted(folder.iterdir()):
            if not sample_dir.is_dir():
                continue

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

            rows.append({
                "sample_name": sample_dir.name,
                "library_type": folder_name,
                "percent_mapped": round(percent_mapped, 2),
            })

    return rows


def main():
    rows = collect_mapping_rates(BASE_DIR)

    if not rows:
        print("No samples found. Check BASE_DIR.")
        return

    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample_name", "library_type", "percent_mapped"])
        writer.writeheader()
        writer.writerows(rows)

    print(f"\nWrote {len(rows)} samples to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
