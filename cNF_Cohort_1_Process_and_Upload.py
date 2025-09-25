#!/usr/bin/env python3
"""
cNF Proteomics Pipeline (Cohort 1)

Runs the following 3 steps:
  1) Download All data and Metadata from Synapse
  2) Propose annotations & perform per-sample splits. Get sex/age.
  3) Package under "Proteomic Data (Cohort 1)" and upload to Synapse with annotations
"""

import os
import re
import math
import shutil
import hashlib
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pandas as pd
import synapseclient
from synapseclient import Folder,File

# =============================================================================
# CONFIG — edit here if needed (no CLI)
# =============================================================================

# Verbosity
VERBOSE = False  # set True to see progress logs

# --- Step 1: Inventory / downloads
PARENT_SYN_ID = "syn51301417"       # Source Synapse folder for raw/normalized files
DOWNLOAD_FILES = True                # If False, inventory only (no file contents)
DOWNLOAD_DIR = Path("step1_downloads")
INVENTORY_CSV = Path("step1_inventory_existing_files.csv")

# --- Step 2: inputs
PROTEOMICS_FILES = ["step1_downloads/DIA_Global_DiaNN_Results.tsv", "step1_downloads/DIA_Global_FP_Results.tsv"]
PHOSPHO_FILES = ["step1_downloads/DIA_Phospho_FP_Results.tsv"]
PHOSPHO_SITEID_FILE = "step1_downloads/DIA_Phospho_FP_Results_SiteID.txt"

# Metadata
META_FILES = [
    {"path": "metadata_cNF_Exp1.xlsx", "sheet": 0},
    {"path": "metadata_cNF_Exp2.xlsx", "sheet": 1},
]

# Patient info
PATIENT_INFO_TSV = "patient_info.tsv"
ADD_AGE_SEX_TO_FILES = False   # If True, write Age/Sex into per-sample TSVs (not needed for these)

# Outputs for step 2
OUTDIR_ANNOT = Path("step2_outputs")
OUTDIR_SPLIT = Path("step2_split")

# --- Step 3: packaging/upload
PACKAGE_DIR = Path("step3_upload_package")
DATA_ROOT_NAME = "Proteomics Data by Patient"
UPLOAD_PARENT_SYN_ID = "syn51301417" # destination in Synapse
ENABLE_SYNAPSE_UPLOAD = True        # If False, do a dry-run (no uploads, use for debugging)
SYNAPSE_PAT_ENVVAR = "SYNAPSE_AUTH_TOKEN"  # Env var name for Synapse PAT

# Study-level annotations (includes your studyId change)
ANNOT_DEFAULTS = {
    "species": "Homo sapiens",
    "studyId": "syn51301409",
    "studyName": "Leveraging patient-derived cutaneous neurofibroma organoid models to identify biomarkers of drug response",
    "tumorType": "cutaneous neurofibroma (cNF)",
    "initiative": "Biology and Therapeutic Development for Cutaneous Neurofibromas",
    "fundingAgency": "NTAP",
    "disease": "Neurofibromatosis type 1",
    "diagnosis": "Neurofibromatosis type 1",
}

DATATYPE_BY_KIND = {
    "proteomics": "DIA Global Proteomics",
    "phosphoproteomics": "DIA Phosphoproteomics",
    "phospho-siteid": "DIA Phospho SiteID",
}

# Data download paths
SYN_FIXED_DOWNLOAD_DIR = Path("syn_inputs")
SYN_ID_METADATA_EXP1   = "syn69920463"
SYN_ID_METADATA_EXP2   = "syn69920464"
SYN_ID_PATIENT_INFO    = "syn69931351"

# =============================================================================
# Logging
# =============================================================================

def log(msg: str):
    if VERBOSE:
        print(msg)

def warn(msg: str):
    print(msg)

# =============================================================================
# Utilities
# =============================================================================

def syn_login():
    syn = synapseclient.Synapse()
    token = os.environ.get(SYNAPSE_PAT_ENVVAR)
    if token:
        syn.login(authToken=token, silent=True)
    else:
        syn.login(silent=True)
    return syn

def md5_file(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def safe_mkdir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def _is_nan(x):
    try:
        return x is None or (isinstance(x, float) and math.isnan(x))
    except Exception:
        return False

# =============================================================================
# Step 1: Inventory / download
# =============================================================================

def fetch_all_files(syn, syn_id: str, download: bool, outdir: Optional[Path]) -> List[dict]:
    records = []
    for item in syn.getChildren(syn_id):
        item_type = item.get('type')
        if item_type in ('org.sagebionetworks.repo.model.FileEntity',
                         'org.sagebionetworks.repo.model.FileHandleAssociate'):
            file_syn_id = item['id']
            entity = syn.get(file_syn_id, downloadFile=download, downloadLocation=str(outdir) if outdir else None)
            ann = syn.get_annotations(file_syn_id)
            rec = {
                "synId": file_syn_id,
                "name": item.get("name"),
                "parentId": item.get("parentId"),
                "versionNumber": getattr(entity, 'versionNumber', None),
                "fileSizeBytes": getattr(entity, 'fileSize', None),
                "md5": getattr(entity, 'md5', None),
            }
            for key, value in (ann.items() if hasattr(ann, "items") else []):
                if isinstance(value, (list, tuple)):
                    rec[f"ann.{key}"] = ";".join(str(v) for v in value)
                else:
                    rec[f"ann.{key}"] = str(value)
            records.append(rec)
        elif item_type == 'org.sagebionetworks.repo.model.Folder':
            log(f"[INFO] Descending into folder {item.get('name')} ({item['id']})")
            records.extend(fetch_all_files(syn, item['id'], download, outdir))
    return records

def run_step1_inventory():
    syn = syn_login()
    if DOWNLOAD_FILES:
        safe_mkdir(DOWNLOAD_DIR)
    log(f"[STEP1] Starting fetch under {PARENT_SYN_ID} (download={DOWNLOAD_FILES})")
    recs = fetch_all_files(syn, PARENT_SYN_ID, DOWNLOAD_FILES, DOWNLOAD_DIR if DOWNLOAD_FILES else None)
    log(f"[STEP1] Collected {len(recs)} file records")
    pd.DataFrame(recs).to_csv(INVENTORY_CSV, index=False)
    log(f"[STEP1] Inventory saved -> {INVENTORY_CSV}")

def download_fixed_inputs_from_synapse():
    """
    Download metadata_cNF_Exp1.xlsx, metadata_cNF_Exp2.xlsx, and patient_info.tsv
    by Synapse ID into SYN_FIXED_DOWNLOAD_DIR, and update the pipeline paths so
    Step 2 uses these downloaded copies.
    """
    id_map = {
        "metadata_exp1": SYN_ID_METADATA_EXP1,
        "metadata_exp2": SYN_ID_METADATA_EXP2,
        "patient_info":  SYN_ID_PATIENT_INFO,
    }
    if not any(id_map.values()):
        return

    safe_mkdir(SYN_FIXED_DOWNLOAD_DIR)
    syn = syn_login()

    def _fetch(syn_id: str, suggested_name: str) -> Optional[Path]:
        if not syn_id:
            return None
        ent = syn.get(syn_id, downloadFile=True, downloadLocation=str(SYN_FIXED_DOWNLOAD_DIR))
        local = Path(getattr(ent, "path", "") or (SYN_FIXED_DOWNLOAD_DIR / suggested_name))
        return local

    local_meta1 = _fetch(SYN_ID_METADATA_EXP1, "metadata_cNF_Exp1.xlsx")
    local_meta2 = _fetch(SYN_ID_METADATA_EXP2, "metadata_cNF_Exp2.xlsx")
    local_pinfo = _fetch(SYN_ID_PATIENT_INFO, "patient_info.tsv")

    # Update META_FILES and PATIENT_INFO_TSV in-place
    global META_FILES, PATIENT_INFO_TSV
    if local_meta1 and local_meta1.exists():
        # update matching entry if present, else append
        found = False
        for entry in META_FILES:
            if Path(entry["path"]).name.lower().startswith("metadata_cnf_exp1"):
                entry["path"] = str(local_meta1); found = True; break
        if not found:
            META_FILES.append({"path": str(local_meta1), "sheet": 0})

    if local_meta2 and local_meta2.exists():
        found = False
        for entry in META_FILES:
            if Path(entry["path"]).name.lower().startswith("metadata_cnf_exp2"):
                entry["path"] = str(local_meta2); found = True; break
        if not found:
            META_FILES.append({"path": str(local_meta2), "sheet": 1})

    if local_pinfo and local_pinfo.exists():
        PATIENT_INFO_TSV = str(local_pinfo)

# =============================================================================
# Helper for part 2
# =============================================================================

def canonicalize_dataset_name(token: str) -> str:
    token = token.replace(".raw", "")
    token = re.sub(r"(BEHCoA)[\.-_](\d{2})[\.-_](\d{2})[\.-_](\d{2})", r"BEHCoA-\2-\3-\4", token, flags=re.I)
    return token.strip("._- ")

def extract_dataset_name_loose(colname: str):
    if not isinstance(colname, str):
        return None
    m = re.search(r"(cNF[^\\/\t\n\r]*?)(?:\.raw\b|$)", colname, flags=re.IGNORECASE)
    if not m:
        return None
    token = m.group(1)
    m2 = re.search(
        r"(cNF[^\\/\t\n\r]*?DIA_[GP]_\d+_[^_]*?_BEHCoA[.\-_]\d{2}[.\-_]\d{2}[.\-_]\d{2})",
        token,
        flags=re.IGNORECASE,
    )
    token = m2.group(1) if m2 else token
    return canonicalize_dataset_name(token)

def load_metadata(meta_specs: List[dict]) -> Tuple[pd.DataFrame, Dict[str, int], Dict[str, int]]:
    meta_list = []
    for mf in meta_specs:
        path = mf["path"]; sheet = mf.get("sheet", 0)
        _m = pd.read_excel(path, sheet_name=sheet)
        _m.columns = _m.columns.str.strip()
        for col in ["Description", "Specimen", "Patient", "SampleAlias", "Date"]:
            if col not in _m.columns:
                _m[col] = pd.NA
        if "DatasetNameGlobal" not in _m.columns or "DatasetNamePhospho" not in _m.columns:
            raise ValueError(f"{path} (sheet {sheet}) missing DatasetNameGlobal/Phospho columns")
        _m["DatasetNameGlobal"] = _m["DatasetNameGlobal"].astype(str).str.strip().map(canonicalize_dataset_name)
        _m["DatasetNamePhospho"] = _m["DatasetNamePhospho"].astype(str).str.strip().map(canonicalize_dataset_name)
        _m["Specimen"] = _m["Specimen"].astype(str).str.strip()
        _m["Patient"] = _m["Patient"].astype(str).str.strip()
        _m["SampleAlias"] = _m["SampleAlias"].astype(str)
        _m["SampleAlias"] = _m["SampleAlias"].str.strip()
        meta_list.append(_m)

    meta = pd.concat(meta_list, ignore_index=True).drop_duplicates()
    global_map = dict(zip(meta["DatasetNameGlobal"], meta.index))
    phospho_map = dict(zip(meta["DatasetNamePhospho"], meta.index))
    return meta, global_map, phospho_map

def derive_sample_and_tumor(meta_row: pd.Series):
    sample = (meta_row.get("Patient") or "").strip() or None
    specimen = (meta_row.get("Specimen") or "").strip()
    alias = (meta_row.get("SampleAlias") or "").strip()
    tumor = None
    m = re.search(r"(?:^|[_\-\s])([Tt]\d+)\b", specimen)
    if m:
        tumor = m.group(1).upper()
    elif re.fullmatch(r"[Tt]\d+", alias or ""):
        tumor = alias.upper()
    return sample, tumor

def find_sample_columns(df: pd.DataFrame, treat_first_as_site: bool = False):
    cols = list(df.columns)
    if treat_first_as_site:
        site_col = cols[0]
        sample_cols = [c for c in cols[1:] if extract_dataset_name_loose(c)]
        if not sample_cols:
            sample_cols = cols[1:]
        return site_col, sample_cols
    sample_cols = [c for c in cols if extract_dataset_name_loose(c)]
    non_sample  = [c for c in cols if c not in sample_cols]
    site_col = non_sample[0] if len(non_sample) == 1 else None
    return site_col, sample_cols

# =============================================================================
# Step 2
# =============================================================================

def annotate_file(file_paths: List[str], lookup_map: Dict[str, int], meta: pd.DataFrame,
                  kind: str, siteid: bool = False) -> pd.DataFrame:
    rows = []
    for path in file_paths:
        if not os.path.exists(path):
            warn(f"[{kind}] Missing file: {path}")
            continue
        dfh = pd.read_csv(path, sep="\t", nrows=1)
        _, sample_cols = find_sample_columns(dfh, treat_first_as_site=siteid)
        log(f"[{kind}] {os.path.basename(path)}: matched {len(sample_cols)} columns")
        for col in sample_cols:
            ds = extract_dataset_name_loose(col)
            if not ds or ds not in lookup_map:
                continue
            row = meta.loc[lookup_map[ds]].copy()
            sample, tumor = derive_sample_and_tumor(row)
            info = row.to_dict()
            info.update({
                "dataset_name": ds,
                "column_header": col,
                "file": os.path.basename(path),
                "data_type": kind,
                "sample": sample,
                "tumor": tumor,
            })
            rows.append(info)
    cols = list(meta.columns) + ["dataset_name", "column_header", "file", "data_type", "sample", "tumor"]
    return pd.DataFrame(rows, columns=cols).drop_duplicates()

def write_per_sample_matrix(path, kind, outdir, id_priority, meta, lookup_map, manifest_rows):
    if not os.path.exists(path):
        warn(f"[{kind}] Missing file: {path}")
        return
    df = pd.read_csv(path, sep="\t")
    site_col, sample_cols = find_sample_columns(df, treat_first_as_site=False)
    if not sample_cols:
        warn(f"[{kind}] No sample columns in {os.path.basename(path)} — skipping.")
        return
    id_col = None
    for cand in (id_priority or []):
        if cand in df.columns:
            id_col = cand; break
    if id_col is None:
        id_col = site_col if (site_col and site_col in df.columns) else df.columns[0]

    bname = os.path.splitext(os.path.basename(path))[0]
    odir = outdir / kind / bname
    safe_mkdir(odir)

    for col in sample_cols:
        ds = extract_dataset_name_loose(col) or "unknown_dataset"
        out = odir / f"{ds}.{kind}.tsv"
        if ds in lookup_map:
            meta_row = meta.loc[lookup_map[ds]]
        else:
            meta_row = pd.Series({})
        sample, tumor = derive_sample_and_tumor(meta_row) if not meta_row.empty else (None, None)

        sub = df[[id_col, col]].rename(columns={id_col: "feature_id", col: "value"})
        sub["sample"] = sample
        sub["tumor"]  = tumor
        sub.to_csv(out, sep="\t", index=False)

        if ds in lookup_map:
            manifest_rows.append({
                "kind": kind,
                "file_path": str(out),
                "dataset_name": ds,
                "sample": sample,
                "tumor": tumor,
                "Specimen": meta_row.get("Specimen"),
                "Patient": meta_row.get("Patient"),
                "SampleAlias": meta_row.get("SampleAlias"),
                "Date": meta_row.get("Date"),
                "Description": meta_row.get("Description"),
            })
    log(f"[{kind}] Wrote {len(sample_cols)} files -> {odir}")

def write_per_sample_siteid(path, outdir, meta, lookup_map, manifest_rows):
    if not os.path.exists(path):
        warn(f"[phospho-siteid] Missing file: {path}")
        return
    df = pd.read_csv(path, sep="\t")
    site_col, sample_cols = find_sample_columns(df, treat_first_as_site=True)
    if not sample_cols:
        warn(f"[phospho-siteid] No sample columns in {os.path.basename(path)} — skipping.")
        return
    bname = os.path.splitext(os.path.basename(path))[0]
    odir = outdir / "phospho-siteid" / bname
    safe_mkdir(odir)
    for col in sample_cols:
        ds = extract_dataset_name_loose(col) or "unknown_dataset"
        out = odir / f"{ds}.phospho_siteid.tsv"
        if ds in lookup_map:
            meta_row = meta.loc[lookup_map[ds]]
        else:
            meta_row = pd.Series({})
        sample, tumor = derive_sample_and_tumor(meta_row) if not meta_row.empty else (None, None)
        sub = df[[site_col, col]].rename(columns={site_col: "site_id", col: "intensity"})
        sub["sample"] = sample
        sub["tumor"]  = tumor
        sub.to_csv(out, sep="\t", index=False)
        if ds in lookup_map:
            manifest_rows.append({
                "kind": "phospho-siteid",
                "file_path": str(out),
                "dataset_name": ds,
                "sample": sample,
                "tumor": tumor,
                "Specimen": meta_row.get("Specimen"),
                "Patient": meta_row.get("Patient"),
                "SampleAlias": meta_row.get("SampleAlias"),
                "Date": meta_row.get("Date"),
                "Description": meta_row.get("Description"),
            })
    log(f"[phospho-siteid] Wrote {len(sample_cols)} files -> {odir}")

def load_patient_info(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df.columns = df.columns.str.strip()
    if "Patient" not in df.columns:
        if "Patient ID" in df.columns:
            df = df.rename(columns={"Patient ID": "Patient"})
        else:
            raise ValueError("patient_info.tsv must contain 'Patient' or 'Patient ID'.")
    df["Patient"] = df["Patient"].astype(str).str.strip()
    keep = ["Patient"]
    if "Age" in df.columns: keep.append("Age")
    if "Sex" in df.columns: keep.append("Sex")
    return df[keep].drop_duplicates()

def merge_patient_info(in_csv: Path, out_csv: Path, pinfo: pd.DataFrame) -> int:
    if not in_csv.exists():
        return 0
    df = pd.read_csv(in_csv)
    df.columns = df.columns.str.strip()
    if "Patient" not in df.columns:
        alt = [c for c in df.columns if c.strip() == "Patient"]
        if alt: df = df.rename(columns={alt[0]: "Patient"})
        else:
            df.to_csv(out_csv, index=False)
            return len(df)
    df["Patient"] = df["Patient"].astype(str).str.strip()
    merged = df.merge(pinfo, on="Patient", how="left")

    cols = list(merged.columns)
    if "Age" in merged.columns or "Sex" in merged.columns:
        new_order, seen = [], set()
        for c in cols:
            new_order.append(c)
            if c == "Patient":
                if "Age" in merged.columns: new_order.append("Age")
                if "Sex" in merged.columns: new_order.append("Sex")
        ordered = []
        for c in new_order:
            if c in merged.columns and c not in seen:
                seen.add(c); ordered.append(c)
        merged = merged[ordered]

    safe_mkdir(out_csv.parent)
    merged.to_csv(out_csv, index=False)
    return len(merged)

def stamp_age_sex_into_files(manifest_csv: Path, pinfo: pd.DataFrame) -> int:
    if not manifest_csv.exists():
        return 0
    man = pd.read_csv(manifest_csv)
    man.columns = man.columns.str.strip()
    if "Patient" not in man.columns:
        alt = [c for c in man.columns if c.strip() == "Patient"]
        if alt: man = man.rename(columns={alt[0]: "Patient"})
        else:
            return 0
    man["Patient"] = man["Patient"].astype(str).str.strip()
    merged = man.merge(pinfo, on="Patient", how="left")
    updated = 0
    for _, row in merged.iterrows():
        fpath = row.get("file_path", None)
        if not isinstance(fpath, str) or not os.path.exists(fpath):
            continue
        try:
            df = pd.read_csv(fpath, sep="\t")
            if "Age" in merged.columns: df["Age"] = row.get("Age")
            if "Sex" in merged.columns: df["Sex"] = row.get("Sex")
            df.to_csv(fpath, sep="\t", index=False)
            updated += 1
        except Exception:
            pass
    return updated

def run_step2_all():
    """Run Step 2: Split by sample, get metadata, create annotations."""
    safe_mkdir(OUTDIR_ANNOT); safe_mkdir(OUTDIR_SPLIT)

    meta, global_map, phospho_map = load_metadata(META_FILES)

    # Proposed annotations
    df_prot = annotate_file(PROTEOMICS_FILES, global_map,  meta, kind="proteomics",        siteid=False)
    df_phos = annotate_file(PHOSPHO_FILES,    phospho_map, meta, kind="phosphoproteomics", siteid=False)
    df_site = annotate_file([PHOSPHO_SITEID_FILE], phospho_map, meta, kind="phospho-siteid", siteid=True)

    df_prot.to_csv(OUTDIR_ANNOT / "step2_proteomics_proposed_annotations.csv", index=False)
    df_phos.to_csv(OUTDIR_ANNOT / "step2_phospho_proposed_annotations.csv", index=False)
    df_site.to_csv(OUTDIR_ANNOT / "step2_phospho_siteid_proposed_annotations.csv", index=False)

    # Missing report
    missing_prot = set(global_map)  - set(df_prot["dataset_name"]) if not df_prot.empty else set(global_map)
    missing_phos = set(phospho_map) - set(df_phos["dataset_name"]) if not df_phos.empty else set(phospho_map)
    missing_site = set(phospho_map) - set(df_site["dataset_name"]) if not df_site.empty else set(phospho_map)
    with open(OUTDIR_ANNOT / "step2_missing_fields_summary.csv", "w") as f:
        f.write("missing_type,dataset_name\n")
        for ds in sorted(missing_prot): f.write(f"proteomics,{ds}\n")
        for ds in sorted(missing_phos): f.write(f"phosphoproteomics,{ds}\n")
        for ds in sorted(missing_site): f.write(f"phospho-siteid,{ds}\n")

    # Split per-sample + manifest
    manifest_rows = []
    for p in PROTEOMICS_FILES:
        write_per_sample_matrix(
            p, kind="proteomics", outdir=OUTDIR_SPLIT,
            id_priority=["Genes", "Protein.Group"], meta=meta, lookup_map=global_map,
            manifest_rows=manifest_rows
        )
    for p in PHOSPHO_FILES:
        write_per_sample_matrix(
            p, kind="phosphoproteomics", outdir=OUTDIR_SPLIT,
            id_priority=["Modified.Sequence", "Stripped.Sequence", "Precursor.Id"],
            meta=meta, lookup_map=phospho_map,
            manifest_rows=manifest_rows
        )
    if os.path.exists(PHOSPHO_SITEID_FILE):
        write_per_sample_siteid(PHOSPHO_SITEID_FILE, OUTDIR_SPLIT, meta, phospho_map, manifest_rows=manifest_rows)

    manifest_path = OUTDIR_ANNOT / "step2_per_sample_manifest.csv"
    pd.DataFrame(manifest_rows).to_csv(manifest_path, index=False)

    # Merge Age/Sex into CSVs (+ optional stamping)
    pinfo = load_patient_info(PATIENT_INFO_TSV)

    man_in  = OUTDIR_ANNOT / "step2_per_sample_manifest.csv"
    man_out = OUTDIR_ANNOT / "step2_per_sample_manifest_with_patient.csv"
    prot_in = OUTDIR_ANNOT / "step2_proteomics_proposed_annotations.csv"
    prot_out= OUTDIR_ANNOT / "step2_proteomics_proposed_annotations_with_patient.csv"
    phos_in = OUTDIR_ANNOT / "step2_phospho_proposed_annotations.csv"
    phos_out= OUTDIR_ANNOT / "step2_phospho_proposed_annotations_with_patient.csv"
    site_in = OUTDIR_ANNOT / "step2_phospho_siteid_proposed_annotations.csv"
    site_out= OUTDIR_ANNOT / "step2_phospho_siteid_proposed_annotations_with_patient.csv"

    merge_patient_info(man_in,  man_out,  pinfo)
    merge_patient_info(prot_in, prot_out, pinfo)
    merge_patient_info(phos_in, phos_out, pinfo)
    if site_in.exists():
        merge_patient_info(site_in, site_out, pinfo)

    if ADD_AGE_SEX_TO_FILES:
        stamp_age_sex_into_files(man_in, pinfo)

# =============================================================================
# Step 3: Package & Upload
# =============================================================================

def read_csv_if_exists(p: Path, usecols=None) -> pd.DataFrame:
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p)
    if usecols:
        keep = [c for c in usecols if c in df.columns]
        if keep: df = df[keep]
    return df

def read_per_sample_manifest_base() -> pd.DataFrame:
    a = OUTDIR_ANNOT / "step2_per_sample_manifest_with_patient.csv"
    b = OUTDIR_ANNOT / "step2_per_sample_manifest.csv"
    if a.exists(): return pd.read_csv(a)
    if b.exists(): return pd.read_csv(b)
    return pd.DataFrame()

def load_age_sex_maps() -> Dict[str, dict]:
    """
    Prefer *_with_patient.csv; then base CSVs; then fallback via patient_info.tsv + manifest.
    Ensures we can annotate files even when TSVs contain no Age/Sex columns.
    """
    def _norm_age_sex(df: pd.DataFrame) -> pd.DataFrame:
        age_cols = [c for c in df.columns if re.fullmatch(r"Age(_[xy])?", c, flags=re.I)]
        sex_cols = [c for c in df.columns if re.fullmatch(r"Sex(_[xy])?", c, flags=re.I)]
        def first_nonnull(row, cols):
            for c in cols:
                if c in row and pd.notna(row[c]): return row[c]
            return pd.NA
        if age_cols: df["Age"] = df.apply(lambda r: first_nonnull(r, age_cols), axis=1)
        if sex_cols: df["Sex"] = df.apply(lambda r: first_nonnull(r, sex_cols), axis=1)
        return df

    maps: Dict[str, dict] = {}

    # 1) merged-with-patient sources (data with age and sex already merged in)
    merged_candidates = [
        OUTDIR_ANNOT / "step2_proteomics_proposed_annotations_with_patient.csv",
        OUTDIR_ANNOT / "step2_phospho_proposed_annotations_with_patient.csv",
        OUTDIR_ANNOT / "step2_phospho_siteid_proposed_annotations_with_patient.csv",
    ]
    for p in merged_candidates:
        if p.exists():
            df = pd.read_csv(p)
            if "dataset_name" in df.columns:
                df = _norm_age_sex(df)
                for _, r in df.iterrows():
                    ds = str(r.get("dataset_name", "")).strip()
                    if ds:
                        maps[ds] = {"Age": r.get("Age", pd.NA), "Sex": r.get("Sex", pd.NA)}

    # 2) fill gaps from base CSVs
    base_candidates = [
        OUTDIR_ANNOT / "step2_proteomics_proposed_annotations.csv",
        OUTDIR_ANNOT / "step2_phospho_proposed_annotations.csv",
        OUTDIR_ANNOT / "step2_phospho_siteid_proposed_annotations.csv",
    ]
    for p in base_candidates:
        if p.exists():
            df = pd.read_csv(p)
            if "dataset_name" in df.columns:
                df = _norm_age_sex(df)
                for _, r in df.iterrows():
                    ds = str(r.get("dataset_name", "")).strip()
                    if ds and ds not in maps:
                        maps[ds] = {"Age": r.get("Age", pd.NA), "Sex": r.get("Sex", pd.NA)}

    # 3) Get patient_info via manifest Patient IDs
    pi = pd.read_csv(PATIENT_INFO_TSV, sep=None, engine="python") if Path(PATIENT_INFO_TSV).exists() else pd.DataFrame()
    pcol = next((c for c in ["Patient ID", "Patient", "patient", "patient_id"] if c in pi.columns), None)
    pi_map = {}
    if not pi.empty and pcol:
        pi["Patient_norm"] = pi[pcol].astype(str).str.strip()
        if "Age" not in pi.columns and "age" in pi.columns:
            pi = pi.rename(columns={"age": "Age"})
        if "Sex" not in pi.columns and "sex" in pi.columns:
            pi = pi.rename(columns={"sex": "Sex"})
        keep = ["Patient_norm"] + [c for c in ["Age", "Sex"] if c in pi.columns]
        pi_map = pi[keep].set_index("Patient_norm").to_dict(orient="index")

    man = read_per_sample_manifest_base()
    if not man.empty and {"dataset_name", "Patient"}.issubset(man.columns):
        man["dataset_name"] = man["dataset_name"].astype(str).str.strip()
        man["Patient"] = man["Patient"].astype(str).str.strip()
        for _, r in man.iterrows():
            ds, pt = r["dataset_name"], r["Patient"]
            if ds and pt and (ds not in maps or pd.isna(maps[ds].get("Age")) or pd.isna(maps[ds].get("Sex"))):
                rec = pi_map.get(pt, {}) if pi_map else {}
                if rec:
                    current = maps.get(ds, {})
                    if "Age" in rec and (not current or pd.isna(current.get("Age"))):
                        current["Age"] = rec.get("Age", pd.NA)
                    if "Sex" in rec and (not current or pd.isna(current.get("Sex"))):
                        current["Sex"] = rec.get("Sex", pd.NA)
                    maps[ds] = current

    return maps

def read_per_sample_manifest() -> pd.DataFrame:
    df = read_per_sample_manifest_base()
    if df.empty:
        raise FileNotFoundError("Missing step2 per-sample manifest.")
    needed = ["dataset_name","sample","tumor","Patient","Description","kind","file_path"]
    for c in needed:
        if c not in df.columns:
            df[c] = pd.NA
    df["dataset_name"] = df["dataset_name"].astype(str).str.strip()
    return df

def clean_and_reorder_data_tsv(tsv_path: Path):
    """Drop Age/Sex columns and move 'sample' to be the first column (in-place)."""
    try:
        df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    except Exception:
        return
    cols_lower = {c.lower(): c for c in df.columns}
    drop_cols = []
    for key in ["age","sex"]:
        if key in cols_lower: drop_cols.append(cols_lower[key])
    if drop_cols:
        df = df.drop(columns=drop_cols, errors="ignore")
    if "sample" in [c.lower() for c in df.columns]:
        sample_col = next(c for c in df.columns if c.lower() == "sample")
        new_order = [sample_col] + [c for c in df.columns if c != sample_col]
        df = df[new_order]
    df.to_csv(tsv_path, sep="\t", index=False)

def ensure_folder(syn, name: str, parent_id: str, cache: dict) -> str:
    key = (parent_id, name)
    if key in cache: return cache[key]
    f = Folder(name=name, parent=parent_id)
    f = syn.store(f)
    cache[key] = f["id"]
    return f["id"]

def ensure_path_folders(syn, parent_id: str, rel_dir: Path, cache: dict) -> str:
    cur = parent_id
    for part in rel_dir.parts:
        if part in ("","."): continue
        cur = ensure_folder(syn, part, cur, cache)
    return cur

def _normalize_for_annotations(d):
    out = {}
    for k, v in d.items():
        if _is_nan(v): continue
        if isinstance(v, (list, tuple)):
            vv = []
            for e in v:
                if _is_nan(e): continue
                if isinstance(e, (bool,int,float,str)): vv.append(e)
                else: vv.append(str(e))
            if vv: out[k] = vv
        elif isinstance(v, (bool,int,float,str)):
            out[k] = v
        else:
            out[k] = str(v)
    return out

def apply_annotations(syn, entity, ann_dict: dict):
    ann = _normalize_for_annotations(ann_dict)
    try:
        anns = syn.get_annotations(entity)
        for k, v in ann.items():
            anns[k] = v
        syn.set_annotations(anns)
    except TypeError:
        syn.set_annotations(entity, annotations=ann)

def build_package() -> pd.DataFrame:
    log("[PKG] Building package structure …")
    package_data = PACKAGE_DIR / DATA_ROOT_NAME
    package_meta = PACKAGE_DIR / "metadata"
    safe_mkdir(package_data); safe_mkdir(package_meta)

    # Copy helpful metadata
    for fn in ["step2_per_sample_manifest_with_patient.csv",
               "step2_per_sample_manifest.csv",
               "step2_proteomics_proposed_annotations.csv",
               "step2_phospho_proposed_annotations.csv",
               "step2_phospho_siteid_proposed_annotations.csv"]:
        src = OUTDIR_ANNOT / fn
        if src.exists(): shutil.copy2(src, package_meta / src.name)
    if Path(PATIENT_INFO_TSV).exists():
        shutil.copy2(PATIENT_INFO_TSV, package_meta / Path(PATIENT_INFO_TSV).name)

    df_ps = read_per_sample_manifest()
    age_sex_map = load_age_sex_maps()

    discovered = []
    targets = [
        ("proteomics", OUTDIR_SPLIT / "proteomics"),
        ("phosphoproteomics", OUTDIR_SPLIT / "phosphoproteomics"),
        ("phospho-siteid", OUTDIR_SPLIT / "phospho-siteid"),
    ]

    for kind, base in targets:
        if not base.exists(): continue
        for source_group in sorted([p for p in base.iterdir() if p.is_dir()]):
            for f in sorted(source_group.glob("*.tsv")):
                ds_name = f.name.split(".", 1)[0]
                rel = Path(DATA_ROOT_NAME) / kind / source_group.name / f.name
                dst = PACKAGE_DIR / rel
                safe_mkdir(dst.parent)
                shutil.copy2(f, dst)

                # Drop Age/Sex, move 'sample' first
                clean_and_reorder_data_tsv(dst)

                size_bytes = dst.stat().st_size
                md5 = md5_file(dst)

                row = df_ps[df_ps["dataset_name"] == ds_name]
                if row.empty:
                    sample = tumor = patient = desc = None
                else:
                    r = row.iloc[0]
                    sample  = r.get("sample", pd.NA)
                    tumor   = r.get("tumor", pd.NA)
                    patient = r.get("Patient", pd.NA)
                    desc    = r.get("Description", pd.NA)

                age = age_sex_map.get(ds_name, {}).get("Age")
                sex = age_sex_map.get(ds_name, {}).get("Sex")

                discovered.append({
                    "rel_path": str(rel).replace("\\","/"),
                    "kind": kind,
                    "dataset_name": ds_name,
                    "sample": sample if pd.notna(sample) else None,
                    "tumor":  tumor  if pd.notna(tumor)  else None,
                    "Patient": patient if pd.notna(patient) else None,
                    "Age": int(age) if (age is not None and not pd.isna(age)) else None,
                    "Sex": sex if (sex is not None and not pd.isna(sex)) else None,
                    "Description": desc if pd.notna(desc) else None,
                    "size_bytes": int(size_bytes),
                    "md5": md5,
                })

    if not discovered:
        raise RuntimeError("No files discovered in step2_split/* to package.")

    df_upload = pd.DataFrame(discovered)
    df_upload["dataType"] = df_upload["kind"].map(DATATYPE_BY_KIND).fillna(df_upload["kind"])
    for k, v in ANNOT_DEFAULTS.items():
        df_upload[k] = v

    up_manifest = (PACKAGE_DIR / "metadata" / "upload_manifest.csv")
    df_upload.to_csv(up_manifest, index=False)
    log(f"[PKG] Wrote {up_manifest} ({len(df_upload)} rows)")
    return df_upload

def upload_to_synapse(df_upload: pd.DataFrame):
    out_idx_path = PACKAGE_DIR / "synapse_upload_index.csv"

    if not ENABLE_SYNAPSE_UPLOAD:
        df_upload.assign(synId=None, parentId="(dry-run)").to_csv(out_idx_path, index=False)
        return

    syn = syn_login()
    folder_cache = {}
    uploaded = []

    for _, row in df_upload.iterrows():
        rel_path = Path(row["rel_path"])
        local_path = (PACKAGE_DIR / rel_path).resolve()
        if not local_path.exists():
            raise FileNotFoundError(f"Missing packaged file: {local_path}")

        syn_folder_id = ensure_path_folders(syn, UPLOAD_PARENT_SYN_ID, rel_path.parent, folder_cache)
        entity = File(path=str(local_path), parent=syn_folder_id)
        entity = syn.store(entity)

        # All Annotation Data
        ann = {
            "sample": row.get("sample"),
            "tumor": row.get("tumor"),
            "Age": row.get("Age"),
            "Sex": row.get("Sex"),
            "Description": row.get("Description"),
            "species": row.get("species"),
            "studyId": row.get("studyId"),
            "dataType": row.get("dataType"),
            "studyName": row.get("studyName"),
            "tumorType": row.get("tumorType"),
            "initiative": row.get("initiative"),
            "fundingAgency": row.get("fundingAgency"),
            "disease": row.get("disease"),
            "diagnosis": row.get("diagnosis"),
        }
        if ann.get("Age") is not None:
            try: ann["Age"] = int(ann["Age"])
            except Exception: ann["Age"] = str(ann["Age"])
        if ann.get("Sex") is not None:
            ann["Sex"] = str(ann["Sex"]).strip()

        apply_annotations(syn, entity, ann)

        uploaded.append({
            "synId": entity["id"],
            "name": entity["name"],
            "parentId": entity["parentId"],
            **row.to_dict()
        })

    pd.DataFrame(uploaded).to_csv(out_idx_path, index=False)

# =============================================================================
# Main — Run it all
# =============================================================================


def main():
    """This runs it all. Comment out steps you do not wish to run."""
    # Step 1: inventory/download from your source Synapse folder
    run_step1_inventory()
    download_fixed_inputs_from_synapse()

    # Step 2: process using the freshly downloaded metadata/patient_info
    run_step2_all()

    # Step 3: package & upload
    df_upload = build_package()
    upload_to_synapse(df_upload)

if __name__ == "__main__":
    main()
