"""
AMR Prediction API  —  v2.0.0
Supports: Salmonella enterica  |  Escherichia coli
Pathogen identity is supplied by the caller (Kraken2 output parsed upstream).
"""

from __future__ import annotations

import base64
import json
import os
import pickle
import re
import shutil
import subprocess
import traceback
import uuid
from datetime import datetime
from io import BytesIO
from pathlib import Path
from typing import Any, Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

# ══════════════════════════════════════════════════════════════════
# APP SETUP
# ══════════════════════════════════════════════════════════════════
app = FastAPI(
    title="AMR Prediction API",
    version="2.0.0",
    description=(
        "Predict antibiotic resistance in Salmonella enterica and "
        "Escherichia coli genomes using ML with SHAP explainability. "
        "Intended for research use; clinical decisions require confirmatory AST."
    ),
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ══════════════════════════════════════════════════════════════════
# DIRECTORY CONSTANTS
# ══════════════════════════════════════════════════════════════════
WORK_DIR   = Path(os.getenv("WORK_DIR",   "/app/work"))
MODELS_DIR = Path(os.getenv("MODELS_DIR", "/app/models"))
SCRIPTS_DIR = Path(os.getenv("SCRIPTS_DIR", "/app/scripts"))
TEMPLATES_DIR = Path(os.getenv("FEATURE_TEMPLATES_DIR", "/app/feature_templates"))

# ══════════════════════════════════════════════════════════════════
# PATHOGEN CONFIGURATION
# Each entry fully describes one pathogen's pipeline.
# ══════════════════════════════════════════════════════════════════

# Kraken2 outputs a taxonomy string like:
#   "Salmonella enterica"  /  "Escherichia coli"
# We normalise to lowercase and do substring matching so minor
# version strings ("Salmonella enterica subsp. enterica") still match.

PATHOGEN_CONFIG: Dict[str, Dict[str, Any]] = {

    # ── Salmonella ──────────────────────────────────────────────
    "salmonella": {
        "display_name": "Salmonella enterica",
        # Strings Kraken2 may produce — any substring match triggers this path
        "kraken_aliases": ["salmonella"],
        "scripts_dir": SCRIPTS_DIR / "salmonella",
        "templates_dir": TEMPLATES_DIR / "salmonella",
        "models_dir": MODELS_DIR / "salmonella",
        "reference": Path(os.getenv("SALMONELLA_REFERENCE", "/app/reference/salmonella_LT2.gbff")),
        "pipeline_steps": [
            {"cmd": "bash",    "script": "01_extract_genes.sh",  "timeout": 300},
            {"cmd": "python3", "script": "01b_process_genes.py", "timeout": 60},
            {"cmd": "bash",    "script": "02_create_blast_db.sh","timeout": 60},
            {"cmd": "python3", "script": "03_extract_kmers.py",  "timeout": 600},
            {"cmd": "bash",    "script": "04_extract_snps.sh",   "timeout": 1200},
            {"cmd": "python3", "script": "04b_process_snps.py",  "timeout": 60},
            {"cmd": "python3", "script": "06_align_features.py", "timeout": 60},
        ],
        # Single shared feature file; models load aligned_full.csv
        "prediction_mode": "single_aligned",
        "antibiotics": {
            "pefoxacin":       {"model": "pefoxacin_salmonella_model.pkl",       "aligned": "aligned_full.csv", "who_aware": "Watch"},
            "trimethoprim":    {"model": "trimethoprim_salmonella_model.pkl",     "aligned": "aligned_full.csv", "who_aware": "Access"},
            "sulfamethoxazole":{"model": "sulfamethoxazole_salmonella_model.pkl", "aligned": "aligned_full.csv", "who_aware": "Access"},
            "nalidixic_acid":  {"model": "nalidixic_salmonella_model.pkl",        "aligned": "aligned_full.csv", "who_aware": "Watch"},
            "ampicillin":      {"model": "ampicillin_salmonella_model.pkl",       "aligned": "aligned_full.csv", "who_aware": "Access"},
            "chloramphenicol": {"model": "chloramphenicol_salmonella_model.pkl",  "aligned": "aligned_full.csv", "who_aware": "Access"},
            "streptomycin":    {"model": "streptomycin_salmonella_model.pkl",     "aligned": "aligned_full.csv", "who_aware": "Access"},
            "sulfisoxazole":   {"model": "sulfisoxazole_salmonella_model.pkl",    "aligned": "aligned_full.csv", "who_aware": "Watch"},
        },
        "metadata": {
            "reference_genome": "Salmonella_Typhimurium_LT2",
            "training_samples": 398,
        },
    },

    # ── E. coli ─────────────────────────────────────────────────
    "ecoli": {
        "display_name": "Escherichia coli",
        "kraken_aliases": ["escherichia coli", "e. coli", "e.coli"],
        "scripts_dir": SCRIPTS_DIR / "e_coli",
        "templates_dir": TEMPLATES_DIR / "e_coli",
        "models_dir": MODELS_DIR / "e_coli",
        "reference": Path(os.getenv("ECOLI_REFERENCE", "/app/reference/escherichia_coli.gbff")),
        "pipeline_steps": [
            {"cmd": "bash",    "script": "01_extract_genes.sh",  "timeout": 300},
            {"cmd": "python3", "script": "01b_process_genes.py", "timeout": 60},
            {"cmd": "bash",    "script": "02_create_blast_db.sh","timeout": 60},
            {"cmd": "python3", "script": "03_extract_kmers.py",  "timeout": 600},
            {"cmd": "bash",    "script": "04_extract_snps.sh",   "timeout": 1200},
            {"cmd": "python3", "script": "04b_process_snps.py",  "timeout": 60},
            {"cmd": "python3", "script": "06_align_features.py", "timeout": 60},
        ],
        # Per-antibiotic aligned CSVs
        "prediction_mode": "per_antibiotic_aligned",
        "antibiotics": {
            "amoxicillin_clavulanic_acid": {"model": "amoxicillin_ecoli_model.pkl",   "aligned": "aligned_amoxicillin.csv", "who_aware": "Access"},
            "cefuroxime":                  {"model": "cefuroxime_ecoli_model.pkl",    "aligned": "aligned_coxime.csv",      "who_aware": "Watch"},
            "doxycycline":                 {"model": "doxycycline_ecoli_model.pkl",   "aligned": "aligned_doxy.csv",        "who_aware": "Access"},
            "levofloxacin":                {"model": "levofloxacin_ecoli_model.pkl",  "aligned": "aligned_levo.csv",        "who_aware": "Watch"},
            "nalidixic_acid":              {"model": "nalidixic_acid_ecoli_model.pkl","aligned": "aligned_nali.csv",        "who_aware": "Watch"},
            "norfloxacin":                 {"model": "norfloxacin_ecoli_model.pkl",   "aligned": "aligned_norflo.csv",      "who_aware": "Watch"},
            "streptomycin":                {"model": "streptomycin_ecoli_model.pkl",  "aligned": "aligned_strep.csv",       "who_aware": "Access"},
            "tetracycline":                {"model": "tetracycline_ecoli_model.pkl",  "aligned": "aligned_tetra.csv",       "who_aware": "Access"},
        },
        "metadata": {
            "reference_genome": "Escherichia_coli_K12_MG1655",
            "training_samples": 571,   
        },
    },
}

# ══════════════════════════════════════════════════════════════════
# PATHOGEN DETECTION  —  Kraken2 output normalisation
# ══════════════════════════════════════════════════════════════════

def resolve_pathogen(kraken_species: str) -> Optional[str]:
    """
    Map a free-text Kraken2 species string to one of our supported
    pathogen keys ('salmonella' | 'ecoli').

    Kraken2 may return strings like:
        "Salmonella enterica"
        "Salmonella enterica subsp. enterica serovar Typhimurium"
        "Escherichia coli O157:H7"
        "Escherichia coli"

    We lowercase and scan for any registered alias as a substring.
    Returns None if no match, allowing the caller to 400-reject.
    """
    normalised = kraken_species.lower().strip()
    # Remove strain / serovar decorators — keep genus+species
    normalised = re.sub(r'\s+(subsp|serovar|str|strain|serotype).*$', '', normalised)

    for key, cfg in PATHOGEN_CONFIG.items():
        for alias in cfg["kraken_aliases"]:
            if alias in normalised:
                return key
    return None


# ══════════════════════════════════════════════════════════════════
# PIPELINE RUNNER
# ══════════════════════════════════════════════════════════════════

def run_pipeline(work_path: Path, pathogen_key: str) -> None:
    """
    Execute the ordered pre-processing steps for one pathogen.
    Each script receives the work directory via an environment variable
    so scripts never need hardcoded paths.
    """
    cfg = PATHOGEN_CONFIG[pathogen_key]
    scripts_dir = cfg["scripts_dir"]

    env = os.environ.copy()
    # Inject runtime paths so every script can read them
    env["WORK_DIR"]       = str(work_path)
    env["SCRIPTS_DIR"]    = str(scripts_dir)
    env["TEMPLATES_DIR"]  = str(cfg["templates_dir"])
    env["MODELS_DIR"]     = str(cfg["models_dir"])
    env["REFERENCE"]      = str(cfg["reference"])

    for step in cfg["pipeline_steps"]:
        script_path = scripts_dir / step["script"]
        cmd = [step["cmd"], str(script_path)]
        try:
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                timeout=step["timeout"],
                cwd=str(work_path),
                env=env,
            )
        except subprocess.TimeoutExpired as exc:
            raise HTTPException(
                status_code=504,
                detail=f"Pipeline timeout in {step['script']} "
                       f"(limit {step['timeout']}s). "
                       "Consider splitting the genome or increasing timeout.",
            ) from exc
        except subprocess.CalledProcessError as exc:
            stderr = exc.stderr.decode(errors="replace") if exc.stderr else "no stderr"
            raise HTTPException(
                status_code=500,
                detail=f"Pipeline failed at {step['script']}: {stderr[:300]}",
            ) from exc


# ══════════════════════════════════════════════════════════════════
# SHAP UTILITIES
# ══════════════════════════════════════════════════════════════════

def _get_shap_values(model, X: pd.DataFrame):
    """Return 1-D SHAP value array for the Resistant class."""
    explainer = shap.TreeExplainer(model)
    sv = explainer.shap_values(X)
    if isinstance(sv, list):          # multi-class tree models
        return sv[1], explainer.expected_value[1]
    return sv[0], explainer.expected_value   # single-output


def get_shap_evidence(model, X: pd.DataFrame, top_n: int = 5) -> List[Dict]:
    """Return top-N SHAP-ranked features with metadata."""
    try:
        sv, _ = _get_shap_values(model, X)
        impacts = pd.DataFrame({"feature": X.columns, "shap_value": sv})
        impacts = impacts.reindex(
            impacts["shap_value"].abs().sort_values(ascending=False).index
        )
        evidence = []
        for _, row in impacts.head(top_n).iterrows():
            feat = row["feature"]
            impact = float(row["shap_value"])
            # Feature-type heuristics preserved from v1
            if re.search(r'NC_\d+.*>', feat):
                ftype = "SNP"
            elif len(feat) == 10 and feat.isalpha() and feat.isupper():
                ftype = "K-mer"
            else:
                ftype = "Gene"
            evidence.append({
                "feature":      feat,
                "type":         ftype,
                "impact_score": round(impact, 4),
                "effect":       "promotes_resistance" if impact > 0 else "promotes_susceptibility",
            })
        return evidence
    except Exception as exc:
        raise HTTPException(500, f"SHAP evidence extraction failed: {exc}") from exc


def create_force_plot(model, X: pd.DataFrame, label: str) -> str:
    """Return base64-encoded PNG SHAP force plot."""
    try:
        sv, ev = _get_shap_values(model, X)
        plt.figure(figsize=(14, 3))
        shap.force_plot(ev, sv, X.iloc[0], matplotlib=True, show=False, text_rotation=10)
        plt.title(f"{label} — Feature Contributions", fontsize=13, fontweight="bold")
        buf = BytesIO()
        plt.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
        plt.close()
        buf.seek(0)
        return "data:image/png;base64," + base64.b64encode(buf.read()).decode()
    except Exception as exc:
        raise HTTPException(500, f"SHAP force plot failed: {exc}") from exc


# ══════════════════════════════════════════════════════════════════
# CONFIDENCE + ACTION HELPERS
# ══════════════════════════════════════════════════════════════════

def _confidence_and_action(prob: float) -> tuple[str, str]:
    if prob >= 0.85:
        return "High", "REPORT_FINAL"
    if prob >= 0.65:
        return "Medium", "CONSIDER_CONFIRMATION"
    return "Low", "CONFIRMATORY_AST_REQUIRED"


# ══════════════════════════════════════════════════════════════════
# PREDICTION ENGINE
# ══════════════════════════════════════════════════════════════════

def _load_model(path: Path):
    with open(path, "rb") as fh:
        saved = pickle.load(fh)
    # Some models are wrapped in a dict {"model": <clf>, ...}
    return saved["model"] if isinstance(saved, dict) else saved


def _reorder_features(model, X: pd.DataFrame) -> pd.DataFrame:
    """Align DataFrame columns to model's expected feature order."""
    if hasattr(model, "get_booster"):          # XGBoost
        expected = model.get_booster().feature_names
    elif hasattr(model, "feature_names_in_"):  # sklearn ≥ 1.0
        expected = list(model.feature_names_in_)
    else:
        return X                               # trust the CSV order

    missing = [f for f in expected if f not in X.columns]
    if missing:
        X = pd.concat(
            [X, pd.DataFrame(0, index=X.index, columns=missing)], axis=1
        )
    return X[expected]


def predict_single_antibiotic(
    antibiotic_key: str,
    cfg: Dict,
    models_dir: Path,
    work_path: Path,
) -> Dict:
    """
    Run prediction for one antibiotic and return a fully structured result dict.
    Produces SHAP evidence + force plot for the winning model.
    """
    model_path = models_dir / cfg["model"]
    aligned_path = work_path / cfg["aligned"]

    if not model_path.exists():
        raise HTTPException(500, f"Model not found: {model_path}")
    if not aligned_path.exists():
        raise HTTPException(500, f"Aligned feature file not found: {aligned_path}")

    model = _load_model(model_path)
    df    = pd.read_csv(aligned_path)
    X     = _reorder_features(model, df.drop(columns=["Genome_ID"]))

    proba        = model.predict_proba(X)[0]
    prob_res     = float(proba[1])
    prob_sus     = float(proba[0])
    prediction   = "Resistant" if prob_res >= 0.5 else "Susceptible"
    final_prob   = prob_res if prediction == "Resistant" else prob_sus
    confidence, action = _confidence_and_action(final_prob)

    evidence   = get_shap_evidence(model, X)
    force_plot = create_force_plot(model, X, antibiotic_key.replace("_", " ").title())

    return {
        "who_aware_category":  cfg.get("who_aware", "Unknown"),
        "phenotype":           prediction,
        "probability_resistant": round(prob_res, 4),
        "probability_susceptible": round(prob_sus, 4),
        "probability_score":   round(final_prob, 4),
        "confidence_category": confidence,
        "action_required":     action,
        "clinical_note":       _clinical_note(prediction, confidence, action),
        "evidence":            evidence,
        "shap_force_plot":     force_plot,
    }


def _clinical_note(phenotype: str, confidence: str, action: str) -> str:
    """
    Plain-language note suitable for a clinical dashboard.
    This is research-grade; always append the standard disclaimer.
    """
    if action == "REPORT_FINAL":
        return (
            f"Model predicts {phenotype.upper()} with HIGH confidence. "
            "Results may be reported pending clinical review."
        )
    if action == "CONSIDER_CONFIRMATION":
        return (
            f"Model predicts {phenotype.upper()} with MEDIUM confidence. "
            "Confirmatory phenotypic AST is advisable before clinical use."
        )
    return (
        f"Model predicts {phenotype.upper()} with LOW confidence. "
        "Phenotypic AST is REQUIRED before any clinical decision."
    )


# ══════════════════════════════════════════════════════════════════
# QUALITY METRICS
# ══════════════════════════════════════════════════════════════════

def collect_quality_metrics(work_path: Path, file_size: int) -> Dict:
    metrics: Dict[str, Any] = {"genome_size_mb": round(file_size / (1024 ** 2), 2)}
    for label, filename in [
        ("genes_detected",  "gene_presence_production.csv"),
        ("kmers_matched",   "kmer_production.csv"),
        ("snps_detected",   "snp_production.csv"),
    ]:
        fp = work_path / filename
        if fp.exists():
            df = pd.read_csv(fp)
            metrics[label] = max(0, df.shape[1] - 1)
        else:
            metrics[label] = 0

    total = metrics.get("genes_detected", 0) + metrics.get("kmers_matched", 0) + metrics.get("snps_detected", 0)
    metrics["total_features_extracted"] = total
    metrics["passed_qc"] = total > 10     # at least some features must be found
    return metrics


# ══════════════════════════════════════════════════════════════════
# GLOBAL EXCEPTION HANDLER
# ══════════════════════════════════════════════════════════════════

@app.exception_handler(Exception)
async def global_exception_handler(request, exc):
    return JSONResponse(
        status_code=500,
        content={
            "status": "error",
            "message": str(exc),
            "type": type(exc).__name__,
            "disclaimer": DISCLAIMER,
            "timestamp": datetime.utcnow().isoformat() + "Z",
        },
    )


# ══════════════════════════════════════════════════════════════════
# STANDARD DISCLAIMER  (injected into every response)
# ══════════════════════════════════════════════════════════════════

DISCLAIMER = (
    "This tool is intended for research and epidemiological surveillance. "
    "Predictions must NOT be used as the sole basis for clinical treatment "
    "decisions. All Resistant predictions at Medium or Low confidence require "
    "confirmatory phenotypic antimicrobial susceptibility testing (AST) "
    "before clinical use. Results should be interpreted by a qualified "
    "clinical microbiologist."
)


# ══════════════════════════════════════════════════════════════════
# ENDPOINTS
# ══════════════════════════════════════════════════════════════════

@app.get("/health")
def health_check():
    """Verify that models for both pathogens are present on disk."""
    missing: List[str] = []
    for pk, pcfg in PATHOGEN_CONFIG.items():
        for ab, acfg in pcfg["antibiotics"].items():
            model_file = pcfg["models_dir"] / acfg["model"]
            if not model_file.exists():
                missing.append(str(model_file))

    ok = len(missing) == 0
    return {
        "status": "healthy" if ok else "degraded",
        "version": "2.0.0",
        "pathogens_configured": list(PATHOGEN_CONFIG.keys()),
        "missing_models": missing,
        "timestamp": datetime.utcnow().isoformat() + "Z",
    }


@app.post("/predict")
async def predict_resistance(
    genome: UploadFile = File(...),
    kraken_species: str = Form(...),
):
    """
    Predict antibiotic resistance from an assembled genome.

    Parameters
    ----------
    genome          : FASTA / FNA / FA file (1 KB – 15 MB)
    kraken_species  : Species string from Kraken2 classification
                      e.g. "Salmonella enterica" or "Escherichia coli"

    Returns
    -------
    Structured JSON with per-antibiotic predictions, SHAP evidence,
    quality metrics, and the standard clinical disclaimer.
    """
    job_id    = str(uuid.uuid4())[:8]
    work_path = WORK_DIR / job_id
    start_time = datetime.utcnow()

    try:
        # ── 1. Resolve pathogen ────────────────────────────────
        pathogen_key = resolve_pathogen(kraken_species)
        if pathogen_key is None:
            raise HTTPException(
                status_code=422,
                detail=(
                    f"Unsupported or unrecognised species: '{kraken_species}'. "
                    f"Supported: Salmonella enterica, Escherichia coli."
                ),
            )

        pcfg = PATHOGEN_CONFIG[pathogen_key]

        # ── 2. Validate genome file ───────────────────────────
        if not genome.filename.endswith((".fna", ".fasta", ".fa")):
            raise HTTPException(400, "Invalid file format. Upload .fna, .fasta, or .fa")

        content   = await genome.read()
        file_size = len(content)

        if file_size > 15 * 1024 * 1024:
            raise HTTPException(400, "File too large. Maximum 15 MB.")
        if file_size < 1000:
            raise HTTPException(400, "File too small for a valid genome assembly.")

        # Quick FASTA sanity check
        first_line = content[:200].decode(errors="replace").strip().splitlines()[0] if content else ""
        if not first_line.startswith(">"):
            raise HTTPException(400, "File does not appear to be a valid FASTA file (missing '>' header).")

        # ── 3. Write genome to work directory ─────────────────
        work_path.mkdir(parents=True, exist_ok=True)
        genome_path = work_path / "query_genome.fna"
        genome_path.write_bytes(content)

        # ── 4. Run bioinformatics pipeline ────────────────────
        run_pipeline(work_path, pathogen_key)

        # ── 5. Collect quality metrics ────────────────────────
        qc = collect_quality_metrics(work_path, file_size)

        if not qc["passed_qc"]:
            raise HTTPException(
                422,
                detail=(
                    "Quality control failed: no genomic features could be extracted. "
                    "Verify the genome assembly is complete and corresponds to "
                    f"{pcfg['display_name']}."
                ),
            )

        # ── 6. Get genome ID from first aligned file ──────────
        first_aligned_name = next(iter(pcfg["antibiotics"].values()))["aligned"]
        first_aligned_path = work_path / first_aligned_name
        genome_id = pd.read_csv(first_aligned_path, usecols=["Genome_ID"])["Genome_ID"].iloc[0]

        # ── 7. Run predictions for all antibiotics ────────────
        susceptibility: Dict[str, Any] = {}
        shap_plots: Dict[str, str] = {}

        for ab_key, ab_cfg in pcfg["antibiotics"].items():
            result = predict_single_antibiotic(
                ab_key, ab_cfg, pcfg["models_dir"], work_path
            )
            # Separate SHAP plots from main body for clean output
            shap_plots[ab_key] = result.pop("shap_force_plot")
            susceptibility[ab_key] = result

        end_time = datetime.utcnow()

        # ── 8. Build response ─────────────────────────────────
        output = {
            "job_id": job_id,
            "status": "completed",
            "disclaimer": DISCLAIMER,
            "sample_info": {
                "sample_id": genome_id,
                "original_filename": genome.filename,
            },
            "organism_identification": {
                "detected_species":      pcfg["display_name"],
                "kraken_input":          kraken_species,
                "resolved_pathogen_key": pathogen_key,
            },
            "timestamps": {
                "submitted_at":           start_time.isoformat() + "Z",
                "completed_at":           end_time.isoformat() + "Z",
                "processing_time_seconds": int((end_time - start_time).total_seconds()),
            },
            "quality_control": {
                "passed_qc": qc["passed_qc"],
                "metrics": {
                    "genome_size_mb":          qc["genome_size_mb"],
                    "genes_detected":          qc["genes_detected"],
                    "kmers_matched":           qc["kmers_matched"],
                    "snps_detected":           qc["snps_detected"],
                    "total_features_extracted": qc["total_features_extracted"],
                },
            },
            "antimicrobial_susceptibility": susceptibility,
            "visualizations": {
                "shap_force_plots": shap_plots,
            },
            "model_metadata": {
                "pipeline_version": "2.0.0",
                "pathogen":         pathogen_key,
                "reference_genome": pcfg["metadata"]["reference_genome"],
                "training_samples": pcfg["metadata"]["training_samples"],
                "card_version":     "2024.01",
                "confidence_thresholds": {
                    "High":   "probability ≥ 0.85  →  REPORT_FINAL",
                    "Medium": "probability 0.65–0.84  →  CONSIDER_CONFIRMATION",
                    "Low":    "probability < 0.65  →  CONFIRMATORY_AST_REQUIRED",
                },
            },
        }

        return output

    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(
            500,
            detail=f"Unexpected error: {exc}\n{traceback.format_exc()[:600]}",
        ) from exc
    finally:
        # Always clean up, even on error
        if work_path.exists():
            shutil.rmtree(work_path, ignore_errors=True)


# ── Local development entry-point ─────────────────────────────
if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)