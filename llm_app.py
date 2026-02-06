import io
import html
import json
import re
import subprocess
import os
from dataclasses import dataclass
from typing import List, Optional, Tuple, Dict

import requests
try:
    from groq import Groq
except Exception:  # pragma: no cover
    Groq = None
import streamlit as st
import streamlit.components.v1 as components
try:
    import vl_convert as vlc  # optional export for PNG/SVG/PDF
except Exception:  # pragma: no cover
    vlc = None
try:
    import cairosvg  # optional PDF export
except Exception:  # pragma: no cover
    cairosvg = None

try:
    import pdfplumber  # optional fallback
except Exception:  # pragma: no cover
    pdfplumber = None
try:
    import pytesseract  # optional OCR
except Exception:  # pragma: no cover
    pytesseract = None


APP_TITLE = "ClinPDF-Report Analyzer"
CBIO_BASE = "https://www.cbioportal.org"
DEFAULT_OLLAMA_URL = os.environ.get("OLLAMA_HOST", "http://127.0.0.1:11434")
DEFAULT_GROQ_MODEL = os.environ.get("GROQ_MODEL", "llama-3.1-70b-versatile")
DEFAULT_MODEL = "llama3.1:70b"
DEFAULT_MAX_CHARS = 16000



def _rows_to_csv(rows: List[Dict[str, str]]) -> str:
    if not rows:
        return ""
    fieldnames = sorted({k for r in rows for k in r.keys()})
    output = io.StringIO()
    output.write(",".join(fieldnames) + "\n")
    for r in rows:
        output.write(",".join(str(r.get(k, "")).replace(",", " ") for k in fieldnames) + "\n")
    return output.getvalue()


def _rows_to_tsv(rows: List[Dict[str, str]]) -> str:
    if not rows:
        return ""
    fieldnames = list(rows[0].keys())
    output = io.StringIO()
    output.write("\t".join(fieldnames) + "\n")
    for r in rows:
        output.write("\t".join(str(r.get(k, "")) for k in fieldnames) + "\n")
    return output.getvalue()


def _parse_vaf(v):
    if v is None:
        return None
    s = str(v).strip()
    if not s:
        return None
    s = s.replace("%", "")
    try:
        return float(s)
    except ValueError:
        return None


def _build_heatmap_data(rows: List[Dict[str, str]]) -> List[Dict[str, object]]:
    bucket: Dict[tuple, float] = {}
    for r in rows:
        patient = r.get("Patient Name") or "Unknown"
        gene = r.get("Gene") or ""
        prot = r.get("Protein change") or ""
        if not gene:
            continue
        gv = f"{gene}:{prot}".strip()
        vaf = _parse_vaf(r.get("VAF"))
        if vaf is None:
            continue
        key = (patient, gv)
        if key not in bucket or vaf > bucket[key]:
            bucket[key] = vaf
    data = [{"Patient": k[0], "GeneVariant": k[1], "VAF": v} for k, v in bucket.items()]
    return data


def _truncate(text: str, n: int = 160) -> str:
    if not text:
        return ""
    s = str(text)
    if len(s) <= n:
        return s
    return s[:n].rstrip() + "..."


def _fill_missing(row: Dict[str, object], fill_value: str = ".") -> Dict[str, object]:
    out = {}
    for k, v in row.items():
        if v is None:
            out[k] = fill_value
        elif isinstance(v, str) and v.strip() == "":
            out[k] = fill_value
        else:
            out[k] = v
    return out


def _render_proteinpaint(gene: str, genome: str, dataset: str, host: str) -> None:
    link = "https://proteinpaint.stjude.org/?gene=FBXW7&genome=hg38&dataset=pediatric"
    st.markdown("**ProteinPaint**")
    st.markdown(
        f"""
        <a href="{link}" target="_blank" rel="noreferrer"
           style="display:inline-block;padding:8px 14px;background:#2563eb;color:#fff;
                  border-radius:6px;text-decoration:none;font-weight:600;">
          Open ProteinPaint (FBXW7)
        </a>
        """,
        unsafe_allow_html=True,
    )


def _build_cbioportal_input(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    data = []
    for r in rows:
        data.append(
            {
                "Hugo_Symbol": r.get("Gene", ""),
                "Sample_ID": r.get("Patient Name", ""),
                "Protein_Change": r.get("Protein change", ""),
            }
        )
    return data


@dataclass
class VariantHit:
    gene: str
    protein_change: str
    cdna_change: Optional[str]
    genomic_change: Optional[str]
    vaf: Optional[str]
    raw_line: str


GENE_LINE_RE = re.compile(r"^\s*([A-Z0-9]{2,})\b.*?p\.([A-Za-z]{3}\d+[A-Za-z0-9]+)")
CNV_GENE_RE = re.compile(r"^\s*([A-Z0-9]{2,})\s*\(CNV\)", re.IGNORECASE)
PROTEIN_RE = re.compile(r"\bp\.([A-Za-z]{3}\d+[A-Za-z0-9]+)\b", re.IGNORECASE)
CDNA_RE = re.compile(r"(?:cDNA\s*change\s*:?\s*)?c\.\s*([0-9_]+[A-Za-z0-9>_]+)", re.IGNORECASE)
CDNA_INLINE_RE = re.compile(r"c\.\s*([0-9_]+[A-Za-z0-9>_]+)", re.IGNORECASE)
GENOMIC_RE = re.compile(r"chr\w+:g\.[0-9_]+(?:[ACGT]>[ACGT]|ins[ACGT]+|del[ACGT]+|dup[ACGT]+)", re.IGNORECASE)
VAF_RE = re.compile(r"(?:Variant Allele Frequency|VAF)\s*[-:]*\s*([0-9.]+%?)", re.IGNORECASE)

NAME_RE = re.compile(r"^([A-Za-z][A-Za-z\s.'-]{1,})$")
AGE_SEX_RE = re.compile(r"AGE:\s*([0-9]+)\s*years\s*\|\s*Gender:\s*([A-Za-z]+)", re.IGNORECASE)
DATE_RECEIVED_RE = re.compile(r"Date Received:\s*([0-9]{1,2}\s+[A-Za-z]+\s+[0-9]{4})", re.IGNORECASE)
DATE_LINE_RE = re.compile(r"\b[0-9]{1,2}\s+[A-Za-z]+\s+[0-9]{4}\b")
DATE_LINE_RE_2 = re.compile(r"\b[0-9]{1,2}[-/][0-9]{1,2}[-/][0-9]{2,4}\b")
DATE_LINE_RE_3 = re.compile(r"\b[0-9]{1,2}\s+[A-Za-z]{3}\s+[0-9]{4}\b", re.IGNORECASE)
CLIN_BG_RE = re.compile(r"CLINICAL BACKGROUND\s*:?\s*(.*)", re.IGNORECASE)
INTERPRETATION_RE = re.compile(r"Interpretation\s*(.+)", re.IGNORECASE)
EXON_RE = re.compile(r"Exon:\s*([0-9]+)", re.IGNORECASE)
NUC_CHANGE_RE = re.compile(r"Nucleotide change:\s*(chr\w+:g\.[0-9_]+[A-Za-z0-9>_]+)", re.IGNORECASE)
GENOMIC_RE = re.compile(r"chr\w+:g\.[0-9_]+[A-Za-z0-9>_]+", re.IGNORECASE)
TRANSCRIPT_RE = re.compile(r"Transcript ID:\s*([A-Za-z0-9_.-]+)", re.IGNORECASE)
VAF_RE = re.compile(r"(?:Variant Allele Frequency|VAF)\s*[-:]*\s*([0-9.]+%?)", re.IGNORECASE)
VARIANT_DEPTH_RE = re.compile(r"Variant Allele Depth/Total depth:\s*([0-9]+/[0-9]+x?)", re.IGNORECASE)
POP_MAF_RE = re.compile(r"Population MAF:\s*([^;\n]+(?:;[^\n]+)?)", re.IGNORECASE)
INSILICO_RE = re.compile(r"In-silico Predictions:\s*([^\n]+)", re.IGNORECASE)
GENE_FUNCTION_RE = re.compile(r"Gene Function:\s*([A-Za-z0-9/ \-]+)", re.IGNORECASE)
GENE_SUMMARY_RE = re.compile(r"Gene Summary:\s*(.*)", re.IGNORECASE)
CNV_TYPE_RE = re.compile(r"CNV type:\s*([^\n]+)", re.IGNORECASE)
CLINVAR_RE = re.compile(r"\b(?:RCV\d+|SCV\d+|VCV\d+|ClinVar\s*ID\s*:\s*\d+)\b", re.IGNORECASE)
PMID_RE = re.compile(r"PMID\s*:\s*(\d+)", re.IGNORECASE)
CLIN_RELEVANCE_RE = re.compile(r"Clinical relevance\s*:?\s*(.*)", re.IGNORECASE)
CLIN_THER_RE = re.compile(r"Clinical and Therapeutic Relevance\s*:?\s*(.*)", re.IGNORECASE)
PUBMED_REFS_RE = re.compile(r"PubMed References", re.IGNORECASE)


@st.cache_data(show_spinner=False)
def pdf_to_text_bytes(pdf_bytes: bytes) -> str:
    # Primary: pdftotext if available
    try:
        proc = subprocess.run(
            ["pdftotext", "-", "-"],
            input=pdf_bytes,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
        )
        text = proc.stdout.decode("utf-8", errors="ignore")
        if text.strip():
            return text
    except Exception:
        pass

    # Fallback: pdfplumber if installed
    if pdfplumber is None:
        return ""

    text_parts = []
    with pdfplumber.open(io.BytesIO(pdf_bytes)) as pdf:
        for page in pdf.pages:
            t = page.extract_text()
            if t:
                text_parts.append(t)
            elif pytesseract is not None:
                img = page.to_image(resolution=300).original
                text_parts.append(pytesseract.image_to_string(img))
    return "\n".join(text_parts)


@st.cache_data(show_spinner=False)
def extract_variants(text: str) -> List[VariantHit]:
    hits: List[VariantHit] = []
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    bad_gene_tokens = {
        "PROTEIN",
        "GENE",
        "VARIANT",
        "CLINICAL",
        "ACTIONABLE",
        "ADDITIONAL",
        "INTERPRETATION",
        "DATE",
        "REPORT",
        "SPECIMEN",
        "SAMPLE",
        "TRANSCRIPT",
        "NUCLEOTIDE",
        "THE",
    }

    for ln in lines:
        m = GENE_LINE_RE.search(ln)
        if not m:
            continue
        if m.group(1).upper() in bad_gene_tokens:
            continue
        gene = m.group(1)
        protein_change = f"p.{m.group(2)}"
        cdna_match = CDNA_RE.search(ln)
        genomic_match = GENOMIC_RE.search(ln)
        vaf_match = VAF_RE.search(ln)
        hits.append(
            VariantHit(
                gene=gene,
                protein_change=protein_change,
                cdna_change=f"c.{cdna_match.group(1)}" if cdna_match else None,
                genomic_change=genomic_match.group(0) if genomic_match else None,
                vaf=vaf_match.group(1) if vaf_match else None,
                raw_line=ln,
            )
        )

    # If nothing matched the strict line pattern, try to assemble from nearby fields.
    if not hits:
        gene_candidates = set(re.findall(r"\b([A-Z0-9]{2,})\b", text))
        # Heuristic: only keep genes that appear near "Variant" or "p."
        for ln in lines:
            if "p." not in ln:
                continue
            p_match = PROTEIN_RE.search(ln)
            if not p_match:
                continue
            protein_change = f"p.{p_match.group(1)}"
            cdna_match = CDNA_RE.search(ln)
            genomic_match = GENOMIC_RE.search(ln)
            vaf_match = VAF_RE.search(ln)

            # pick the first gene-like token on the line if present
            gene = None
            for token in re.findall(r"\b([A-Z0-9]{2,})\b", ln):
                if token in gene_candidates:
                    gene = token
                    break
            if gene:
                hits.append(
                    VariantHit(
                        gene=gene,
                        protein_change=protein_change,
                        cdna_change=f"c.{cdna_match.group(1)}" if cdna_match else None,
                        genomic_change=genomic_match.group(0) if genomic_match else None,
                        vaf=vaf_match.group(1) if vaf_match else None,
                        raw_line=ln,
                    )
                )

    return hits


def _extract_variant_blocks(text: str) -> List[Dict[str, str]]:
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    blocks: List[List[str]] = []
    current: List[str] = []

    # Prefer ACTIONABLE BIOMARKER DETAILS section if present
    actionable_start = None
    actionable_end = None
    for i, ln in enumerate(lines):
        if "ACTIONABLE BIOMARKER DETAILS" in ln.upper():
            actionable_start = i + 1
            break
    if actionable_start is not None:
        for j in range(actionable_start, len(lines)):
            up = lines[j].upper()
            if "ADDITIONAL BIOMARKERS DETECTED" in up or "GLOSSARY" in up or "DISCLAIMER" in up:
                actionable_end = j
                break
    if actionable_start is not None:
        lines = lines[actionable_start:actionable_end]

    def _skip_noise(ln: str) -> bool:
        up = ln.upper()
        if up.startswith("MEDGENOME LABS") or up.startswith("TEL") or up.startswith("WWW."):
            return True
        if up.startswith("PAGE ") or up.startswith("AGE:") or up.startswith("COUPON/UID:"):
            return True
        if up.startswith("PROTEIN CHANGE:"):
            return True
        return False

    bad_gene_tokens = {
        "PROTEIN",
        "GENE",
        "VARIANT",
        "CLINICAL",
        "ACTIONABLE",
        "ADDITIONAL",
        "INTERPRETATION",
        "DATE",
        "REPORT",
        "SPECIMEN",
        "SAMPLE",
        "TRANSCRIPT",
        "NUCLEOTIDE",
        "THE",
    }

    for ln in lines:
        if _skip_noise(ln):
            continue
        m = GENE_LINE_RE.search(ln)
        m_cnv = CNV_GENE_RE.search(ln)
        if m_cnv and m_cnv.group(1).upper() not in bad_gene_tokens:
            if current:
                blocks.append(current)
            current = [ln]
            continue
        if m and m.group(1).upper() not in bad_gene_tokens:
            if current:
                blocks.append(current)
            current = [ln]
        elif current:
            current.append(ln)
    if current:
        blocks.append(current)

    results: List[Dict[str, str]] = []
    for block in blocks:
        joined = " ".join(block)
        # Normalize missing spaces before common labels (e.g., "KRASExon: 2")
        joined = re.sub(r"([A-Za-z])([A-Z][a-z]+:)", r"\1 \2", joined)
        joined = re.sub(r"(Gene:)([A-Z0-9])", r"\1 \2", joined)
        joined = re.sub(r"(Exon:)(\\d)", r"\1 \2", joined)
        joined = re.sub(r"(Variant Allele Frequency:)([0-9])", r"\1 \2", joined)
        joined = re.sub(r"(Nucleotide change:)(chr)", r"\1 \2", joined, flags=re.IGNORECASE)
        joined = re.sub(r"(cDNA change:)(c\\.)", r"\1 \2", joined, flags=re.IGNORECASE)
        m = GENE_LINE_RE.search(block[0])
        m_cnv = CNV_GENE_RE.search(block[0])
        if not m and not m_cnv:
            continue
        if m_cnv:
            gene = m_cnv.group(1)
            protein_change = "CNV"
            variant_type = "CNV"
        else:
            gene = m.group(1)
            protein_change = f"p.{m.group(2)}"
            variant_type = "SNV/Indel"
        cdna_match = CDNA_RE.search(joined)
        if cdna_match is None:
            cdna_match = CDNA_INLINE_RE.search(joined)
        genomic_match = GENOMIC_RE.search(joined)
        vaf_match = VAF_RE.search(joined)

        exon = None
        nucleotide_change = None
        transcript_id = None
        variant_depth = None
        vaf = vaf_match.group(1) if vaf_match else ""
        population_maf = None
        insilico = None
        gene_function = None
        gene_summary = None
        cnv_type = None
        cnv_raw_full = None
        interpretation = None
        clinical_relevance = None
        clinvar_ids = set()
        pubmed_ids = set()

        in_gene_summary = False
        gene_summary_parts: List[str] = []
        in_clin_relevance = False
        clin_rel_parts: List[str] = []
        in_pubmed_refs = False
        pubmed_numbers: List[str] = []

        for idx, ln in enumerate(block):
            if _skip_noise(ln):
                continue
            # handle split headings
            next_ln = block[idx + 1] if idx + 1 < len(block) else ""
            if ln.strip().upper() == "GENE SUMMARY" and gene_summary is None:
                in_gene_summary = True
                continue
            if ln.strip().upper() == "CLINICAL AND" and "THERAPEUTIC RELEVANCE" in next_ln.upper():
                in_clin_relevance = True
                continue

            if exon is None:
                em = EXON_RE.search(ln)
                if em:
                    exon = em.group(1)
            if nucleotide_change is None:
                nm = NUC_CHANGE_RE.search(ln)
                if nm:
                    nucleotide_change = nm.group(1).strip()
            if nucleotide_change is None:
                gm = GENOMIC_RE.search(ln)
                if gm:
                    nucleotide_change = gm.group(0).strip()
            if transcript_id is None:
                tm = TRANSCRIPT_RE.search(ln)
                if tm:
                    transcript_id = tm.group(1).strip()
            if variant_depth is None:
                dm = VARIANT_DEPTH_RE.search(ln)
                if dm:
                    variant_depth = dm.group(1).strip()
            if population_maf is None:
                pm = POP_MAF_RE.search(ln)
                if pm:
                    population_maf = pm.group(1).strip()
            if insilico is None:
                im2 = INSILICO_RE.search(ln)
                if im2:
                    insilico = im2.group(1).strip()
            if gene_function is None:
                gf = GENE_FUNCTION_RE.search(ln)
                if gf:
                    gene_function = gf.group(1).strip()
            if cnv_type is None:
                ct = CNV_TYPE_RE.search(ln)
                if ct:
                    cnv_raw = ct.group(1).strip()
                    next_ln = block[idx + 1].strip() if idx + 1 < len(block) else ""
                    # handle split: "Copy" + "Number Gain/Loss"
                    if cnv_raw.lower() == "copy" and next_ln.lower().startswith("number"):
                        cnv_raw = f"{cnv_raw} {next_ln}"
                    elif cnv_raw.lower() == "copy number" and (next_ln.lower().startswith("gain") or next_ln.lower().startswith("loss")):
                        cnv_raw = f"{cnv_raw} {next_ln}"
                    cnv_raw_full = cnv_raw
                    if any(k in cnv_raw.lower() for k in ["gain", "dup", "ampl"]):
                        cnv_type = "Gain"
                    elif any(k in cnv_raw.lower() for k in ["loss", "del"]):
                        cnv_type = "Loss"
                    else:
                        cnv_type = cnv_raw
                else:
                    low = ln.lower()
                    if any(k in low for k in ["duplication", "amplification", "copy number gain"]):
                        cnv_type = "Gain"
                    elif any(k in low for k in ["deletion", "loss", "copy number loss"]):
                        cnv_type = "Loss"
            if gene_summary is None:
                gs = GENE_SUMMARY_RE.search(ln)
                if gs:
                    val = gs.group(1).strip()
                    if val:
                        gene_summary_parts.append(val)
                        in_gene_summary = True
                        continue
                    else:
                        in_gene_summary = True
                        continue
            if in_gene_summary:
                if GENE_LINE_RE.search(ln) or CLIN_THER_RE.search(ln) or CLIN_RELEVANCE_RE.search(ln):
                    in_gene_summary = False
                else:
                    gene_summary_parts.append(ln)
            if interpretation is None:
                im = INTERPRETATION_RE.search(ln)
                if im:
                    val = im.group(1).strip()
                    if val:
                        interpretation = val
            if clinical_relevance is None:
                cm = CLIN_RELEVANCE_RE.search(ln) or CLIN_THER_RE.search(ln)
                if cm:
                    val = cm.group(1).strip()
                    if val:
                        clin_rel_parts.append(val)
                        in_clin_relevance = True
                        continue
                    else:
                        in_clin_relevance = True
                        continue
            if in_clin_relevance:
                if GENE_LINE_RE.search(ln):
                    in_clin_relevance = False
                else:
                    clin_rel_parts.append(ln)
            for cid in CLINVAR_RE.findall(ln):
                clinvar_ids.add(cid.strip())
            for pmid in PMID_RE.findall(ln):
                pubmed_ids.add(pmid.strip())
            if PUBMED_REFS_RE.search(ln):
                in_pubmed_refs = True
                nums = re.findall(r"\b\d{7,8}\b", ln)
                pubmed_numbers.extend(nums)
                continue
            if in_pubmed_refs:
                if GENE_LINE_RE.search(ln) or "MedGenome" in ln or "Page " in ln:
                    in_pubmed_refs = False
                else:
                    nums = re.findall(r"\b\d{7,8}\b", ln)
                    pubmed_numbers.extend(nums)

        if gene_summary_parts:
            gene_summary = "\n".join(gene_summary_parts).strip()
        if clin_rel_parts:
            clinical_relevance = "\n".join(clin_rel_parts).strip()

        if not nucleotide_change:
            gmatch = GENOMIC_RE.search(joined)
            if gmatch:
                nucleotide_change = gmatch.group(0)

        if not genomic_match and nucleotide_change and nucleotide_change.lower().startswith("chr"):
            genomic_match = re.match(r".+", nucleotide_change)

        if variant_type == "CNV":
            if cnv_type is None or (isinstance(cnv_type, str) and cnv_type.lower() in {"copy", "copy number", "."}):
                joined_low = joined.lower()
                if any(k in joined_low for k in ["gain", "dup", "amplification", "copy number gain"]):
                    cnv_type = "Gain"
                elif any(k in joined_low for k in ["loss", "deletion", "copy number loss"]):
                    cnv_type = "Loss"

        if pubmed_numbers:
            for n in pubmed_numbers:
                pubmed_ids.add(n)

        results.append(
            {
                "gene": gene,
                "protein_change": protein_change,
                "cdna_change": f"c.{cdna_match.group(1)}" if cdna_match else "",
                "genomic_change": genomic_match.group(0) if genomic_match else "",
                "vaf": vaf,
                "interpretation": interpretation or "",
                "clinical_relevance": clinical_relevance or "",
                "exon": exon or "",
                "nucleotide_change": nucleotide_change or "",
                "transcript_id": transcript_id or "",
                "variant_depth": variant_depth or "",
                "population_maf": population_maf or "",
                "insilico_predictions": insilico or "",
                "gene_function": gene_function or "",
                "gene_summary": gene_summary or "",
                "variant_type": (cnv_raw_full or cnv_type or ".") if variant_type == "CNV" else variant_type,
                "clinvar_ids": ", ".join(sorted(clinvar_ids)) if clinvar_ids else "",
                "pubmed_ids": ", ".join(sorted(pubmed_ids)) if pubmed_ids else "",
                "evidence": block[0],
            }
        )

    # Deduplicate by gene + protein change, preferring the most complete record
    merged: Dict[tuple, Dict[str, str]] = {}
    for v in results:
        key = (v["gene"], v["protein_change"])
        if key not in merged:
            merged[key] = v
            continue
        cur = merged[key]
        # merge fields
        for field in [
            "cdna_change",
            "genomic_change",
            "vaf",
            "interpretation",
            "clinical_relevance",
            "exon",
            "nucleotide_change",
            "transcript_id",
            "variant_depth",
            "population_maf",
            "insilico_predictions",
            "gene_function",
            "gene_summary",
            "variant_type",
        ]:
            cur_val = cur.get(field)
            if (cur_val is None or cur_val == "" or cur_val == ".") and v.get(field):
                cur[field] = v[field]
        # merge ids
        for field in ["clinvar_ids", "pubmed_ids"]:
            if v.get(field):
                cur_set = set([x.strip() for x in cur.get(field, "").split(",") if x.strip()])
                new_set = set([x.strip() for x in v.get(field, "").split(",") if x.strip()])
                combined = sorted(cur_set | new_set)
                cur[field] = ", ".join(combined)
        # prefer longer evidence if current is short
        if len(v.get("evidence", "")) > len(cur.get("evidence", "")):
            cur["evidence"] = v.get("evidence", "")

    return list(merged.values())


@st.cache_data(show_spinner=False)
def extract_patient_info(text: str) -> Dict[str, Optional[str]]:
    info: Dict[str, Optional[str]] = {
        "patient_name": None,
        "age": None,
        "sex": None,
        "disease_name": None,
        "clinical_background_full": None,
        "date_received": None,
        "interpretation": None,
        "exon": None,
        "nucleotide_change": None,
        "clinvar_ids": None,
        "pubmed_ids": None,
        "clinical_relevance": None,
        "_sources": {},
    }

    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]

    # Patient name and age/sex are usually the first two lines
    if len(lines) >= 1:
        # Prefer the first line even if it fails strict regex
        if "AGE:" not in lines[0].upper():
            info["patient_name"] = lines[0]
            info["_sources"]["patient_name"] = lines[0]
    if len(lines) >= 2:
        m = AGE_SEX_RE.search(lines[1])
        if m:
            info["age"] = m.group(1)
            info["sex"] = m.group(2)
            info["_sources"]["age_sex"] = lines[1]

    # Age / Sex
    for ln in lines:
        m = AGE_SEX_RE.search(ln)
        if m:
            info["age"] = m.group(1)
            info["sex"] = m.group(2)
            info["_sources"]["age_sex"] = ln
            break

    # Disease name from Clinical Background (same line or next)
    for i, ln in enumerate(lines):
        m = CLIN_BG_RE.search(ln)
        if m:
            value = m.group(1).strip()
            if value:
                info["disease_name"] = value
                info["_sources"]["disease_name"] = ln
            else:
                for j in range(i + 1, min(i + 5, len(lines))):
                    if lines[j].strip():
                        info["disease_name"] = lines[j].strip()
                        info["_sources"]["disease_name"] = lines[j].strip()
                        break
            break

    # Full Clinical Background block
    for i, ln in enumerate(lines):
        if ln.strip().upper() == "CLINICAL BACKGROUND":
            block = []
            for j in range(i + 1, len(lines)):
                nxt = lines[j].strip()
                if not nxt:
                    continue
                upper = nxt.upper()
                if upper in {
                    "TEST RESULT SUMMARY",
                    "ACTIONABLE BIOMARKER DETAILS",
                    "ADDITIONAL BIOMARKERS DETECTED",
                    "GLOSSARY",
                    "REPORT DETAILS",
                }:
                    break
                if re.fullmatch(r"[A-Z0-9\\s/()\\-]+", upper) and len(upper.split()) <= 4:
                    break
                block.append(nxt)
            if block:
                info["clinical_background_full"] = " ".join(block).strip()
                info["_sources"]["clinical_background_full"] = " | ".join(block)
            break

    # Date Received (same line, nearby next/previous date-like line)
    for i, ln in enumerate(lines):
        m = DATE_RECEIVED_RE.search(ln)
        if m:
            info["date_received"] = m.group(1).strip()
            info["_sources"]["date_received"] = ln
            break
        if ln.lower().startswith("date received"):
            # look ahead
            for j in range(i + 1, min(i + 30, len(lines))):
                if "report date" in lines[j].lower():
                    continue
                m2 = DATE_LINE_RE.search(lines[j]) or DATE_LINE_RE_3.search(lines[j]) or DATE_LINE_RE_2.search(lines[j])
                if m2:
                    info["date_received"] = m2.group(0)
                    info["_sources"]["date_received"] = lines[j].strip()
                    break
            # look backward if still empty
            if not info["date_received"]:
                for j in range(max(0, i - 10), i):
                    if "report date" in lines[j].lower():
                        continue
                    m2 = DATE_LINE_RE.search(lines[j]) or DATE_LINE_RE_3.search(lines[j]) or DATE_LINE_RE_2.search(lines[j])
                    if m2:
                        info["date_received"] = m2.group(0)
                        info["_sources"]["date_received"] = lines[j].strip()
                        break
            if info["date_received"]:
                break

    # Interpretation (same line or next)
    for i, ln in enumerate(lines):
        m = INTERPRETATION_RE.search(ln)
        if m:
            value = m.group(1).strip()
            if value:
                info["interpretation"] = value
                info["_sources"]["interpretation"] = ln
            else:
                for j in range(i + 1, min(i + 5, len(lines))):
                    if lines[j].strip():
                        info["interpretation"] = lines[j].strip()
                        info["_sources"]["interpretation"] = lines[j].strip()
                        break
            break

    # Exon and Nucleotide change (genomic)
    for ln in lines:
        m = EXON_RE.search(ln)
        if m:
            info["exon"] = m.group(1)
            info["_sources"]["exon"] = ln
            break
    for ln in lines:
        m = NUC_CHANGE_RE.search(ln)
        if m:
            info["nucleotide_change"] = m.group(1).strip()
            info["_sources"]["nucleotide_change"] = ln
            break

    # Clinical relevance (normalized)
    for i, ln in enumerate(lines):
        m = CLIN_RELEVANCE_RE.search(ln)
        if m:
            value = m.group(1).strip()
            if value:
                info["clinical_relevance"] = value
                info["_sources"]["clinical_relevance"] = ln
            else:
                for j in range(i + 1, min(i + 5, len(lines))):
                    if lines[j].strip():
                        info["clinical_relevance"] = lines[j].strip()
                        info["_sources"]["clinical_relevance"] = lines[j].strip()
                        break
            break

    # ClinVar IDs and PubMed IDs
    clinvar_ids = set()
    pubmed_ids = set()
    for ln in lines:
        for cm in CLINVAR_RE.findall(ln):
            clinvar_ids.add(cm.strip())
            info["_sources"].setdefault("clinvar_ids", []).append(ln)
        for pm in PMID_RE.findall(ln):
            pubmed_ids.add(pm.strip())
            info["_sources"].setdefault("pubmed_ids", []).append(ln)

    info["clinvar_ids"] = ", ".join(sorted(clinvar_ids)) if clinvar_ids else None
    info["pubmed_ids"] = ", ".join(sorted(pubmed_ids)) if pubmed_ids else None

    return info


@st.cache_data(show_spinner=False)
def pubmed_search(query: str, retmax: int = 20) -> List[dict]:
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    esearch = f"{base}esearch.fcgi"
    esummary = f"{base}esummary.fcgi"

    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": retmax,
    }
    r = requests.get(esearch, params=params, timeout=20)
    r.raise_for_status()
    data = r.json()
    id_list = data.get("esearchresult", {}).get("idlist", [])
    if not id_list:
        return []

    s_params = {
        "db": "pubmed",
        "id": ",".join(id_list),
        "retmode": "json",
    }
    sr = requests.get(esummary, params=s_params, timeout=20)
    sr.raise_for_status()
    sdata = sr.json().get("result", {})
    results = []
    for pmid in id_list:
        item = sdata.get(pmid)
        if not item:
            continue
        results.append(
            {
                "pmid": pmid,
                "title": item.get("title", ""),
                "authors": ", ".join(a.get("name", "") for a in item.get("authors", []) if a.get("name")),
                "journal": item.get("fulljournalname", ""),
                "year": item.get("pubdate", ""),
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            }
        )
    return results


@st.cache_data(show_spinner=False)
def clinvar_search(query: str, retmax: int = 20) -> List[dict]:
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    esearch = f"{base}esearch.fcgi"
    esummary = f"{base}esummary.fcgi"

    params = {
        "db": "clinvar",
        "term": query,
        "retmode": "json",
        "retmax": retmax,
    }
    r = requests.get(esearch, params=params, timeout=20)
    r.raise_for_status()
    data = r.json()
    id_list = data.get("esearchresult", {}).get("idlist", [])
    if not id_list:
        return []

    s_params = {
        "db": "clinvar",
        "id": ",".join(id_list),
        "retmode": "json",
    }
    sr = requests.get(esummary, params=s_params, timeout=20)
    sr.raise_for_status()
    sdata = sr.json().get("result", {})
    results = []
    for cid in id_list:
        item = sdata.get(cid)
        if not item:
            continue
        title = item.get("title", "")
        results.append(
            {
                "clinvar_id": cid,
                "title": title,
                "url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{cid}/",
            }
        )
    return results


@st.cache_data(show_spinner=False)
def ollama_list_models(base_url: str) -> List[str]:
    try:
        r = requests.get(f"{base_url}/api/tags", timeout=10)
        r.raise_for_status()
        data = r.json()
        return [m.get("name", "") for m in data.get("models", []) if m.get("name")]
    except Exception:
        return []


@st.cache_data(show_spinner=False)
def ollama_health(base_url: str) -> Tuple[bool, str]:
    try:
        r = requests.get(f"{base_url}/api/tags", timeout=5)
        r.raise_for_status()
        return True, "OK"
    except Exception as e:
        return False, str(e)


def _select_quality_model(models: List[str]) -> Optional[str]:
    if not models:
        return None
    # Prefer largest parameter sizes for quality.
    preferred_sizes = ["70b", "65b", "34b", "33b", "30b", "13b", "8b", "7b", "3b"]
    lowered = [m.lower() for m in models]
    for size in preferred_sizes:
        for i, name in enumerate(lowered):
            if size in name:
                return models[i]
    return models[0]


def groq_chat(
    api_key: str,
    model: str,
    messages: List[Dict[str, str]],
    temperature: float = 0.1,
) -> str:
    if Groq is None:
        raise RuntimeError("groq package is not installed")
    client = Groq(api_key=api_key)
    resp = client.chat.completions.create(
        model=model,
        messages=messages,
        temperature=temperature,
    )
    return resp.choices[0].message.content or ""


@st.cache_data(show_spinner=False)
def ollama_chat(
    model: str,
    messages: List[Dict[str, str]],
    base_url: str,
    temperature: float = 0.1,
) -> str:
    payload = {
        "model": model,
        "messages": messages,
        "stream": False,
        "format": "json",
        "options": {"temperature": temperature},
    }
    chat_url = f"{base_url}/api/chat"
    r = requests.post(chat_url, json=payload, timeout=180)
    if r.status_code == 404:
        # Older Ollama versions may not support /api/chat. Fall back to /api/generate.
        system = next((m.get("content", "") for m in messages if m.get("role") == "system"), "")
        user = next((m.get("content", "") for m in messages if m.get("role") == "user"), "")
        prompt = (system + "\n\n" + user).strip()
        gen_payload = {
            "model": model,
            "prompt": prompt,
            "stream": False,
            "format": "json",
            "options": {"temperature": temperature},
        }
        gr = requests.post(f"{base_url}/api/generate", json=gen_payload, timeout=180)
        gr.raise_for_status()
        gdata = gr.json()
        return gdata.get("response", "")
    r.raise_for_status()
    data = r.json()
    return data.get("message", {}).get("content", "")


def _find_json_object(text: str) -> Optional[str]:
    start = text.find("{")
    if start == -1:
        return None
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[start : i + 1]
    return None


def _prepare_llm_input(text: str, max_chars: int = DEFAULT_MAX_CHARS) -> str:
    if len(text) <= max_chars:
        return text
    # Try to pull the actionable section if present.
    up = text.upper()
    marker = "ACTIONABLE BIOMARKER DETAILS"
    idx = up.find(marker)
    if idx != -1:
        sub = text[idx : idx + max_chars]
        return sub
    # Fallback: head + tail
    head = text[: int(max_chars * 0.75)]
    tail = text[-int(max_chars * 0.25) :]
    return head + "\n...\n" + tail


def _merge_patient_info(llm_info: Dict[str, str], regex_info: Dict[str, Optional[str]]) -> Dict[str, str]:
    out = {k: (llm_info.get(k) if llm_info.get(k) not in [None, "", "."] else None) for k in regex_info.keys()}
    for k, v in regex_info.items():
        if k == "_sources":
            continue
        if out.get(k) in [None, "", "."] and v not in [None, "", "."]:
            out[k] = v
    return out


def _normalize_variant(v: Dict[str, object]) -> Dict[str, str]:
    def _s(key: str) -> str:
        val = v.get(key)
        if val is None:
            return ""
        return str(val)

    return {
        "gene": _s("gene"),
        "protein_change": _s("protein_change"),
        "cdna_change": _s("cdna_change"),
        "genomic_change": _s("genomic_change"),
        "vaf": _s("vaf"),
        "interpretation": _s("interpretation"),
        "clinical_relevance": _s("clinical_relevance"),
        "exon": _s("exon"),
        "nucleotide_change": _s("nucleotide_change"),
        "transcript_id": _s("transcript_id"),
        "variant_depth": _s("variant_depth"),
        "population_maf": _s("population_maf"),
        "insilico_predictions": _s("insilico_predictions"),
        "gene_function": _s("gene_function"),
        "gene_summary": _s("gene_summary"),
        "variant_type": _s("variant_type"),
        "clinvar_ids": _s("clinvar_ids"),
        "pubmed_ids": _s("pubmed_ids"),
    }


def _merge_variants(llm_variants: List[Dict[str, str]], regex_variants: List[Dict[str, str]]) -> List[Dict[str, str]]:
    if not llm_variants:
        return regex_variants
    if not regex_variants:
        return llm_variants

    merged: Dict[tuple, Dict[str, str]] = {}
    for v in llm_variants:
        key = (v.get("gene", ""), v.get("protein_change", ""))
        merged[key] = dict(v)

    for rv in regex_variants:
        key = (rv.get("gene", ""), rv.get("protein_change", ""))
        if key not in merged:
            merged[key] = dict(rv)
            continue
        cur = merged[key]
        for field in [
            "cdna_change",
            "genomic_change",
            "vaf",
            "interpretation",
            "clinical_relevance",
            "exon",
            "nucleotide_change",
            "transcript_id",
            "variant_depth",
            "population_maf",
            "insilico_predictions",
            "gene_function",
            "gene_summary",
            "variant_type",
            "clinvar_ids",
            "pubmed_ids",
        ]:
            if (cur.get(field) in [None, "", "."]) and rv.get(field) not in [None, "", "."]:
                cur[field] = rv.get(field, "")
        merged[key] = cur
    return list(merged.values())


def _parse_llm_json(raw: str) -> Dict[str, object]:
    try:
        return json.loads(raw)
    except Exception:
        obj = _find_json_object(raw)
        if not obj:
            raise ValueError("LLM returned non-JSON output")
        return json.loads(obj)


@st.cache_data(show_spinner=False)
def llm_extract(
    text: str,
    model: str,
    base_url: str,
    max_chars: int = DEFAULT_MAX_CHARS,
    retries: int = 1,
    temperature: float = 0.1,
    backend: str = "ollama",
    groq_api_key: Optional[str] = None,
) -> Dict[str, object]:
    llm_text = _prepare_llm_input(text, max_chars=max_chars)
    system = (
        "You are a medical report extraction engine. "
        "Return ONLY valid JSON with the exact schema. "
        "Never include markdown, comments, or extra keys. "
        "If a field is unknown, use an empty string."
    )
    user = (
        "Extract patient and variant information from the report text. "
        "Schema:\n"
        "{\n"
        "  \"patient_info\": {\n"
        "    \"patient_name\": \"\",\n"
        "    \"age\": \"\",\n"
        "    \"sex\": \"\",\n"
        "    \"disease_name\": \"\",\n"
        "    \"clinical_background_full\": \"\",\n"
        "    \"date_received\": \"\",\n"
        "    \"interpretation\": \"\",\n"
        "    \"exon\": \"\",\n"
        "    \"nucleotide_change\": \"\",\n"
        "    \"clinvar_ids\": \"\",\n"
        "    \"pubmed_ids\": \"\",\n"
        "    \"clinical_relevance\": \"\"\n"
        "  },\n"
        "  \"variants\": [\n"
        "    {\n"
        "      \"gene\": \"\",\n"
        "      \"protein_change\": \"\",\n"
        "      \"cdna_change\": \"\",\n"
        "      \"genomic_change\": \"\",\n"
        "      \"vaf\": \"\",\n"
        "      \"interpretation\": \"\",\n"
        "      \"clinical_relevance\": \"\",\n"
        "      \"exon\": \"\",\n"
        "      \"nucleotide_change\": \"\",\n"
        "      \"transcript_id\": \"\",\n"
        "      \"variant_depth\": \"\",\n"
        "      \"population_maf\": \"\",\n"
        "      \"insilico_predictions\": \"\",\n"
        "      \"gene_function\": \"\",\n"
        "      \"gene_summary\": \"\",\n"
        "      \"variant_type\": \"\",\n"
        "      \"clinvar_ids\": \"\",\n"
        "      \"pubmed_ids\": \"\"\n"
        "    }\n"
        "  ]\n"
        "}\n\n"
        "Report text:\n" + llm_text
    )

    last_err = None
    for attempt in range(retries + 1):
        extra = ""
        if attempt > 0:
            extra = (
                "\n\nIMPORTANT: Your previous response was invalid. "
                "Return ONLY valid JSON with the exact schema and no extra text."
            )
        messages = [
            {"role": "system", "content": system},
            {"role": "user", "content": user + extra},
        ]
        if backend == "groq":
            if not groq_api_key:
                raise ValueError("GROQ_API_KEY is not set")
            content = groq_chat(
                api_key=groq_api_key,
                model=model,
                messages=messages,
                temperature=temperature,
            )
        else:
            content = ollama_chat(
                model=model,
                messages=messages,
                base_url=base_url,
                temperature=temperature,
            )
        raw = content.strip()
        try:
            return _parse_llm_json(raw)
        except Exception as e:
            last_err = e
            continue
    raise last_err if last_err else ValueError("LLM extraction failed")


st.set_page_config(page_title=APP_TITLE, layout="wide")

st.title(APP_TITLE)
st.caption("Authors: Mukesh Kumar, Surabhi Seth")
st.markdown(
    "Upload clinical report PDFs to extract gene/variant details using a local LLM, "
    "then optionally search PubMed and ClinVar for evidence."
)

with st.sidebar:
    st.header("Inputs")
    default_backend = "Groq (cloud)" if os.environ.get("STREAMLIT_CLOUD") else "Ollama (local)"
    backend = st.selectbox(
        "LLM backend",
        ["Ollama (local)", "Groq (cloud)"],
        index=0 if default_backend == "Ollama (local)" else 1,
    )
    use_groq = backend.startswith("Groq")
    if use_groq:
        groq_api_key = st.text_input("GROQ_API_KEY", value=os.environ.get("GROQ_API_KEY", ""), type="password")
        model = st.text_input("Groq model", value=DEFAULT_GROQ_MODEL)
        base_url = DEFAULT_OLLAMA_URL
        if Groq is None:
            st.error("Groq SDK not installed. Add `groq` to requirements.txt.")
        elif groq_api_key:
            st.success("Groq key set")
        else:
            st.error("Groq key missing")
    else:
        base_url = st.text_input("Ollama base URL", value=DEFAULT_OLLAMA_URL)
        st.caption("If Ollama is running elsewhere, set OLLAMA_HOST or change this URL.")
        ollama_ok, _ = ollama_health(base_url)
        if ollama_ok:
            st.success("Ollama: reachable")
        else:
            st.error("Ollama: not reachable")
            st.caption("Start it with `ollama serve`, then refresh this page.")
        models = ollama_list_models(base_url)
        if models:
            quality_model = _select_quality_model(models)
            default_index = models.index(quality_model) if quality_model in models else 0
            model = st.selectbox("Model", models, index=default_index)
        else:
            model = st.text_input("Model", value=DEFAULT_MODEL)
        groq_api_key = None
    if use_groq and not groq_api_key:
        st.info("Add your own Groq API key above to enable cloud LLM extraction.")
    use_llm_default = bool(groq_api_key) if use_groq else True
    use_llm = st.checkbox("Use LLM extraction", value=use_llm_default)
    max_chars = st.slider("LLM input max chars", 3000, 30000, DEFAULT_MAX_CHARS, 1000)
    llm_retries = st.slider("LLM JSON retries", 0, 2, 1, 1)
    temperature = st.slider("LLM temperature", 0.0, 0.5, 0.0, 0.1)
    enable_external = st.checkbox(
        "Enable external lookups (PubMed/ClinVar/ProteinPaint/cBioPortal)",
        value=True,
    )
    pubmed_retmax = st.slider("Max PubMed results", 5, 50, 20, 5)
    st.markdown("---")
    st.subheader("ProteinPaint")
    pp_host = st.text_input("ProteinPaint host", value="https://proteinpaint.stjude.org")
    pp_genome = st.text_input("Genome", value="hg38")
    pp_dataset = st.text_input("Dataset", value="pediatric")
    st.markdown("---")
    st.caption("Tip: use a specific variant query like `PTPN11 p.Gly503Ala` for precision.")

st.info(
    "LLM inference can run locally via Ollama. External lookups "
    "(PubMed/ClinVar/ProteinPaint/cBioPortal) only use gene/variant queries "
    "and can be disabled in the sidebar."
)
if use_groq:
    st.warning("Groq backend sends report text to a cloud API.")
if not enable_external:
    st.warning("External lookups are disabled. PubMed/ClinVar/ProteinPaint/cBioPortal sections are hidden.")

tab_single, tab_multi = st.tabs(["Single File", "Multiple Files"])


def _extract_with_fallback(text: str) -> Tuple[Dict[str, str], List[Dict[str, str]]]:
    regex_info = extract_patient_info(text)
    regex_variants = _extract_variant_blocks(text)
    if not regex_variants:
        hits = extract_variants(text)
        regex_variants = [
            {
                "gene": h.gene,
                "protein_change": h.protein_change,
                "cdna_change": h.cdna_change or "",
                "genomic_change": h.genomic_change or "",
                "vaf": h.vaf or "",
                "interpretation": regex_info.get("interpretation") or "",
                "clinical_relevance": regex_info.get("clinical_relevance") or "",
                "exon": regex_info.get("exon") or "",
                "nucleotide_change": regex_info.get("nucleotide_change") or "",
                "transcript_id": "",
                "variant_depth": "",
                "population_maf": "",
                "insilico_predictions": "",
                "gene_function": "",
                "gene_summary": "",
                "variant_type": "SNV/Indel",
                "clinvar_ids": regex_info.get("clinvar_ids") or "",
                "pubmed_ids": regex_info.get("pubmed_ids") or "",
            }
            for h in hits
        ]

    if not use_llm:
        return regex_info, regex_variants

    if use_groq and not groq_api_key:
        return regex_info, regex_variants
    try:
        llm_data = llm_extract(
            text,
            model=model,
            base_url=base_url,
            max_chars=max_chars,
            retries=llm_retries,
            temperature=temperature,
            backend="groq" if use_groq else "ollama",
            groq_api_key=groq_api_key,
        )
    except Exception as e:
        st.warning(f"LLM extraction failed, using regex fallback. Error: {e}")
        return regex_info, regex_variants

    llm_info = llm_data.get("patient_info", {}) if isinstance(llm_data, dict) else {}
    llm_variants_raw = llm_data.get("variants", []) if isinstance(llm_data, dict) else []
    llm_variants = [_normalize_variant(v) for v in llm_variants_raw if isinstance(v, dict)]

    merged_info = _merge_patient_info(llm_info, regex_info)
    merged_variants = _merge_variants(llm_variants, regex_variants)
    return merged_info, merged_variants


with tab_single:
    uploaded = st.file_uploader("Upload a PDF", type=["pdf"], key="single_upload")

    pdf_bytes = None
    pdf_name = None

    if uploaded is not None:
        pdf_bytes = uploaded.read()
        pdf_name = uploaded.name

    if pdf_bytes:
        with st.spinner("Extracting text from PDF..."):
            text = pdf_to_text_bytes(pdf_bytes)

        if not text.strip():
            st.warning("No text could be extracted. If this is a scanned PDF, OCR is required.")
        else:
            with st.spinner("Extracting structured data..."):
                info, variant_blocks = _extract_with_fallback(text)

            st.subheader("Patient Short Summary")
            summary_row = {
                "Patient Name": info.get("patient_name") or "",
                "Age": info.get("age") or "",
                "Sex": info.get("sex") or "",
                "Clinical Background (Full)": info.get("clinical_background_full") or "",
                "Date Received": info.get("date_received") or "",
            }
            st.dataframe([summary_row], use_container_width=True)

            if enable_external and info.get("pubmed_ids"):
                st.subheader("PubMed Links (from PDF)")
                for pmid in str(info["pubmed_ids"]).split(","):
                    pmid = pmid.strip()
                    if not pmid:
                        continue
                    st.markdown(f"PMID: {pmid}  \\\nhttps://pubmed.ncbi.nlm.nih.gov/{pmid}/")

            st.subheader("Extracted Variants")
            if not variant_blocks:
                st.info("No variants detected with the current parser.")
            else:
                rows = [
                    {
                        "Patient Name": info.get("patient_name") or "",
                        "Gene": v.get("gene", ""),
                        "Protein change": v.get("protein_change", ""),
                        "cDNA": v.get("cdna_change", ""),
                        "Genomic": v.get("genomic_change", ""),
                        "VAF": v.get("vaf", ""),
                        "Variant Allele Depth/Total depth": v.get("variant_depth", ""),
                        "Population MAF": v.get("population_maf", ""),
                        "In-silico Predictions": v.get("insilico_predictions", ""),
                        "Gene Function": v.get("gene_function", ""),
                        "Gene Summary": v.get("gene_summary", ""),
                        "Variant Type": v.get("variant_type", ""),
                        "Interpretation": v.get("interpretation", ""),
                        "Clinical Relevance": v.get("clinical_relevance", ""),
                        "Exon": v.get("exon", ""),
                        "Nucleotide Change (Genomic)": v.get("nucleotide_change", ""),
                        "Transcript ID": v.get("transcript_id", ""),
                        "ClinVar IDs (from PDF)": v.get("clinvar_ids", ""),
                        "PubMed IDs (from PDF)": v.get("pubmed_ids", ""),
                    }
                    for v in variant_blocks
                ]
                rows = [_fill_missing(r) for r in rows]
                display_rows = []
                for r in rows:
                    r2 = dict(r)
                    r2["Gene Summary"] = _truncate(r.get("Gene Summary", ""))
                    r2["Clinical Relevance"] = _truncate(r.get("Clinical Relevance", ""))
                    display_rows.append(r2)
                st.dataframe(display_rows, use_container_width=True)
                with st.expander("Read more (Gene Summary / Clinical Relevance)"):
                    patient_ids = sorted({r.get("Patient Name", "") for r in rows if r.get("Patient Name")})
                    selected = st.selectbox("Select Patient", ["All"] + patient_ids, key="readmore_select_single")
                    for r in rows:
                        if selected != "All" and r.get("Patient Name") != selected:
                            continue
                        gs = r.get("Gene Summary", "")
                        cr = r.get("Clinical Relevance", "")
                        if not gs and not cr:
                            continue
                        st.markdown(
                            f"**{r.get('Patient Name','')} | {r.get('Gene','')} {r.get('Protein change','')}**"
                        )
                        if gs:
                            st.markdown("Gene Summary:")
                            st.write(gs)
                        if cr:
                            st.markdown("Clinical and Therapeutic Relevance:")
                            st.write(cr)
                pmid_items = []
                for r in rows:
                    pmids = str(r.get("PubMed IDs (from PDF)", "")).strip()
                    if pmids:
                        pmid_items.append(
                            {
                                "Patient Name": r.get("Patient Name", ""),
                                "Gene": r.get("Gene", ""),
                                "Protein change": r.get("Protein change", ""),
                                "PMIDs": pmids,
                            }
                        )
                if enable_external and pmid_items:
                    with st.expander("PMID Links"):
                        for item in pmid_items:
                            pmid_links = []
                            for p in item["PMIDs"].split(","):
                                p = p.strip()
                                if not p:
                                    continue
                                pmid_links.append(
                                    f'<a href="https://pubmed.ncbi.nlm.nih.gov/{p}/" target="_blank" rel="noreferrer">{p}</a>'
                                )
                            st.markdown(
                                f'**{item["Patient Name"]} | {item["Gene"]} {item["Protein change"]}**: '
                                + ", ".join(pmid_links),
                                unsafe_allow_html=True,
                            )
                genes = sorted({r.get("Gene") for r in rows if r.get("Gene")})
                if enable_external and genes:
                    st.subheader("ProteinPaint Lollipop")
                    gene_sel = st.selectbox("Select gene", genes, key="pp_gene_single")
                    _render_proteinpaint(gene_sel, pp_genome, pp_dataset, pp_host)

                first_gene = rows[0].get("Gene", "")
                first_prot = rows[0].get("Protein change", "")
                default_query = " ".join(x for x in [first_gene, first_prot] if x)
                first_cdna = rows[0].get("cDNA", "")
                first_genomic = rows[0].get("Genomic", "")
                if first_cdna:
                    clinvar_default_query = f"{first_gene} {first_cdna}"
                elif first_genomic:
                    clinvar_default_query = f"{first_gene} {first_genomic}"
                else:
                    clinvar_default_query = default_query

                download_payload = {
                    "patient_info": info,
                    "variants": rows,
                }
                st.download_button(
                    label="Download JSON",
                    data=json.dumps(download_payload, indent=2),
                    file_name="report_extraction.json",
                    mime="application/json",
                )

                if enable_external:
                    st.subheader("PubMed Search")
                    manual_query = st.text_input("PubMed query", value=default_query, key="pubmed_query_single")
                    if st.button("Search PubMed", key="pubmed_btn_single"):
                        with st.spinner("Querying PubMed..."):
                            try:
                                results = pubmed_search(manual_query, retmax=pubmed_retmax)
                            except Exception as e:
                                st.error(f"PubMed query failed: {e}")
                                results = []

                        if not results:
                            st.warning("No PubMed results found.")
                        else:
                            for r in results:
                                st.markdown(
                                    f"**{r['title']}**  \\\n                                {r['authors']}  \\\n                                {r['journal']} ({r['year']})  \\\n                                PMID: {r['pmid']}  \\\n                                {r['url']}"
                                )

                    st.subheader("ClinVar Search")
                    clinvar_query = st.text_input(
                        "ClinVar query", value=clinvar_default_query, key="clinvar_query_single"
                    )
                    if st.button("Search ClinVar", key="clinvar_btn_single"):
                        with st.spinner("Querying ClinVar..."):
                            try:
                                cresults = clinvar_search(clinvar_query, retmax=pubmed_retmax)
                            except Exception as e:
                                st.error(f"ClinVar query failed: {e}")
                                cresults = []

                        if not cresults:
                            st.warning("No ClinVar results found.")
                        else:
                            for r in cresults:
                                st.markdown(f"**{r['title']}**  \\\nClinVar ID: {r['clinvar_id']}  \\\n{r['url']}")
    else:
        st.info("Provide a PDF to begin.")

with tab_multi:
    uploaded_files = st.file_uploader(
        "Upload multiple PDFs", type=["pdf"], accept_multiple_files=True, key="multi_upload"
    )
    if uploaded_files:
        all_summary_rows = []
        all_variant_rows = []
        for up in uploaded_files:
            with st.spinner(f"Extracting {up.name}..."):
                text = pdf_to_text_bytes(up.read())
            if not text.strip():
                st.warning(f"No text extracted from {up.name}")
                continue
            info, variant_blocks = _extract_with_fallback(text)
            patient_name = info.get("patient_name") or up.name
            all_summary_rows.append(
                {
                    "Patient Name": patient_name,
                    "Age": info.get("age") or "",
                    "Sex": info.get("sex") or "",
                    "Clinical Background (Full)": info.get("clinical_background_full") or "",
                    "Date Received": info.get("date_received") or "",
                }
            )

            for v in variant_blocks:
                all_variant_rows.append(
                    {
                        "Patient Name": patient_name,
                        "Gene": v.get("gene", ""),
                        "Protein change": v.get("protein_change", ""),
                        "cDNA": v.get("cdna_change", ""),
                        "Genomic": v.get("genomic_change", ""),
                        "VAF": v.get("vaf", ""),
                        "Variant Allele Depth/Total depth": v.get("variant_depth", ""),
                        "Population MAF": v.get("population_maf", ""),
                        "In-silico Predictions": v.get("insilico_predictions", ""),
                        "Gene Function": v.get("gene_function", ""),
                        "Gene Summary": v.get("gene_summary", ""),
                        "Variant Type": v.get("variant_type", ""),
                        "Interpretation": v.get("interpretation", ""),
                        "Clinical Relevance": v.get("clinical_relevance", ""),
                        "Exon": v.get("exon", ""),
                        "Nucleotide Change (Genomic)": v.get("nucleotide_change", ""),
                        "Transcript ID": v.get("transcript_id", ""),
                        "ClinVar IDs (from PDF)": v.get("clinvar_ids", ""),
                        "PubMed IDs (from PDF)": v.get("pubmed_ids", ""),
                    }
                )

        st.subheader("Patient Short Summary")
        if all_summary_rows:
            st.dataframe(all_summary_rows, use_container_width=True)
            st.download_button(
                label="Download Summary CSV",
                data=_rows_to_csv(all_summary_rows),
                file_name="patient_summary.csv",
                mime="text/csv",
            )
        else:
            st.info("No summaries to display.")

        st.subheader("Extracted Variants")
        if all_variant_rows:
            all_variant_rows = [_fill_missing(r) for r in all_variant_rows]
            display_rows = []
            for r in all_variant_rows:
                r2 = dict(r)
                r2["Gene Summary"] = _truncate(r.get("Gene Summary", ""))
                r2["Clinical Relevance"] = _truncate(r.get("Clinical Relevance", ""))
                display_rows.append(r2)
            st.dataframe(display_rows, use_container_width=True)
            with st.expander("Read more (Gene Summary / Clinical Relevance)"):
                patient_ids = sorted({r.get("Patient Name", "") for r in all_variant_rows if r.get("Patient Name")})
                selected = st.selectbox("Select Patient", ["All"] + patient_ids, key="readmore_select_multi")
                for r in all_variant_rows:
                    if selected != "All" and r.get("Patient Name") != selected:
                        continue
                    gs = r.get("Gene Summary", "")
                    cr = r.get("Clinical Relevance", "")
                    if not gs and not cr:
                        continue
                    st.markdown(
                        f"**{r.get('Patient Name','')} | {r.get('Gene','')} {r.get('Protein change','')}**"
                    )
                    if gs:
                        st.markdown("Gene Summary:")
                        st.write(gs)
                    if cr:
                        st.markdown("Clinical and Therapeutic Relevance:")
                        st.write(cr)
            pmid_items = []
            for r in all_variant_rows:
                pmids = str(r.get("PubMed IDs (from PDF)", "")).strip()
                if pmids:
                    pmid_items.append(
                        {
                            "Patient Name": r.get("Patient Name", ""),
                            "Gene": r.get("Gene", ""),
                            "Protein change": r.get("Protein change", ""),
                            "PMIDs": pmids,
                        }
                    )
            if enable_external and pmid_items:
                with st.expander("PMID Links"):
                    for item in pmid_items:
                        pmid_links = []
                        for p in item["PMIDs"].split(","):
                            p = p.strip()
                            if not p:
                                continue
                            pmid_links.append(
                                f'<a href="https://pubmed.ncbi.nlm.nih.gov/{p}/" target="_blank" rel="noreferrer">{p}</a>'
                            )
                        st.markdown(
                            f'**{item["Patient Name"]} | {item["Gene"]} {item["Protein change"]}**: '
                            + ", ".join(pmid_links),
                            unsafe_allow_html=True,
                        )
            st.download_button(
                label="Download Variants CSV",
                data=_rows_to_csv(all_variant_rows),
                file_name="variants.csv",
                mime="text/csv",
            )
            genes = sorted({r.get("Gene") for r in all_variant_rows if r.get("Gene")})
            if enable_external and genes:
                st.subheader("ProteinPaint Lollipop")
                gene_sel = st.selectbox("Select gene", genes, key="pp_gene_multi")
                _render_proteinpaint(gene_sel, pp_genome, pp_dataset, pp_host)

                st.subheader("cBioPortal")
                cbio_link = f"{CBIO_BASE}/mutation_mapper?gene={gene_sel}"
                st.markdown(
                    "**Steps**\n"
                    "1. Copy the Input Data below\n"
                    "2. Click the button to open cBioPortal\n"
                    "3. Paste in the input and visualize"
                )
                st.markdown(
                    f"""
                    <a href="{cbio_link}" target="_blank" rel="noreferrer"
                       style="display:inline-block;padding:8px 14px;background:#2563eb;color:#fff;
                              border-radius:6px;text-decoration:none;font-weight:600;">
                      Open cBioPortal Mutation Mapper
                    </a>
                    """,
                    unsafe_allow_html=True,
                )
                cbio_rows = _build_cbioportal_input(all_variant_rows)
                if cbio_rows:
                    tsv = _rows_to_tsv(cbio_rows)
                    st.markdown("**Input Data (copy this):**")
                    components.html(
                        f"""
                        <div style="display:flex; justify-content:flex-start; margin-bottom:6px;">
                          <button onclick="navigator.clipboard.writeText(document.getElementById('cbio_text_multi').value)">
                            Copy
                          </button>
                        </div>
                        <textarea id="cbio_text_multi" style="width:100%;height:200px;">{html.escape(tsv)}</textarea>
                        """,
                        height=260,
                    )
                else:
                    st.info("No variants available for cBioPortal input.")

            st.subheader("VAF Heatmap")
            heatmap_data = _build_heatmap_data(all_variant_rows)
            if heatmap_data:
                heatmap_height = st.slider("Heatmap height", 200, 1000, 400, 50, key="heatmap_height")
                heatmap_width = st.slider("Heatmap width", 300, 2000, 900, 50, key="heatmap_width")
                st.caption("Interactive heatmap. Hover to see values. Scroll to pan.")
                spec = {
                    "mark": "rect",
                    "selection": {"grid": {"type": "interval", "bind": "scales"}},
                    "encoding": {
                        "x": {"field": "GeneVariant", "type": "nominal", "axis": {"labelAngle": -45}},
                        "y": {"field": "Patient", "type": "nominal"},
                        "color": {"field": "VAF", "type": "quantitative", "scale": {"scheme": "viridis"}},
                        "tooltip": [
                            {"field": "Patient", "type": "nominal"},
                            {"field": "GeneVariant", "type": "nominal"},
                            {"field": "VAF", "type": "quantitative"},
                        ],
                    },
                    "width": heatmap_width,
                    "height": heatmap_height,
                }
                st.vega_lite_chart(heatmap_data, spec, use_container_width=True)
                vega_spec = dict(spec)
                vega_spec["data"] = {"values": heatmap_data}
                html_payload = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8"/>
  <script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
</head>
<body>
  <div id="vis"></div>
  <script>
    const spec = {json.dumps(vega_spec)};
    vegaEmbed('#vis', spec, {{actions: false}});
  </script>
</body>
</html>"""
                st.download_button(
                    label="Download VAF Heatmap (HTML)",
                    data=html_payload,
                    file_name="vaf_heatmap.html",
                    mime="text/html",
                )
                if vlc is not None:
                    try:
                        png_bytes = vlc.vegalite_to_png(vega_spec)
                        st.download_button(
                            label="Download VAF Heatmap (PNG)",
                            data=png_bytes,
                            file_name="vaf_heatmap.png",
                            mime="image/png",
                        )
                    except Exception:
                        st.info("PNG export unavailable.")
                    try:
                        svg_text = vlc.vegalite_to_svg(vega_spec)
                        st.download_button(
                            label="Download VAF Heatmap (SVG)",
                            data=svg_text,
                            file_name="vaf_heatmap.svg",
                            mime="image/svg+xml",
                        )
                    except Exception:
                        st.info("SVG export unavailable.")
                    if cairosvg is not None:
                        try:
                            pdf_bytes = cairosvg.svg2pdf(bytestring=svg_text.encode("utf-8"))
                            st.download_button(
                                label="Download VAF Heatmap (PDF)",
                                data=pdf_bytes,
                                file_name="vaf_heatmap.pdf",
                                mime="application/pdf",
                            )
                        except Exception:
                            st.info("PDF export unavailable.")
                else:
                    st.info("PNG/PDF export requires `vl-convert-python` (and `cairosvg` for PDF).")
            else:
                st.info("Heatmap unavailable (no VAF values).")
        else:
            st.info("No variants to display.")
    else:
        st.info("Upload one or more PDFs to begin.")

st.markdown("---")
st.caption("This tool is made for report interpretation.")
