import io
import json
import re
import subprocess
from dataclasses import dataclass
from typing import List, Optional, Tuple, Dict

import requests
import streamlit as st
import streamlit.components.v1 as components
 

try:
    import pdfplumber  # optional fallback
except Exception:  # pragma: no cover
    pdfplumber = None
try:
    import pytesseract  # optional OCR
except Exception:  # pragma: no cover
    pytesseract = None


APP_TITLE = "Clinical Report Variant Finder + PubMed"

def _rows_to_csv(rows: List[Dict[str, str]]) -> str:
    if not rows:
        return ""
    fieldnames = sorted({k for r in rows for k in r.keys()})
    output = io.StringIO()
    output.write(",".join(fieldnames) + "\\n")
    for r in rows:
        output.write(",".join(str(r.get(k, "")).replace(",", " ") for k in fieldnames) + "\\n")
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


def _render_proteinpaint(gene: str, genome: str, dataset: str, host: str) -> None:
    if not gene:
        return
    host = host.rstrip("/")
    link = f"{host}/?gene={gene}&genome={genome}"
    if dataset:
        link += f"&dataset={dataset}"
    st.markdown(
        f'ProteinPaint link: <a href="{link}" target="_blank" rel="noreferrer">{link}</a>',
        unsafe_allow_html=True,
    )


 


@dataclass
class VariantHit:
    gene: str
    protein_change: str
    cdna_change: Optional[str]
    genomic_change: Optional[str]
    vaf: Optional[str]
    raw_line: str


GENE_LINE_RE = re.compile(r"^\s*([A-Z0-9]{2,})\b.*?p\.([A-Za-z]{3}\d+[A-Za-z]{3}(?:fsTer\d+)?)")
PROTEIN_RE = re.compile(r"\bp\.([A-Za-z]{3}\d+[A-Za-z]{3}(?:fsTer\d+)?)\b", re.IGNORECASE)
CDNA_RE = re.compile(r"c\.\s*([0-9_]+(?:[ACGT]>[ACGT]|ins[ACGT]+|del[ACGT]+|dup[ACGT]+))", re.IGNORECASE)
GENOMIC_RE = re.compile(r"chr\w+:g\.[0-9_]+(?:[ACGT]>[ACGT]|ins[ACGT]+|del[ACGT]+|dup[ACGT]+)", re.IGNORECASE)
VAF_RE = re.compile(r"(?:Variant Allele Frequency|VAF)\s*[-:]*\s*([0-9.]+%?)", re.IGNORECASE)

NAME_RE = re.compile(r"^([A-Za-z][A-Za-z\s.'-]{1,})$")
AGE_SEX_RE = re.compile(r"AGE:\s*([0-9]+)\s*years\s*\|\s*Gender:\s*([A-Za-z]+)", re.IGNORECASE)
DATE_RECEIVED_RE = re.compile(r"Date Received:\s*([0-9]{1,2}\s+[A-Za-z]+\s+[0-9]{4})", re.IGNORECASE)
DATE_LINE_RE = re.compile(r"\b[0-9]{1,2}\s+[A-Za-z]+\s+[0-9]{4}\b")
CLIN_BG_RE = re.compile(r"CLINICAL BACKGROUND\s*:?\s*(.*)", re.IGNORECASE)
INTERPRETATION_RE = re.compile(r"Interpretation\s*(.+)", re.IGNORECASE)
EXON_RE = re.compile(r"Exon:\s*([0-9]+)", re.IGNORECASE)
NUC_CHANGE_RE = re.compile(r"Nucleotide change:\s*(chr\w+:g\.[0-9_]+[A-Za-z0-9>_]+)", re.IGNORECASE)
TRANSCRIPT_RE = re.compile(r"Transcript ID:\s*([A-Za-z0-9_.-]+)", re.IGNORECASE)
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
            if "ADDITIONAL BIOMARKERS DETECTED" in lines[j].upper() or "GLOSSARY" in lines[j].upper():
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
        if not m:
            continue
        gene = m.group(1)
        protein_change = f"p.{m.group(2)}"
        cdna_match = CDNA_RE.search(joined)
        genomic_match = GENOMIC_RE.search(joined)
        vaf_match = VAF_RE.search(joined)

        exon = None
        nucleotide_change = None
        transcript_id = None
        interpretation = None
        clinical_relevance = None
        clinvar_ids = set()
        pubmed_ids = set()

        in_clin_relevance = False
        clin_rel_parts: List[str] = []
        in_pubmed_refs = False
        pubmed_numbers: List[str] = []

        for ln in block:
            if _skip_noise(ln):
                continue

            if exon is None:
                em = EXON_RE.search(ln)
                if em:
                    exon = em.group(1)
            if nucleotide_change is None:
                nm = NUC_CHANGE_RE.search(ln)
                if nm:
                    nucleotide_change = nm.group(1).strip()
            if transcript_id is None:
                tm = TRANSCRIPT_RE.search(ln)
                if tm:
                    transcript_id = tm.group(1).strip()
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
                        clinical_relevance = val
                        in_clin_relevance = False
                    else:
                        in_clin_relevance = True
                        continue
            if in_clin_relevance:
                if PUBMED_REFS_RE.search(ln) or GENE_LINE_RE.search(ln) or "Gene Summary" in ln:
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

        if clinical_relevance is None and clin_rel_parts:
            clinical_relevance = " ".join(clin_rel_parts).strip()

        if not nucleotide_change:
            gmatch = GENOMIC_RE.search(joined)
            if gmatch:
                nucleotide_change = gmatch.group(0)

        if not genomic_match and nucleotide_change and nucleotide_change.lower().startswith("chr"):
            genomic_match = re.match(r".+", nucleotide_change)

        if pubmed_numbers:
            for n in pubmed_numbers:
                pubmed_ids.add(n)

        results.append(
            {
                "gene": gene,
                "protein_change": protein_change,
                "cdna_change": f"c.{cdna_match.group(1)}" if cdna_match else "",
                "genomic_change": genomic_match.group(0) if genomic_match else "",
                "vaf": vaf_match.group(1) if vaf_match else "",
                "interpretation": interpretation or "",
                "clinical_relevance": clinical_relevance or "",
                "exon": exon or "",
                "nucleotide_change": nucleotide_change or "",
                "transcript_id": transcript_id or "",
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
        ]:
            if not cur.get(field) and v.get(field):
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

    # Date Received (same line or next date-like line)
    for i, ln in enumerate(lines):
        m = DATE_RECEIVED_RE.search(ln)
        if m:
            info["date_received"] = m.group(1).strip()
            info["_sources"]["date_received"] = ln
            break
        if ln.lower().startswith("date received"):
            for j in range(i + 1, min(i + 30, len(lines))):
                if "report date" in lines[j].lower():
                    continue
                m2 = DATE_LINE_RE.search(lines[j])
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


st.set_page_config(page_title=APP_TITLE, layout="wide")

st.title(APP_TITLE)
st.caption("Authors: Mukesh Kumar, Surabhi Seth")

st.markdown(
    "Upload clinical report PDFs to extract gene/variant details, then search PubMed and ClinVar for evidence."
)

with st.sidebar:
    st.header("Inputs")
    pubmed_retmax = st.slider("Max PubMed results", 5, 50, 20, 5)
    st.markdown("---")
    st.subheader("ProteinPaint")
    pp_host = st.text_input("ProteinPaint host", value="https://proteinpaint.stjude.org")
    pp_genome = st.text_input("Genome", value="hg38")
    pp_dataset = st.text_input("Dataset", value="pediatric")
    st.markdown("---")
    st.caption("Tip: use a specific variant query like `PTPN11 p.Gly503Ala` for precision.")

tab_single, tab_multi = st.tabs(["Single File", "Multiple Files"])

with tab_single:
    pdf_source = st.radio("PDF Source", ["Upload", "Local Path"], horizontal=True, key="single_source")
    default_path = "/home/adminsb/disk2/report/TRN4475861_9068316_1257283_clinical_report_1_1745247720810608.pdf"
    local_path = st.text_input("Local PDF path", value=default_path, key="single_path")

    pdf_bytes = None
    pdf_name = None

    if pdf_source == "Upload":
        uploaded = st.file_uploader("Upload a PDF", type=["pdf"], key="single_upload")
        if uploaded is not None:
            pdf_bytes = uploaded.read()
            pdf_name = uploaded.name
    else:
        if local_path:
            try:
                with open(local_path, "rb") as f:
                    pdf_bytes = f.read()
                    pdf_name = local_path.split("/")[-1]
            except Exception as e:
                st.error(f"Failed to read PDF: {e}")

    if pdf_bytes:
        with st.spinner("Extracting text from PDF..."):
            text = pdf_to_text_bytes(pdf_bytes)

        if not text.strip():
            st.warning("No text could be extracted. If this is a scanned PDF, OCR is required.")
        else:
            info = extract_patient_info(text)

            st.subheader("Patient / Report Summary")
            summary_rows = [
                {"Field": "Patient Name", "Info": info.get("patient_name") or ""},
                {"Field": "Age", "Info": info.get("age") or ""},
                {"Field": "Sex", "Info": info.get("sex") or ""},
                {"Field": "Disease Name (Clinical Background)", "Info": info.get("disease_name") or ""},
                {"Field": "Date Received", "Info": info.get("date_received") or ""},
                {"Field": "Interpretation", "Info": info.get("interpretation") or ""},
                {"Field": "Clinical Relevance", "Info": info.get("clinical_relevance") or ""},
                {"Field": "Exon", "Info": info.get("exon") or ""},
                {"Field": "Nucleotide Change (Genomic)", "Info": info.get("nucleotide_change") or ""},
                {"Field": "ClinVar IDs (from PDF)", "Info": info.get("clinvar_ids") or ""},
                {"Field": "PubMed IDs (from PDF)", "Info": info.get("pubmed_ids") or ""},
            ]
            st.table(summary_rows)

            if info.get("pubmed_ids"):
                st.subheader("PubMed Links (from PDF)")
                for pmid in info["pubmed_ids"].split(", "):
                    st.markdown(f"PMID: {pmid}  \\\nhttps://pubmed.ncbi.nlm.nih.gov/{pmid}/")

            st.subheader("Extracted Variants")
            hits = extract_variants(text)
            variant_blocks = _extract_variant_blocks(text)
            if not hits:
                st.info("No variants detected with the current parser.")
            else:
                if variant_blocks:
                    rows = [
                        {
                            "Patient Name": info.get("patient_name") or "",
                            "Gene": v["gene"],
                            "Protein change": v["protein_change"],
                            "cDNA": v["cdna_change"],
                            "Genomic": v["genomic_change"],
                            "VAF": v["vaf"],
                            "Interpretation": v["interpretation"],
                            "Clinical Relevance": v["clinical_relevance"],
                            "Exon": v["exon"],
                            "Nucleotide Change (Genomic)": v["nucleotide_change"],
                            "Transcript ID": v["transcript_id"],
                            "ClinVar IDs (from PDF)": v["clinvar_ids"],
                            "PubMed IDs (from PDF)": v["pubmed_ids"],
                            "Evidence": v["evidence"],
                        }
                        for v in variant_blocks
                    ]
                else:
                    rows = [
                        {
                            "Patient Name": info.get("patient_name") or "",
                            "Gene": h.gene,
                            "Protein change": h.protein_change,
                            "cDNA": h.cdna_change or "",
                            "Genomic": h.genomic_change or "",
                            "VAF": h.vaf or "",
                            "Interpretation": info.get("interpretation") or "",
                            "Clinical Relevance": info.get("clinical_relevance") or "",
                            "Exon": info.get("exon") or "",
                            "Nucleotide Change (Genomic)": info.get("nucleotide_change") or "",
                            "Transcript ID": "",
                            "ClinVar IDs (from PDF)": info.get("clinvar_ids") or "",
                            "PubMed IDs (from PDF)": info.get("pubmed_ids") or "",
                            "Evidence": h.raw_line,
                        }
                        for h in hits
                    ]
                st.dataframe(rows, use_container_width=True)
                genes = sorted({r.get("Gene") for r in rows if r.get("Gene")})
                if genes:
                    st.subheader("ProteinPaint Lollipop")
                    gene_sel = st.selectbox("Select gene", genes, key="pp_gene_single")
                    _render_proteinpaint(gene_sel, pp_genome, pp_dataset, pp_host)

                first = hits[0]
                default_query = f"{first.gene} {first.protein_change}"
                if first.cdna_change:
                    clinvar_default_query = f"{first.gene} {first.cdna_change}"
                elif first.genomic_change:
                    clinvar_default_query = f"{first.gene} {first.genomic_change}"
                else:
                    clinvar_default_query = default_query

                download_payload = {
                    "patient_info": {k: v for k, v in info.items() if k != "_sources"},
                    "patient_info_evidence": info.get("_sources", {}),
                    "variants": rows,
                }
                st.download_button(
                    label="Download JSON",
                    data=json.dumps(download_payload, indent=2),
                    file_name="report_extraction.json",
                    mime="application/json",
                )

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
                                f"**{r['title']}**  \
                                {r['authors']}  \
                                {r['journal']} ({r['year']})  \
                                PMID: {r['pmid']}  \
                                {r['url']}"
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
            info = extract_patient_info(text)
            patient_name = info.get("patient_name") or up.name
            all_summary_rows.append(
                {
                    "Patient Name": patient_name,
                    "Age": info.get("age") or "",
                    "Sex": info.get("sex") or "",
                    "Disease Name (Clinical Background)": info.get("disease_name") or "",
                    "Date Received": info.get("date_received") or "",
                    "Interpretation": info.get("interpretation") or "",
                    "Clinical Relevance": info.get("clinical_relevance") or "",
                    "Exon": info.get("exon") or "",
                    "Nucleotide Change (Genomic)": info.get("nucleotide_change") or "",
                    "ClinVar IDs (from PDF)": info.get("clinvar_ids") or "",
                    "PubMed IDs (from PDF)": info.get("pubmed_ids") or "",
                }
            )

            variant_blocks = _extract_variant_blocks(text)
            hits = extract_variants(text)
            if variant_blocks:
                for v in variant_blocks:
                    all_variant_rows.append(
                        {
                            "Patient Name": patient_name,
                            "Gene": v["gene"],
                            "Protein change": v["protein_change"],
                            "cDNA": v["cdna_change"],
                            "Genomic": v["genomic_change"],
                            "VAF": v["vaf"],
                            "Interpretation": v["interpretation"],
                            "Clinical Relevance": v["clinical_relevance"],
                            "Exon": v["exon"],
                            "Nucleotide Change (Genomic)": v["nucleotide_change"],
                            "Transcript ID": v["transcript_id"],
                            "ClinVar IDs (from PDF)": v["clinvar_ids"],
                            "PubMed IDs (from PDF)": v["pubmed_ids"],
                            "Evidence": v["evidence"],
                        }
                    )
            else:
                for h in hits:
                    all_variant_rows.append(
                        {
                            "Patient Name": patient_name,
                            "Gene": h.gene,
                            "Protein change": h.protein_change,
                            "cDNA": h.cdna_change or "",
                            "Genomic": h.genomic_change or "",
                            "VAF": h.vaf or "",
                            "Interpretation": info.get("interpretation") or "",
                            "Clinical Relevance": info.get("clinical_relevance") or "",
                            "Exon": info.get("exon") or "",
                            "Nucleotide Change (Genomic)": info.get("nucleotide_change") or "",
                            "Transcript ID": "",
                            "ClinVar IDs (from PDF)": info.get("clinvar_ids") or "",
                            "PubMed IDs (from PDF)": info.get("pubmed_ids") or "",
                            "Evidence": h.raw_line,
                        }
                    )

        st.subheader("Patient / Report Summary")
        if all_summary_rows:
            fields = [
                "Age",
                "Sex",
                "Disease Name (Clinical Background)",
                "Date Received",
                "Interpretation",
                "Clinical Relevance",
                "Exon",
                "Nucleotide Change (Genomic)",
                "ClinVar IDs (from PDF)",
                "PubMed IDs (from PDF)",
            ]
            wide_rows = []
            for field in fields:
                row = {"Field": field}
                for s in all_summary_rows:
                    row[s["Patient Name"]] = s.get(field, "")
                wide_rows.append(row)
            st.dataframe(wide_rows, use_container_width=True)
            st.download_button(
                label="Download Summary CSV",
                data=_rows_to_csv(wide_rows),
                file_name="patient_summary.csv",
                mime="text/csv",
            )
        else:
            st.info("No summaries to display.")

        st.subheader("Extracted Variants")
        if all_variant_rows:
            st.dataframe(all_variant_rows, use_container_width=True)
            st.download_button(
                label="Download Variants CSV",
                data=_rows_to_csv(all_variant_rows),
                file_name="variants.csv",
                mime="text/csv",
            )
            genes = sorted({r.get("Gene") for r in all_variant_rows if r.get("Gene")})
            if genes:
                st.subheader("ProteinPaint Lollipop")
                gene_sel = st.selectbox("Select gene", genes, key="pp_gene_multi")
                _render_proteinpaint(gene_sel, pp_genome, pp_dataset, pp_host)

            st.subheader("VAF Heatmap")
            heatmap_data = _build_heatmap_data(all_variant_rows)
            if heatmap_data:
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
                    "width": "container",
                    "height": 400,
                }
                st.vega_lite_chart(heatmap_data, spec, use_container_width=True)
            else:
                st.info("Heatmap unavailable (no VAF values).")
        else:
            st.info("No variants to display.")
    else:
        st.info("Upload one or more PDFs to begin.")

st.markdown("---")
st.caption("This tool extracts variants heuristically and is not a clinical diagnostic device.")
