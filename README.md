# Local LLaMA PDF Extractor (Streamlit + Ollama)

This app extracts structured data from patient PDF/DOCX/TXT reports using a **local** LLM (Ollama) or a **cloud** LLM (Groq) and displays it in Streamlit. It is designed for **on‑prem / local use** to keep report data private, with an optional cloud mode for hosted deployments.

## Requirements

- **Ubuntu 24.04** (or similar Linux)
- Python 3.10+
- Ollama installed locally (for local LLM mode)
- Groq API key (for Streamlit Cloud / hosted mode)
- Optional system tools for better PDF/OCR:
  - `poppler-utils` (for `pdftotext`)
  - `tesseract-ocr` (for OCR on scanned PDFs)

## 1) Install Ollama (local LLM server)

```bash
curl -fsSL https://ollama.com/install.sh | sh
```

Start Ollama (keep this running):

```bash
ollama serve
```

Verify Ollama:

```bash
ollama -v
curl -i http://127.0.0.1:11434/api/tags
```

## 2) Pull a model (local)

Quality model (large):

```bash
ollama pull llama3.1:70b
```

If your machine is smaller, use a lighter model:

```bash
ollama pull llama3.1:8b
```

## 3) Create and activate conda env

```bash
conda create -n llm_pdf python=3.10 -y
conda activate llm_pdf
```

## 4) Install Python dependencies

```bash
pip install -r llm_requirements.txt
```

Required Python packages (from `llm_requirements.txt`):
- `streamlit`
- `requests`
- `pdfplumber`
- `pytesseract`
- `pillow`
- `vl-convert-python`
- `cairosvg`
- `groq`
- `python-docx`

## 5) Run the Streamlit app

```bash
streamlit run llm_app.py
```

Then open the local URL shown by Streamlit.

## Streamlit Cloud (hosted)

Streamlit Community Cloud cannot reach your local Ollama. To run there, use Groq.

1. Create a Groq account and API key.
2. Add a Streamlit secret or environment variable:

```bash
GROQ_API_KEY=your_key_here
```

3. In the app sidebar, select `Groq (cloud)` as the LLM backend.
4. (Optional) Add `APP_URL` in secrets to show share buttons.
5. (Optional) Add `GITHUB_URL` in secrets to override the GitHub button.

## Local bot mode (chat)

- Use the **Chat** tab after uploading a report.
- The chatbot will answer questions using the latest extracted report context.

## Notes

- **Ollama models are not installed inside the Python venv.** They are managed by Ollama and stored system‑wide.
- If Ollama is on another host or port, set it in the app sidebar or set:

```bash
export OLLAMA_HOST="http://127.0.0.1:11434"
```

## Privacy & Hosting

This app is intended to run **locally** with Ollama on the same machine for privacy.

- **Streamlit Community Cloud** cannot access your local Ollama.
- For cloud hosting, use **Groq** or deploy **both** Ollama and Streamlit on the same server.

## Troubleshooting

**Ollama: not reachable**

- Ensure `ollama serve` is running.
- Confirm the URL in the sidebar is correct.
- Check:

```bash
curl -i http://127.0.0.1:11434/api/tags
```

If it returns `{"models":[]}`, pull a model using `ollama pull <model>`.

## Optional: Install system tools

```bash
sudo apt-get update
sudo apt-get install -y poppler-utils tesseract-ocr
```

## Supported File Types

- PDF (with optional OCR for scanned pages)
- DOCX (requires `python-docx`)
- TXT (UTF‑8)
