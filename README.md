# Local LLaMA PDF Extractor (Streamlit + Ollama)

This app extracts structured data from patient PDF reports using a **local** LLM (Ollama) and displays it in Streamlit. It is designed for **on‑prem / local use** to keep report data private.

## Requirements

- **Ubuntu 24.04** (or similar Linux)
- Python 3.10+
- Ollama installed locally

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

## 5) Run the Streamlit app

```bash
streamlit run llm_app.py
```

Then open the local URL shown by Streamlit.

## Notes

- **Ollama models are not installed inside the Python venv.** They are managed by Ollama and stored system‑wide.
- If Ollama is on another host or port, set it in the app sidebar or set:

```bash
export OLLAMA_HOST="http://127.0.0.1:11434"
```

## Privacy & Hosting

This app is intended to run **locally** with Ollama on the same machine.

- **Streamlit Community Cloud** cannot access your local Ollama, so LLM extraction will fail there.
- For remote hosting, deploy **both** Ollama and Streamlit on the **same server** and point the app to that server’s Ollama URL.

## Troubleshooting

**Ollama: not reachable**

- Ensure `ollama serve` is running.
- Confirm the URL in the sidebar is correct.
- Check:

```bash
curl -i http://127.0.0.1:11434/api/tags
```

If it returns `{"models":[]}`, pull a model using `ollama pull <model>`.
