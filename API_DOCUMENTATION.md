
# AMR Prediction API - Documentation (v2.0.0)

## Overview
ML-powered API for predicting antibiotic resistance in **Salmonella enterica** and **Escherichia coli** genomes. This version introduces automated pathogen resolution based on Kraken2 taxonomy and WHO AWaRe classification.

**Version:** 2.0.0  
**Docker Image:** `ejiaborrita/requence_amr_project:v2.0.0`

---

## Deployment

### Pull Image
```bash
docker pull ejiaborrita/requence_amr_project:v2.0.0
```

### Run Container
Ensure you mount your local models and work directories if they are not baked into the image.
```bash
docker run -d -p 8000:8000 --name amr-api-v2 ejiaborrita/requence_amr_project:v2.0.0
```

---

## API Endpoints

### 1. Health Check
**GET** `/health`  
Checks if all required ML models for both pathogens are loaded and ready.

**Response (200 OK):**
```json
{
  "status": "healthy",
  "version": "2.0.0",
  "pathogens_configured": ["salmonella", "ecoli"]
}
```

---

### 2. Predict Resistance
**POST** `/predict`

Processes a genome assembly and returns resistance predictions.

**Request:**
- **Content-Type:** `multipart/form-data`
- **Body:**
  - `genome` (file): FASTA/FNA format (1KB - 15MB).
  - `kraken_species` (string): The species name from your Kraken2 output (e.g., "Salmonella enterica").

**cURL Example:**
```bash
curl -X POST http://localhost:8000/predict \
  -F "genome=@my_genome.fna" \
  -F "kraken_species=Escherichia coli"
```

**Response (200 OK) Highlights:**
```json
{
  "job_id": "8a2b3c4d",
  "organism_identification": {
    "detected_species": "Escherichia coli",
    "resolved_pathogen_key": "ecoli"
  },
  "antimicrobial_susceptibility": {
    "levofloxacin": {
      "who_aware_category": "Watch",
      "phenotype": "Resistant",
      "probability_score": 0.89,
      "confidence_category": "High",
      "clinical_note": "Model predicts RESISTANT with HIGH confidence...",
      "evidence": [...]
    }
  },
  "visualizations": {
    "shap_force_plots": {
       "levofloxacin": "data:image/png;base64,..."
    }
  }
}
```

---

## New Features in v2.0.0

### Pathogen Support
| Pathogen | Reference Genome | Training Samples |
| :--- | :--- | :--- |
| **Salmonella enterica** | LT2 (NC_003197) | 398 |
| **Escherichia coli** | K12 MG1655 | 571 |

### WHO AWaRe Integration
Each prediction now includes a `who_aware_category`:
* **Access:** First or second-choice antibiotics.
* **Watch:** Antibiotics with higher resistance potential.
* **Reserve:** Last-resort options.

### Quality Control (QC)
The API now performs a strict QC check. If the pipeline extracts fewer than 10 total features (Genes + SNPs + K-mers), the request will fail with a `422 Unprocessable Entity` to prevent false "Susceptible" calls due to poor sequence quality.

---

## Response Fields (Updated)

| Field | Description |
| :--- | :--- |
| `who_aware_category` | WHO classification for stewardship. |
| `clinical_note` | Plain-language interpretation of the result. |
| `shap_force_plot` | Base64-encoded PNG found under `visualizations`. |
| `passed_qc` | Boolean indicating if the genome met feature density requirements. |

---

## Integration Example (Python)

```python
import requests

payload = {'kraken_species': 'Salmonella enterica'}
files = [('genome', ('genome.fna', open('genome.fna', 'rb'), 'application/octet-stream'))]

response = requests.post('http://localhost:8000/predict', data=payload, files=files)
data = response.json()

if data['status'] == 'completed':
    species = data['organism_identification']['detected_species']
    print(f"Results for {species}:")
    for ab, res in data['antimicrobial_susceptibility'].items():
        print(f"- {ab}: {res['phenotype']} ({res['confidence_category']} confidence)")
```

---

## Disclaimer
**Research Use Only.** Predictions must NOT be used as the sole basis for clinical treatment. Medium/Low confidence results require confirmatory phenotypic AST.

---
