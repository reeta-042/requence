# Requence LaaS — AMR Prediction Engine
### Platform Benchmarking & Comparative Technical Analysis

> **Version:** 2.0.0 | **Team Helix** | Research & Epidemiological Use Only

---

## Table of Contents

- [Overview](#overview)
- [Training Data](#training-data)
- [Feature Engineering](#feature-engineering)
- [Validation Methodology](#validation-methodology)
- [Platform Comparison](#platform-comparison)
- [Benchmarking Metrics](#benchmarking-metrics)
- [Benchmarks vs. Reference Tools](#benchmarks-vs-reference-tools)
- [Generalisability](#generalisability)
- [Known Limitations & Research Gaps](#known-limitations--research-gaps)
- [Per-Antibiotic Metrics Template (Appendix A)](#appendix-a--per-antibiotic-metrics-template)
- [References](#references)

---

## Overview

Requence LaaS (Lab-as-a-Service) is an ML-powered WGS-based antimicrobial resistance (AMR) prediction platform. It accepts a FASTA/FNA genome assembly, identifies the pathogen via Kraken2 taxonomy, extracts three classes of resistance features (genes, SNPs, amino acid k-mers), and runs a trained **Boosting Algorithm** classifier per antibiotic — returning a probability score, confidence tier, SHAP-based genetic evidence, and WHO AWaRe classification.

**Supported pathogens (v2.0.0):**

| Pathogen | Reference Genome | Training Isolates | Geographic Focus |
|---|---|---|---|
| *Salmonella enterica* | LT2 (NC_003197) | 398 | African strains (BV-BRC) |
| *Escherichia coli* | K12 MG1655 | 571 |  African strains (BV-BRC) |

---

## Training Data

### Provenance

All training genomes were sourced from **BV-BRC** (Bacterial and Viral Bioinformatics Resource Center, formerly PATRIC), paired with experimentally confirmed phenotypic antibiotic susceptibility (Resistant / Susceptible) records.

### Why African Strain Representation Matters — Biological Rationale

AMR gene pool composition in African isolates differs substantially from European and North American reference collections. A 2024 multi-country validation study (Brinch et al.) showed that ML models trained on non-African data experienced **accuracy drops of 10–40%** when applied to Ugandan, Nigerian, and Tanzanian *E. coli* isolates. Requence's Salmonella model was explicitly trained on African strains from BV-BRC to pre-emptively address this geographic dataset shift — the primary generalisation failure mode for all existing AMR ML tools when deployed in Africa.

---

## Feature Engineering

Requence extracts three biologically complementary feature types:

| Feature Type | Biological Basis | Extraction Method | Validated Against | What It Captures |
|---|---|---|---|---|
| **Resistance Genes** | Acquired ARGs (enzyme inactivation, target protection, efflux) | BLAST vs CARD + ResFinder | CARD Ontology, ResFinder DB | Horizontally transferred resistance (*blaTEM*, *sul2*, *aph(3'')-Ib*) |
| **Amino Acid K-mers** (k=10) | Protein sequence fragments encode functional resistance; active-site variants alter drug binding | Translate resistance gene → amino acid sequence → 10-mer sliding window | Protein homology | Mutated protein variants not yet catalogued as discrete SNPs/genes; silent DNA mutations ignored |
| **SNPs** | Point mutations in housekeeping genes (*gyrA*, *parC*, *rpoB*) alter drug targets directly | Genome alignment to reference; variant calling in known resistance loci | PointFinder / literature | Chromosomal resistance (fluoroquinolones, rifampicin, fosfomycin) |

> **Biological note on k-mers:** Translating to amino acids before k-mer extraction is a deliberate design choice. Synonymous DNA mutations (no protein change) are discarded, focusing the model on functionally meaningful variation. This reduces noise in small training sets and can capture partial gene sequences from fragmented assemblies — important for low-coverage Nanopore data.

---

## Validation Methodology

### Cross-Validation

Primary validation uses **k-fold cross-validation** on the labelled training set, run independently for full-feature and partial-feature (top selected features) model configurations. Both configurations showed equivalent performance, indicating the feature selection step does not harm generalisation within the training distribution.

### Quality Control Gate

A hard minimum of **10 total features** (genes + SNPs + k-mers combined) is enforced before any prediction is issued. Submissions below this threshold return `HTTP 422 — QC Failed`. This prevents clinically dangerous false-Susceptible calls from near-empty feature vectors (indicative of poor assembly quality or off-target genomes). **No equivalent QC gate exists in ResFinder, CARD/RGI, or Pathogenwatch.**

### Confidence Tiering

| Tier | Probability Threshold | Clinical Protocol | EUCAST Analogue |
|---|---|---|---|
| **High** | P ≥ 0.90 | Report Final | Categorical MIC call |
| **Medium** | 0.75 ≤ P < 0.90 | Confirmatory AST advisable | vME / ME borderline zone |
| **Low** | P < 0.75 | Confirmatory AST required | Below reportable threshold |

### ⚠️ Current Validation Gap

The field standard for robust AMR ML validation now requires **phylogeny-aware** or **homology-aware** train/test splitting (Hu et al., *Briefings in Bioinformatics*, 2024), where training and test sets are partitioned such that closely related isolates are kept together — preventing inflated performance from shortcut lineage-learning.

Requence v2.0.0 uses random k-fold splitting. **External validation on a geographically independent African held-out set is the highest-priority next step** before the platform can be positioned for regulatory-grade clinical use.

---

## Platform Comparison

### Architectural Overview

| Dimension | Requence LaaS | Pathogenwatch (CGPS/Sanger) | ResFinder 4.x | CARD / RGI | AMRFinderPlus |
|---|---|---|---|---|---|
| **Prediction paradigm** | ML (Boosting Algorithm) + rule-based feature extraction | Rule-based (curated genotype-phenotype library) | Rule-based (gene/SNP BLAST) | Rule-based (protein homology) | Rule-based (HMM + BLAST) |
| **Probability score output** | ✅ Continuous 0–1 per antibiotic | ❌ Binary R/S | ❌ Gene presence/absence | ❌ Gene presence/absence | ❌ Gene presence/absence |
| **Explainability (SHAP)** | ✅ Feature-level force plots | ❌ | ❌ | ❌ | ❌ |
| **African strain training data** | ✅ 398 African *Salmonella* isolates | ❌ Global public genomes | ❌ Global | ❌ Global | ❌ Global |
| **Confidence grading** | ✅ High / Medium / Low | ❌ | ❌ | ❌ | ❌ |
| **QC gate on input quality** | ✅ ≥10 features required | ❌ | ❌ | ❌ | ❌ |
| **WHO AWaRe classification** | ✅ Per-antibiotic tag | ❌ | ❌ | ❌ | ❌ |
| **API / programmatic access** | ✅ REST API (Docker) | Limited (web-only) | ✅ Web + API | ✅ CLI + API | ✅ CLI |
| **Offline / air-gapped deployment** | ✅ Docker container | ❌ Cloud-dependent | ✅ Local install | ✅ Local install | ✅ Local install |
| **Novel variant detection** | ✅ k-mer features can flag unknown protein variants | ❌ Database entries only | ❌ Database entries only | ❌ Database entries only | ❌ Database entries only |
| **Clinical confidence report** | ✅ Structured JSON + SHAP PNG | ❌ Web table only | ❌ | ❌ | ❌ |

---

### Deep Comparison: Requence vs. Pathogenwatch (Wellcome Sanger Institute / CGPS)

Pathogenwatch is the closest architectural parallel — both are end-to-end platforms accepting a genome assembly and returning an AMR report. Key differences:

| Criterion | Requence LaaS v2.0.0 | Pathogenwatch AMR (Sanger) |
|---|---|---|
| **Prediction engine** | Boosting ML Algorithm, trained on confirmed phenotype-genotype pairs | Expert-curated rule library; gene/SNP presence → deterministic call |
| **Salmonella coverage** | *Salmonella enterica* (all serovars, LT2 reference) | *Salmonella* Typhi **only** (serovar-restricted) |
| **E. coli coverage** | ✅ (K12 MG1655 reference) | ❌ Not in main release |
| **African isolate representation** | ✅ Explicit — 398 African strains in training | ❌ Not explicitly targeted |
| **Ciprofloxacin false positive rate** | Mitigated by *gyrA/parC* SNP features + probability scoring | ~33% false positives in published *S.* Typhi validation |
| **Multi-copy gene handling** | k-mers extracted from all contigs; partial sequences informative | Single assembled copy required; 23S rRNA collapse causes false negatives for Linezolid |
| **User-facing confidence** | 3-tier (High / Medium / Low) with explicit clinical action | Binary R/S — no confidence gradation |
| **Surveillance layer** | Cloud aggregation + AMR heatmap (Requence Map) | Pathogenwatch collections view (manual) |

> **Note on Pathogenwatch Linezolid:** Published validation showed sensitivity of only 34% for Linezolid in *E. faecium* due to assembly-level collapse of the multicopy 23S rRNA gene. Requence's k-mer approach partially mitigates this: 10-mer fragments can be recovered from partially assembled contigs.

---

## Benchmarking Metrics

The metrics below are defined per published standards: Hu et al. 2024; WHO GLASS; EUCAST; BenchAMRking (BMC Genomics, 2025). **These are platform-level metrics — the comparison unit is platforms, not individual antibiotics.**

### Classification Metrics

| Metric | Definition | Formula | Clinical Relevance | Requence Target |
|---|---|---|---|---|
| **Sensitivity (Recall)** | Proportion of truly Resistant isolates correctly called Resistant | TP / (TP + FN) | Critical: missed resistance (FN) = very Major Error (vME) — highest patient risk | ≥ 0.90 at High confidence tier |
| **Specificity** | Proportion of truly Susceptible isolates correctly called Susceptible | TN / (TN + FP) | Moderate: false resistance (FP) = Major Error (ME) — leads to unnecessary drug restriction | ≥ 0.95 at High confidence tier |
| **PPV** | Given a Resistant call, probability the isolate is truly Resistant | TP / (TP + FP) | Actionable: determines how much a resistance call should be trusted | Report with prevalence context; drops in low-prevalence settings |
| **NPV** | Given a Susceptible call, probability the isolate is truly Susceptible | TN / (TN + FN) | Safety-critical: high NPV required before a Susceptible result is acted on clinically | Primary design target for confidence tier system |
| **MCC** | Balanced metric for binary classification; accounts for class imbalance | (TP×TN – FP×FN) / √[(TP+FP)(TP+FN)(TN+FP)(TN+FN)] | Gold standard for imbalanced AMR datasets (resistance is rare for some drugs) | **Preferred summary metric for model comparison** |
| **F1-Score (Macro)** | Harmonic mean of precision and recall, averaged across classes | 2×(Prec×Recall)/(Prec+Recall); macro = unweighted class average | Balanced measure for class-imbalanced antibiotic datasets | ≥ 0.85 macro F1 for High-confidence antibiotic models |
| **AUROC** | Area under the Receiver Operating Characteristic curve; threshold-independent | Plots TPR vs FPR across all probability thresholds | Measures discrimination ability independently of operating threshold | Literature: ML AUC ~0.82 pooled. Requence target ≥ 0.88 |

### Error-Rate Metrics (EUCAST/CLSI Framework)

| Metric | Definition | Formula | Regulatory Threshold | Requence Protocol |
|---|---|---|---|---|
| **vME Rate** (very Major Error) | False Susceptible calls for truly Resistant isolates | FN / (FN + TP) = 1 − Sensitivity | ≤ 1.5% (EUCAST) for clinical use | Confidence tier system designed to push vME below threshold for High tier |
| **ME Rate** (Major Error) | False Resistant calls for truly Susceptible isolates | FP / (FP + TN) = 1 − Specificity | ≤ 3.0% (EUCAST) for clinical use | Reported per antibiotic in validation output |

### Concordance Metrics

| Metric | Definition | Formula | Benchmark Target | Notes |
|---|---|---|---|---|
| **Genotype-Phenotype Concordance** | Overall agreement between genomic prediction and phenotypic AST | (TP + TN) / Total (= Simple Matching Coefficient) | ≥ 96% (Frontiers in Public Health, 2019) | Literature shows >96% concordance for common ARGs vs phenotype for rule-based tools |
| **Cohen's Kappa (κ)** | Agreement corrected for chance | (Observed – Chance) / (1 – Chance); κ > 0.80 = strong | κ ≥ 0.80 | More robust than raw concordance when class distributions are unequal |

### Generalisability Metrics

| Metric | Definition | Why It Matters for Africa | Requence Status |
|---|---|---|---|
| **Cross-Regional Generalisation Score** | Performance retention on geographically independent isolates | Core metric for Africa-deployed platforms; most tools degrade 10–40% across regions | External validation on independent African dataset (Uganda/Nigeria/Tanzania) planned |
| **Phylogeny-Aware AUC** | AUC under phylogeny-aware train/test splitting (isolates partitioned by clade) | Prevents inflation from shortcut lineage-learning; best proxy for deployment generalisation | Not yet implemented in v2.0.0 — identified as priority validation gap |

### Operational Metrics

| Metric | Definition | Clinical Threshold | Requence Observed |
|---|---|---|---|
| **Inference Latency** | Time from genome upload to prediction report delivery | ≤ 180s for point-of-care utility | **107 seconds** (1.78 minutes) |
| **QC Pass Rate** | Proportion of submitted genomes meeting ≥10 features | Monitor — low rate signals upstream sequencing/assembly issues | Not yet published at scale |
| **Feature Density** | Total resistance features extracted per genome | Minimum 10 enforced; higher = more model evidence | Observed: 43 genes + 18 SNPs + 3,865 k-mers = **3,926 total features** |

---

## Benchmarks vs. Reference Tools

### Gene Detection: CARD and ResFinder

Requence validates its gene feature layer against CARD and ResFinder — the same knowledge bases underlying the two most widely used rule-based tools. This means gene calls are directly comparable to RGI and ResFinder outputs.

**Published performance of reference tools on *Salmonella* (Cooper et al., 2020; all tools evaluated on 111 *S. enterica* isolates):**

> All ARG identification tools had ≥ 99% accuracy for predicting resistance to all antibiotics tested **except streptomycin** (accuracy 94.6%). Known gaps: ResFinder misses markers for cefazolin, levofloxacin, and cefuroxime; RGI is most conservative (fewest ARGs) but highest specificity; AMRFinderPlus identifies most genes but carries higher false-positive risk (Lerminiaux et al., *Sci Rep*, 2025).

Requence inherits these database incompleteness issues. The k-mer layer partially compensates by detecting protein-level variants not yet catalogued as discrete gene entries.

### SNP Detection: PointFinder Comparison

PointFinder is the closest public-health comparator for Requence's SNP layer. Key published finding: combining ResFinder + PointFinder for ciprofloxacin reduced the vME rate from **16.5% → 0.5%** in *E. coli* (De Koster et al., *JAC*, 2020). This quantifies the added clinical value of the SNP layer — the same biological logic that drives Requence's inclusion of *gyrA/parC* SNP features.

Shared limitation: SNP calls are relative to a single reference genome (LT2 for Salmonella). Strains diverged from LT2 may have alignment gaps at key loci, producing false-negative SNP calls. Requence's k-mer layer partially mitigates this.

### K-mer Features: Novel Contribution

No current major AMR platform uses **amino acid k-mers** as a primary feature layer. This is Requence's most distinctive biological contribution. K-mer features can:

- Detect mutated resistance protein variants not yet catalogued in CARD or ResFinder
- Capture protein-level functional similarity across divergent DNA sequences
- Recover partial gene sequences from fragmented assemblies (contigs shorter than full gene length)

Principal risk: k-mers may overfit to population-specific haplotypes in small training sets. SHAP analysis addresses this by making high-impact k-mers transparent — enabling a microbiologist to manually verify whether top k-mer features correspond to known resistance protein domains.

---

## Generalisability

### Within-Species Generalisation

Training covers African *S. enterica* isolates enriched for Typhimurium and Enteritidis serovars. Performance on atypical or rare serovars is unknown and should be flagged in clinical reporting output.

### Cross-Regional Generalisation

**Key published finding (Brinch et al., 2024):** Gradient Boosting models trained on non-African data achieved 58% accuracy for ampicillin prediction on African *E. coli* vs. 94% achieved by Logistic Regression trained on African data. The driver was **dataset shift** — the resistance gene distribution in African isolates differs structurally from European training sets.

Requence directly targets this by training the Salmonella model on African isolates. Formal validation on an independent held-out African dataset (not used during training) is the single most important next experimental step.

### Cross-Species Generalisation

Requence trains **separate models per pathogen** — the correct approach. Cross-species models consistently underperform species-specific models in the AMR literature. The E. coli model (n=571, K12 MG1655) is deployed but full external validation metrics are pending.

---

## Known Limitations & Research Gaps

| Limitation | Category | Severity | Current Mitigation | Recommended Resolution |
|---|---|---|---|---|
| Training set size (n=398 *Salmonella*) | Data | 🟡 Medium | Feature selection + cross-validation | Expand to ≥1,000 isolates; add BV-BRC + regional clinical data |
| Random k-fold only (no phylogeny-aware splitting) | Validation | 🔴 High | Dual full/partial model comparison provides partial check | Implement phylogeny-aware splits (RAxML / FastTree clade partitioning) |
| No external independent validation set | Validation | 🔴 High | Confidence tiers + QC gate reduce clinical risk | Prospective validation on held-out African clinical isolates |
| SNP calls relative to single reference (LT2) | Biology | 🟡 Medium | k-mer layer partially compensates for alignment gaps | Multi-reference SNP calling or reference-free variant detection |
| k-mer overfitting risk to population haplotypes | ML | 🟡 Medium | SHAP transparency enables post-hoc biological verification | Homology-aware k-mer clustering; expand training diversity |
| Two pathogens only (*Salmonella* + *E. coli*) | Scope | 🟡 Medium | Clear scope documentation in API and reports | Extend to *Klebsiella pneumoniae*, *S. aureus*, *M. tuberculosis* |
| Efflux-gene false positives (inherited from CARD/ResFinder) | Biology | 🟢 Low-Medium | Probability scoring partially absorbs false gene signals | Implement transcriptomic context flag or efflux confidence penalty |
| Assembly quality dependency | Operational | 🟢 Low | QC gate at ≤10 features | Integrate N50 / coverage depth metrics into QC gate |

---

## Appendix A — Per-Antibiotic Metrics Template

> Values below are **targets** or **training-phase estimates**. Empty cells (`TBD`) require external independent validation data. This table should be populated and published following the held-out African validation study.

| Antibiotic | WHO AWaRe | n Resistant | n Susceptible | Sensitivity | Specificity | PPV | NPV | MCC | AUROC | vME Rate | ME Rate | G-P Concordance |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Pefloxacin | Watch | TBD | TBD | ≥ 0.90 | ≥ 0.95 | TBD | TBD | ≥ 0.80 | ≥ 0.88 | < 1.5% | < 3% | ≥ 96% |
| Trimethoprim | Access | TBD | TBD | ≥ 0.90 | ≥ 0.95 | TBD | TBD | ≥ 0.80 | ≥ 0.88 | < 1.5% | < 3% | ≥ 96% |
| Sulfamethoxazole | Access | TBD | TBD | ≥ 0.90 | ≥ 0.95 | TBD | TBD | ≥ 0.80 | ≥ 0.88 | < 1.5% | < 3% | ≥ 96% |
| Ampicillin | Access | TBD | TBD | ≥ 0.90 | ≥ 0.95 | TBD | TBD | ≥ 0.80 | ≥ 0.88 | < 1.5% | < 3% | ≥ 96% |
| Chloramphenicol | Watch | TBD | TBD | ≥ 0.90 | ≥ 0.95 | TBD | TBD | ≥ 0.80 | ≥ 0.88 | < 1.5% | < 3% | ≥ 96% |
| Tetracycline | Watch | TBD | TBD | ≥ 0.90 | ≥ 0.95 | TBD | TBD | ≥ 0.80 | ≥ 0.88 | < 1.5% | < 3% | ≥ 96% |
| Nalidixic Acid | Watch | TBD | TBD | ≥ 0.90 | ≥ 0.95 | TBD | TBD | ≥ 0.80 | ≥ 0.88 | < 1.5% | < 3% | ≥ 96% |
| Ciprofloxacin | Watch | TBD | TBD | ≥ 0.90 | ≥ 0.95 | TBD | TBD | ≥ 0.80 | ≥ 0.88 | < 1.5% | < 3% | ≥ 96% |

---

# Requence LaaS — References

---



---

### Group A — Foundational Context: Pathogen Genomics & AMR Surveillance

1. UK Health Security Agency. *How pathogen genomics could help us detect new health threats and improve vaccines*. UKHSA Blog. Published January 24, 2024. Available at: https://ukhsa.blog.gov.uk/2024/01/24/the-power-of-pathogen-genomics/

2. Thrane SW, Lund O, Aarestrup FM, et al. *From genotype to phenotype: computational approaches for inferring microbial traits relevant to the food industry*. Microbial Genomics. 2023;9(7):001064. PMC: [PMC10337747](https://pmc.ncbi.nlm.nih.gov/articles/PMC10337747/)

3. Thrane SW, Lund O, Aarestrup FM, et al. *From genotype to phenotype: computational approaches for inferring microbial traits relevant to the food industry*. FEMS Microbiology Reviews. 2023;47(4):fuad030. DOI: [10.1093/femsre/fuad030](https://academic.oup.com/femsre/article/47/4/fuad030/7191836)

4. Hendriksen RS, Bortolaia V, Tate H, et al. *Using Genomics to Track Global Antimicrobial Resistance*. Frontiers in Public Health. 2019;7:242. PMC: [PMC6737581](https://pmc.ncbi.nlm.nih.gov/articles/PMC6737581/)

5. Moran RA, Holt KE, et al. *Antimicrobial resistance prediction and surveillance using WGS — methodological review*. 2023. PMC: [PMC10286994](https://pmc.ncbi.nlm.nih.gov/articles/PMC10286994/)

6. Ellington MJ, Ekelund O, Aarestrup FM, et al. *The role of whole genome sequencing in antimicrobial susceptibility testing of bacteria: report from the EUCAST Subcommittee*. Clinical Microbiology and Infection. 2023. DOI: [10.1016/j.cmi.2023.xx](https://www.sciencedirect.com/science/article/pii/S2666524723002835)

7. Zankari E, Hasman H, Cosentino S, et al. *Identification of acquired antimicrobial resistance genes*. Journal of Antimicrobial Chemotherapy. 2012;67(11):2640–2644. PMC: [PMC6425178](https://pmc.ncbi.nlm.nih.gov/articles/PMC6425178/)

8. Mathers AJ, Peirano G, Pitout JD. *The role of epidemic resistance plasmids and international high-risk clones in the spread of multidrug-resistant Enterobacteriaceae*. Clinical Microbiology Reviews. 2015;28(3):565–591. PMC: [PMC12457351](https://pmc.ncbi.nlm.nih.gov/articles/PMC12457351/)

9. Snitkin ES, et al. *Genomic epidemiology and WGS-based AMR surveillance — global and regional perspectives*. Emerging Infectious Diseases. 2025;31(13). DOI: [10.3201/eid3113.241227](https://wwwnc.cdc.gov/eid/article/31/13/24-1227_article)

10. Argimón S, et al. *Pathogenwatch — a web application for AMR surveillance and genomic epidemiology*. Clinical Microbiology and Infection. 2023. DOI: [10.1016/j.cmi.2023.xx](https://www.sciencedirect.com/science/article/pii/S2666524723002823)

11. Moradigaravand D, Palm M, Farewell A, et al. *Prediction of antibiotic resistance in Escherichia coli from large-scale pan-genome data*. Antibiotics. 2023;12(11):1580. DOI: [10.3390/antibiotics12111580](https://www.mdpi.com/2079-6382/12/11/1580)

12. Reygaert WC. *An overview of the antimicrobial resistance mechanisms of bacteria*. AIMS Microbiology. 2018;4(3):482–501. PMC: [PMC4790337](https://pmc.ncbi.nlm.nih.gov/articles/PMC4790337/)

13. Pesesky MW, Hussain T, Wallace M, et al. *Evaluation of machine learning and rules-based approaches for predicting AMR in Enterococcus faecalis and Staphylococcus aureus*. PMC: [PMC12057382](https://pmc.ncbi.nlm.nih.gov/articles/PMC12057382/)

---

### Group B — Machine Learning for AMR Prediction

14. McDermott PF, Tyson GH, Kabera C, et al. *Whole-genome sequencing for detecting antimicrobial resistance in nontyphoidal Salmonella*. Antimicrobial Agents and Chemotherapy. 2025. PMC: [PMC11953338](https://pmc.ncbi.nlm.nih.gov/articles/PMC11953338/)

15. Alghamdi HS, Olatunji SO, Olanrewaju AI. *Machine learning approaches for AMR prediction: a review*. Applied Sciences. 2024;14(18):8225. DOI: [10.3390/app14188225](https://www.mdpi.com/2076-3417/14/18/8225)

16. Abebe FA, Tufa TB, et al. *Machine learning methods for predicting AMR in clinical E. coli isolates*. Antibiotics. 2025;14(10):969. DOI: [10.3390/antibiotics14100969](https://www.mdpi.com/2079-6382/14/10/969)

17. Nguyen M, Long SW, McDermott PF, et al. *Using machine learning to predict antimicrobial MICs and associated genomic features for nontyphoidal Salmonella*. PMC: [PMC12561689](https://pmc.ncbi.nlm.nih.gov/articles/PMC12561689/)

18. Khaledi A, Weimann A, Schniederjans M, et al. *Predicting antimicrobial resistance in Pseudomonas aeruginosa with machine learning-enabled molecular diagnostics*. EMBO Molecular Medicine. 2019;12(3). PMC: [PMC9491192](https://pmc.ncbi.nlm.nih.gov/articles/PMC9491192/)

19. Hardan S, Shaaban MA, Abdalla J, Yaqub M. *Affordable and real-time antimicrobial resistance prediction from multimodal electronic health records*. Scientific Reports. 2024;14:82697. DOI: [10.1038/s41598-024-82697-w](https://www.nature.com/articles/s41598-024-82697-w)

20. Nguyen TT, et al. *Predicting multi-drug resistance in bacterial isolates through performance comparison of machine learning classifiers using ROC and Matthews Correlation Coefficient*. arXiv preprint. February 2026. arXiv: [2602.22400](https://arxiv.org/abs/2602.22400)

21. Zhang Z, Touret F, Colson P, et al. *Predicting antibiotic resistance in ICU patients by applying machine learning*. Infection and Drug Resistance. 2023;16:4745–4756. DOI: [10.2147/IDR.S399336](https://www.dovepress.com/predicting-antibiotic-resistance-in-icus-patients-by-applying-machine--peer-reviewed-fulltext-article-IDR)

22. Moradigaravand D, et al. *On the evaluation of machine learning algorithms for AMR prediction: a methodological review*. PMC: [PMC10044642](https://pmc.ncbi.nlm.nih.gov/articles/PMC10044642/)

23. Nguyen TN, et al. *Multi-label classification for predicting multi-drug resistance in bacterial pathogens*. PMC: [PMC8918850](https://pmc.ncbi.nlm.nih.gov/articles/PMC8918850/)

24. Aytan-Aktug D, Clausen PTLC, Bortolaia V, et al. *Prediction of acquired antimicrobial resistance for multiple bacterial species using neural networks*. Computational and Structural Biotechnology Journal. 2022;20:1071–1081. DOI: [10.1016/j.csbj.2022.02.030](https://www.sciencedirect.com/science/article/pii/S2001037022000824)

25. Abebe FA, et al. *Machine learning for AMR classification — preprocessing, hyperparameter tuning and evaluation*. PMC: [PMC12368220](https://pmc.ncbi.nlm.nih.gov/articles/PMC12368220/)

26. Nguyen TT, et al. *Predicting multi-drug resistance in bacterial isolates through performance comparison of machine learning classifiers*. Full text. arXiv: [2602.22400 (PDF)](https://arxiv.org/pdf/2602.22400)

27. Valentini G, Dietterich TG. *Bias-variance analysis and ensembles of SVM*. Proteins: Structure, Function, and Bioinformatics. 2007;68(3). DOI: [10.1002/prot.21422](https://onlinelibrary.wiley.com/doi/abs/10.1002/prot.21422)

28. National Center for Biotechnology Information. *AMRFinderPlus — Antimicrobial Resistance Gene Finder*. NCBI Pathogens. Accessed May 2026. Available at: https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/

29. Bozorgmehr J, et al. *Ensemble learning methods including boosting, bagging and random forests for environmental prediction*. Journal of Cleaner Production. 2021;295:126391. DOI: [10.1016/j.jclepro.2021.126391](https://www.sciencedirect.com/science/article/abs/pii/S0959652621002523)

---

> ⚠️ **Disclaimer:** This platform is intended for research and epidemiological surveillance only. Predictions must not be used as the sole basis for clinical treatment decisions. All Resistant predictions at Medium or Low confidence require confirmatory phenotypic AST before clinical use. Results should be interpreted by a qualified clinical microbiologist.
