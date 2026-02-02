---
name: peptidase-cleavage-analysis
description: Automates the process of characterizing peptidase specificity through literature review and predicting cleavage sites in protein sequences. Integrates literature search to identify cleavage patterns and applies them to find potential cleavage sites.
---

# Peptidase Cleavage Analysis

## Overview

This skill provides a workflow to analyze peptidases (proteases) by combining automated literature review with sequence analysis. It helps researchers determine the cleavage specificity of a given peptidase from scientific literature and then uses that information to predict cleavage sites in a target protein sequence.

## When to Use This Skill

Use this skill when you:
- Need to find the cleavage sites of a specific peptidase in a protein sequence.
- Don't know the exact cleavage pattern of a peptidase and need to search the literature for it.
- Want to batch process multiple peptidases or substrates.
- Need to summarize known specificity data for a less common enzyme.

## Core Capabilities

### 1. Literature Review for Specificity

Uses AI-powered search (e.g., Perplexity or simple web search) to find consensus cleavage patterns from literature.

```bash
python scripts/analyze_peptidase.py --peptidase "Caspase-3" --mode literature_only
```

### 2. Cleavage Site Prediction

Predicts cleavage sites based on a provided or discovered pattern.

```bash
python scripts/analyze_peptidase.py --peptidase "Trypsin" --sequence "MKWVTFISLLLLFSSAYSRGVFRR" --pattern "(?<=[KR])(?!P)"
```

### 3. End-to-End Analysis

Combines literature search and prediction.

```bash
python scripts/analyze_peptidase.py --peptidase "Factor Xa" --sequence "MGRLVGG"
```

## Dependencies

- `biopython`: For sequence handling.
- `perplexity-search` (optional but recommended): For literature review capability.
- `regex`: For advanced pattern matching.

## Installation

```bash
uv pip install biopython regex
```

## Usage

### Basic Usage

**Find cleavage sites for a known enzyme:**

```bash
python scripts/analyze_peptidase.py --peptidase "Trypsin" --sequence "MKWVTFISLLLLFSSAYSRGVFRRDAH KSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLKECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYARRHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFEQLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTLSEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLVAASQAALGL"
```

**Find specificity from literature:**

```bash
python scripts/analyze_peptidase.py --peptidase "TEV protease" --mode literature_only
```
