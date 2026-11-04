# Ref-Guided-Assembly-Scaffolder

[![License: MIT](https://shields.io)](https://opensource.org)

An automated, memory-efficient Python pipeline for chromosome-scale reference-guided genome assembly and scaffolding. This framework bridges the gap between raw sequence data parsing and ordered chromosomal stitching.

## 🧬 Architectural Features

- **Low-RAM Stream Processing:** Implements line-by-line streaming file decompression and sequence parsing. This memory-optimized architecture prevents RAM exhaustion, allowing multi-gigabyte genomic FASTA datasets to be safely processed on consumer-grade hardware.
- **Strict Spatial Filtering:** Utilizes an integrated dual-filtration layer using Map Quality thresholds (MAPQ >= 40) and query alignment coverage to clean repetitive sequence noise and sequence insertion anomalies.
- **Interval-Overlap Resolution:** Features custom coordinate interval tracking to dynamically identify nested scaffold fragments, adjust adjacent edge duplication arrays, and correct orientation mismatches via reverse-complement parsing.
- **Standardized Exporting:** Re-constructs oriented, trimmed layout consensus blocks linked by custom structural gap placeholders (`100-N` strings) and formats the final export file into industry-standard 80-character wrapped lines.

## 🦬 Validation Case Study: *Bubalus bubalis*

To evaluate and validate the architectural workflow of the pipeline, a case study was executed to construct a chromosome-scale sequence map for the water buffalo (*Bubalus bubalis*):
1. **Reference Chromosome:** Water buffalo Chromosome 6 sequence template (NCBI Accession: `CM034276.1`).
2. **Input Sequences:** Unoriented, fragmented scaffolds from the *Bubalus bubalis* genomic draft pool.
3. **Execution Pipeline:** Successfully streamed input FASTA data, executed strict alignment mapping via `minimap2 -x asm5`, filtered repetitive element fragments, adjusted spatial boundaries, and stitched a unified chromosome assembly.

## 💻 Prerequisites & Environment

The core script operates as a standalone Python engine but coordinates with open-source command-line biocomputing utilities. Ensure the following are available in your environment path:

- Python 3.x
- `minimap2` 
- `samtools`

## 🛠️ Usage

1. Clone this repository to your local system:
   ```bash
   git clone https://github.com
   cd ref-guided-assembly-scaffolder
   ```
2. Place your `reference.fasta` and `scaffolds.fasta` files in the directory.
3. Run the primary scaffolding engine:
   ```bash
   python assemble.py
   ```

## 📄 License & Intellectual Property

This software framework is released under the **MIT License**. You are free to view, download, modify, and distribute this pipeline for personal, commercial, or academic purposes, provided that original authorship credit is permanently preserved and referenced.


