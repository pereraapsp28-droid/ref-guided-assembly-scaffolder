!apt-get update && apt-get install -y minimap2 samtools



import urllib.request
import os

# 1. PASTE YOUR URLS HERE
REFERENCE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=CM034276.1&rettype=fasta&retmode=text"
SCAFFOLDS_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/794/615/GCA_004794615.1_Bubbub1.0/GCA_004794615.1_Bubbub1.0_genomic.fna.gz"

# Local filenames to save as
REF_FILE = "reference.fasta"
SCAFFOLDS_FILE = "scaffolds.fasta"

print("Downloading reference genome...")
urllib.request.urlretrieve(REFERENCE_URL, REF_FILE)
print("Reference download complete.")

print("Downloading 2GB scaffold file (this may take a few minutes)...")
urllib.request.urlretrieve(SCAFFOLDS_URL, SCAFFOLDS_FILE)
print("Scaffolds download complete.")





import gzip
import shutil
import os

# Your downloaded file names
REF_FILE = "reference.fasta"
SCAFFOLDS_FILE = "scaffolds.fasta"

def decompress_if_gzipped(filepath):
    # Check for the unique gzip magic number signature (0x1f 0x8b)
    if os.path.exists(filepath):
        with open(filepath, "rb") as f:
            is_gzip = f.read(2) == b"\x1f\x8b"
            
        if is_gzip:
            print(f"Detected compressed gzip format for {filepath}. Extracting...")
            temp_unzipped = filepath + ".tmp"
            
            # Decompress line-by-line to avoid loading 2GB into RAM
            with gzip.open(filepath, "rb") as f_in:
                with open(temp_unzipped, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            
            # Replace the old compressed file with the clean uncompressed file
            os.remove(filepath)
            os.rename(temp_unzipped, filepath)
            print(f"Successfully unzipped: {filepath}")
        else:
            print(f"{filepath} is already uncompressed text. No unzipping needed.")
    else:
        print(f"Error: {filepath} does not exist.")

# Run the decompression on both files
print("Checking reference file...")
decompress_if_gzipped(REF_FILE)

print("\nChecking 2GB scaffold file...")
decompress_if_gzipped(SCAFFOLDS_FILE)





import os
for f in ["reference.fasta", "scaffolds.fasta"]:
    if os.path.exists(f):
        print(f"{f}: {os.path.getsize(f) / (1024*1024):.2f} MB")
    else:
        print(f"⚠️ {f} is missing completely!")






# Print the headers inside your reference file to find the exact name
with open(REF_FILE, "r") as f:
    for _ in range(50): # Check first 50 lines
        line = f.readline()
        if not line:
            break
        if line.startswith(">"):
            print(line.strip())




with open("reference.fasta", "r") as f:
    for line in f:
        if line.startswith(">"):
            # Splits the header at the first space to get the exact clean name
            print(f"Your True Reference Name Is: '{line.strip()[1:].split()[0]}'")
            break






# -x asm5: Keeps the optimization for under 1% divergence (same species)
# -c: Forces minimap2 to generate CIGAR strings (calculates exact alignment identity)
# --min-occ 200: Ignores highly repetitive, low-complexity genomic repeats
!minimap2 -t 32 -x asm5 -c --min-occ 200 --secondary=no reference.fasta scaffolds.fasta > alignments.paf







# --- CONFIGURATION ---
TRUE_REF_NAME = "CM034276.1"  
PAF_FILE = "alignments.paf"
SCAFFOLDS_FILE = "scaffolds.fasta"
CLEAN_OUTPUT_FASTA = "filtered_target_scaffolds.fasta"

# MIN_COVERAGE: At least 50% of the scaffold length must match Chromosome 6
MIN_COVERAGE = 0.50  
# ---------------------

print(f"Analyzing alignment coverage fractions against {TRUE_REF_NAME}...")

scaffold_match_bases = {}
scaffold_lengths = {}

with open(PAF_FILE, "r") as paf:
    for line in paf:
        parts = line.strip().split("\t")
        if len(parts) < 12:
            continue
            
        scaffold_id = parts[0]
        scaffold_len = int(parts[1])
        target_name = parts[5]
        alignment_len = int(parts[10]) # Number of residue matches/aligned bases
        mapq = int(parts[11])
        
        if target_name == TRUE_REF_NAME and mapq >= 40:
            # Track scaffold properties
            scaffold_lengths[scaffold_id] = scaffold_len
            # Accumulate matching blocks for each scaffold
            scaffold_match_bases[scaffold_id] = scaffold_match_bases.get(scaffold_id, 0) + alignment_len

# Apply structural coverage logic
valid_scaffolds = set()
for scaf_id, matched_bp in scaffold_match_bases.items():
    total_len = scaffold_lengths[scaf_id]
    coverage_fraction = matched_bp / total_len
    
    # Keep only if it meets our physical sequence similarity coverage threshold
    if coverage_fraction >= MIN_COVERAGE:
        valid_scaffolds.add(scaf_id)

print(f"\nFiltered out the noise!")
print(f"Scaffolds dropped as repetitive elements: {len(scaffold_match_bases) - len(valid_scaffolds)}")
print(f"Scaffolds retained with high structural similarity: {len(valid_scaffolds)}")

print("\nStreaming 2GB file to save true structural sequences...")
scaffolds_saved = 0
with open(SCAFFOLDS_FILE, "r") as infile, open(CLEAN_OUTPUT_FASTA, "w") as outfile:
    keep_sequence = False
    for line in infile:
        if line.startswith(">"):
            clean_id = line[1:].strip().split()[0]
            if clean_id in valid_scaffolds:
                keep_sequence = True
                outfile.write(line)
                scaffolds_saved += 1
            else:
                keep_sequence = False
        else:
            if keep_sequence:
                outfile.write(line)

print(f"Done! Successfully saved {scaffolds_saved} high-similarity scaffolds to {CLEAN_OUTPUT_FASTA}")






import os

# --- CONFIGURE STRINGS ---
PAF_FILE = "alignments.paf"
FILTERED_SCAFFOLDS = "filtered_target_scaffolds.fasta"
TARGET_ACCESSION = "CM034276.1" 
# -------------------------

print("1. Loading your strictly filtered sequences into memory...")
scaffold_sequences = {}
current_id, current_seq = None, []

with open(FILTERED_SCAFFOLDS, "r") as f:
    for line in f:
        if line.startswith(">"):
            if current_id: 
                scaffold_sequences[current_id] = "".join(current_seq)
            current_id = line[1:].strip().split()[0]
            current_seq = []
        else:
            current_seq.append(line.strip())
    if current_id: 
        scaffold_sequences[current_id] = "".join(current_seq)

# This confirms it matches your first cell perfectly
print(f"Loaded exactly {len(scaffold_sequences)} unique whitelisted scaffolds.")

print("\n2. Parsing STRICT alignment coordinates (MAPQ >= 40)...")
layout = []
with open(PAF_FILE, "r") as paf:
    for line in paf:
        parts = line.strip().split("\t")
        if len(parts) < 12: continue
        
        scaffold_id = parts[0]
        target_name = parts[5]
        mapq = int(parts[11])
        
        # CRITICAL FIX 1: Only look at the 1,519 scaffolds that passed your strict coverage test
        if scaffold_id not in scaffold_sequences:
            continue
            
        # CRITICAL FIX 2: Only allow maximum confidence alignments to your target chromosome
        if TARGET_ACCESSION.lower() in target_name.lower() and mapq >= 40:
            layout.append({
                'id': scaffold_id,
                'ref_start': int(parts[7]), 
                'ref_end': int(parts[8]), 
                'strand': parts[4]
            })

# Sort all matching premium fragments by their starting position on the chromosome
layout.sort(key=lambda x: x['ref_start'])
print(f"Found and sorted {len(layout)} ultra-strict layout intervals matching your chromosome.")







print("Cleaning layout artifacts and preparing unique puzzle pieces (Optimized)...")

ordered_scaffold_ids = []
last_end = -1
duplicates_removed = 0
overlaps_adjusted = 0

# Track scaffolds we already committed to avoid processing them multiple times
seen_scaffolds = set()

for item in layout:
    scaf_id = item['id']
    if scaf_id not in scaffold_sequences: 
        continue
    
    # 1. Skip if we already used this exact scaffold piece (drops internal duplicates immediately)
    if scaf_id in seen_scaffolds:
        duplicates_removed += 1
        continue

    # 2. FIX NESTED DUPLICATES: Skip if it falls inside a larger, previous scaffold coordinate footprint
    if item['ref_end'] <= last_end:
        duplicates_removed += 1
        continue
        
    # 3. Track overlapping adjustments without modifying the text strings yet
    if item['ref_start'] < last_end:
        overlaps_adjusted += 1

    # Commit the clean layout piece
    ordered_scaffold_ids.append({
        'id': scaf_id,
        'strand': item['strand'],
        'overlap_from_start': max(0, last_end - item['ref_start'])
    })
    
    seen_scaffolds.add(scaf_id)
    last_end = item['ref_end']

print("\nProcessing final physical sequences and building unified strand...")

# Helper function to flip backward sequences safely
def reverse_complement(seq):
    comp = {'A':'T','T':'A','C':'G','G':'C','N':'N','a':'t','t':'a','c':'g','g':'c','n':'n'}
    return "".join(comp.get(base, base) for base in reversed(seq))

ordered_sequences = []
for entry in ordered_scaffold_ids:
    seq = scaffold_sequences[entry['id']]
    
    # Flip if orientation is minus
    if entry['strand'] == "-":
        seq = reverse_complement(seq)
        
    # Trim the edge overlap if needed
    if entry['overlap_from_start'] > 0 and entry['overlap_from_start'] < len(seq):
        seq = seq[entry['overlap_from_start']:]
        
    if len(seq) > 0:
        ordered_sequences.append(seq)

print(f"\n✅ Duplication filtration complete!")
print(f"Dropped completely duplicated/nested sequences: {duplicates_removed:,}")
print(f"Trimmed and resolved overlapping boundaries: {overlaps_adjusted:,}")
print(f"Retained {len(ordered_sequences):,} clean unique pieces for final merging.")








# Ensure TARGET_ACCESSION is set correctly from your previous setup
TARGET_ACCESSION = "CM034276.1" # <-- Fixed to match your exact chromosome ID!
FINAL_CHROMOSOME_FASTA = "chromosome_6_assembled_clean.fasta"
gap_filler = "N" * 100  # Standard genomic placeholder gap

# len(ordered_sequences) dynamically prints your actual sweet-spot clean count
print(f"Stitching your {len(ordered_sequences)} clean fragments with 100-N structural gaps...")
full_sequence = gap_filler.join(ordered_sequences)

print(f"Saving final single sequence to your workspace as: {FINAL_CHROMOSOME_FASTA}...")
with open(FINAL_CHROMOSOME_FASTA, "w") as out:
    out.write(f">{TARGET_ACCESSION}_assembled_chromosome_scale\n")
    # Save in standard 80-character line width format
    for i in range(0, len(full_sequence), 80):
        out.write(full_sequence[i:i+80] + "\n")

print("\n--- 🏁 PIPELINE COMPLETED SUCCESSFULLY ---")
print(f"Your final single clean sequence file is saved as: {FINAL_CHROMOSOME_FASTA}")
print(f"Final Assembled Chromosome Length: {len(full_sequence):,} base pairs.")


"""
STANDALONE CORRECTED ASSEMBLY SCRIPT
Independent of your original notebook cells — doesn't read or rely on any
variables from your original code. Only needs these files to already exist:
  - alignments.paf
  - scaffolds.fasta
  - reference.fasta

Fixes applied (vs. the original approach):
  1. Coverage filter uses matching bases (PAF col 9), not block length (col 10)
  2. Each scaffold is clipped to ONLY its aligned region (q_start:q_end)
     before being used — no more inserting whole unaligned scaffolds
  3. Reverse-complement is applied AFTER clipping, so orientation is
     correct for the actual aligned sub-sequence
  4. Overlap trimming happens AFTER clipping+orientation, so it trims the
     correct end regardless of strand
  5. Dedup is based only on reference-coordinate overlap, not scaffold ID
     — so a scaffold with multiple genuine, separate, non-overlapping
       alignments contributes all of them, not just one
  6. No per-block identity filter (tested separately; it fragments the
     assembly further and lowers overall identity — not worth it)

Output: corrected_assembly.fasta
"""

import os

PAF_FILE = "alignments.paf"
SCAFFOLDS_FILE = "scaffolds.fasta"
TARGET_ACCESSION = "CM034276.1"
MIN_MAPQ = 40
MIN_COVERAGE = 0.50
OUTPUT_FASTA = "chromosome_6_assembled_clean.fasta"


def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
             'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
    return "".join(comp.get(base, base) for base in reversed(seq))


# --- Step 1: parse PAF, compute per-scaffold coverage using MATCHING BASES ---
print("Step 1: parsing alignments.paf for coverage filtering (using true matching bases)...")
scaffold_match_bases = {}
scaffold_lengths = {}

with open(PAF_FILE, "r") as paf:
    for line in paf:
        parts = line.strip().split("\t")
        if len(parts) < 12:
            continue
        scaffold_id = parts[0]
        scaffold_len = int(parts[1])
        target_name = parts[5]
        matches = int(parts[9])   # <-- true matching bases, not block length
        mapq = int(parts[11])

        if target_name == TARGET_ACCESSION and mapq >= MIN_MAPQ:
            scaffold_lengths[scaffold_id] = scaffold_len
            scaffold_match_bases[scaffold_id] = scaffold_match_bases.get(scaffold_id, 0) + matches

valid_scaffolds = set()
for scaf_id, matched_bp in scaffold_match_bases.items():
    total_len = scaffold_lengths[scaf_id]
    if total_len > 0 and (matched_bp / total_len) >= MIN_COVERAGE:
        valid_scaffolds.add(scaf_id)

print(f"  {len(scaffold_match_bases)} scaffolds had qualifying alignments; "
      f"{len(valid_scaffolds)} passed the {MIN_COVERAGE:.0%} coverage threshold.")

# --- Step 2: load only the valid scaffolds' sequences into memory ---
print("\nStep 2: loading valid scaffold sequences from scaffolds.fasta...")
scaffold_sequences = {}
current_id, current_seq, keep = None, [], False

with open(SCAFFOLDS_FILE, "r") as f:
    for line in f:
        if line.startswith(">"):
            if current_id and keep:
                scaffold_sequences[current_id] = "".join(current_seq)
            current_id = line[1:].strip().split()[0]
            current_seq = []
            keep = current_id in valid_scaffolds
        else:
            if keep:
                current_seq.append(line.strip())
    if current_id and keep:
        scaffold_sequences[current_id] = "".join(current_seq)

print(f"  Loaded {len(scaffold_sequences)} scaffold sequences.")

# --- Step 3: re-parse PAF, this time keeping query coordinates too ---
print("\nStep 3: re-parsing alignments.paf for query coordinates (q_start/q_end)...")
layout = []
with open(PAF_FILE, "r") as paf:
    for line in paf:
        parts = line.strip().split("\t")
        if len(parts) < 12:
            continue
        scaffold_id = parts[0]
        q_start, q_end = int(parts[2]), int(parts[3])
        strand = parts[4]
        target_name = parts[5]
        ref_start, ref_end = int(parts[7]), int(parts[8])
        mapq = int(parts[11])

        if scaffold_id not in valid_scaffolds:
            continue
        if target_name != TARGET_ACCESSION or mapq < MIN_MAPQ:
            continue

        layout.append({'id': scaffold_id, 'ref_start': ref_start, 'ref_end': ref_end,
                        'q_start': q_start, 'q_end': q_end, 'strand': strand})

layout.sort(key=lambda x: x['ref_start'])
print(f"  Found {len(layout)} candidate alignment blocks.")

# --- Step 4: dedup by REFERENCE OVERLAP ONLY (not by scaffold ID) ---
print("\nStep 4: removing overlapping/nested blocks (keeping ALL non-overlapping ones per scaffold)...")
last_end = -1
kept = []
dropped = 0
for item in layout:
    if item['ref_end'] <= last_end:
        dropped += 1
        continue
    ref_overlap = max(0, last_end - item['ref_start'])
    kept.append({**item, 'ref_overlap': ref_overlap})
    last_end = item['ref_end']

print(f"  Kept {len(kept)} non-overlapping blocks (dropped {dropped} overlapping/nested).")

# --- Step 5: clip -> orient -> trim -> stitch ---
print("\nStep 5: clipping to aligned region, orienting, trimming, and stitching...")
pieces = []
for entry in kept:
    raw_seq = scaffold_sequences[entry['id']]
    clipped = raw_seq[entry['q_start']:entry['q_end']]          # clip FIRST
    if entry['strand'] == "-":
        clipped = reverse_complement(clipped)                    # orient SECOND
    if entry['ref_overlap'] > 0 and entry['ref_overlap'] < len(clipped):
        clipped = clipped[entry['ref_overlap']:]                 # trim THIRD
    if len(clipped) > 0:
        pieces.append(clipped)

print(f"  Built {len(pieces)} clipped/oriented/trimmed pieces.")

gap_filler = "N" * 100
full_sequence = gap_filler.join(pieces)

with open(OUTPUT_FASTA, "w") as out:
    out.write(f">{TARGET_ACCESSION}_corrected_assembly\n")
    for i in range(0, len(full_sequence), 80):
        out.write(full_sequence[i:i + 80] + "\n")

print(f"\n--- DONE ---")
print(f"Saved: {OUTPUT_FASTA}")
print(f"Assembled length: {len(full_sequence):,} bp")
print(f"File size on disk: {os.path.getsize(OUTPUT_FASTA) / 1e6:.2f} MB")



















"""
ASSEMBLY QC REPORT
===================
Standalone quality-control script for a reference-guided chromosome
assembly. Independent of any other script — does not regenerate, modify,
or depend on any other pipeline step. Reads two existing files and reports
on assembly quality only.

Inputs (must already exist):
    - reference.fasta   : the reference chromosome
    - ASSEMBLY_FILE      : your assembled chromosome FASTA (set below)

Outputs:
    - A printed QC report (coverage, identity, contiguity, gap profile)
    - assembly_qc_report.txt : the same report saved to disk

Usage:
    Edit ASSEMBLY_FILE below to point at whichever assembly you want to
    QC, then run this as a single cell / script.
"""

import os
import subprocess
from datetime import datetime

# ----------------------------- CONFIGURATION -----------------------------
REFERENCE_FILE = "reference.fasta"
ASSEMBLY_FILE = "chromosome_6_assembled_clean.fasta"   # <-- point this at the file you want to QC
MIN_MAPQ = 0                                  # report includes all mapped blocks; adjust if needed
PAF_OUTPUT = "qc_alignment.paf"
REPORT_FILE = "assembly_qc_report.txt"
# ---------------------------------------------------------------------------


def run_minimap2(reference, query, out_paf):
    cmd = ["minimap2", "-x", "asm5", "-c", "--secondary=no", reference, query]
    with open(out_paf, "w") as out:
        result = subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"minimap2 failed: {result.stderr}")


def fasta_length(path):
    total = 0
    with open(path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def parse_paf_blocks(paf_path, min_mapq=0):
    blocks = []
    with open(paf_path, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            mapq = int(parts[11])
            if mapq < min_mapq:
                continue
            blocks.append({
                "ref_start": int(parts[7]),
                "ref_end": int(parts[8]),
                "matches": int(parts[9]),
                "block_len": int(parts[10]),
                "mapq": mapq,
            })
    return blocks


def merged_reference_span(blocks):
    """Non-redundant bp of the reference covered by at least one block."""
    if not blocks:
        return 0
    intervals = sorted((b["ref_start"], b["ref_end"]) for b in blocks)
    span = 0
    cur_start, cur_end = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_end:
            cur_end = max(cur_end, e)
        else:
            span += cur_end - cur_start
            cur_start, cur_end = s, e
    span += cur_end - cur_start
    return span


def n_gap_profile(path, min_gap_len=50):
    seq_parts = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq_parts.append(line.strip())
    seq = "".join(seq_parts).upper()

    gaps = []
    i, n = 0, len(seq)
    while i < n:
        if seq[i] == "N":
            j = i
            while j < n and seq[j] == "N":
                j += 1
            if (j - i) >= min_gap_len:
                gaps.append((i, j, j - i))
            i = j
        else:
            i += 1
    return gaps, len(seq)


def n50(lengths):
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    running = 0
    for length in lengths:
        running += length
        if running >= total / 2:
            return length
    return 0


def build_report():
    lines = []

    def log(msg=""):
        print(msg)
        lines.append(msg)

    log("=" * 78)
    log("ASSEMBLY QUALITY CONTROL REPORT")
    log(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"Reference: {REFERENCE_FILE}")
    log(f"Assembly:  {ASSEMBLY_FILE}")
    log("=" * 78)

    if not os.path.exists(REFERENCE_FILE):
        log(f"ERROR: reference file not found: {REFERENCE_FILE}")
        return lines
    if not os.path.exists(ASSEMBLY_FILE):
        log(f"ERROR: assembly file not found: {ASSEMBLY_FILE}")
        return lines

    ref_len = fasta_length(REFERENCE_FILE)
    asm_len = fasta_length(ASSEMBLY_FILE)

    log(f"\nReference length:  {ref_len:,} bp")
    log(f"Assembly length:   {asm_len:,} bp")

    log("\nAligning assembly to reference (minimap2 -x asm5)...")
    run_minimap2(REFERENCE_FILE, ASSEMBLY_FILE, PAF_OUTPUT)

    blocks = parse_paf_blocks(PAF_OUTPUT, min_mapq=MIN_MAPQ)
    if not blocks:
        log("\nNo alignment blocks found — assembly does not align to reference.")
        return lines

    aligned_span = merged_reference_span(blocks)
    total_matches = sum(b["matches"] for b in blocks)
    total_block_len = sum(b["block_len"] for b in blocks)
    mean_identity = (total_matches / total_block_len * 100) if total_block_len else 0.0
    unaligned_bp = asm_len - aligned_span
    unaligned_pct = (unaligned_bp / asm_len * 100) if asm_len else 0.0

    log("\n" + "-" * 78)
    log("COVERAGE & IDENTITY")
    log("-" * 78)
    log(f"  Reference coverage:        {aligned_span:,} bp  ({aligned_span/ref_len*100:.2f}% of reference)")
    log(f"  Mean alignment identity:   {mean_identity:.2f}%")
    log(f"  Alignment blocks:          {len(blocks):,}")
    log(f"  Sequence not accounted for by alignment: {unaligned_bp:,} bp ({unaligned_pct:.2f}% of assembly)")
    log("    (large values here indicate padding/unaligned sequence in the assembly)")

    # Contiguity: treat each PAF block's length as a "fragment" length for N50
    block_lengths = [b["ref_end"] - b["ref_start"] for b in blocks]
    n50_value = n50(block_lengths)
    log("\n" + "-" * 78)
    log("CONTIGUITY")
    log("-" * 78)
    log(f"  N50 (of aligned blocks):   {n50_value:,} bp")
    log(f"  Largest aligned block:     {max(block_lengths):,} bp")
    log(f"  Smallest aligned block:    {min(block_lengths):,} bp")
    log(f"  Median aligned block:      {sorted(block_lengths)[len(block_lengths)//2]:,} bp")

    # Gap profile
    gaps, seq_len = n_gap_profile(ASSEMBLY_FILE)
    total_gap_bp = sum(g[2] for g in gaps)
    log("\n" + "-" * 78)
    log("GAP PROFILE (N-runs in assembly)")
    log("-" * 78)
    log(f"  Total gap regions (>=50 N's): {len(gaps):,}")
    log(f"  Total gap bp:                 {total_gap_bp:,} bp ({total_gap_bp/seq_len*100:.2f}% of assembly)")
    if gaps:
        top_gaps = sorted(gaps, key=lambda g: -g[2])[:10]
        log("\n  10 largest gaps:")
        for start, end, length in top_gaps:
            log(f"    pos {start:>12,}-{end:<12,}  length {length:>8,} bp")

    # Overall quality summary
    log("\n" + "-" * 78)
    log("SUMMARY")
    log("-" * 78)
    coverage_pct = aligned_span / ref_len * 100
    log(f"  Reference coverage: {coverage_pct:.2f}%")
    log(f"  Mean identity:      {mean_identity:.2f}%")
    log(f"  Padding/unaligned:  {unaligned_pct:.2f}% of file")
    log(f"  Contiguity (N50):   {n50_value:,} bp")
    log("=" * 78)

    return lines


if __name__ == "__main__":
    report_lines = build_report()
    with open(REPORT_FILE, "w") as f:
        f.write("\n".join(report_lines))
    print(f"\nFull report saved to: {REPORT_FILE}")











for plants repeat filter and code before qc

"""
REPEAT REGION FILTER
=====================
Standalone script — does not touch or depend on any other script.

Problem: repetitive DNA (common in plant genomes like coconut) causes many
scaffold regions to align to multiple different reference positions, or
many different scaffolds to pile up on the same reference region. This
inflates block count, causes overlap redundancy, and drags down average
identity (repeats tend to align messily).

What this does:
  1. Reads alignments.paf
  2. Measures how many alignment blocks cover each region of the reference
     (depth of coverage)
  3. Flags reference regions with unusually high depth as likely repetitive
     (default: depth > 3x the median depth)
  4. Writes out a NEW paf file with repetitive-region blocks removed, so you
     can re-run your assembly/stitching pipeline on the cleaned file
  5. Prints a summary so you can see how much of the genome/alignments were
     flagged as repetitive

Input:  alignments.paf
Output: alignments_repeat_filtered.paf
        repeat_regions_report.txt
"""

import os
from collections import defaultdict

PAF_FILE = "alignments.paf"
OUTPUT_PAF = "alignments_repeat_filtered.paf"
REPORT_FILE = "repeat_regions_report.txt"
DEPTH_MULTIPLIER = 3.0   # flag regions with depth > this x the median depth
BIN_SIZE = 1000          # bp resolution for depth calculation


def parse_paf(path):
    records = []
    with open(path, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            records.append({
                "raw": line,
                "scaffold_id": parts[0],
                "target_name": parts[5],
                "ref_start": int(parts[7]),
                "ref_end": int(parts[8]),
            })
    return records


def compute_depth(records, bin_size):
    """Depth per reference bin, per target sequence (in case of multiple
    reference sequences in the PAF)."""
    depth = defaultdict(lambda: defaultdict(int))  # target -> bin -> count
    for r in records:
        target = r["target_name"]
        start_bin = r["ref_start"] // bin_size
        end_bin = r["ref_end"] // bin_size
        for b in range(start_bin, end_bin + 1):
            depth[target][b] += 1
    return depth


def median(values):
    s = sorted(values)
    n = len(s)
    if n == 0:
        return 0
    mid = n // 2
    if n % 2 == 0:
        return (s[mid - 1] + s[mid]) / 2
    return s[mid]


def flag_repetitive_bins(depth, multiplier):
    """Return set of (target, bin) tuples considered repetitive."""
    flagged = set()
    thresholds = {}
    for target, bins in depth.items():
        values = list(bins.values())
        med = median(values)
        threshold = max(med * multiplier, 1)  # avoid zero threshold
        thresholds[target] = (med, threshold)
        for b, count in bins.items():
            if count > threshold:
                flagged.add((target, b))
    return flagged, thresholds


def block_is_repetitive(record, flagged_bins, bin_size, overlap_fraction=0.5):
    """A block is considered repetitive if more than `overlap_fraction` of
    its bins fall in flagged (high-depth) regions."""
    target = record["target_name"]
    start_bin = record["ref_start"] // bin_size
    end_bin = record["ref_end"] // bin_size
    total_bins = end_bin - start_bin + 1
    flagged_count = sum(1 for b in range(start_bin, end_bin + 1) if (target, b) in flagged_bins)
    return (flagged_count / total_bins) > overlap_fraction


if __name__ == "__main__":
    lines = []

    def log(msg=""):
        print(msg)
        lines.append(msg)

    log("=" * 74)
    log("REPEAT REGION DETECTION & FILTERING")
    log("=" * 74)

    if not os.path.exists(PAF_FILE):
        raise SystemExit(f"Missing input file: {PAF_FILE}")

    records = parse_paf(PAF_FILE)
    log(f"Loaded {len(records):,} alignment blocks from {PAF_FILE}")

    depth = compute_depth(records, BIN_SIZE)
    flagged_bins, thresholds = flag_repetitive_bins(depth, DEPTH_MULTIPLIER)

    log(f"\nDepth thresholds (bin size = {BIN_SIZE:,} bp, flag if depth > {DEPTH_MULTIPLIER}x median):")
    for target, (med, thr) in thresholds.items():
        log(f"  {target}: median depth = {med:.1f}, threshold = {thr:.1f}")

    total_bins = sum(len(b) for b in depth.values())
    log(f"\nTotal reference bins covered: {total_bins:,}")
    log(f"Bins flagged as repetitive:   {len(flagged_bins):,} "
        f"({len(flagged_bins)/total_bins*100:.2f}% of covered regions)")

    kept_records = []
    dropped_records = []
    for r in records:
        if block_is_repetitive(r, flagged_bins, BIN_SIZE):
            dropped_records.append(r)
        else:
            kept_records.append(r)

    log(f"\nAlignment blocks kept (non-repetitive):    {len(kept_records):,}")
    log(f"Alignment blocks dropped (repetitive):     {len(dropped_records):,} "
        f"({len(dropped_records)/len(records)*100:.2f}% of all blocks)")

    with open(OUTPUT_PAF, "w") as out:
        for r in kept_records:
            out.write(r["raw"])

    log(f"\nSaved filtered alignments to: {OUTPUT_PAF}")
    log("Use this file in place of alignments.paf in your assembly pipeline")
    log("to exclude repeat-driven noise from scaffold placement and stitching.")
    log("=" * 74)

    with open(REPORT_FILE, "w") as f:
        f.write("\n".join(lines))
    print(f"\nReport saved to: {REPORT_FILE}")
























