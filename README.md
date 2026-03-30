# Alpaqa: Assembly-Level Profiling And Quality Assessment
<img src="https://github.com/user-attachments/assets/8ff14b77-4ac3-459a-89e0-be984b87b6cc" align="right" width="170" style="margin-left: 15px;">

Alpaqa analyzes base-level quality scores in bacterial genome assemblies generated from Oxford Nanopore sequencing data.

The tool masks noisy regions and quantifies the number of Low-Quality Bases (LQBs) per megabase, a metric we found to correlate with overall base-level assembly accuracy. In addition, Alpaqa performs binomial tests to identify 4-, 5-, and 6-mer motifs that are significantly associated with LQBs. These motifs often represent systematic, error-prone sequence contexts and may correspond to targets of DNA modification systems.

Alpaqa operates on polished assemblies generated with Dorado polish or Medaka2. By detecting systematic error signatures, Alpaqa allows users to flag assemblies that may appear complete but contain motif-specific inaccuracies. This is particularly valuable in applications requiring high base-level accuracy, such as outbreak investigation and transmission analysis.

An automated bacterial genome assembly pipeline for ONT data which includes Medaka2 polishing, quality assessment with alpaqa, and masking of low-quality bases is available at [boap](https://github.com/MBiggel/boap).

## Installation

### Via pip
Alpaqa can be installed directly from PyPI:

```bash
pip install alpaqa-bio
```

### Via uv (recommended)
If you use [uv](https://github.com/astral-sh/uv), you can run it without installation:

```bash
uvx --from alpaqa-bio alpaqa --help
```

Or install it:

```bash
uv tool install alpaqa-bio
```

## Generating input files for alpaqa
Alpaqa relies on fastq assembly files with phred quality scores generated with [dorado polish](https://github.com/nanoporetech/dorado) or [medaka2](https://github.com/nanoporetech/medaka) (-q flag). The tool has been tested with data generated with SUP@v5.0, SUP@v5.2, HAC@v5.0, and HAC@v5.2 basecalling models. SUP data is strongly recommended.

#### Using medaka2

```bash
READS=reads.fastq.gz
ASSEMBLY=draft_assembly.fasta
THREADS=16

# medaka2
medaka_consensus -i $READS -d $ASSEMBLY -o polished_assembly -t $THREADS --bacteria -q
```
#### Using dorado polish

```bash
READS=reads.fastq.gz
ASSEMBLY=draft_assembly.fasta
THREADS=16

dorado aligner $ASSEMBLY $READS | samtools sort -@ $THREADS -o aligned.bam
samtools index aligned.bam
dorado polish aligned.bam $ASSEMBLY -t $THREADS -o polished_assembly --bacteria -q
```

## Alpaqa usage

Once installed, run `alpaqa` by passing your FASTQ assembly files:

```bash
alpaqa -i *.fastq --threads 16
```

## Options

```text
usage: alpaqa -i assembly.fastq -o output.tsv --threads 16 [options]

ALPAQA v0.1.1

options:
  -h, --help            show this help message and exit
  -i INPUT [INPUT ...], --input INPUT [INPUT ...]
                        Input fastq assembly file(s).
  -o OUTPUT, --output OUTPUT
                        Output TSV filename. (default: alpaqa_report.tsv)
  --report              Generate detailed reports on kmer analysis.
  -t THREADS, --threads THREADS
                        Number of threads. (default: 1)
  -v, --version         show program's version number and exit

Advanced:
  --all-contigs         Analyze ALL contigs > min_contig_len (Default is longest contig only).
  --no-mask             Disable masking of LQB dense regions.
  --lqb-threshold LQB_THRESHOLD
                        Phred Q-score threshold (inclusive). (default: 5)
  --min-contig-len MIN_CONTIG_LEN
                        Minimum contig length to include. (default: 50000)
  --window-size WINDOW_SIZE
                        Window size (bp) for scanning LQB dense regions. (default: 5000)
  --density-multiplier DENSITY_MULTIPLIER
                        Multiplier for baseline density to trigger masking. (default: 5.0)
  --mask-floor MASK_FLOOR
                        Absolute minimum LQB density (0.001 = 0.1%). (default: 0.001)
```

## Output

Alpaqa generates a tab-separated (TSV) report. Each row represents one input assembly file with the following fields:

| Field                | Description                                                                                                                                                                                                                    |
| :------------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Filename**         | Name of the input FASTQ file.                                                                                                                                                                                                  |
| **AvgQ**             | Average Phred quality score across all bases in the file.                                                                                                                                                                      |
| **LQB_raw/Mbp**      | Low-Quality Bases per Megabase (Raw): Total count of bases with a Phred score between 1 and 5 in the entire file, normalized per million bases.                                                                                |
| **LQB/Mbp**          | **Low-Quality Bases per Megabase (Filtered):** The density of low-quality bases in the longst contig after excluding (masking) unreliable LQB-dense regions.                                                                   |
| **MaskThresh**       | Dynamic Masking Threshold: The density threshold that triggered masking. Calculated as 5x the baseline density, but is at least 0.1% (0.001) to prevent over-masking assemblies.                                               |
| **Masked_Bases**           | Total number of bases removed from the  analysis.                                                                                                                                                                              |
| **Contigs_Analyzed** | Shown as X//Y, where X is the number of contigs analyzed (default: longest contig only) and Y is the total number of contigs in the file.                                                                                      |
| **Bases_Analyzed**   | Total number of bases analyzed.                                                                                                                                                                                                |
| **Sig4m / 5m / 6m**  | **Significant Motifs:** The top 3 DNA patterns (4, 5, and 6-mers) most significantly linked to quality drops (Q1-5) based on binomial testing. Includes the motif and the percentage of its occurrences that were low quality. |

## Interpretation
For assemblies generated with SUP@v5.2 data, following thresholds were established:
* **<5 LQBs/Mbp**: Highly accurate. These assemblies typically yield identical or near-identical cgMLST profiles when compared to Illumina-polished references.
* **5 to 10 LQBs/Mbp**: Potentially reliable but should be interpreted with caution. Masking low-quality bases with fastq2a.py (see below) is recommended for downstream cgMLST analyses.
* **>10 LQBs/Mbp**: Usually unsuitable for high-resolution genotyping due to excessive base-level errors.

Beyond systematic sequencing errors, elevated LQB counts may indicate sample contamination, insufficient sequencing depth, or low read quality, which can be further investigated using tools such as [nanoq](https://github.com/esteinig/nanoq) and [CheckM](https://github.com/Ecogenomics/CheckM).

![LQB_vs_errors](https://github.com/user-attachments/assets/62f0f612-9a29-48c0-b8d7-c59f0e229238)

## Masking low-quality bases for downstream genotyping
The helper script `fastq2a.py` can be used to convert FASTQ assembly files into FASTA format while masking low-quality bases. Replacing bases below a user-defined Q-score threshold with "N" helps maintain robustness in downstream bacterial cgMLST analyses. Tools such as Ridom SeqSphere exclude alleles containing these ambiguous bases to prevent false distance calculations. However, the final count of remaining target loci should be monitored to ensure sufficient resolution for high-resolution typing. For assemblies with ~5 to 10 LQBs/Mpb, masking bases with qscores ≤10 provides a good balance between accuracy and genomic resolution.

Note: if you installed via pip/uv, you can run this helper script using:
```bash
python -m alpaqa.fastq2a -i assembly.fastq -o assembly.masked.fasta -q 10
```

## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
