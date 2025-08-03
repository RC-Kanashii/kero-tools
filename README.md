# KERO Tools

This repository contains a suite of command-line tools for manipulating KERO files, a format for storing k-mers and associated data efficiently. `kero-tools` is a single executable that provides several sub-commands for different file operations.

## Install

To build the project, you will need a C++ compiler and CMake (version 3.19 or newer).

```bash
git clone --recursive <your-repository-url>/kero-tools.git
cd kero-tools
mkdir build && cd build
cmake ..
make -j 16
```

The executable `kero-tools` will be located in the `kero-tools/bin/<BuildType>/` directory (e.g., `kero-tools/bin/Release/`).

## Commands

### `kero-tools instr`

Convert a text file (containing sequences or k-mers) into the KERO format.

The input file can contain one long sequence per line, or one k-mer per line. If data-size is non-zero, the tool expects data values (integers) to follow the sequence on the same line.

**Example with data (k=3):**

```
ACGTACGT 10,12,15,18,20,5
TTGCA 8,9,2
```

**Parameters:**

- `-i, --infile <input.txt>` [required]: Input text file with sequences or k-mers.
- `-o, --outfile <output.kero>` [required]: Output KERO file.
- `-k, --kmer-size <int>` [required]: The size of the k-mers.
- `-d, --data-size <int>`: Size of data associated with each k-mer in bytes (Default: 0).
- `-m, --max-kmer-seq <int>`: Maximum number of k-mers to store per sequence block. Longer sequences will be split (Default: 255).
- `--delimiter <char>`: Delimiter between a sequence and its data (Default: ' ').
- `--data-delimiter <char>`: Delimiter between data values for different k-mers (Default: ',').

**Usage:**

```bash
# Convert a file of sequences to KERO format for k=31
kero-tools instr -i sequences.txt -o sequences.kero -k 31

# Convert a file of k-mers with 1-byte counts
kero-tools instr -i counts.txt -o counts.kero -k 21 -d 1
```

### `kero-tools outstr`

Read a KERO file and print its k-mer sequences and data to standard output as strings.

**Parameters:**

- `-i, --infile <input.kero>` [required]: The KERO file to read.
- `-c, --reverse-complement`: Print the lexicographically canonical k-mer (between the forward and reverse-complement strand).

**Usage:**

```bash
kero-tools outstr -i my_kmers.kero > my_kmers.txt
```

### `kero-tools query`

Query a KERO file to find the count or associated data for one or more k-mers. This tool requires an indexed KERO file (created with the bucket and compact tools).

**Parameters:**

- `-i, --infile <input.kero>` [required]: The indexed KERO file to query.
- `-o, --output <output.txt>` [required]: The file to write query results to.
- `-s, --kmer <string>`: A single k-mer or a longer sequence from which to query all k-mers.
- `-S, --kmer-file <filename>`: A file containing k-mers or sequences to query (one per line).
- `--kmer-per-line`: Flag to indicate that the --kmer-file contains one k-mer per line, not long sequences.
- `--no-index`: Force a sequential scan of the file instead of using the index (very slow).

**Usage:**

```bash
# Query a single k-mer
kero-tools query -i db.kero -s AGCT... -o result.txt

# Query all k-mers from a file of sequences
kero-tools query -i db.kero -S reads.fasta -o results.txt
```

### `kero-tools bucket`

Partition k-mers from raw sections into minimizer-based sections. All k-mers that share the same smallest m-mer (minimizer) will be grouped into the same section, which is a prerequisite for compaction.

**Parameters:**

- `-i, --infile <input.kero>` [required]: The KERO file to bucketize.
- `-o, --outfile <output.kero>` [required]: The output file containing minimizer sections.
- `-m, --minimizer-size <int>` [required]: The size of the minimizer (m).

**Usage:**

```bash
kero-tools bucket -i raw_kmers.kero -o bucketed.kero -m 15
```

### `kero-tools compact`

Compact overlapping k-mers within each minimizer section into longer "super-k-mers". This significantly reduces file size and is essential for efficient querying.

**Parameters:**

- `-i, --infile <input.kero>` [required]: A bucketed KERO file.
- `-o, --outfile <output.kero>` [required]: The output compacted file.
- `-s, --sorted`: Sort super-k-mers to allow binary search, which may slightly reduce compaction ratio.

**Usage:**

```bash
kero-tools compact -i bucketed.kero -o compacted_db.kero
```

### `kero-tools disjoin`

The reverse of compact. It splits every super-k-mer block containing n k-mers back into n separate blocks, each with a single k-mer.

**Parameters:**

- `-i, --infile <input.kero>` [required]: The compacted KERO file.
- `-o, --outfile <output.kero>` [required]: The output file with disjoint k-mers.

**Usage:**

```bash
kero-tools disjoin -i compacted.kero -o disjoint.kero
```

### `kero-tools merge`

Merge multiple KERO files into a single output file. The files must share the same nucleotide encoding.

**Parameters:**

- `-i, --inputs <file1.kero> <file2.kero> ...`: A space-separated list of input files.
- `-f, --input-filelist <list.txt>`: A file containing a list of input file paths, one per line.
- `-o, --outfile <output.kero>` [required]: The name for the merged output file.

**Usage:**

```bash
kero-tools merge -i part1.kero part2.kero -o merged.kero
```

### `kero-tools split`

Split a single KERO file into multiple files, with one new file created for each section in the original.

**Parameters:**

- `-i, --infile <input.kero>` [required]: The KERO file to split.
- `-o, --outdir <path>`: The directory where output files will be saved (Default: ./).

**Usage:**

```bash
kero-tools split -i my_file.kero -o my_file_sections/
```

### `kero-tools translate`

Rewrite a KERO file with a new nucleotide encoding scheme.

**Parameters:**

- `-i, --infile <input.kero>` [required]: The file to translate.
- `-o, --outfile <output.kero>` [required]: The translated output file.
- `-e, --encoding <string>` [required]: A 4-character string defining the new encoding. For example, AGTC sets A=0, G=1, T=2, C=3.

**Usage:**

```bash
kero-tools translate -i file.kero -o translated.kero -e AGTC
```

### `kero-tools data-rm`

Remove all associated data (e.g., counts) from the k-mers in a file, setting the data size to 0.

**Parameters:**

- `-i, --infile <input.kero>` [required]: The KERO file with data.
- `-o, --outfile <output.kero>` [required]: The output file without data.

**Usage:**

```bash
kero-tools data-rm -i with_counts.kero -o no_counts.kero
```

### `kero-tools validate`

Read a KERO file and check for structural integrity and corruption. The program will exit with an error if a problem is found.

**Parameters:**

- `-i, --infile <input.kero>` [required]: The file to validate.
- `-o, --outfile <output.txt>` [required]: The file to write the validation report to.
- `-v, --verbose`: Print detailed information about the file's structure.
- `--index-only`: Restrict validation to the index sections only.

**Usage:**

```bash
kero-tools validate -i my_file.kero -o validation_report.txt -v
```