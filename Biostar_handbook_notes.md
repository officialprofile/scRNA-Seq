**Assay** is an investigative procedure used to assess or measure the presence, amount, or function of some target (like a DNA fragment).

**“Forward”** and **“reverse”** are just labels. The choice of labeling is arbitrary and does not depend on any inherent property of the DNA. The forward strand is not “special.” Scientists decide which to call “forward” and which to call “reverse” when they first analyze the DNA of an organism. Even though the decision is arbitrary, it’s important to maintain consistency with that decision for the sake of clear communication.

When a process occurs in the expected direction, its directionality may be called **sense**; if it is going against the normal direction, its directionality may be called **anti-sense**. It is very important not to collate the concepts of **forward/reverse** with **sense/anti-sense**, as these concepts are completely unrelated. The sense/anti-sense is relative to a sequence’s direction; the sequence, in turn, may come from a forward or reverse strand.

The forward and reverse strands may also be denoted with different terms. For example, in some datasets you may find them labeled as + and -. They might also be called **top and bottom strands**, or even Watson and Crick strands.

The “official” definition of the term **gene** in the Sequence Ontology is: *A region (or regions) that includes all of the sequence elements necessary to encode a functional transcript. A gene may include regulatory regions, transcribed regions and other functional sequence regions.*

The region of the mRNA before the start codon (or the corresponding genomic region) is called the **5' UTR** (5 prime UTR), or untranslated region; the portion from the stop codon to the start of the poly-A tail is the 3' UTR (three prime UTR).

The genomic region just before the 5' UTR may contain patterns of nucleotides, called the **promoter**, that are used to position the molecular machinery that performs the transcription to RNA. For about 60% of human genes, the **promoter occurs in a CpG island**: a region of DNA where C and G are far more frequent than A and T. Other patterns in the DNA tell the cell when (how frequently and in what tissues) to transcribe the gene; that is, they regulate transcription. A pattern that increases the frequency of transcription operations is an **enhancer**, while one that decreases the rate is a **silencer**.

**Enhancers** help a gene turn on and regulate where, when, and how much the gene should be expressed for the cell to function properly. One enhancer can regulate many genes, and one gene can have many enhancers. Since DNA is folded and crumpled up in a cell, an enhancer can be far away from a gene it regulates in terms of distance measured along the DNA strand, but physically very close to the gene promoter. The best way to identify a sequence as an enhancer is to experimentally disrupt it and observe gene expression changes. Since this is cumbersome in practice, we use proxies like histone modifications (chemical decorations on the structural backbone upon which DNA is coiled) to identify enhancers. Silencers work similarly – they can be far from a gene, many-to-many relationships, hard to find – but cause a gene to be less, rather than more, expressed.

Two regions of DNA that evolved from the same sequence (through processes of duplication of genomic regions and separation of two species) are **homologous**, or **homologs** of one another.

**Homology is not a synonym of sequence similarity!** Homologous sequences are usually similar to one another, but similarity of sequences does not indicate homology.

Cut out the 3rd column of a tab-delimited text file and sort it to only show unique lines (i.e. remove duplicates):
```bash
cut -f 3 file.txt | sort -u
```

Count how many lines in a file contain the words ‘cat’ or ‘bat’ (-c option of grep counts lines):
```bash
grep -c '[bc]at' file.txt
```

Turn lower-case text into upper-case (using tr command to ‘transliterate’):
```bash
cat file.txt | tr 'a-z' 'A-Z'
```

Change all occurrences of ‘Chr1’ to ‘Chromosome 1’ and write changed output to a new file (using sed command):
```bash
cat file.txt | sed 's/Chr1/Chromosome 1/' > file2.txt
```

To find out the most frequent unique elements, we need to sort the output of uniq.
```bash
cat types.txt | sort | uniq -c | sort -rn | head
```

The **GFF/GTF/BED** formats are the so-called “interval” formats that retain only the coordinate positions for a region in a genome. Each is tab delimited and will contain information on the chromosomal coordinate, start, end, strand, value and other attributes, though the order of the columns will depend on the feature.

Up to 12 columns may be specified for BED formats. The coordinates start at 0 (that is the first valid coordinate). The GFF/GTF formats1 are 9 column tab-delimited formats. The coordinates start at value 1.

The **SAM/BAM** formats are so-called Sequence Alignment Maps. These files typically represent the results of aligning a FASTQ file to a reference FASTA file and describe the individual, pairwise alignments that were found. Different algorithms may create different alignments (and hence BAM files).

The **VCF (Variant Call Format)** describes the variation of alignments relative to a reference. Extensions: .vcf, vcf.gz, .bcf. Think of the VCF file as a file that captures the differences of for each of the sequences in the FASTQ file relative to the genome in the FASTA file.

**Entrez** is the interface into the data stored at NCBI. NCBI and Entrez (the interface into GenBank) are the best resources to obtain any sequence that humanity has ever found and stored. **Ensembl** is the interface into the data store at EBI.

An **accession number**, for example, AF086833, applies to the complete database record and remains stable even if updates/revisions are made to the record. The version number is formed by adding a dotted number such as .2 to the accession number to form: AF086833.2. This number is a unique identifier for the sequence data within its record.

As it happens, the **FASTA format** is not “officially” defined - even though it carries the majority of data information on living systems. Its origins go back to a software tool called Fasta written by David Lipman (a scientist that later became, and still is, the director of NCBI) and William R. Pearson of the University of Virginia. The tool itself has (to some extent) been superceded by the BLAST suite but the format (named after the tool itself) not only survived, but has become the de facto standard.

The numbers represent the error probabilities in a somewhat convoluted manner via the formula Error=10ˆ(-P/10) basically summarized as:
 - P=0 means 1/1 (100% error)
 - P=10 means 1/10 (10% error)
 - P=20 means 1/100 (1% error)
 - P=30 means 1/1000 (0.1% error)
 - P=40 means 1/10000 (0.01% error(

**When qualities look like “swearing”** in a comic book !"#$%&'()*+,-. it means that the data is of low quality. Imagine the sequencing instrument telling you: Hey, look at all this !$!!$*#@! data! When the quality scores are readable, perhaps jumbled letters, ABAFEFGGII your data has great quality. Imagine the sequencing instrument telling you: Hey, look at all this IAIAIAIAIA data!

Most confusingly, in this encoding the characters that indicate the highest base quality in Sanger will indicate one of the lowest qualities. It is easy to recognize this, though. When our FASTQ quality contains lower case letters abcdefghi it means that your data is in this older format.

Find Phred quality score by the Sanger encoding of letter A (prints 32)”
```bash
python -c 'print (ord("A")-33)'
```

Find Sanger encoded letter that corresponds to Phred quality score 32 (prints A):
```bash
python -c 'print (chr(32+33))'
```

Compute the error probability represented by the Pred quality score 32 (prints 0.00063095734448):
```bash
python -c 'from math import*; print (10**-(32/10.0))'
```

```bash
seqkit stat *.gz
```

### How to extract a subset of sequences with name/ID list file?
This is a frequently used manipulation. Let’s create a sample ID list file, which may also come from some other method like a mapping result.
```bash
seqkit sample --proportion 0.001 duplicated-reads.fq.gz | seqkit seq --name --only-id > id.
```
ID list file:
`head id.txt` shows:
```
SRR1972739.2996
SRR1972739.3044
SRR1972739.3562
```

Searching by ID list file:
```bash
seqkit grep --pattern-file id.txt duplicated-reads.fq.gz > duplicated-reads.subset.fq.gz
```

### How do I split FASTA sequences according to information in the header?
Overview of FASTA Headers:
`seqkit head -n 3 viral.2.protein.faa.gz | seqkit seq --name` generates:
```
gi|526245011|ref|YP_008320337.1| terminase small subunit [Paenibacillus phage phiIBB_Pl23]
gi|526245012|ref|YP_008320338.1| terminase large subunit [Paenibacillus phage phiIBB_Pl23]
gi|526245013|ref|YP_008320339.1| portal protein [Paenibacillus phage phiIBB_Pl23]
```

`seqkit split` can split FASTA/Q files according to ID, number of parts, size of each part, or sequence region. In this case, we’ll split according to sequence ID (species name) which can be specified by flag `--id-regexp`.

Default ID:
`seqkit head -n 3 viral.2.protein.faa.gz | seqkit seq --name --only-id` outputs:
```
gi|526245011|ref|YP_008320337.1|
gi|526245012|ref|YP_008320338.1|
gi|526245013|ref|YP_008320339.1|
```

New ID:
```bash
seqkit head -n 3 viral.2.protein.faa.gz | seqkit seq --name --only-id --id-regexp "\[(.+)\]"
```
prints:
```
Paenibacillus phage phiIBB_Pl23
Paenibacillus phage phiIBB_Pl23
Paenibacillus phage phiIBB_Pl23
```

Split:
```bash
seqkit split --by-id --id-regexp "\[(.+)\]" viral.1.protein.faa.gz
```

### How to concatenate two FASTA sequences in to one?

Data (not in same order):
```bash
cat 1.fa
```
>seq1
aaaaa
>seq2
ccccc
>seq3
ggggg
```bash
cat 2.fa
```
>seq3
TTTTT
>seq2
GGGGG
>seq1
CCCCC

Step 1. Convert FASTA to tab-delimited (3 columns, the 3rd column is blank (no quality for FASTA)) file:
```bash
seqkit fx2tab 1.fa > 1.fa.tsv
seqkit fx2tab 2.fa > 2.fa.tsv
```
```
cat -A 1.fa.tsv
seq1^Iaaaaa^I$
seq2^Iccccc^I$
seq3^Iggggg^I$
```
Step 2. Merge two table files:
```bash
csvtk join -H -t 1.fa.tsv 2.fa.tsv | cat -A
```
```
seq1^Iaaaaa^I^ICCCCC^I$
seq2^Iccccc^I^IGGGGG^I$
seq3^Iggggg^I^ITTTTT^I$
```
Step 3. Note that there are two TAB between the two sequences, so we can remove them to join the sequences
```bash
csvtk join -H -t 1.fa.tsv 2.fa.tsv | sed 's/\t\t//'
```
```
seq1
aaaaaCCCCC
seq2
cccccGGGGG
seq3
gggggTTTTT
```

Step 4. Convert tab-delimited file back to FASTA file:
```bash
csvtk join -H -t 1.fa.tsv 2.fa.tsv | sed 's/\t\t//' | seqkit tab2fx
```
```
>seq1
aaaaaCCCCC
>seq2
cccccGGGGG
>seq3
gggggTTTTT
```

All in one command:
```bash
csvtk join -H -t <(seqkit fx2tab 1.fa) <(seqkit fx2tab 2.fa) | sed 's/\t\t//' | seqkit tab2fx
```


Most genome browsers draw linear tracks that represent the forward strand of the genome from its 5’ (left) towards the 3’ (right) direction. 

Glyphs are a visual representation of some genomic features:
 - Horizontal intervals: directions, genes, alignments
 - Values over intervals: coverages, probabilities
 - Attributes at locations: mutations, deletions, junctions

An ontology is a controlled vocabulary of terms or concepts and a restricted set of relationships between those terms.
 - Where does the product exhibit its effect? -> Cellular Component (CC)
 - How does it work? -> Molecular Function (MF)
 - What is the purpose of the gene product? –> Biological Process (BP)

In the most common interpretation, functional analysis means interpreting data in the con- text of all known information gathered up present time.

ORA analysis does not account for the magnitude of expression levels. The gene is either in the list or not.

Gene set enrichment analysis” refers to the process of discovering the common characteristics potentially present in a list of genes. When these characteristics are GO terms, the process is called “functional enrichment

Reproducibility is not a binary classification: reproducible vs irreproducible. Reproducibility is a scale. Reproducibility is not a measure of scientific validity. It is a measure of transparency.

Paired-end sequencing is more expensive than single end sequencing but not radically so (perhaps 20% more).

We recommend that genomic variation and genome assembly analyses use paired-end sequencing as much as possible. On the downside, paired-end sequencing measures the same fragment twice, hence at the same genomic coverage, it will use half as many unique fragments. For analysis methods that use sequencing for quantification such as RNA-Seq and ChIP-Seq, this may be a disadvantage.

High-throughput sequencing data broadly fall into two categories:
• Short reads (300 bp or smaller).
• Long reads (>350 bp to several kilobases).
What goes into a high-throughput sequencing experiment may be familiar to many wet-lab scientists.

Do NOT use spaces or the following characters (? ( ) [ ] / \ = + < > : ; " ' , * ˆ | & .) in sample identifiers. These characters can cause problems with Illumina (and/or other technology) data pre-processing software and should be avoided.

The simple definition of **coverage** is

`C = number of sequenced bases / total genome size`

The only slightly confounding factor is that the number of sequenced bases come in units of “read lengths” rather than one base at a time.

Typically the coverage is expressed with the x symbol. For example, 10x indicates that on average, each base of the genome will be covered by 10 measurements.

For a sequencing instrument that produces variable reads lengths we have:

`C = SUM(Li) / G`

Where the sum goes over N the read count:

`Li = length of read i
G = size of the genome`

For an instrument with fixed read length L, that produces N reads the formula simplifies to:

`C = N * L / G`

At a coverage C the probability of a base not being sequenced is:

`P = exp(-C)`

**Example**
How many bases will not be sequenced when we cover a 20KB virus at 10x coverage.

`P = exp(-10)`

From command line:
```
python -c "import math; print (math.exp(-10))"
```
Prints: 4.53999297625e-05

To figure out how many bases will not be sequenced, we need to multiply this result by the genome length:
```
python -c "import math; print (math.exp(-10) * 20000)"
```
Prints: 0.90799859525

So at 10x coverage, we will miss about one base. When we cover the human genome at the same coverage:
```
python -c "import math; print (math.exp(-10) * 3E9)"
```
136199.789287

We will miss 136,199 bases. Note, though, that these gaps will be theoretically dispersed across the entire genome.

**FastQC** generates its reports by evaluating a small subset of the data and extrapolating those findings to the entirety of the dataset. Many of the metrics are only computed on the first 200,000 measurements then are being tracked through the rest of the data.

When FastQC runs, it generates “stoplight” icons for each analysis having “pass,” “warning,” and “error” symbols. Most of the time, these symbols are not meaningful. They were developed for only a particular class of samples and library preparation methods and just for certain types of instruments.

In very simple terms, current sequencing technology begins by breaking up long pieces of DNA into lots more short pieces of DNA. The resultant set of DNA is called a “library” and the short pieces are called “fragments”. Each of the
fragments in the library are then sequenced individually and in parallel. There are two ways of sequencing a fragment - either just from one end, or from both ends of a fragment. If only one end is sequenced, you get a single read. If your technology can sequence both ends, you get a “pair” of reads for each fragment. These “paired-end” reads are standard practice on Illumina instruments like the GAIIx, HiSeq, and MiSeq.

Now, for single-end reads, you need to make sure your read length (L) is shorter than your fragment length (F) or otherwise the sequence will run out of DNA to read! Typical Illumina fragment libraries would use F ~ 450bp, but this is variable. For paired-end reads, you want to make sure that F is long enough to fit two reads. This means you need F to be at least 2L. As L=100 or 150bp these days for most people, using F~450bp is fine, there is a still a safety margin in the middle.

However, some things have changed in the Illumina ecosystem this year. Firstly, read lengths are now moving to >150bp on the HiSeq (and have already been on the GAIIx) and to >250bp on the MiSeq, with possibilities of longer ones coming soon! This means that the standard library size F~450bp has become too small, and paired-end reads will overlap. Secondly, the new enzymatic Nextera library preparation system produces a wide spread of F sizes compared to the previous TruSeq system. With Nextera, we see F ranging from 100bp to 900bp in the same library. So some reads will overlap, and others won’t. It’s starting to get messy.

The whole point of paired-end reads is to get the benefit of longer reads without actually being able to sequence reads that long. A paired-end read (two reads of length L) from a fragment of length F, is a bit like a single-read of length F, except a bunch of bases in the middle of it are unknown, and how many of them there are is only roughly known (as libraries are only nominally of length F, each read will vary). This gives the reads a more extended context, and this mainly helps in de novo assembly and in aligning more reads unambiguously to a reference genome. However, many software tools will get confused if you give them overlapping pairs, and if we could overlap them and turn them into longer single-end reads, many tools will produce better results, and faster.

Suppose this is the full path.
```
FULL=/data/foo/genome.fasta.tar.gz
echo $FULL
```

To make it: genome.fasta.tar.gz
```
NAME=$(basename ${FULL})
echo $NAME
```
To make it: fasta.tar.gz
```
EXT1=${FULL#*.}
echo $EXT1
```
To get only the extension: gz
```
EXT2=${FULL##*.}
```

### How does Awk work?

An awk program works on a line by line basis and for each line of a file it attempts the following:
```
awk 'CONDITION { ACTIONS }'
```
For each line it awk tries to match the CONDITION, and if that condition matches it performs the ACTIONS. Do note the curly brackets and the quotes. Note that multiple conditions and actions may be present:
```
awk 'CONDITION1 { ACTIONS1 } CONDITION2 { ACTIONS2 } CONDITION3 { ACTIONS3 }'
```

In general, you always want to specify the splitting character when using awk:
```
awk -F '\t'
```
The flag above will set the field separator to be the “tab” character.

### What are some special awk variables I should know about?

When awk runs, many variables are set beyond the column variables $1, $2 etc.
 - $0 is the original line.
 - NF number of fields in the current line (number of columns that awk recognizes)
 - NR number of records, the number of lines processed (line number)
 - OFS output fields separator, the character placed between items when printed

As an example usage:
```
cat SGD_features.tab | awk '{ print NR, NF }' | head
```
This prints:
```
1 30
2 9
3 40
4 9
5 14
6 40
...
```
The first number is the line number; the second number is how many columns does awk think that the file has. Look what happens if you split by tab characters:
```
cat SGD_features.tab | awk -F '\t' '{ print NR, NF }' | head
```
you will get:
```
1 16
2 16
3 16
4 16
5 16
6 16
...
```

Make sure to use double quotes within patterns as to not conflict with the single quotes that the whole awk program is enclosed within.

Do not use **--color=always** when saving data in a new file for other downstream processing. When you use coloring, additional color-related data is added to the file, and that will affect downstream processing.

**grep** and **extended** grepegrep are tools that operate on one line at a time and can be used to identify lines with patterns. When sequences wrap over many lines we need dedicated tools like dreg1 or fuzznuc2 (Emboss Tools).

This will miss patterns that wrap around new lines.
```
cat KU182908.fa | grep AAAAAA
```
The dreg tool matches and reports the locations.
```
cat KU182908.fa | dreg -filter -pattern AAAAAA
```
Search a pattern with ambiguous N bases.
```
cat KU182908.fa | fuzznuc -filter -pattern 'AANAA
```

### Example regular expression searches:
Get a FASTQ dataset.
```
fastq-dump --split-files SRR519926
```
Find an ATG anchored at the start of the line
```
cat SRR519926_1.fastq | egrep "^ATG" --color=always | head
```
Find an ATG anchored at the end of the line
```
cat SRR519926_1.fastq | egrep "ATG\$" --color=always | head
```
Find TAATA or TATTA patterns, this is a range of characters
```
cat SRR519926_1.fastq | egrep "TA[A,T]TA" --color=always | head
```
Find TAAATA or TACCTA, these are groups of words
```
cat SRR519926_1.fastq | egrep "TA(AA|CC)TA" --color=always | head
```
Find TA followed by zero or or more A followed by TA
```
cat SRR519926_1.fastq | egrep "TA(A*)TA" --color=always | head
```
Find TA followed by one or or more A followed by TA
```
cat SRR519926_1.fastq | egrep "TA(A+)TA" --color=always | head
```
Find TA followed by two to five As followed by TA
```
cat SRR519926_1.fastq | egrep "TAA{2,5}TA" --color=always | head
```
Match Ilumina adaptors at the end of the reads
Match AGATCGG anywhere followed by any number of bases
```
cat SRR519926_1.fastq | egrep "AGATCGG.*" --color=always | head
```

When the sequences are very similar (nearly identical) the choice of scoring or even algorithms may not matter at all. The results will be robust - producing identical alignments across different algorithms and parameter settings.

When the sequences are dissimilar the choice of algorithm and scoring will typically have radical impacts on the results. In general, the more different the sequences, the more sensitive the alignment becomes to parameter choices.

The **CIGAR** string (Compact Idiosyncratic Gapped Alignment Report) is an alignment format used within the Sequence Alignment Map (SAM) files that form the backbone of most bioinformatics high-throughput analyses. For the same alignment from above:
```
ATGCAAATGACAAATAC
||||   |||.||.|
ATGC---TGATAACT--
```
the reported CIGAR format would be:
```
4M3D3M1X2M1X1M2D
```
We read out this “idiosyncratic” construct like so `4M + 3D + 3M + 1X + 2M + 1X + 1M + 2D`:
 - 4 matches followed by
 - 3 deletions,
 - 3 matches,
 - 1 mismatch,
 - 2 matches,
 - 1 mismatch,
 - 1 match,
 - 2 deletions.

The format above is in a so-called “Extended CIGAR,” meaning that it employs the X symbol for mismatches.

## Short read aligners

How are short reads different from long reads? Modern sequencing instruments brought about two significant changes:
 1. The read lengths became very short (50-300bp)
 2. The number of reads grew extremely high (hundreds of millions)

The rationale behind the new paradigm was to treat sequencing reads as small, independent measurements to be subsequently matched either against a known genome (called resequencing) or against one another (called assembly)

### How do short read aligners work?
First, let us make something abundantly clear: short read aligners are marvels of modern scientific software engineering. A well-optimized short read aligner can match over ten thousand sequences per second (!) against the 3 billion bases of the human genome. Their performance, implementation, and level of execution is nothing short of astounding.
But there is another side to the story. The rational expectation is that the published aligners would produce reasonably similar alignments and that the difference between various implementation would manifest primarily in their performance, hardware requirements, and small fluctuations in accuracy and precision.

In reality, the results that two aligners produce can be substantially different

While there is no shortage of scientific publications that claim to compare the accuracy and performance of various tools, most of these papers fail to capture the essential differences between the methods.

The minimum length for the read is algorithm- and implementation-dependent; it is rarely (if ever) mentioned in the documentation. Many tools stop working correctly when the read lengths drop under 30 bases. When studying small RNAs, for example, microRNAs, we have to look for tools that can align very short reads; the selection is surprisingly sparse.

The purpose of a mapper tool is locating a region in a genome, not producing an optimal alignment to that region.

Mapping
 - A mapping is a region where a read sequence is placed.
 - A mapping is regarded to be correct if it overlaps the true region.

Alignment
 - An alignment is the detailed placement of each base in a read.
 - An alignment is regarded to be correct if each base is placed correctly.

It is essential, however, to remember that the distinction between mapping and alignment does exist, and to recognize further that different studies have different requirements in regards to the use of these concepts.

For example, studies examining SNPs and variations in a genome would be primarily alignment-oriented. By contrast, studies using ChIP-Seq data would be essentially mapping-oriented.

So how do we pick the best aligner? There are good starting choices, such as bwa and bowtie2, that will perform well in all domains of application.

There is also a cultural and historical element to it. Scientists at the Broad Institute will likely use bwa because it was developed there.

In general, all short read aligners operate on the same principles:
 1. First build an “index” from a reference genome (this only needs to be done once).
 2. The reads in FASTA/FASTQ files are then aligned against the index created in step

Index building consists of pre-processing the reference genome so that the program can search it efficiently. Each program will build a different type of index; sometimes it may produce multiple files with odd-looking names or extensions. For this reason, it is best to place the reference genome in separate directories.

A SAM file encompasses all known information about the sample and its alignment; typically, we never look at the FastQ file again, since the SAM format contains all (well almost all) information that was also present in the FastQ measurements.

The first version of **bowtie aligner** was the first implementation of the Burrows-Wheeler algorithm for short read alignment, and with that it has opened a new era in processing high-throughput data.

The latest version of the algorithm bowtie2, is almost always preferred. In this book when we talk about the bowtie program we mean bowtie version 2.

The **SAM format** is a TAB-delimited, line-oriented text format consisting of a
 1. Header section, where each line contains some metadata
 2. Alignment section where each line provides information on an alignment

A **BAM file** is a binary, compressed (and almost always sorted) representation of the SAM information. Generally, BAM files are sorted, usually by the alignment coordinate and more rarely by the read names.
 - Sorting by coordinate allows fast query on the information by location.
 - Sorting by read name is required when the identically named read pairs need to be accessed quickly as read pairs will then be stored in adjacent lines.

**CRAM files** are conceptually similar to BAM files. CRAM files represent a more efficient version of the information in BAM files, where the gain is made by storing the reference separately. The essential difference between BAM and CRAM format is that a CRAM file can only be read if the reference sequence is present at a given location on the target computer. In contrast, the information in the BAM file is complete.

The downside of this approach is that losing the reference sequence means losing the data (that is, being unable to decode it).

Even practicing bioinformaticians frequently have misconceptions about SAM files. Among them, the most prevalent is the assumption that BAM files produced by different tools will have a comparable amount of information and the differences are mostly in the accuracy or performance characteristics of the alignments.

The **FLAG is the field of the SAM spec** where the designers of the format have strayed the furthest from what one might call common sense. It was a way to “cram” (no pun intended) as much information as possible into a single number. To do so, they came up with a contorted rationalization that no biologist will ever manage to remember, hence ensuring the longterm employment of bioinformaticians across the globe. Perhaps that was the plan all along, in which case I must admit, the FLAG is not such a bad thing after all.

**Column 5-6: Mapping Quality (MAPQ) and Compact Idiosyncratic Gapped Alignment Representation (CIGAR)**

These columns reflect the Mapping Quality (MAPQ) and the so-called Compact Idiosyncratic Gapped Alignment Representation (CIGAR).
```
samtools view SRR1972739.bwa.bam | cut -f 5,6 | head -5
```
to produce:
```
60 101M
60 101M
60 101M
60 101M
60 101M
```
The values in the MAPQ column here are all 60. This column was designed to indicate the likelihood of the alignment being placed incorrectly. It is the same Phred score that we encountered in the FASTQ files. And we read it the same way, 60/10 = 6 so the chance of seeing this alignment being wrong is 10ˆ-6 or 1/1,000,000 one in a million.

The CIGAR string is a different beast altogether. It is meant to represent the alignment via numbers followed by letters:
 - M match of mismatch
 - I insertion
 - D deletion
 - S soft clip
 - H hard clip
 - N skipping

These are also meant to be “readable”; the 18S83M says that 18 bases soft clipped and the next 83 are a match or mismatch.

First, let’s clarify some wording:
 1. Selecting means to keep alignments that match a condition.
 2. Filtering means to remove alignments that match a condition.

### How can I find out the depth of coverage?
There are several tools that can compute the number of reads that overlap each base. With
samtools you can do:

Compute the depth of coverage.
```
samtools depth SRR1972739.bwa.bam | head
```
This produces a file with name, coordinate and coverage on that coordinate:
```
AF086833 46 1
AF086833 47 1
AF086833 48 1
...
```
Run bamtools coverage:
```
bamtools coverage -in SRR1972739.bwa.bam | head
```
The output seems identical
```
AF086833 45 1
AF086833 46 1
AF086833 47 1
...
```
Except note how it is “one off”! The index starts at 0 rather than 1 in the latter case. These type of inconsistencies are extraordinarily common in bioinformatics and are the source of endless errors. You could ask for all positions to be shown with samtools and you can verify that it starts at 1
```
samtools depth -a SRR1972739.bwa.bam | head
```
it produces:
```
AF086833 1 0
AF086833 2 0
AF086833 3 0
```

Multiple mappings are caused primarily by repeats. They are less frequent given longer reads. If a read has multiple mappings, and all these mappings are almost entirely overlapping with each other; except the single-best optimal mapping, all the other mappings get mapping quality <Q3 and are ignored by most SNP/INDEL callers.

**Genomic variations** are typically categorized into different classes and are often denoted with a shortened acronym:
 - SNP, a single nucleotide polymorphism - A change of a single base.
 - INDEL, an insertion or a deletion - A single base added or removed.
 - SNV, a single nucleotide variant. A SNPs or an INDEL. The change is still single basepair long.
 - MNP, a multi-nucleotide polymorphism - Several consecutive SNPs in a block.
 - MNV, a multi-nucleotide variant - Multiples SNPs and INDELs as a block.
 - Short variations. Variations (MNVs) that are typically less than 50bp in length.
 - Large-scale variations. Variations that are larger 50bp.
 - Structural variants. Variations on the kilobase scale that involve (one or more) chromosomes.

SNPs are in a sense the holy grail of genomics. Their “popularity” perhaps is a manifestation of our hopes and wishful thinking. After all wouldn’t it be great if every disease or human condition could be distilled down to a few changes and mutations in the genome? We would need to find that one SNP, and with that, we solve the problem, collect our reward and move onto the next challenge.

And sometimes (rarely though) things do work out this way. Some individual mutations do indeed cause disease(s). The OMIM database (Online Mendelian Inheritance in Man) is an online catalog of human genes and genetic disorders collect this data.

In the context of genomic variations, the genotyping typically refers to identifying one, more than one, or all mutations of the sequence under investigation.

Specifically, the word “genotyping” seems to be used very liberally in publications and training materials. Most of the time, when reading these it is not clear what the authors mean when they state that they have “genotype” data. What they probably and usually mean is that they have established with some confidence that a variant is present.

A haplotype is a group of genes in an organism that is inherited together from a single parent.

In the context of genomic variations, a haplotype is a group of mutations that are inherited together. But just as with the “genotyping,” the word “haplotyping” may be used in a broader and narrower sense.

There are simulators available for the whole genome, exome, RNA-seq and chip-seq data simulation. Some tools run from command line; others have graphical user interfaces.
 - wgsim/dwgsim simulates obtaining sequencing reads from whole genomes.
 - msbar EMBOSS tool can simulate mutations in a single sequence.
 - biosed EMBOSS tool provides control on what mutations we introduce.
 - ReadSim is a sequencing read simulator to target long reads such as PacBio or Nanopore.
 - Art is the most sophisticated and up to date simulator that attempts to mimic the errors that sequencing instruments produce.
 - Metasim simulates reads from metagenomes (from a diverse taxonomical composition present in different abundances).
 - Polyester an RNA-Seq1 measurement simulator

It is also very likely that the human reference genome is not even a viable - meaning that if the sequence were “implanted” into a human cell, it would not function. Now, knowing that, think about the repercussions of this for a second. Scientists are using a non-existing, potentially non-viable human reference to characterize all other “real” human genomes

**Variant calling** is the process of identifying and cataloging the differences between the observed sequencing reads and a reference genome.

The **ploidy** is the number of complete sets of chromosomes in a cell and hence represent the number of possible alleles (variant forms) of the DNA of a single cell.

A **haplotype** (haploid genotype) is a group of alleles in an organism that are (or may be) inherited together from a single parent. Variants that occur on the same DNA molecule form a haplotype.

Typically a **VCF file** is a result of running a “variant caller” (or “SNP caller” as some people call it) on one or more BAM alignment files. The result of running the variant caller will be a VCF file that contains a column for each sample. Samples may be stored in different BAM files or may a single BAM file may include multiple samples tagged via read-groups.

The **BCF format** is a binary representation of the VCF format. It has a similar relationship as the SAM/BAM but with the exception that the BCF files are not sorted and indexed and cannot be queried the same way as typical BAM files do.

**Variant annotation** means predicting the effects of genetic variants (SNPs, insertions, deletions, copy number variations (CNV) or structural variations (SV)) on the function of genes, transcripts, and protein sequence, as well as regulatory regions.

**Gene isoforms** are mRNAs that are produced from the same locus but are different in their transcription start sites (TSSs), protein-coding DNA sequences (CDSs) and untranslated regions (UTRs), potentially altering gene function.

The sequencing coverage drops towards the ends of DNA molecules due to the lower chances of producing fragments that include the end of the molecule.

### What is the RPKM?
If we wanted to compare the number of reads mapped to one given transcript to another transcript of the same sample, we have to account for the fact that longer transcripts will produce more DNA fragments merely because they are longer.

FPKM is an extension of the already flawed concept of RPKM to paired-end reads. Whereas RPKM refers to reads, FPKM computes the same values over read pair fragments. Conceptually is even worse as the word “fragment” only adds another level of ambiguity to an already questionable concept

Trimmed mean of M values (TMM) normalization estimates sequencing depth after excluding genes for which the ratio of counts between a pair of experiments is too extreme or for which the average expression is too extreme. The edgeR software implements a TMM normalization.

The DESeq normalization method (implemented in the DEseq R package) estimates sequencing depth based on the count of the gene with the median count ratio across all genes.

## Introduction to ChIP-Seq analysis
ChIP-Seq stands for chromatin immunoprecipitation followed by sequencing. In a nutshell, the process consists of a laboratory protocol (abbreviated as ChIP) by the end of which the full DNA content of a cell is reduced to a much smaller subset of it. This subset of DNA is then sequenced (abbreviated as Seq) and is mapped against a known reference genome.

You will see how ChIP-Seq analysis is all about peaks - a term, that in our opinion is adopted widely yet lacks any exact definition.

The DNA fragments coming from a ChIP-Seq study are much shorter than a transcriptstudied in RNA-Seq analysis. Instead of counting measurements distributed over transcripts a single read may cover most of the fragment (or may even be longer than the original DNA fragment).

As a matter of fact, for ChIP seq results we need to compare to control (baseline) to verify that the protocol worked at all. Whereas in RNA-Seq we can immediately see the success of it by having reads cover transcript, in ChIP-Seq study no such certainty exists. Also when comparing peaks across conditions each peak has to first pass the comparison to a background, and then only those peaks that pass the first stage are compared to one another in a second stage.

Chromatin immunoprecipitation (ChIP) followed by high-throughput DNA sequencing (ChIP-seq) is a technique to map genome-wide transcription factor binding sites and histone-modification enriched regions.

Briefly, DNA bounding proteins and DNA (Chromatin) are cross-linked by formaldehyde and the chromatin is sheared by sonication into small fragments (typically 200 ~ 600 bp). The protein-DNA complex is immnuoprecipitated by an antibody specific to the protein of interest. Then the DNA is purified and made to a library for sequencing. After aligning the sequencing reads to a reference genome, the genomic regions with many reads enriched are where the protein of interest bounds. ChIP-seq is a critical method for dissecting the regulatory mechanisms of gene expression.

**HOMER** is a very nice tool for finding enriched peaks, doing motif analysis and many more. It requires making a tag directory first. ROSE22 is the tool from Richard Young’s lab.

For ChIP-seq, it is important to filter artifact regions that has abnomral high signals. From the website of Anshul Kundaje in Stanford University:

These regions are often found at specific types of repeats such as centromeres, telomeres and satellite repeats. It is especially important to remove these regions that computing measures of similarity such as Pearson correlation between genome-wide tracks that are especially affected by outliers.

## Introduction to metagenomics
Metagenomics is the study of genetic material recovered directly from environmental samples. The broad field may also be referred to as environmental genomics, ecogenomics or community genomic

From a bioinformatician’s point of view, metagenomics is the application of high-throughput genomic techniques to the study of communities of microbial organisms.

Traditionally genetics has focused on **vertical transfer** of genetic material: from parent to offspring. There is another mechanism of evolution, the so called **horizontal gene transfer** that allows some types of organisms to share genetic material across different species. Among single-celled organisms, it is perhaps the dominant form of genetic transfer.

Important: It is essential to recognize from the very beginning that most microorganisms cannot exist independently and their survival requires the presence of diverse communities of other microbes. Thus we need to keep in mind that organisms need to be treated in the context of their communities rather than individual species.

Current practices in metagenomics fall into the categories of:
 1. 16S rRNA sequencing
 2. Whole-metagenome sequencing

In a nutshell the 16S gene forms a peculiar part of a genome that has both very highly conserved regions (identical across just about all bacteria) as well as highly variable regions (hypervariable) that may be specific to a taxonomical level such as genus or species. You can use the conserved regions to isolate the variable regions that may be used to identify the bacteria. Since only a tiny section of the entire bacteria is genome will be sequenced, the 16S method is very “economical” requires a fraction of the coverage and produces far lower quantities of data.

Whereas during 16S sequencing each read corresponds to a bacterial cell, in the case of whole genome sequencing the length of the genome also matters. Longer bacteria will produce more DNA and will be seen in more sequencing reads. Another way to say this is that ten read counts of one bacteria vs. one read count of another bacteria does not mean that there were more of the first than the second.
