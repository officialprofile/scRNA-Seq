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
