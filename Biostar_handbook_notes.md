**Assay** is an investigative procedure used to assess or measure the presence, amount, or function of some target (like a DNA fragment).

**“Forward”** and **“reverse”** are just labels. The choice of labeling is arbitrary and does not depend on any inherent property of the DNA. The forward strand is not “special.” Scientists decide which to call “forward” and which to call “reverse” when they first analyze the DNA of an organism. Even though the decision is arbitrary, it’s important to maintain consistency with that decision for the sake of clear communication.

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
