<img src="https://github.com/chasewnelson/trevolver/blob/master/trevolver_logo.png?raw=true" title="Trevolver logo by Mitch Lin" alt="Trevolver logo by Mitch Lin" align="middle">

**Trevolver** is a Perl program for simulating non-reversible DNA sequence evolution on a fixed bifurcating tree using trinucleotide context. It relies on no external dependencies, facilitating maximum portability. Just download and run.

To test the simulation with the example data, execute the program at the Unix command line or Mac Terminal as follows:

### FORMAT:

	trevolver.pl --tree=<newick>.txt --seed_sequence=<seed>.fa \
	--rate_matrix=<64x4>.txt  --branch_unit=<#> --track_mutations \
	--tracked_motif=<ACGT> --verbose > <output_name>.txt

Find some real [examples](#examples) below. Check out our <a target="_blank" rel="noopener noreferrer" href="https://www.biorxiv.org/content/10.1101/672717v1">preprint on bioRxiv</a>.

## <a name="contents"></a>Contents

* [Description](#description)
* [How it Works](#how-it-works)
* [Options](#options)
* [Examples](#examples)
* [Output](#output)
	* [VCF Output](#vcf-output)
* [Troubleshooting](#troubleshooting)
* [Acknowledgments](#acknowledgments)
* [Citation](#citation)
* [Contact](#contact)
* [References](#references)

## <a name="description"></a>Description

New mutation data, including *de novo* mutations detected in whole genome sequencing of father/mother-child 'trios', can be used to empirically estimate context-dependent mutation rates. The most commonly used context is the trinucleotide, where the flanking nucleotide on either side of a position are considered (one 5', one 3'). Unfortunately, a tool is lacking for simulation of non-reversible (*i.e.*, asymmetric rates) DNA evolution on a fixed bifurcating (binary; fully resolved) gene tree, such as that produced by coalescence simulations. **Trevolver** was made to fill this gap. The user must provide a bifurcating tree, a seed sequence, a 64 × 4 (trinucleotide × nucleotide) rate matrix, and a branch unit (see [options](#options)). **Trevolver** outputs a <a target="_blank" rel="noopener noreferrer" href="https://github.com/samtools/hts-specs">Variant Call Format</a> (VCF) SNP report for all tree tips (leaves; taxa), as well as the history of mutations and motifs for each lineage (see [output](#output)).

## <a name="how-it-works"></a>How it Works

**Trevolver** begins by placing the seed sequence on the root of the tree. It then proceeds to recursively traverse the tree from root to tip as an ordered binary tree structure. In other words, at each internal node, it processes the left child (and its left child, and so on) before processing the right child, where the left (first) descendant node is considered to be that which comes first using the alphabetic (not numeric) `sort()` function. Unique node ID's are assigned as the letter 'n' followed by an integer (*e.g.*, n5), where the root is considered n=0.

When a single branch (as opposed to a cluster) is encountered, **Trevolver** begins to evolve the sequence along the branch. If the parent node of the branch also happens to be the root, the seed (ancestral) sequence is used as a starting point. Otherwise, a new and temporary copy of the seed sequence is created, into which the most recent allele at each mutated site is imputed. This saves computer memory, as each node need only inherit its mutational history, not an entire sequence, from its predeccesor.

After a new sequence is initiated on a branch, it begins to evolve under trinucleotide-context-dependent mutation rates specified by the user input. The mutation rate matrix in **Trevolver** follows the format of the forward-time simulation <a target="_blank" href="https://messerlab.org/slim/">**SLiM**</a> (Haller and Messer 2019): the 4<sup>3</sup> = 64 rows correspond to the initial states of the 64 alphabetically-ordered trinucleotides, while the 4 columns correspond to the possible derived states of the central nucleotide. For example, the third column of the first row should contain the AAA➞AGA mutation rate. Mutation rates for identities (*e.g.*, AAA➞AAA) should be 0.

When evolution begins on a branch, **Trevolver** first calculates the overall mutation rate of the sequence of length *L* by tallying the *L* - 2 trinucleotides in the sequence (*i.e.*, the first and last nucleotides, lacking trinucleotide context, are ignored). The overall mutation rate of the sequence (*u*) is then used to calculate a random exponentially-distributed waiting time to the next mutation as *g*<sub>w</sub> = -(1 / *u*) × ln(*x*<sub>1</sub>) generations, where *x* is a random number between 0 (inclusive) and 1 (exclusive) (Yang 2014). If the expected waiting time (*g*<sub>w</sub>) is less than the length of the branch, a mutation occurs. The sequence is then examined one trinucleotide at a time, until the cumulative mutation rate along the sequence exceeds a second randomly chosen value between 0 (inclusive) and *u* (exclusive). A mutation then occurs at the central position of the first trinucleotide for which this condition holds and for which the mutation rate is greater than 0. Because of this context-dependence, it is possible for highly mutable trinucleotides to 'erode' over time, and the overall mutation rate of the sequence is an emergent (rather than pre-specified) value depending on the rate matrix. Indeed, and the overall mutation rate may decreased (or increase) until an equilibrium triuncleotide composition is reached.

Happy Trevolving!

## <a name="options"></a>Options

Call **Trevolver** using the following options:

* `--tree` (**REQUIRED**): file containing a bifurcating evolutionary tree in newick format with branch lengths. **NO NODE NAMES OR SUPPORT VALUES AT THIS TIME**. Only the first encountered tree is used, *i.e.*, embarassingly parallel analyses must be executed one level up.
* `--seed_sequence` (**REQUIRED**): FASTA file containing starting (seed) sequence at tree root, to be evolved. Only the first sequence encountered is used.
* `--rate_matrix` (**REQUIRED**): file containing 64 × 4 tab-delimited trinucleotide rate matrix where rows and columns are in alphabetical order. For example, values in the first row correspond to AAA➞AAA, AAA➞ACA, AAA➞AGA, and AAA➞ATA.
* `--branch_unit` (**REQUIRED**): a scaling factor, equal to 2*N*<sub>0</sub> or 4*N*<sub>0</sub> in most coalescent simulations (*e.g.*, <a target="_blank" href="https://home.uchicago.edu/rhudson1/source/mksamples.html">**ms**</a>; Hudson 2002), which can be multiplied by a branch length to obtain the length of the lineage in generations. Branch lengths will be multiplied by this value and rounded up to the nearest integer to determine number of generations. For example, given the value 144740, a branch length of 0.228826612 in the phylogenetic tree would correspond to 144740 × 0.228826612 = 33,120.364 generations.
* `--random_seed` (*OPTIONAL*): integer with which to seed the random number generator. If not supplied, one will be chosen and reported. More specifically, we implement the checksum approach given in <a target="_blank" href="https://en.wikipedia.org/wiki/Programming_Perl">*Programming Perl*</a>: 

		srand(time ^ $$ ^ unpack "%32L*", `ps wwaxl | gzip`)
		
* `--tracked_motif` (*OPTIONAL*): a motif to track after each mutation. For example, to report the number of CpG sites over the course of a run, specify CG like so: `--tracked_motif=CG`.
* `--track_mutations` (*OPTIONAL*): reports the mutation rate and count over time.
* `--excluded_taxa` (*OPTIONAL*): path of **file containing a list of taxa names** (comma-separated) to be treated as outgroups, *i.e.*, excluded from the VCF file and the consensus sequence(s). This might be desirable if a small number of taxa represent outgroups, to which polymorphism in an ingroup is being compared.
* `--outgroups` (*OPTIONAL*): **number of outgroups** to be excluded for identification of the ingroup most recent common ancestor (MRCA) and for calculation of ingroup variant frequencies in the VCF file. Outgroups are considered to be the most deeply-branching terminal taxa (external nodes). For example, if `--outgroups=2` is specified, the fixed tree is navigated starting at the root. At each internal node, the branch containing the fewest terminal taxa is considered to contain the outgroup(s). Once the user-specified number of outgroups is identified, the most recent common ancestor (MRCA) node of the remaining (ingroup) taxa is identified and reported. If a set of non-arbitrary outgroup taxa does not exist for the user-specified number (*e.g.*, if `--outgroups=2` is called, but the two deepest splits contain three rather than two taxa), a warning is printed and no outgroups are used.
* `--burn_in` (*OPTIONAL*): number of generations to 'burn in' sequence evolution before initiating evolution at the tree root. Note that this is measured in absolute number of generations, **not** in branch units. Standard practice is 10*N* or 20*N* generations, where *N* is the (effective) population size.
* `--vcf_output` (*OPTIONAL*): name of a [VCF format output file](#vcf-output) to be generated in the working directory, unless a full path name is given. If not specified, a file will be printed in the working directory with a `.vcf` extension using the name of the tree file as a prefix.
* `--suppress_seed_seq` (*OPTIONAL*): suppress printing the ancestral (seed) sequence in the output. This might be desirable if the seed sequence is very large and its inclusion in the output consumes too much disk space.
* `--suppress_MRCA_seq` (*OPTIONAL*): prevents the MRCA (ingroup most recent common ancestor) sequence from being printed.
* `--suppress_consensus_seq` (*OPTIONAL*): prevents the consensus sequence (containing the `REF` allele, here defined as the major allele, at each site) from being printed.
* `--suppress_MUTATION` (*OPTIONAL*): prevents the `//MUTATION` data from being printed.
* `--verbose` (*OPTIONAL*): tell Trevolver to report EVERYTHING that happens. Not recommended except for development and debugging purposes.

## <a name="examples"></a>EXAMPLES

Example input and output files are available in the `EXAMPLE_INPUT` and `EXAMPLE_OUTPUT` directories at this GitHub page, where reproducible examples are numbered (*e.g.*, **output_example1.txt**). When the random seed has not been specified, exact results can be reproduced by using the same random number seed reported in the output file present in `EXAMPLE_OUTPUT`. Note that, if your input file(s) (*e.g.*, **tree\_6taxa.txt**) are not in the working directory (*i.e.*, where your Terminal is currently operating), you will need to specify the full path of the file name (*e.g.*, **/Users/ohta/Desktop/trevolver\_practice/tree\_6taxa.txt**). Also note that, in the examples below, a `\` is used simply to continue the previous command on the line.

### EXAMPLE 1: A SIMPLE SIMULATION

	trevolver.pl --tree=tree_6taxa.txt --seed_sequence=seed_sequence.fa \
	--rate_matrix=mutation_CpGx20.txt --vcf_output=example1.vcf --branch_unit=10000 \
	> example1.txt

### EXAMPLE 2: A SIMPLE SIMULATION WITH 1,000-GENERATION BURN-IN

	trevolver.pl --tree=tree_6taxa.txt --seed_sequence=seed_sequence.fa \
	--rate_matrix=mutation_CpGx20.txt --vcf_output=example2.vcf --branch_unit=10000 \
	--burn_in=1000 > example2.txt

### EXAMPLE 3: TYPICAL USAGE (program decides random seed; not verbose)

	trevolver.pl --tree=tree_6taxa.txt --seed_sequence=seed_sequence.fa \
	--rate_matrix=mutation_CpGx20.txt --branch_unit=144740 --track_mutations \
	--tracked_motif=CG --vcf_output=example3.vcf > example3.txt

### EXAMPLE 4: SIMULATION WITH TWO OUTGROUPS

	trevolver.pl --tree=tree_10taxa.txt --seed_sequence=seed_sequence.fa \
	--rate_matrix=mutation_CpGx20.txt --branch_unit=144740 --track_mutations \
	--tracked_motif=CG --vcf_output=example4.vcf --outgroups=2 > example4.txt
	
### EXAMPLE 5: EVOLVE A SINGLE SEQUENCE FOR 1 MILLION GENERATIONS

To simulate the evolution of a single sequence, provide a tree with only one taxon and branch length. For example, to simulate the evolution of a single sequence named "my_creature" for 1 million generations, the tree file would contain, simply, `(my_creature:1000000);`. Provide a scaling factor (`--branch_unit`) of 1, and you're good to go:

	trevolver.pl --tree=tree_1taxon.txt --seed_sequence=seed_sequence.fa \
	--rate_matrix=mutation_equal.txt --vcf_output=example5.vcf --branch_unit=1 \
	> example5.txt
	
### EXAMPLE 6: MANY OPTIONS USED

	trevolver.pl --tree=tree_7taxa.txt --seed_sequence=seed_sequence.fa \
	--rate_matrix=mutation_equal.txt --branch_unit=144740 --random_seed=123456789 \
	--tracked_motif=CG --track_mutations --vcf_output=example6.vcf --outgroups=2 \
	--burn_in=500 --suppress_seed_seq --suppress_consensus_seq --verbose \
	> example6.txt

### EXAMPLE 7: MINIMUM OPTIONS WITH OUTPUT TO SCREEN

This will automatically generate a VCF file named **tree\_7taxa\_trevolver.vcf** in the working directory.

	trevolver.pl --tree=tree_7taxa.txt --seed_sequence=seed_sequence.fa \
	--rate_matrix=mutation_CpGx20.txt --branch_unit=1447

## <a name="output"></a>Output

Depending on the options specified, **Trevolver** will output the following data:

### <a name="standard-output"></a>Standard Output

The beginning of the output will report:

* `COMMAND`: how Trevolver was called.
* `RANDOM_SEED`: the random seed used.
* `INPT_SEQUENCE`: the nucleotide sequence used as input.
* `SEED_SEQUENCE`: the nucleotide sequence used to seed the simulation. If there is no burn-in time, this is identical to the `INPT_SEQUENCE`.
* `MRCA_SEQUENCE`: nucleotide sequence of the ingroup MRCA (most recent common ancestor).
* `CONS_SEQUENCE`: consensus nucleotide sequence of the ingroup at the end of the simulation (*i.e.*, tree tips).
* `TREE`: the fixed tree on which the simulation took place.
* `MRCA_GENERATION`: the generation (from time 0 at the root) in which the most recent common ancestor (MRCA) of the ingroup lived.
* `MRCA_SUBTREE`: the subtree for which the MRCA is the root.
* `OUTGROUPS`: names of the outgroups, if applicable.
* `MRCA_NODE_ID`: node ID of the MRCA.
* `VCF_OUTPUT_FILE`: file name of the VCF output.

Additionally, following a brief **SUMMARY OF RESULTS**, the following flags indicate separate sections of more detailed output:

* `//MUTATION`: the following lines contain a full mutation history with three columns: taxon, site, and mutation. The mutation column data is in the format `[generation]-[ancestral allele]>[derived allele]`. For example, a C>T mutation which occurred in generation 1,988 would be listed as `1988-C>T`.
* `//TRACKED`: the following lines contain tracked mutation rates and/or motif data with five columns: lineage, generation, mutation\_rate, mutation\_count, and motif\_count.

### <a name="vcf-output"></a>VCF Output

The <a target="_blank" href="https://github.com/samtools/hts-specs">Variant Call Format</a> (VCF) output conforms to format VCFv4.1, such as used by the 1000 Genomes Project GRCh38/hg38 release, with some notable exceptions. For convenience and reproducibility, additional metadata headers (lines beginning with `##`) are used to indicate arguments used as **Trevolver** input, key results, and experiment specifications, including those described in [Standard Output](#standard-output) as well as:

* `burn_in_mutations`: total number of mutations during burn-in period.
* `total_mutations`: total number of mutations on all branches (excluding burn-in).
* `tree_length`: total branch length (generations).
* `experiment_length`: root-to-tip length (generations). It is currently required that all root-to-tip lengths, measured in generations, are equal.
* `simulation_time`: length of run (seconds).

Headers that are irrelevant (*e.g.*, non-single nucleotide variant descriptors) have been removed. The input, seed, consensus (ingroup), and MRCA (ingroup) sequences are printed in the header metadata for convenience unless `--suppress_input_seq`, `--suppress_seed_seq`, `--suppress_consensus_seq`, or `--suppress_MRCA_seq` are called. Additionally, numerous data types (many unknowable in real-life evolutionary analyses) have been added to the `INFO` column:

* `REF`/`REF_OG`: the consensus (major) allele of the ingroup/outgroup(s), which may or may not match the `AA` (ancestral allele).
* `AA`: ancestral allele of whole tree (the allele of the seed sequence).
* `MRCA`: Most Recent Common Ancestor (MRCA) allele with respect to the ingroup, if `--outgroups` is used.
* `MUTATIONS`/`MUTATIONS_OG`: all unique mutations that have occurred at this site, *e.g.*, `G>A`. Multiple mutations are comma-separated (*e.g.*, `G>A,A>G`) in **chronological order**, for convenience in downstream analyses. If outgroups or excluded taxa are specified, `MUTATIONS` refers to only the ingroup, while `MUTATIONS_OG` refers to only the outgroup(s).
* `GENERATIONS`/`GENERATIONS_OG`: time (generation) at which the unique mutation(s) occurred, comma-separated in the same order (**chronological**). If outgroups or excluded taxa are specified, `GENERATIONS ` refers to only the ingroup, while `GENERATIONS_OG` refers to only the outgroup(s).
* `TAXA`/`TAXA_OG`: number of taxa which share the unique mutation(s), comma-separated in the same order (**chronological**). If outgroups or excluded taxa are specified, `TAXA` refers to only the ingroup, while `TAXA_OG` refers to only the outgroup(s).
* `ARBITRARY_REF`: flag indicating there was a tie for the major/consensus/most common allele at this site. If the highest allele frequency is shared by two alleles, the one that comes first alphabetically is reported.
* `MULTIHIT`/`MULTIH_OG`: flags indicating a site has experienced more than one mutation, in either the same or a distinct lineage, in the ingroup/outgroup(s).
* `MULTIALLELIC`/`MULTIA_OG`: flags indicating a site has multiple minor (non-reference) alleles (*i.e.*, >2 alleles) in the ingroup/outgroup(s). All multiallelic sites are multihit, but the reverse is not true.
* `BACK_MUTATION`/`BACK_M_OG`: flag indicating a site has experienced back mutation, returning to a previous state/allele, in the ingroup/outgroup(s). All sites with back mutation have experienced multiple hits, but the reverse is not true.
* `RECURRENT_MUTATION`/`RECURRENT_M_OG`: flag indicating a site has experienced recurrent mutation, i.e., the same change occurring multiple times independently, in the ingroup/outgroup(s).
* `INVARIANT_ANCESTRAL`: flag indicating a site has no polymorphism in the ingroup leaves (extant taxa), and that the fixed state matches the ancestral allele (AA). Thus, even if a mutation occurred in the history of the site, the derived allele was lost via back mutation.
* `INVARIANT_DERIVED`: flag indicating a site has no polymorphism in the ingroup leaves (extant taxa), and that the fixed state matches a derived allele that resulted from mutation. This implies that all extant ingroup sequences are descended from a mutated ancestor.
* `NO_ANCESTRAL`: flag indicating that no ancestral (seed) alleles remain in the extant individuals of the ingroup. Note that this is compatible with the presence of polymorphism, so long as none of the alleles match the ancestral allele.
* `REF_OG`: the consensus (major) allele of outgroup(s), which may or may not match the AA.
* `REF_OG_COUNT`: count of outgroup allele matching the outgroup consensus.
* `REF_OG_AF`: frequency of outgroup allele matching the outgroup consensus.
* `ALLELES_OG`: comma-separated list of all alleles present in the outgroup(s).
* `ALLELE_COUNTS_OG`: comma-separated listed of all allele counts for alleles present in the outgroup(s), in the same order as `ALLELES_OG`.
* `OG_FIXED`: flag indicating a site is fixed for one allele (*i.e.*, has no variation) in the outgroup(s). This will be true by definition if there is only only outgroup. Note that frequences of outgroup alleles can be retrieved from the `ALLELE_COUNTS_OG` data.
* `OG_DIVERGED`: flag indicating one or more outgroup alleles *differs from* one or more ingroup alleles. The outgroups may be diverged from the ingroup but not fixed.
* `OG_SHARE`: flag indicating one or more outgroup alleles *matches* one or more ingroup alleles. The outgroups may share alleles with the ingroup but not be fixed.

## <a name="troubleshooting"></a>Troubleshooting

If you have questions about **Trevolver**, please click on the <a target="_blank" href="https://github.com/chasewnelson/trevolver/issues">Issues</a> tab at the top of this page and begin a new thread, so that others might benefit from the discussion. Common questions will be addressed in this section.

## <a name="acknowledgments"></a>Acknowledgments
**Trevolver** was written with support from a Gerstner Scholars Fellowship from the Gerstner Family Foundation at the American Museum of Natural History to C.W.N. (2016-2019), and is maintained with support from the same. The logo image was designed by Mitch Lin (2019); copyright-free DNA helix obtained from Pixabay. Thanks to Reed A. Cartwright, Michael Dean, Dan Graur, Ming-Hsueh Lin, Lisa Mirabello, Sergios Orestis-Kolokotronis, Michael Tessler, and Meredith Yeager for discussion.

## <a name="citation"></a>Citation

When using this software, please refer to and cite:

>Nelson CW, Fu Y, Li W-H. Trevolver: simulating non-reversible DNA sequence evolution in trinucleotide context on a bifurcating tree. Submitted to *Bioinformatics*. bioRxiv doi: <a target="_blank" rel="noopener noreferrer" href="https://www.biorxiv.org/content/10.1101/672717v1">https://doi.org/10.1101/672717</a>

and this page:

>https://github.com/chasewnelson/trevolver

## <a name="contact"></a>Contact
If you have questions about **Trevolver**, please click on the <a target="_blank" href="https://github.com/chasewnelson/trevolver/issues">Issues</a> tab at the top of this page and begin a new thread, so that others might benefit from the discussion.

Other correspondence should be addressed to Chase W. Nelson: 

* cnelson <**AT**> amnh <**DOT**> org

## <a name="references"></a>References

* Haller BC, Messer PW. 2019. <a target="_blank" href="https://academic.oup.com/mbe/article/36/3/632/5229931">SLiM 3: forward genetic simulations beyond the Wright–Fisher model</a>. *Molecular Biology and Evolution* **36**:632–637.
* Hudson RR. 2002. <a target="_blank" href="https://academic.oup.com/bioinformatics/article/18/2/337/225783">Generating samples under a Wright-Fisher neutral model of genetic variation</a>. *Bioinformatics* **18**:337–338.
* Yang Z. 2014. <a target="_blank" href="https://www.oxfordscholarship.com/view/10.1093/acprof:oso/9780199602605.001.0001/acprof-9780199602605">*Molecular Evolution: A Statistical Approach*</a>. New York, NY: Oxford University Press.