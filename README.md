<img src="https://github.com/chasewnelson/trivolver/blob/master/trivolver_logo.png?raw=true" title="trivolver logo by Mitch Lin" alt="trivolver logo by Mitch Lin" align="middle">

**trivolver** is a Perl program for simulating non-reversible DNA sequence evolution on a fixed bifurcating tree using trinucleotide context. It relies on no external dependencies, facilitating maximum portability.

To test the simulation with the example data, run the following at the Unix command line or Mac Terminal:

`trivolver.pl --tree=test_tree.txt --seed_sequence=HsGgAncestor_chr9_4600001_4610000.fa --rate_matrix=mutation_CpGx20.txt --branch_unit=144740 > output.txt`

Find more [examples](#examples) below.

## <a name="description"></a>Description

New mutation data, including *de novo* mutations detected in whole genome sequencing of father/mother-child 'trios', can be used to empirically estimate context-dependent mutation rates. The most commonly used context is the trinucleotide, where the flanking nucleotide on either side of a position are considered (one 5', one 3'). Unfortunately, a tool is lacking for simulation of non-reversible (*i.e.*, asymmetric rates) DNA evolution on a fixed bifurcating (binary; fully resolved) gene tree, such as that produced by coalescence simulations. **trivolver** was made to fill this gap. The user must provide a bifurcating tree, a seed sequence, a 64 × 4 (trinucleotide × nucleotide) rate matrix, and a branch unit (see [options](#options)). **trivolver** outputs a <a target="_blank" href="https://github.com/samtools/hts-specs">Variant Call Format</a> (VCF) SNP report for all tree tips (leaves; taxa), as well as the history of mutations and motifs for each lineage (see [output](#output)).

## <a name="how-it-works"></a>How it Works

**trivolver** begins by placing the seed sequence on the root of the tree. It then proceeds to recursively traverse the tree from root to tip as an ordered binary tree structure. In other words, at each internal node, it processes the left child (and its left child, and so on) before processing the right child, where the left (first) descendant node is considered to be that which comes first using the alphabetic (not numeric) `sort()` function. Unique node ID's are assigned as the letter 'n' followed by an integer (*e.g.*, n5), where the root is considered n=1.

When a single branch (as opposed to a cluster) is encountered, **trivolver** begins to evolve the sequence along the branch. If the parent node of the branch also happens to be the root, the seed (ancestral) sequence is used as a starting point. Otherwise, a new and temporary copy of the seed sequence is created, into which the most recent allele at each mutated site is imputed. This saves computer memory, as each node need only inherit its mutational history, not an entire sequence, from its predeccesor.

After a new sequence is initiated on a branch, it begins to evolve under trinucleotide-context-dependent mutation rates specified by the user input. The mutation rate matrix in **trivolver** follows the format of the forward-time simulation <a target="_blank" href="https://messerlab.org/slim/">**SLiM**</a> (Haller and Messer 2019): the 4<sup>3</sup> = 64 rows correspond to the initial states of the 64 alphabetically-ordered trinucleotides, while the 4 columns correspond to the possible derived states of the central nucleotide. For example, the third column of the first row should contain the AAA➞AGA mutation rate. Mutation rates for identities (*e.g.*, AAA➞AAA) should be 0.

When evolution begins on a branch, trivolver first calculates the overall mutation rate of the sequence of length *L* by tallying the *L* - 2 trinucleotides in the sequence (*i.e.*, the first and last nucleotides, lacking trinucleotide context, are ignored). The overall mutation rate of the sequence (*u*) is then used to calculate a random exponentially-distributed waiting time to the next mutation as *g*<sub>w</sub> = -(1 / *u*) × ln(*x*<sub>1</sub>) generations, where *x* is a random number between 0 (inclusive) and 1 (exclusive) (Yang 2014). If the expected waiting time (*g*<sub>w</sub>) is less than the length of the branch, a mutation occurs. The sequence is then examined one trinucleotide at a time, until the cumulative mutation rate along the sequence exceeds a second randomly chosen value between 0 (inclusive) and *u* (exclusive). A mutation then occurs at the central position of the first trinucleotide for which this condition holds and for which the mutation rate is greater than 0. Because of this context-dependence, it is possible for highly mutable trinucleotides to 'erode' over time, and the overall mutation rate of the sequence is an emergent (rather than pre-specified) value depending on the rate matrix. Indeed, and the overall mutation rate may decreased (or increase) until an equilibrium triuncleotide composition is reached.

Happy trivolving!

## <a name="options"></a>Options

Call trivolver using the following options:

* `--tree` (**REQUIRED**): file containing a bifurcating evolutionary tree in newick format with branch lengths. NO NODE NAMES OR SUPPORT VALUES AT THIS TIME. Only the first encountered tree is used, *i.e.*, embarassingly parallel analyses must be executed one level up.
* `--seed_sequence` (**REQUIRED**): FASTA file containing starting (seed) sequence at tree root, to be evolved. Only the first sequence encountered is used.
* `--rate_matrix` (**REQUIRED**): file containing 64 × 4 tab-delimited trinucleotide rate matrix where rows and columns are in alphabetical order. For example, values in the first row correspond to AAA➞AAA, AAA➞ACA, AAA➞AGA, and AAA➞ATA.
* `--branch_unit` (**REQUIRED**): a scaling factor, equal to 2*N*<sub>0</sub> or 4*N*<sub>0</sub> in most coalescent simulations (*e.g.*, <a target="_blank" href="https://home.uchicago.edu/rhudson1/source/mksamples.html">**ms**</a>; Hudson 2002), which can be multiplied by a branch length to obtain the length of the lineage in generations. Branch lengths will be multiplied by this value and rounded up to the nearest integer to determine number of generations. For example, given the value 144740, a branch length of 0.228826612 in the phylogenetic tree would correspond to 144740 × 0.228826612 = 33,120.364 generations.
* `--random_seed` (OPTIONAL): integer with which to seed the random number generator. If not supplied, one will be chosen and reported. More specifically, we implement the checksum approach given in <a target="_blank" href="https://en.wikipedia.org/wiki/Programming_Perl">*Programming Perl*</a>: 

		srand(time ^ $$ ^ unpack "%32L*", `ps wwaxl | gzip`)
		
* `--tracked_motif` (OPTIONAL): a motif to track after each mutation. For example, to report the number of CpG sites over the course of a run, specify CG:

		--tracked_motif=CG

* `--track_mutations` (OPTIONAL): reports the mutation rate and count over time.
* `--verbose` (OPTIONAL): tell trivolver to tell you EVERYTHING that happens.

## <a name="examples"></a>EXAMPLES

### FORMAT:

	trivolver.pl --tree=<newick>.txt --seed_sequence=<seed>.fa --rate_matrix=<64x4>.txt \\ 
	--branch_unit=<#> --track_mutations --tracked_motif=<ACGT> --verbose > output.txt;
	
### ALL OPTIONS USED:

	trivolver.pl --tree=my_tree.txt --seed_sequence=my_ancestor.fa --rate_matrix=my_mutations.txt \\
	--branch_unit=144740 --random_seed=123456789 --tracked_motif=CG --track_mutations --verbose > my_output.txt
	
### TYPICAL USAGE (program decides random seed; not verbose):

	trivolver.pl --tree=my_tree.txt --seed_sequence=my_ancestor.fa --rate_matrix=my_mutations.txt \\
	--branch_unit=144740 --tracked_motif=CG --track_mutations > my_output.txt

### EVEN FEWER OPTIONS AND OUTPUT TO SCREEN:

	trivolver.pl --tree=my_tree.txt --seed_sequence=my_ancestor.fa --rate_matrix=my_mutations.txt \\
	--branch_unit=144740;

## <a name="output"></a>Output

Depending on the options specified, **trivolver** will output a mutation history file, a <a target="_blank" href="https://github.com/samtools/hts-specs">Variant Call Format</a> (VCF) SNP report for leaf (terminal) nodes, and/or full sequences for some or all leaf nodes.

## <a name="acknowledgments"></a>Acknowledgments
**trivolver** was written by with support from a Gerstner Scholars Fellowship from the Gerstner Family Foundation at the American Museum of Natural History to C.W.N. (2016-2019), and is maintained with support from the same. The logo image was designed by Mitch Lin (2019); copyright-free DNA helix obtained from Pixabay.

## <a name="contact"></a>Contact
If you have questions about **trivolver**, please click on the <a target="_blank" href="https://github.com/chasewnelson/trivolver/issues">Issues</a> tab at the top of this page and begin a new thread, so that others might benefit from the discussion.

Other correspondence should be addressed to Chase W. Nelson: 

* cnelson <**AT**> amnh <**DOT**> org

## <a name="references"></a>References

* Haller BC, Messer PW. 2019. SLiM 3: forward genetic simulations beyond the Wright–Fisher model. *Molecular Biology and Evolution* **36**:632–637.
* Hudson RR. 2002. Generating samples under a Wright-Fisher neutral model of genetic variation. *Bioinformatics* **18**:337–338.
* Yang Z. 2014. *Molecular Evolution: A Statistical Approach*. New York, NY: Oxford University Press.