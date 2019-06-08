<img src="https://github.com/chasewnelson/trivolver/blob/master/trivolver_logo.png?raw=true" title="trivolver logo by Mitch Lin" alt="trivolver logo by Mitch Lin" align="middle">

**trivolver** is a Perl script for simulating non-reversible DNA sequence evolution on a fixed bifurcating tree using trinucleotide context.

To test the simulation with the example data, run the following at the Unix command line or Mac Terminal:

`trivolver.pl --tree=test_tree.txt --seed_sequence=HsGgAncestor_chr9_4600001_4610000.fa --rate_matrix=mutation_rates_FAST.txt --branch_unit=144740 > output.txt`

Call trivolver using the following options:

* `--tree` (**REQUIRED**): file containing a bifurcating evolutionary tree in newick format with branch lengths. NO NODE NAMES OR SUPPORT VALUES AT THIS TIME. Only the first encountered tree is used.
* `--seed_sequence` (**REQUIRED**): FASTA file containing starting (seed) sequence at tree root, to be evolved. Only the first sequence encountered is used.
* `--rate_matrix` (**REQUIRED**): file containing 64 x 4 tab-delimited trinucleotide rate matrix in alphabetical order. First row values correspond to: AAA>AAA, AAA>ACA, AAA>AGA, and AAA>ATA".
* `--branch_unit` (**REQUIRED**): a scaling factor, equal to 2*N*<sub>0</sub> or 4*N*<sub>0</sub> in most coalescent simulations, which can be multiplied by a branch length to obtain the length of the lineage in generations. Branch lengths will be multiplied by this value and rounded up to the nearest integer to determine number of generations. For example, given the value 144740, a branch length of 0.228826612 in the phylogenetic tree would correspond to 144740 × 0.228826612 = 33,120.364 generations.
* `--random_seed` (OPTIONAL): integer with which to seed the random number generator. If not supplied, one will be chosen and reported. More specifically, we implement the checksum approach given in Programming Perl: 

		srand(time ^ $$ ^ unpack "%32L*", `ps wwaxl | gzip`)
		
* `--tracked_motif` (OPTIONAL): a motif to track after each mutation. For example, to report the number of CpG sites over the course of a run, specify CG:

		--tracked_motif=CG

* `--track_mutations` (OPTIONAL): reports the mutation rate and count over time.
* `--verbose` (OPTIONAL): tell trivolver to tell you EVERYTHING that happens.

## EXAMPLES

### FORMAT:

	trivolver.pl --tree=<newick>.txt --seed_sequence=<seed>.fa --rate_matrix=<64x4>.txt \\ 
	--branch_unit=<#> --track_mutations --tracked_motif=<ACGT> --verbose > output.txt";
	
### EXAMPLE USING ALL OPTIONS:

	trivolver.pl --tree=my_tree.txt --seed_sequence=my_ancestor.fa --rate_matrix=my_mutations.txt \\
	--branch_unit=144740 --random_seed=123456789 --tracked_motif=CG --track_mutations --verbose > my_output.txt
	
### EXAMPLE OF TYPICAL USAGE (program decides random seed; not verbose):

	trivolver.pl --tree=my_tree.txt --seed_sequence=my_ancestor.fa --rate_matrix=my_mutations.txt \\
	--branch_unit=144740 --tracked_motif=CG --track_mutations > my_output.txt

### EXAMPLE WITH EVEN FEWER OPTIONS AND OUTPUT TO SCREEN:

	trivolver.pl --tree=my_tree.txt --seed_sequence=my_ancestor.fa --rate_matrix=my_mutations.txt \\
	--branch_unit=144740";

**trivolver** uses an ordered binary tree structure, where the left (first) descendant node is considered to be that which comes first using the alphabetic (not numeric) sort() function.

## <a name="acknowledgments"></a>Acknowledgments
**trivolver** was written with support from a Gerstner Scholars Fellowship from the Gerstner Family Foundation at the American Museum of Natural History (2016-2019), and is maintained with support from the same. The logo image was designed by Mitch Lin (2019); copyright-free DNA helix obtained from Pixabay.