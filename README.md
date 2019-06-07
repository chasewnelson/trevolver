<img src="https://github.com/chasewnelson/trivolver/blob/master/trivolver_logo.png?raw=true" title="trivolver logo by Mitch Lin" alt="trivolver logo by Mitch Lin" align="middle">

**trivolver** is a Perl script for simulating non-reversible DNA sequence evolution on a fixed bifurcating tree using trinucleotide context.

To test the simulation with the example data, run the following at the Unix command line or Mac Terminal:

`trivolver.pl --tree=test_tree.txt --seed_sequence=HsGgAncestor_chr9_4600001_4610000.fa --rate_matrix=mutation_rates_FAST.txt --branch_unit=144740 > output.txt`

Here, `$branch_unit` is a scaling factor, equal to 2*N*<sub>0</sub> or 4*N*<sub>0</sub> in most coalescent simulations, which can be multiplied by a branch length to obtain the length of the lineage in generations. For example, given the value 144740, a branch length of 0.228826612 in the phylogenetic tree would correspond to 144740 × 0.228826612 = 33,120.364 generations.

## <a name="acknowledgments"></a>Acknowledgments
**trivolver** was written with support from a Gerstner Scholars Fellowship from the Gerstner Family Foundation at the American Museum of Natural History (2016-2019), and is maintained with support from the same. The logo image was designed by Mitch Lin (2019); copyright-free DNA helix obtained from Pixabay.