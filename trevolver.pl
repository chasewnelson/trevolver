#! /usr/bin/perl

# Perl script to simulate sequence evolution on a BIFURCATING TREE provided in Newick
#	format with a user-provided TRINUCLEOTIDE (64 x 4) rate matrix, such as that described
#	in SLiM. This allows non-reversible context-dependent evolution with back mutation.
# OUTPUTS: tree with node labels, accumulated mutations, and extant sequences.

#########################################################################################
# EXAMPLE CALL:
#########################################################################################
# trevolver.pl --tree=<treefile>.newick --seed_sequence=<sequence>.fasta --rate_matrix=<k x 4 table>.txt --branch_unit=<4N0 integer> --random_seed=<optional integer>
#########################################################################################

# Copyright (C) 2019 Chase W. Nelson
# DATE CREATED: June 2019
# AUTHOR: Chase W. Nelson
# CONTACT1: cnelson@amnh.org
# AFFILIATION: Sackler Institute for Comparative Genomics, American Museum of Natural History, 
#	New York, NY 10024, USA

# ACKNOWLEDGMENTS: written by C.W.N. with support from a Gerstner Scholars Fellowship from
# 	the Gerstner Family Foundation at the American Museum of Natural History, New York.

use strict;
use List::Util qw(max sum);
use Getopt::Long;

STDOUT->autoflush(1);

# Get the time
my $time1 = time;
my $local_time1 = localtime;

my @commands = @ARGV;


#########################################################################################
# INITIALIZE INPUT VARIABLES
my $tree; # file containing bifurcating evolutionary tree with branch lengths
my $seed_sequence; # file containing starting (seed) sequence at tree root, to be evolved
my $rate_matrix; # file containing 64 x 4 tab-delimited rate matrix in alphabetical order. First row values for: AAA>AAA\tAAA>ACA\tAAA>AGA\tAAA>ATA\n
my $branch_unit; # branch lengths will be multiplied by this value and rounded up to the nearest integer to determine number of generations
my $random_seed; # integer with which to seed the random number generator
my $tracked_motif;
my $track_mutations;
my $vcf_output;
my $excluded_taxa;
my $suppress_ancestral_seq;
my $verbose;

# Get user input, if given. If a Boolean argument is passed, its value is 1; else undef
GetOptions( "tree=s" => \$tree,
			"seed_sequence=s" => \$seed_sequence,
			"rate_matrix=s" => \$rate_matrix,
			"branch_unit=f" => \$branch_unit,
			"random_seed=i" => \$random_seed,
			"tracked_motif=s" => \$tracked_motif,
			"track_mutations" => \$track_mutations,
			"vcf_output=s" => \$vcf_output,
			"excluded_taxa=s" => \$excluded_taxa,
			"suppress_ancestral_seq" => \$suppress_ancestral_seq,
			"verbose" => \$verbose
			)
			
			or print_usage_message("### WARNING: Error in command line arguments (option misspelled?). trevolver terminated.");

#print "\nvcf_output=$vcf_output\n"; # the value of a called flag is 1

unless(-f "$tree") {
	my $specific_warning = "### WARNING: A valid --tree option must be provided.";
	print_usage_message($specific_warning);
}

unless(-f "$seed_sequence") {
	my $specific_warning = "### WARNING: A valid --seed_sequence option must be provided.";
	print_usage_message($specific_warning);
}

unless(-f "$rate_matrix") {
	my $specific_warning = "### WARNING: A valid --rate_matrix option must be provided.";
	print_usage_message($specific_warning);
}

unless($branch_unit =~ /\d/) {
	my $specific_warning = "### WARNING: A valid --branch_unit option must be provided.";
	print_usage_message($specific_warning);
}

if (-f "$vcf_output" || $vcf_output =~ /^\-/) {
	my $specific_warning = "### WARNING: A valid --vcf_output file name must be provided (or already exists).";
	print_usage_message($specific_warning);
}

my $working_directory = `pwd`;
chomp($working_directory);
#print "\nworking_directory=$working_directory\n"; # does NOT include a trailing forward slash, but does end in a newline


##########################################################################################
# Extract tree file prefix
my $file_prefix;
if($tree =~/\/([^\/]+)\..+/) { 
	$file_prefix = $1;
} else {
	$file_prefix = 'trevolver_input';
}


##########################################################################################
# Store the tree

open(IN_TREE, "$tree") or die "Could not open file $tree\n";

my $tree = '';

# Store first tree in file
while(<IN_TREE>) {
	chomp($_);
	my $line = $_;
	if($line =~ /^\(/) { # line starts with an opening parentheses
		if($line =~ /\);$/) { # line ends with a closing parentheses and semicolon
			$tree = $line;
			last;
		}
	}
}

close IN_TREE;

if($tree eq '') {
	my $specific_warning = "### WARNING: NO TREE IN TREE FILE. Must begin with '(' and end with ');'.";
	print_usage_message($specific_warning);
}

# Trim end of tree if ';' or whitespace
#print "\nTREE:\n$tree\n\n";
while(substr($tree, -1) =~ /[\s\;]/) {
	chop($tree);
}
#print "\nPROCESSED TREE:\n$tree\n\n";

if($tree =~ /(\,[a-zA-Z0-9\.\-\|\:\']+\,)/) {
	my $specific_warning = "### WARNING: TREE MUST BE STRICTLY BIFURCATING BUT CONTAINS A POLYTOMY: $1";
	print_usage_message($specific_warning);
}

if($tree =~ /(\(\,[^\(]*\,\()/) {
	my $specific_warning = "### WARNING: TREE MUST BE STRICTLY BIFURCATING BUT CONTAINS A POLYTOMY: $1";
	print_usage_message($specific_warning);
}

if($tree =~ /(\)\d+)/) {
	my $specific_warning = "### WARNING: CURRENTLY, SUPPORT VALUES ARE NOT ALLOWED BUT ARE PRESENT: $1";
	print_usage_message($specific_warning);
}

if($tree =~ /(\)\[)/) {
	my $specific_warning = "### WARNING: CURRENTLY, SUPPORT VALUES ARE NOT ALLOWED BUT ARE PRESENT: $1";
	print_usage_message($specific_warning);
}

if($tree =~ /(\],)/) {
	my $specific_warning = "### WARNING: CURRENTLY, SUPPORT VALUES ARE NOT ALLOWED BUT ARE PRESENT: $1";
	print_usage_message($specific_warning);
}

my $tree_opening_paren_count = ($tree =~ s/\(/\(/g);
my $tree_closing_paren_count = ($tree =~ s/\)/\)/g);

unless($tree_opening_paren_count > 0 && $tree_opening_paren_count == $tree_closing_paren_count) {
	my $specific_warning = "### WARNING: TREE MUST CONTAIN EQUAL NUMBERS OF OPENING (now $tree_opening_paren_count\) AND CLOSING " . 
		"(now $tree_closing_paren_count\) PARENTHESES.";
	print_usage_message($specific_warning);
}


##########################################################################################
# Read in the seed sequence from the fasta file
my $seed_seq = '';
my $seq_num = 0;

open(IN_FASTA, "$seed_sequence") or die "Could not open file $seed_sequence\n";

if($verbose) {
	print "\n################################################################################";
	print "\nRecording seed sequence data from $seed_sequence...\n";
}

my $fasta_header = '';

while(<IN_FASTA>) {
	chomp;
	if(/>([\S]+)\s*/) {
		$fasta_header = $1;
		if($seq_num == 0) {
			$seq_num ++;
		} else {
			last;
		}
	} else {
		$seed_seq .= $_;
	}
}

close IN_FASTA;

$seed_seq = uc($seed_seq);
$seed_seq =~ tr/U/T/;

unless($seed_seq =~ /[ACGT]/ && length($seed_seq) >= 3) {
	my $specific_warning = "### WARNING: Sequence does not contain more than two nucleotides (A, C, G, T).";
	print_usage_message($specific_warning);
}

my $seed_seq_length = length($seed_seq);


##########################################################################################
# Read in the rate matrix
my %rate_matrix;
my @ordered_trinucleotides = qw/AAA AAC AAG AAT
								ACA ACC ACG ACT
								AGA AGC AGG AGT
								ATA ATC ATG ATT
								CAA CAC CAG CAT
								CCA CCC CCG CCT
								CGA CGC CGG CGT
								CTA CTC CTG CTT
								GAA GAC GAG GAT
								GCA GCC GCG GCT
								GGA GGC GGG GGT
								GTA GTC GTG GTT
								TAA TAC TAG TAT
								TCA TCC TCG TCT
								TGA TGC TGG TGT
								TTA TTC TTG TTT/;
# ^ each of those corresponds to ONE ROW of the rate matrix

open(IN_RATE_MATRIX, "$rate_matrix") or die "Could not open file $rate_matrix\n";

if ($verbose) {
	print "\n################################################################################";
	print "\nRecording rate matrix from $rate_matrix...\n";
}

my $row_index = 0;
my $max_whitespaces = 0;

while(<IN_RATE_MATRIX>) {
	chomp;
	my $line = $_;
	my $num_whitespaces = ($line =~ s/(\s+)/$1/g);
	if($num_whitespaces > $max_whitespaces) { $max_whitespaces = $num_whitespaces }
	
	if (/([0-9\.eE\-]+)\s+([0-9\.eE\-]+)\s+([0-9\.eE\-]+)\s+([0-9\.eE\-]+)$/) { # in case there's a header or a column of names in front
		$rate_matrix{$ordered_trinucleotides[$row_index]}->{'A'} = $1;
		$rate_matrix{$ordered_trinucleotides[$row_index]}->{'C'} = $2;
		$rate_matrix{$ordered_trinucleotides[$row_index]}->{'G'} = $3;
		$rate_matrix{$ordered_trinucleotides[$row_index]}->{'T'} = $4;
		
		$row_index++;
	}
}

close IN_RATE_MATRIX;

unless($row_index == 64 || $max_whitespaces <= 5) { # i.e., one more than the last real index
	my $specific_warning = "### WARNING: There must be 64 rows X 4 columns of data in the rate matrix.";
	print_usage_message($specific_warning);
}

if ($verbose) { print "\n" }

print "################################################################################".
	"\n##                                                                            ##".
	"\n##               Evolution On Tree Using Custom Rates Initiated!              ##".
	"\n##                                                                            ##".
	"\n################################################################################\n";

print "\nAnalysis initiated at local time $local_time1\n";

print "\nCOMMAND: trevolver.pl @commands\n";

##########################################################################################
# Generate or assign random seed value
if($random_seed) {
	print "\nRANDOM SEED: $random_seed\n";
	srand($random_seed);
} else {
	$random_seed = srand(time ^ $$ ^ unpack "%32L*", `ps wwaxl | gzip`); # (Programming Perl, p. 955)
	print "\nRANDOM SEED: $random_seed\n";
}


##########################################################################################
# Read in the excluded taxa names
my %excluded_taxa_names;

if ($excluded_taxa) {
	if (-f "$excluded_taxa") {
		open(IN_EXCLUDED_TAXA, "$excluded_taxa");
		
		if($verbose) {
			print "\n################################################################################";
			print "\nRecording excluded taxa from $excluded_taxa...\n";
		}
		
		while(<IN_EXCLUDED_TAXA>) {
			chomp;
			my $line = $_;
			
			if($line =~ /\w/) {
				my @this_line_taxa = split(/\,\s*/, $line);
				
				print "\nEXCLUDED TAXA NAMES: ";
				
				my $excluded_taxa_line = '';
				foreach (@this_line_taxa) {
					$excluded_taxa_names{$_} = 1;
					$excluded_taxa_line .= "$_\,";
				}
				
				chop($excluded_taxa_line);
				
				print "$excluded_taxa_line\n";
			}
		}
		
		close IN_EXCLUDED_TAXA;
		
	} else {
		print "\n### WARNING: could not open file $excluded_taxa\. Excluding no taxa.\n";
	}
}


##########################################################################################
# STORE THE TREE AS A MULTIDMINENSIONAL HASH

if ($verbose) {
	print "\n################################################################################\n";
	print "Recursively trivolving sequences from the root...\n";
}

my %tree;

my $num_mutations = 0;
my $total_branch_length = 0;
my $branch_to_tip_length = 0;
my %taxa_histories;
my %generational_histories;
my $node_id = 1;
my %mutated_sites;

unless($suppress_ancestral_seq) {
	print "\nANCESTRAL SEQUENCE: $seed_seq\n";
}

print "\nTREE: $tree\n";

##########################################################################################
# THE SIMULATION: recursive evolution approach using the subroutine evolve_two_subtrees()
##########################################################################################
evolve_two_subtrees($tree, 0, 0, 'n1=root,');
##########################################################################################
##########################################################################################


##########################################################################################
# PRINT RESULTS
##########################################################################################

print "\n\n#########################################################################################\n";
print "# SUMMARY OF RESULTS:\n";
print "#########################################################################################\n";

print "\ntotal number of mutations on all branches: $num_mutations\n";
print "\ntotal branch length (generations): $total_branch_length\n";
print "\nbranch-to-tip length (generations): $branch_to_tip_length\n";
print "\ntotal number of mutations occurred: $num_mutations\n";
print "\nlength of run in seconds: " . (time - $time1) . "\n\n";

print "\n#########################################################################################\n";
print "# RESULTS:\n";
print "#########################################################################################\n";

print "\n//MUTATION\n";

print "taxon\tsite\tmutations\n";
foreach my $taxon (sort {$a <=> $b} keys %taxa_histories) {
	foreach my $mutated_site (sort {$a <=> $b} keys %{$taxa_histories{$taxon}}) {
		print "$taxon\t$mutated_site\t" . $taxa_histories{$taxon}->{$mutated_site} . "\n";
	}
}

if($track_mutations || $tracked_motif) {
	print "\n//TRACKED\n";

	my $out_line = "lineage\tgeneration\t";
	
	if ($track_mutations) { $out_line .= "mutation_rate\tmutation_count\t" }
	if ($tracked_motif) { $out_line .= "motif_count\t" }
	chop($out_line);
	print "$out_line\n";
	$out_line = '';
	
	
	foreach my $these_nodes (sort keys %generational_histories) {
		foreach my $generations_elapsed (sort {$a <=> $b} keys %{$generational_histories{$these_nodes}}) {
			$out_line .= "$these_nodes";
			chop($out_line);
			$out_line .= "\t$generations_elapsed\t";
			
			if ($generational_histories{$these_nodes}->{$generations_elapsed}->{rate}) {
				$out_line .= $generational_histories{$these_nodes}->{$generations_elapsed}->{rate} . 
				"\t" . $generational_histories{$these_nodes}->{$generations_elapsed}->{count} . "\t";
			}
			
			foreach my $motif (sort keys %{$generational_histories{$these_nodes}->{$generations_elapsed}}) {
				if ($motif ne 'rate' && $motif ne 'nodes' && $motif ne 'count' && $motif ne '') {
					$out_line .= $generational_histories{$these_nodes}->{$generations_elapsed}->{$motif} . ",";
				}
			}
			
			chop($out_line);
			print "$out_line\n";
			$out_line = '';
			
		}
	}
}


##########################################################################################
##########################################################################################
### VCF output
##########################################################################################
##########################################################################################
if ($vcf_output =~ /\w+/) {
	chomp($vcf_output);
	
	my (undef, undef, undef, $day, $month, $year) = localtime;
	$year = $year + 1900;
	$month += 1;
	if (length($month)  == 1) { $month = "0$month" }
	if (length($day) == 1) { $day = "0$day" }
	my $today = "$year$month$day";
	
	if ($vcf_output =~ /\//) { # a path was provided
		my $vcf_output_dir = $vcf_output;
		$vcf_output_dir =~ s/\/[^\/]+$//;
		
		if (-d "$vcf_output_dir") {
			if (-f "$vcf_output") { # file already exist?
				my $new_vcf_output = "$working_directory\/$file_prefix\_$random_seed\.vcf";
				print "\n\/\/VCF\n$vcf_output already exists. VCF output being placed in $new_vcf_output.\n";
				$vcf_output = $new_vcf_output;
			} else { # else we're good to go
				print "\n\/\/VCF\nVCF output being placed in $vcf_output.\n";
			}
		} else {
			my $new_vcf_output = "$working_directory\/$file_prefix\_$random_seed\.vcf";
			print "\n\/\/VCF\n$vcf_output_dir does not exist. VCF output being placed in $new_vcf_output.\n";
			$vcf_output = $new_vcf_output;
		}
		
	} else {
		print "\n\/\/VCF\nVCF output being placed in $working_directory\/$vcf_output.\n";
		$vcf_output = "$working_directory\/$vcf_output";
	}
	
	open(OUT_TREVOLVER_VCF, ">>$vcf_output");
	print OUT_TREVOLVER_VCF "##fileformat=VCFv4.1\n" . 
		"##FILTER=<ID=PASS,Description=\"All filters passed\">\n" .
		"##fileDate=$today\n" . 
		"##reference=https://github.com/chasewnelson/trevolver\n" . 
		"##source=<TREVOLVER,Description=\"trevolver.pl @commands\">\n" . 
		"##tree=$tree\n" . 
		"##seed_sequence=$seed_sequence\n" . 
		"##rate_matrix=$rate_matrix\n" . 
		"##branch_unit=$branch_unit\n" . 
		"##random_seed=$random_seed\n" . 
		"##tracked_motif=$tracked_motif\n" . 
		"##contig=<ID=$file_prefix,assembly=1,length=$seed_seq_length>\n" . 
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" . 
		"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">\n" . 
		"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1)\">\n" . 
		"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data, equivalent to number of extant sequences\">\n" . 
		"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes, equivalent to number of extant sequences\">\n" . 
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth; here, equivalent to number of extant sequences\">\n" . 
		"##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele. AA: Ancestral allele, REF:Reference Allele, ALT:Alternate Allele\">\n" . 
		"##INFO=<ID=VT,Number=.,Type=String,Description=\"indicates what type of variant the line represents\">\n" . 
		"##INFO=<ID=MUTATIONS,Number=.,Type=String,Description=\"unique mutations that occurred at this site\">\n" . 
		"##INFO=<ID=GENERATIONS,Number=.,Type=String,Description=\"the generations at which unique mutations occurred at this site\">\n" . 
		"##INFO=<ID=TAXA,Number=.,Type=String,Description=\"the number of extant taxa (leaves) sharing the unique mutations that occurred at this site\">\n" . 
		"##INFO=<ID=MULTIHIT,Number=0,Type=Flag,Description=\"indicates whether a site has experienced multiple hits, i.e., more than one mutation\">\n" . 
		"##INFO=<ID=MULTIALLELIC,Number=0,Type=Flag,Description=\"indicates whether a site is multi-allelic\">\n" . 
		"##INFO=<ID=BACK_MUTATION,Number=0,Type=Flag,Description=\"indicates whether a site has experienced back mutation, i.e., return to a previous state via MULTIHIT\">\n" .
		"##INFO=<ID=RECURRENT_MUTATION,Number=0,Type=Flag,Description=\"indicates whether a site has experienced recurrent mutation, i.e., the same change occurring multiple times independently\">\n" . 
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
	
	my %site_to_alleles;
	my %site_to_outgroups;
	
	my $num_taxa = 0;
	my $num_outgroups = 0;
	
	# KEYS of %mutated_sites
	foreach my $mutated_site (sort {$a <=> $b} keys %mutated_sites) {
		
		my $AA = substr($seed_seq, $mutated_site - 1, 1);
		# Store history
		my %prev_state_1;
		my %prev_state_2;
		my %prev_change;
		
		foreach my $taxon (sort {$a <=> $b} keys %taxa_histories) {
		
			# OUTGROUP(S), if any
			if ($excluded_taxa_names{$taxon}) {
				
				if ($taxa_histories{$taxon}->{$mutated_site}) { # this taxon HAS the mutated site
					
					my $mutation_history = $taxa_histories{$taxon}->{$mutated_site};
					#print "mutation_history=$mutation_history\n";
					
					my @mutation_history_events = split(/,/, $mutation_history);
	
					my $extant_nt = substr($mutation_history, -1);
					#print "extant_nt=$extant_nt\n";
					
					# Store nucleotide
					$site_to_outgroups{$mutated_site}->{$extant_nt}++;
					
					foreach my $event (@mutation_history_events) { # these are ordered, but won't matter with new approach
						
						if ($event =~ /(\d+)\-\w\>\w/) {
							#91404-A>G
							#868627-G>A
							
							my $generation = $1;
							
							$site_to_outgroups{$mutated_site}->{history}->{$generation}->{$event}++;
						}
					}
					
				} else { # this taxon does NOT have the mutated site, add the ancestral (seed)
					$site_to_outgroups{$mutated_site}->{$AA}++;
				}
			
			# INGROUP (or, if no outgroups, everything)
			} else {
			
				if ($taxa_histories{$taxon}->{$mutated_site}) { # this taxon HAS the mutated site
					
					my $mutation_history = $taxa_histories{$taxon}->{$mutated_site};
					#print "mutation_history=$mutation_history\n";
					
					my @mutation_history_events = split(/,/, $mutation_history);
	
					my $extant_nt = substr($mutation_history, -1);
					#print "extant_nt=$extant_nt\n";
					
					# Store nucleotide
					$site_to_alleles{$mutated_site}->{$extant_nt}++;
					
					foreach my $event (@mutation_history_events) { # these are ordered, but won't matter with new approach
						
						if ($event =~ /(\d+)\-\w\>\w/) {
							#91404-A>G
							#868627-G>A
							
							my $generation = $1;
							
							$site_to_alleles{$mutated_site}->{history}->{$generation}->{$event}++;
						}
					}
					
				} else { # this taxon does NOT have the mutated site, add the ancestral (seed)
					$site_to_alleles{$mutated_site}->{$AA}++;
				}
			}
		}
	}
	
	
	foreach my $taxon (sort {$a <=> $b} keys %taxa_histories) {
		
		# OUTGROUP(S), if any
		if ($excluded_taxa_names{$taxon}) { # exclude outgroups
			$num_outgroups++; # COMEBACK -- we don't appear to be using these anymore. Consider pitching.
		
		# INGROUP (all if no outgroups)
		} else {
			$num_taxa++; # COMEBACK -- will this capture ALL taxa in ALL situations? What if no history? I think yes.
		}	
			
	} # end all taxa
	
	my $out_line = '';
	
	#my %multi_hit;
	my @nts = qw/A C G T/;
	foreach my $mutated_site (sort {$a <=> $b} keys %site_to_alleles) {
		#print "mutated_site=$mutated_site\n";
		
		# AA
		my $AA = substr($seed_seq, $mutated_site - 1, 1);
		
		#print "AA=$AA\n";
		
		# Define REF as the major (consensus) allele
		my $REF = '';
		my $REF_count = 0;
		my $arbitrary_REF = 0;
		
		foreach my $nt (@nts) {
			if ($site_to_alleles{$mutated_site}->{$nt} > $REF_count) {
				$REF = $nt;
				$REF_count = $site_to_alleles{$mutated_site}->{$nt};
				$arbitrary_REF = 0;
			} elsif ($site_to_alleles{$mutated_site}->{$nt} == $REF_count) {
				$arbitrary_REF = 1;
			}
		}
		
		#print "REF=$REF\n";
		
		my $ALT = '';
		my $AC = '';
		my $AN_NS_DP = 0;
		my $num_alleles = 0;
		my @all_nts_present;
		
		foreach my $nt (@nts) {
			#print "nt=$nt\n";
			
			if ($site_to_alleles{$mutated_site}->{$nt} > 0) { # here, possible that back mutation eliminates variation
				#print "nt=$nt\n";
				push(@all_nts_present, $nt);
				
				$AN_NS_DP += $site_to_alleles{$mutated_site}->{$nt};
				$num_alleles++;
				
				if ($nt ne $REF) {
					$ALT .= "$nt\,";
					$AC .= $site_to_alleles{$mutated_site}->{$nt} . ',';
				}
			}
		}
		
		# This will only happen if the sample is 100% REF; print it
		if ($ALT eq '' && $AC eq '') {
			$ALT = "$REF\,";
			$AC = "$AN_NS_DP\,";
		}
		
		chop($ALT);
		chop($AC);
		
		#print "ALT=$ALT\nAC=$AC\n";
		
		# AFs
		my $AF = '';
		foreach my $nt (@nts) {
			if ($nt ne $REF && $site_to_alleles{$mutated_site}->{$nt} > 0) {
				my $this_AF = $site_to_alleles{$mutated_site}->{$nt} / $AN_NS_DP;
				$AF .= "$this_AF\,";
			}
		}
		
		# This will only happen if the sample is 100% REF; print it
		if ($AF eq '') {
			$AF = "1\,";
		}
		
		chop($AF);
		
		$out_line .= "$fasta_header\t$mutated_site\t.\t$REF\t$ALT\t100\tPASS\t";
		$out_line .= "AC=$AC\;AF=$AF\;AN=$AN_NS_DP\;NS=$AN_NS_DP\;DP=$AN_NS_DP\;AA=$AA\;VT=SNP\;"; 
		
		##################################################################################
		# INGROUP (SNP) INFORMATION
		my $mutations = '';
		my $generations = '';
		my $taxa = '';
		
		my $hits = 0;
		my %prev_change;
		my %prev_state_1;
		my %prev_state_2;
		my $back_mutation = 0;
		my $recurrent_mutation = 0;
		
		foreach my $generation (sort {$a <=> $b} keys %{$site_to_alleles{$mutated_site}->{history}}) {
			#$site_to_alleles{$mutated_site}->{history}->{$generation}->{$event}++;
			$hits++;
			
			# SHOULD ONLY BE ONE EVENT PER GENERATION
			foreach my $event (%{$site_to_alleles{$mutated_site}->{history}->{$generation}}) {
				
				my $event_num = $site_to_alleles{$mutated_site}->{history}->{$generation}->{$event};
				
				if($event =~ /(\d+)\-([\>\w+]+)/) {
					my $this_mutation = $2;
					$mutations .= "$2\,";
					$generations .= "$1\,";
					$taxa .= $event_num . ',';
					
					my @two_states = split(/>/, $this_mutation);
					my $state1 = $two_states[0];
					my $state2 = $two_states[1];
					
					# Could be recurrent/parallel
					if ($prev_change{$this_mutation}) {
						$recurrent_mutation++;
					}
					
					$prev_change{$this_mutation}++;
					
					# Could be back mutation
					if ($prev_state_1{$state2}) {
						$back_mutation++;
					}
					
					$prev_state_1{$state1}++;
					$prev_state_2{$state2}++;
				}
			}
		}
		
		chop($mutations);
		chop($generations);
		chop($taxa);
		
		if ($mutations ne '' || $generations ne '' || $taxa ne '') {
			$out_line .= "MUTATIONS=$mutations\;GENERATIONS=$generations\;TAXA=$taxa\;";
		}
		
		##################################################################################
		# OUTGROUP (SUB) INFORMATION
		my $mutations_out = '';
		my $generations_out = '';
		my $taxa_out = '';
		
		my $hits_out = 0;
		my %prev_change_out;
		my %prev_state_1_out;
		my %prev_state_2_out;
		my $back_mutation_out = 0;
		my $recurrent_mutation_out = 0;
		
		foreach my $generation (sort {$a <=> $b} keys %{$site_to_outgroups{$mutated_site}->{history}}) {
		
			$hits_out++;
			
			# SHOULD ONLY BE ONE EVENT PER GENERATION
			foreach my $event (%{$site_to_outgroups{$mutated_site}->{history}->{$generation}}) {
				
				my $event_num = $site_to_outgroups{$mutated_site}->{history}->{$generation}->{$event};
				
				if($event =~ /(\d+)\-([\>\w+]+)/) {
					my $this_mutation = $2;
					$mutations_out .= "$2\,";
					$generations_out .= "$1\,";
					$taxa_out .= $event_num . ',';
					
					my @two_states = split(/>/, $this_mutation);
					my $state1 = $two_states[0];
					my $state2 = $two_states[1];
					
					# Could be recurrent/parallel
					if ($prev_change_out{$this_mutation}) {
						$recurrent_mutation_out++;
					}
					
					$prev_change_out{$this_mutation}++;
					
					# Could be back mutation
					if ($prev_state_1_out{$state2}) {
						$back_mutation_out++;
					}
					
					$prev_state_1_out{$state1}++;
					$prev_state_2_out{$state2}++;
				}
			}
		}
		
		chop($mutations_out);
		chop($generations_out);
		chop($taxa_out);
		
		if ($mutations_out ne '' || $generations_out ne '' || $taxa_out ne '') {
			$out_line .= "MUTATIONS_OUTGROUP=$mutations_out\;GENERATIONS_OUTGROUP=$generations_out\;TAXA_OUTGROUP=$taxa_out\;";
		}
		
		
		##################################################################################
		# INGROUP (SNP) FLAGS
		if ($arbitrary_REF == 1) {
			$out_line .= "ARBITRARY_REF\;";
		}
		
		if ($hits > 1) {
			$out_line .= "MULTIHIT\;";
		}
		
		if ($num_alleles > 2) {
			$out_line .= "MULTIALLELIC\;";
		}
		
		if ($back_mutation > 0) {
			$out_line .= "BACK_MUTATION\;";
		}
		
		if ($recurrent_mutation > 0) {
			$out_line .= "RECURRENT_MUTATION\;";
		}
		
		if ($AF eq '1') {
			if ($REF eq $AA) {
				$out_line .= "INVARIANT_ANCESTRAL\;";
			} else {
				$out_line .= "INVARIANT_DERIVED\;";
			}
		}
		
		if ($REF ne $AA) {
			my $match_AA = 0;
			
			foreach (@all_nts_present) {
				if ($_ eq $AA) {
					$match_AA = 1;
				}
			}
			
			# Unless we have a match to the ancestral allele SOMEWHERE
			unless ($match_AA == 1) {
				$out_line .= "NO_ANCESTRAL\;";
			}
		}
		
		chop($out_line);
		
		print OUT_TREVOLVER_VCF "$out_line\n";
		
		$out_line = '';
	}
	
	close OUT_TREVOLVER_VCF;
} # finish VCF output

exit;


##########################################################################################
##########################################################################################
### SUBROUTINE: recursive evolution approach
##########################################################################################
##########################################################################################
sub evolve_two_subtrees {
	
	my($tree, $generations_elapsed, $mutation_count, $nodes, $mutation_history_ref) = @_;
#	if($verbose) { print "We got the arguments: @_\n" }
	
	my %mutation_history;
	
	if($mutation_history_ref ne '') {
		%mutation_history = %{$mutation_history_ref};
	}
	
	if ($verbose) {
		print "\n##########################################################################################\n";
		print "Analyzing node $node_id\nGenerations elapsed: " . sprintf("%.3f", $generations_elapsed) . "\n";
		print "tree: $tree\n";
	}
	
	# IF WE HAVE MORE THAN ONE TAXON LEFT, WE'LL HAVE A COMMA (,)
	if ($tree =~ /\,/) {
		my $tree_stripped = $tree;
		
		if(substr($tree_stripped, 0, 1) eq '(' && substr($tree_stripped, -1) eq ')') { # first and last are parentheses
			$tree_stripped =~ s/^\(//; # remove first opening parentheses
			$tree_stripped =~ s/\)$//; # remove last closing parentheses
		}
		
		if ($verbose) { print "tree_stripped: $tree_stripped\n" }
		
		# NOW SIX BASIC PATTERN POSSIBILITIES:
		# (1) (-----):0.01					<-- here, first matching closing parentheses coincides with branch length
		# (2) (-----):0.01,(-----):0.01			<-- sub-case of #1 in which more follows, enclosed ()
		# (3) (-----):0.01,A:0.01				<-- sub-case of #1 in which more follows, naked
		# (4) A:0.01,(-----):0.01			<-- here, first opening parentheses happens later (keep track); in what follows matched closing, there's more enclosed
		# (5) A:0.01,B:0.01					<-- here, no parentheses
		# (6) A:0.01						<-- here, terminal taxon
		
		my $comma_count = ($tree_stripped =~ s/,/,/g);
		
		# These two categories are mutually exclusive. Easy to test later.
		my $internal_node = ''; # category 1
		my $branch_length = ''; # category 1
		my $subtree1 = ''; # category 2
		my $subtree2 = ''; # category 2
		
		# ONE COMMA
		if ($comma_count == 1) {
			#die "\none comma\n";
			# TWO POSSIBILITIES:
			# (1) (-----):0.01
			# (5) A:0.01,B:0.01				<-- here, no parentheses
			
			##############################################################################
			### PATTERN (1) (-----):0.01. Ancestral branch underlying internal node. EVOLVE IT!
			##############################################################################
		
			# one opening parentheses
			if(($tree_stripped =~ s/\(/\(/g) == 1) {
				# pattern 1 with a SINGLE internal pair
				# (1) (-----):0.01
				
				if ($verbose) { print "pattern 1a: (A:0.01,B:0.01):0.01\n" }
				
				#$pattern = 1;
				
				if($tree_stripped =~ /^(\(.+\)):([0-9\.eE\-]+)$/) {
					# EXAMPLE: ((303:0.000054266,1177:0.000054266):0.000224311,(163:0.000054489,498:0.000054489):0.000224088):0.000134308
					
					$internal_node = $1;
					$branch_length = $2;
					
				}
			
			##################################################################################
			### PATTERN (5) A:0.01,B:0.01
			##################################################################################
		
			# no opening parentheses
			} elsif (($tree_stripped =~ s/\(/\(/g) == 0) {
				# pattern 5 with a single naked pair
				# (5) A:0.01,B:0.01
				
				if ($verbose) { print "pattern 5: A:0.01,B:0.01\n" }
				
				#$pattern = 5;
				
				if($tree_stripped =~ /^(.+:[0-9\.eE\-]+),(.+:[0-9\.eE\-]+)$/) {
					# EXAMPLE: 
					
					# SUBMIT THE TWO TERMINAL BRANCHES AND RETURN FROM SUBROUTINE
					
					$subtree1 = $1;
					$subtree2 = $2;
					
					($subtree1, $subtree2) = sort($subtree1, $subtree2);
					
					if ($verbose) { print "subtree1: $subtree1\nsubtree2: $subtree2\n" }
					
					# Recursively evolve subtrees
					if ($verbose) { print "Analyzing subtree1...\n" }
					evolve_two_subtrees($subtree1, $generations_elapsed, $mutation_count, $nodes, \%mutation_history);
					
					if ($verbose) { print "Analyzing subtree2...\n" }
					evolve_two_subtrees($subtree2, $generations_elapsed, $mutation_count, $nodes, \%mutation_history);
					
					return;
				}
				
			} else {
				die "\n### TERMINATED: Only one comma in tree but more than one opening parentheses. Format unexpected.\n";
			}
		
		# MORE THAN ONE COMMA
		} elsif ($comma_count > 1) {
			# FOUR POSSIBILITIES:
			# (1) (-----):0.01					<-- here, first matching closing parentheses coincides with branch length
			# (2) (-----):0.01,(-----):0.01			<-- sub-case of #1 in which more follows, enclosed ()
			# (3) (-----):0.01,A:0.01				<-- sub-case of #1 in which more follows, naked
			# (4) A:0.01,(-----):0.01			<-- here, first opening parentheses happens later (keep track); in what follows matched closing, there's more enclosed
			
			# PROCESS THE TREE BY MATCHING PARENTHESES TO IDENTIFY WHICH PATTERN HOLDS
			# This should work: split the tree in half at the first point that the number of 
			# close parentheses ')' is equal to the number of opening parentheses '('. Trees 
			# aren't really rate-limiting, so don't have to be smart or efficient (of which I'm
			# neither, by nature).
			
			my $opening_count = 0;
			my $closing_count = 0;
			my $split_index = 0;
			
			my $first_paren_index = undef;
			
			# IDENTIFY TWO SUBTREES; BUT COULD ALSO BE A SINGLE NODE, TOO (EXAMINE SUBTREE2)
			IDENTIFY_SUBTREES: for (my $tree_index = 0; $tree_index < length($tree_stripped); $tree_index++) {
				# count an opening
				if (substr($tree_stripped, $tree_index, 1) eq '(') {
					$opening_count++;
					
					unless(defined($first_paren_index)) {
						$first_paren_index = $tree_index;
					}
				
				# count a closing
				} elsif (substr($tree_stripped, $tree_index, 1) eq ')') {
					$closing_count++;
				}
				
				# detect if matched parentheses pair has been reached
				if ($opening_count > 0 && $opening_count == $closing_count) {
					$split_index = $tree_index;
					
					#print "opening_count=$opening_count\nclosing_count=$closing_count\nsplit_index=$split_index\n";
					
					# walk ahead to capture only the following branch length
					my $new_char = '';
					while ($split_index < length($tree_stripped) && $new_char ne ',' && $new_char ne '(' && $new_char ne ')') {
						$split_index++;
						$new_char = substr($tree_stripped, $split_index, 1);
						#print "new_char=$new_char\n";
					}
					
					$split_index--;
					
					if($verbose) { print "FINAL string index at which to split tree: $split_index\n" }
					last IDENTIFY_SUBTREES;
				}
			}
			
			# store detected subtrees, skipping the trailing comma at $split_index + 1
			my $tree_portion1;
			my $tree_portion2;
			
			if ($first_paren_index == 0) {
				$tree_portion1 = substr($tree_stripped, 0, $split_index + 1);
				$tree_portion2 = substr($tree_stripped, $split_index + 2);
			} else { # the first parentheses MUST be preceded by a comma
				$tree_portion1 = substr($tree_stripped, 0, $first_paren_index - 1);
				$tree_portion2 = substr($tree_stripped, $first_paren_index);
			}
			
			#print "tree_portion1=$tree_portion1\ntree_portion2=$tree_portion2\n";
					
			# check whether it's two subtrees or one internal branch
			if (length($tree_portion2) > 0) { # two subtrees
				($subtree1, $subtree2) = ($tree_portion1, $tree_portion2);
			} elsif ($tree_portion1 =~ /^(.+):([0-9\.eE\-]+)$/) { 
				($internal_node, $branch_length) = ($1, $2);
			} else {
				die "### TERMINATED: tree contains neither two subtrees nor a single internal node.\n";
			}
		
		# NO COMMA
		} elsif ($comma_count == 0) {
			# ONLY ONE POSSIBILITY: A TERMINAL TAXON
			die "\n### WARNING: No comma but a comma? TERMINATED.\n\n";
		} else {
			die "\n### WARNING: Number of commas not qunatifiable; TERMINATED.\n\n";
		}
		
		# print "subtree1: $subtree1\nsubtree2: $subtree2\n";
		
		# NOW REMIND THE SIX BASIC PATTERN POSSIBILITIES:
		# (1) (-----):0.01					<-- here, first matching closing parentheses coincides with branch length
		# (2) (-----):0.01,(-----):0.01			<-- sub-case of #1 in which more follows, enclosed ()
		# (3) (-----):0.01,A:0.01				<-- sub-case of #1 in which more follows, naked
		# (4) A:0.01,(-----):0.01			<-- here, first opening parentheses happens later (keep track); in what follows matched closing, there's more enclosed
		# (5) A:0.01,B:0.01					<-- here, no parentheses
		# (6) A:0.01						<-- here, terminal taxon
		
		##################################################################################
		### PATTERN (1) (-----):0.01. Ancestral branch underlying internal node. EVOLVE IT!
		##################################################################################
		if (length($subtree1) == 0 && length($subtree2) == 0 && length($internal_node) > 1 && length($branch_length) > 1) {
			
			$node_id++;
			$nodes .= "n$node_id\,";
			
			# Construct the current state of the sequence given the evolutionary history
			# Perhaps not quite as time-efficient, but MUCH MORE MEMORY EFFICIENT.
			my $curr_sequence = $seed_seq;
			#print "seed=$seed_seq\n";
			foreach my $mutated_site (sort {$a <=> $b} keys %mutation_history) {
				my $latest_nt = $mutation_history{$mutated_site};
				chop($latest_nt); # remove ending comma (,)
				$latest_nt = chop($latest_nt); # return new end, which is latest nucleotide
#				print "latest_nt=$latest_nt\n";
				substr($curr_sequence, $mutated_site - 1, 1, $latest_nt);
			}
			#print "curr=$curr_sequence\n";
			
			# Calculate number of generations on the branch
			my $generations = $branch_length * $branch_unit;
			$total_branch_length += $generations;
			#$generations = $generations;
			if($verbose) { print "pattern 1: (-----):0.01\ngenerations to evolve on ancestral branch: " . sprintf("%.3f", $generations) . "\n" }
			
			# Tally all trinucleotides in sequence; first and last nucleotides ignored.
			my %trint_counts;
			my $trint_counts_sum = 0;
			for (my $seq_index = 0; ($seq_index + 2) < length($curr_sequence); $seq_index++) {
				$trint_counts{substr($curr_sequence, $seq_index, 3)}++;
				$trint_counts_sum++;
			}
			
			#print "trinucleotide counts, summing to $trint_counts_sum\:\n";
			#foreach (sort keys %trint_counts) {
			#	print "$_\=" . $trint_counts{$_} . "\n";
			#}
			
			# Calculate weighted mean mutation rate per generation.
			my $mut_rate_overall = 0;
			foreach (keys %trint_counts) {
				my $rate_row_sum = $rate_matrix{$_}->{'A'} + $rate_matrix{$_}->{'C'} + 
										$rate_matrix{$_}->{'G'} + $rate_matrix{$_}->{'T'};
				
				$mut_rate_overall += ($trint_counts{$_} * $rate_row_sum);
			}
			
			# Here's the issue: EACH "generations elapsed" is UNIQUE, so of course we'll only store one value per.
			if ($verbose) { print "starting mutation rate: $mut_rate_overall\n" }
			if ($generations_elapsed == 0 && $track_mutations) {
				$generational_histories{$nodes}->{$generations_elapsed}->{rate} = $mut_rate_overall;
				$generational_histories{$nodes}->{$generations_elapsed}->{count} = $mutation_count;
			}
			if ($generations_elapsed == 0 && $tracked_motif) {
				my @motifs_overlapping = ($curr_sequence =~ /(?=$tracked_motif)/g); # ?= means overlapping matches
				$generational_histories{$nodes}->{$generations_elapsed}->{$tracked_motif} = scalar(@motifs_overlapping);
			}
			
			# Generate random exponential waiting time to next mutation
			my $rand_number = rand();
			my $waiting_time = -(1 / $mut_rate_overall) * log($rand_number);
			#print "waiting_time=$waiting_time\n";
			
			while ($waiting_time < $generations) { # a new mutation!
				$num_mutations++; # global
				$mutation_count++; # lineage-specific
				
				# which site and nucleotide change?
				my $rand_number2 = rand($mut_rate_overall);
				
				my $curr_rate_max = 0;
				my $mutation_complete = 0;
				
				my $mutation_site_index;
				my $trint_AA = '';
				my $trint_DA = '';
				my $prev_trint_AA = '';
				my $prev_trint_DA = '';
				my $next_trint_AA = '';
				my $next_trint_DA = '';
				my $nt_AA = '';
				my $nt_DA = '';
				
				DETERMINE_MUTATION: for (my $seq_index = 0; ($seq_index + 2) < length($curr_sequence); $seq_index++) {
					my $trint = substr($curr_sequence, $seq_index, 3);
					
					foreach my $nt (qw/A C G T/) {
						#print "nt=$nt\n";
						my $this_rate = $rate_matrix{$trint}->{$nt};
						$curr_rate_max += $this_rate;
						
						if ($rand_number2 < $curr_rate_max && $this_rate > 0) { # this one was it!
							
							$mutation_site_index = $seq_index + 1;
							
							if($verbose) { print "MUTATION at generation " . ($generations_elapsed + $waiting_time) . "! $trint" . ($mutation_site_index + 1) }
							
							# Prev and next overlapping trinucleotides
							
							# don't look at prev trint for FIRST trint
							unless($seq_index == 0) {
								$prev_trint_AA = substr($curr_sequence, $seq_index - 1, 3);
								$prev_trint_DA = $prev_trint_AA;
								substr($prev_trint_DA, 2, 1, $nt);
							}
							
							# don't look at next trint for LAST trint
							unless(($seq_index + 3) >= length($curr_sequence)) {
								$next_trint_AA = substr($curr_sequence, $seq_index + 1, 3);
								$next_trint_DA = $next_trint_AA;
								substr($next_trint_DA, 0, 1, $nt);
							}
							
							$trint_AA = $trint;
							$nt_AA = substr($trint, 1, 1, $nt); # returns what was there before replacement
							
							if($verbose) { print "$trint" }
							$trint_DA = $trint;
							$nt_DA = $nt;
							
							# change state of actual sequence
							substr($curr_sequence, $mutation_site_index, 1, $nt);
							
							$mutation_complete = 1;
							last DETERMINE_MUTATION;
						}
					}
				}
				
				# For that EXTREMELY RARE case where it's the last site, last nucleotide (T)
				#my $trint = substr($curr_sequence, length($curr_sequence) - 3, 3);
				#print "VERY LAST TRINT!! $trint\n";
				unless($mutation_complete == 1) { # very last site, nucleotide T
					my $trint = substr($curr_sequence, length($curr_sequence) - 3, 3);
					
					$mutation_site_index = (length($curr_sequence) - 2);
					
					if($verbose) { print "MUTATION at generation " . ($generations_elapsed + $waiting_time) . "! $trint" . ($mutation_site_index + 1) }
					
					# Prev and next overlapping trinucleotides
					
					# don't look at prev trint for FIRST trint
					unless(length($curr_sequence) == 3) {
						$prev_trint_AA = substr($curr_sequence, (length($curr_sequence) - 4), 3);
						$prev_trint_DA = $prev_trint_AA;
						substr($prev_trint_DA, 2, 1, 'T');
					}
					
					$trint_AA = $trint;
					$nt_AA = substr($trint, 1, 1, 'T'); # returns what was there before replacement
					
					#substr($trint, 1, 1, 'T'); redundant
					if($verbose) { print "$trint" }
					$trint_DA = $trint;
					$nt_DA = 'T';
					
					# NO NEXT TRINUCLEOTIDE!
					
					# change state of actual sequence
					substr($curr_sequence, $mutation_site_index, 1, 'T'); # DOUBLE CHECK THAT COMEBACK
					$mutation_complete = 1;
				}
				
				# UPDATE the number of generations REMAINING on the branch
				$generations = $generations - $waiting_time;
				$generations_elapsed += $waiting_time;
				
				# UPDATE mutation history
				#$mutation_history{$mutation_site_index + 1} .= int($generations_elapsed + 1) . "\-$trint_AA\-$trint_DA\,";
				$mutation_history{($mutation_site_index + 1)} .= int($generations_elapsed + 1) . "\-$nt_AA\>$nt_DA\,";
				
				# UPDATE number of trinucleotides in sequence, taking advantage of 1 mutation 
				# at a time. One nucleotide change will affect 3 overlapping trinucleotides.
				# Sum (number of trints, 2 less the number of sites) will remain the same. 
				# Also update the overall mutation rate.
				
				# SUBTRACT obliterated (ancestral) trinucleotide values
				if($prev_trint_AA ne '') {
					$trint_counts{$prev_trint_AA}--;
					my $rate_row_sum_prev_trint_AA = $rate_matrix{$prev_trint_AA}->{'A'} + $rate_matrix{$prev_trint_AA}->{'C'} + 
												$rate_matrix{$prev_trint_AA}->{'G'} + $rate_matrix{$prev_trint_AA}->{'T'};
					$mut_rate_overall -= $rate_row_sum_prev_trint_AA; # implicitly, * 1
				}
				
				$trint_counts{$trint_AA}--;
				my $rate_row_sum_AA = $rate_matrix{$trint_AA}->{'A'} + $rate_matrix{$trint_AA}->{'C'} + 
											$rate_matrix{$trint_AA}->{'G'} + $rate_matrix{$trint_AA}->{'T'};
				$mut_rate_overall -= $rate_row_sum_AA;
				
				if($next_trint_AA ne '') {
				#if ($mutation_site_index < (length($curr_sequence) - 2)) { # if not LAST trinucleotide
					$trint_counts{$next_trint_AA}--;
					my $rate_row_sum_next_trint_AA = $rate_matrix{$next_trint_AA}->{'A'} + $rate_matrix{$next_trint_AA}->{'C'} + 
												$rate_matrix{$next_trint_AA}->{'G'} + $rate_matrix{$next_trint_AA}->{'T'};
					$mut_rate_overall -= $rate_row_sum_next_trint_AA;
				}
				
				# ADD spontaneously generated (derived) trinucleotide values
				if($prev_trint_DA ne '') {
					$trint_counts{$prev_trint_DA}++;
					my $rate_row_sum_prev_trint_DA = $rate_matrix{$prev_trint_DA}->{'A'} + $rate_matrix{$prev_trint_DA}->{'C'} + 
												$rate_matrix{$prev_trint_DA}->{'G'} + $rate_matrix{$prev_trint_DA}->{'T'};
					$mut_rate_overall += $rate_row_sum_prev_trint_DA;
				}
				
				$trint_counts{$trint_DA}++;
				my $rate_row_sum_DA = $rate_matrix{$trint_DA}->{'A'} + $rate_matrix{$trint_DA}->{'C'} + 
											$rate_matrix{$trint_DA}->{'G'} + $rate_matrix{$trint_DA}->{'T'};
				$mut_rate_overall += $rate_row_sum_DA; # for example, if CpG created, may be higher
				
				if($next_trint_DA ne '') {
				#if ($mutation_site_index < (length($curr_sequence) - 2)) { # if not LAST trinucleotide
					$trint_counts{$next_trint_DA}++;
					my $rate_row_sum_next_trint_DA = $rate_matrix{$next_trint_DA}->{'A'} + $rate_matrix{$next_trint_DA}->{'C'} + 
												$rate_matrix{$next_trint_DA}->{'G'} + $rate_matrix{$next_trint_DA}->{'T'};
					$mut_rate_overall += $rate_row_sum_next_trint_DA;
				}
				
				if($verbose) { print "; new mutation rate: $mut_rate_overall\n" }
				
#				if ($track_mutations || $tracked_motif) { 
#					$generational_histories{$generations_elapsed}->{nodes} .= "$node_id\,";
#				}
				if ($track_mutations) { 
					$generational_histories{$nodes}->{$generations_elapsed}->{rate} = $mut_rate_overall;
					$generational_histories{$nodes}->{$generations_elapsed}->{count} = $mutation_count;
				}
				if ($tracked_motif) { # COMEBACK that regex work?
					my @motifs_overlapping = ($curr_sequence =~ /(?=$tracked_motif)/g); # ?= means overlapping matches
					$generational_histories{$nodes}->{$generations_elapsed}->{$tracked_motif} = scalar(@motifs_overlapping);
				}
				
	#			$final_mutation_rate = $mut_rate_overall;
				
	#			print "generations_elapsed=$generations_elapsed\n";
	#			print "remaining on branch, generations=$generations\n";
	#			print "mutation_history:\n";
	#			
	#			foreach my $mutated_site (sort {$a <=> $b} keys %mutation_history) {
	#				print "$mutated_site:" . $mutation_history{$mutated_site} . "\n";
	#			}
	#			
	#			print "mut_rate_overall=$mut_rate_overall\n";
	#			print "\n";
				
	#			print "updated trinucleotide counts\:\n";
	#			foreach (sort keys %trint_counts) {
	#				print "$_\=" . $trint_counts{$_} . "\n";
	#			}
				
	#			print "mut_rate_overall=$mut_rate_overall\n";
				
				# New waiting time
				$rand_number = rand();
				$waiting_time = -(1 / $mut_rate_overall) * log($rand_number);
				#print "waiting_time=$waiting_time\n";
			} # finish a new mutation
			
			# REMOVE $curr_sequence because it is no longer needed!
			undef($curr_sequence);
			
			# Add the remaining time in which mutation DIDN'T occur.
			$generations_elapsed += $generations;
			
			if($verbose) { print "evolution of node complete; generations_elapsed: " . sprintf("%.3f", $generations_elapsed) . "\n" }
##			print "mutation_history:\n";
##			
##			foreach my $mutated_site (sort {$a <=> $b} keys %mutation_history) {
##				print "$mutated_site:" . $mutation_history{$mutated_site} . "\n";
##			}
#			
#		#	print "new sequence:\n";
#		#	print "$curr_sequence\n";
#		##	print "mut_rate_overall=$final_mutation_rate\n"; # wasn't final because wasn't updated since last mutation
#		#	print "\n";
			
			if($verbose) { print "now submitting internal node for evolution:\ninternal_node: $internal_node\n" }
			
			# Recursively evolve subtrees
			#my %mutation_history_copy = %mutation_history;
			
			if($verbose) { print "Analyzing internal_node...\n" }
			evolve_two_subtrees($internal_node, $generations_elapsed, $mutation_count, $nodes, \%mutation_history);
			
		##################################################################################
		### SUBTREE PATTERN
		##################################################################################
		} elsif (length($subtree1) > 1 && length($subtree2) > 1 && length($internal_node) == 0 && length($branch_length) == 0) { # CONTAINS PARENTHESES
		#} elsif ($tree_stripped =~ /(^\(.+\):[0-9\.eE\-]+),(\(.+\):[0-9\.eE\-]+$)/) { # NOT TESTED
			if($verbose) { print "pattern 2: (-----):0.01,(-----):0.01\n" }
			#die "pattern 2\n";
			
			($subtree1, $subtree2) = sort($subtree1, $subtree2);
			
			if($verbose) { print "subtree1: $subtree1\nsubtree2: $subtree2\n" }
			
			# Recursively evolve subtrees
			if($verbose) { print "Analyzing subtree1...\n" }
			evolve_two_subtrees($subtree1, $generations_elapsed, $mutation_count, $nodes, \%mutation_history);
			
			if($verbose) { print "Analyzing subtree2...\n" }
			evolve_two_subtrees($subtree2, $generations_elapsed, $mutation_count, $nodes, \%mutation_history);
		
		##################################################################################
		## NO RECOGNIZABLE PATTERN: ABORT!
		##################################################################################
		} else {
			print "$tree_stripped\n";
			my $specific_warning = "### TERMINATED: TREE ABOVE DIDN'T MATCH ANY EXPECTED TREE FORMAT.";
			print_usage_message($specific_warning);
		}
		
	} else { # BASE CASE: A SINGLE TAXON (NO COMMA). EVOLVE IT AND STORE ITS MUTATION HISTORY.
#		die "Might we have a wen ti? #2.\n";
		#return $mutation_history;
		#%taxa_histories
		
		$node_id++;
		$nodes .= "n$node_id\,";
		
		# Construct the current state of the sequence given the evolutionary history
		# Perhaps not quite as time-efficient, but MUCH MORE MEMORY EFFICIENT.
		my $curr_sequence = $seed_seq;
		#print "seed=$seed_seq\n";
		foreach my $mutated_site (sort {$a <=> $b} keys %mutation_history) {
			my $latest_nt = $mutation_history{$mutated_site};
			chop($latest_nt); # remove ending comma (,)
			$latest_nt = chop($latest_nt); # return new end, which is latest nucleotide
#			print "latest_nt=$latest_nt\n";
			substr($curr_sequence, $mutated_site - 1, 1, $latest_nt);
		}
		#print "curr=$curr_sequence\n";
		
		my $taxon;
		my $branch_length;
		
		if ($tree =~ /^(.+)\:(.+)$/) {
			$taxon = $1;
			$branch_length = $2;
			$taxon =~ s/^\(//; # in case it's a 1-taxon tree, i.e., a single sequence
		} else {
			die "Might we have a wen ti? #3.\n";
		}
		
		chop($nodes);
		$nodes .= "=$taxon\,";
		
#		die "Might we have a wen ti? #4.\n";
		
		if($verbose) { print "ANALYZING A TERMINAL TAXON: $taxon\n" }
		
		# Calculate number of generations on the branch
		my $generations = $branch_length * $branch_unit;
		$total_branch_length += $generations;
		#$generations = ($generations + 1);
		if($verbose) { print "Generations on branch: " . sprintf("%.3f", $generations) . "\n" }
		
		# Tally all trinucleotides in sequence; first and last nucleotides ignored.
		my %trint_counts;
		my $trint_counts_sum = 0;
		for (my $seq_index = 0; ($seq_index + 2) < length($curr_sequence); $seq_index++) {
			$trint_counts{substr($curr_sequence, $seq_index, 3)}++;
			$trint_counts_sum++;
		}
		
		#print "trinucleotide counts, summing to $trint_counts_sum\:\n";
		#foreach (sort keys %trint_counts) {
		#	print "$_\=" . $trint_counts{$_} . "\n";
		#}
		
		# Calculate weighted mean mutation rate per generation.
		my $mut_rate_overall = 0;
		foreach (keys %trint_counts) {
			my $rate_row_sum = $rate_matrix{$_}->{'A'} + $rate_matrix{$_}->{'C'} + 
									$rate_matrix{$_}->{'G'} + $rate_matrix{$_}->{'T'};
			
			$mut_rate_overall += ($trint_counts{$_} * $rate_row_sum);
		}
		
	#		$final_mutation_rate = $mut_rate_overall;
	#		my $mut_rate_mean = $mut_rate_overall / $trint_counts_sum;
		
		if($verbose) { print "starting mutation rate: $mut_rate_overall\n" }
		
		# Generate random exponential waiting time to next mutation
		my $rand_number = rand();
		my $waiting_time = -(1 / $mut_rate_overall) * log($rand_number);
#		print "waiting_time=$waiting_time\n";
		
		while ($waiting_time < $generations) { # a new mutation!
			$num_mutations++; # global
			$mutation_count++; # lineage-specific
			
			# which site and nucleotide change?
			my $rand_number2 = rand($mut_rate_overall);
			
			my $curr_rate_max = 0;
			my $mutation_complete = 0;
			
			my $mutation_site_index;
			my $trint_AA = '';
			my $trint_DA = '';
			my $prev_trint_AA = '';
			my $prev_trint_DA = '';
			my $next_trint_AA = '';
			my $next_trint_DA = '';
			my $nt_AA = '';
			my $nt_DA = '';
			
			DETERMINE_MUTATION: for (my $seq_index = 0; ($seq_index + 2) < length($curr_sequence); $seq_index++) {
				my $trint = substr($curr_sequence, $seq_index, 3);
				
				foreach my $nt (qw/A C G T/) {
					#print "nt=$nt\n";
					$curr_rate_max += $rate_matrix{$trint}->{$nt};
					
					if ($rand_number2 < $curr_rate_max) { # this one was it!
						
						$mutation_site_index = $seq_index + 1;
						
						if($verbose) { print "MUTATION at generation " . ($generations_elapsed + $waiting_time) . "! $trint" . ($mutation_site_index + 1) }
						
						# Prev and next overlapping trinucleotides
						
						# don't look at prev trint for FIRST trint
						unless($seq_index == 0) {
							$prev_trint_AA = substr($curr_sequence, $seq_index - 1, 3);
							$prev_trint_DA = $prev_trint_AA;
							substr($prev_trint_DA, 2, 1, $nt);
						}
						
						# don't look at next trint for LAST trint
						unless(($seq_index + 3) >= length($curr_sequence)) {
							$next_trint_AA = substr($curr_sequence, $seq_index + 1, 3);
							$next_trint_DA = $next_trint_AA;
							substr($next_trint_DA, 0, 1, $nt);
						}
						
						$trint_AA = $trint;
						$nt_AA = substr($trint, 1, 1, $nt); # returns what was there before replacement
						
						if($verbose) { print "$trint" }
						$trint_DA = $trint;
						$nt_DA = $nt;
						
						# change state of actual sequence
						substr($curr_sequence, $mutation_site_index, 1, $nt);
						
						$mutation_complete = 1;
						last DETERMINE_MUTATION;
					}
				}
			}
			
			# For that EXTREMELY RARE case where it's the last site, last nucleotide (T)
			#my $trint = substr($curr_sequence, length($curr_sequence) - 3, 3);
			#print "VERY LAST TRINT!! $trint\n";
			unless($mutation_complete == 1) { # very last site, nucleotide T
				my $trint = substr($curr_sequence, length($curr_sequence) - 3, 3);
				
				$mutation_site_index = (length($curr_sequence) - 2);
				
				if($verbose) { print "MUTATION at generation " . ($generations_elapsed + $waiting_time) . "! $trint" . ($mutation_site_index + 1) }
				
				# Prev and next overlapping trinucleotides
				
				# don't look at prev trint for FIRST trint
				unless(length($curr_sequence) == 3) {
					$prev_trint_AA = substr($curr_sequence, (length($curr_sequence) - 4), 3);
					$prev_trint_DA = $prev_trint_AA;
					substr($prev_trint_DA, 2, 1, 'T');
				}
				
				$trint_AA = $trint;
				$nt_AA = substr($trint, 1, 1, 'T'); # returns what was there before replacement
				
				#substr($trint, 1, 1, 'T'); # redundant
				if($verbose) { print "$trint" }
				$trint_DA = $trint;
				$nt_DA = 'T';
				
				# NO NEXT TRINUCLEOTIDE!
				
				# change state of actual sequence
				substr($curr_sequence, $mutation_site_index, 1, 'T'); # DOUBLE CHECK THAT
				$mutation_complete = 1;
			}
			
			# UPDATE the number of generations REMAINING on the branch
			$generations = $generations - $waiting_time;
			$generations_elapsed += $waiting_time;
			
			# UPDATE mutation history
			#$mutation_history{$mutation_site_index + 1} .= int($generations_elapsed + 1) . "\-$trint_AA\-$trint_DA\,";
			$mutation_history{($mutation_site_index + 1)} .= int($generations_elapsed + 1) . "\-$nt_AA\>$nt_DA\,";
			
			# UPDATE number of trinucleotides in sequence, taking advantage of 1 mutation 
			# at a time. One nucleotide change will affect 3 overlapping trinucleotides.
			# Sum (number of trints, 2 less the number of sites) will remain the same. 
			# Also update the overall mutation rate.
			
			# SUBTRACT obliterated (ancestral) trinucleotide values
			if($prev_trint_AA ne '') {
				$trint_counts{$prev_trint_AA}--;
				my $rate_row_sum_prev_trint_AA = $rate_matrix{$prev_trint_AA}->{'A'} + $rate_matrix{$prev_trint_AA}->{'C'} + 
											$rate_matrix{$prev_trint_AA}->{'G'} + $rate_matrix{$prev_trint_AA}->{'T'};
				$mut_rate_overall -= $rate_row_sum_prev_trint_AA;
			}
			
			$trint_counts{$trint_AA}--;
			my $rate_row_sum_AA = $rate_matrix{$trint_AA}->{'A'} + $rate_matrix{$trint_AA}->{'C'} + 
										$rate_matrix{$trint_AA}->{'G'} + $rate_matrix{$trint_AA}->{'T'};
			$mut_rate_overall -= $rate_row_sum_AA;
			
			if($next_trint_AA ne '') {
			#if ($mutation_site_index < (length($curr_sequence) - 2)) { # if not LAST trinucleotide
				$trint_counts{$next_trint_AA}--;
				my $rate_row_sum_next_trint_AA = $rate_matrix{$next_trint_AA}->{'A'} + $rate_matrix{$next_trint_AA}->{'C'} + 
											$rate_matrix{$next_trint_AA}->{'G'} + $rate_matrix{$next_trint_AA}->{'T'};
				$mut_rate_overall -= $rate_row_sum_next_trint_AA;
			}
			
			# ADD spontaneously generated (derived) trinucleotide values
			if($prev_trint_DA ne '') {
				$trint_counts{$prev_trint_DA}++;
				my $rate_row_sum_prev_trint_DA = $rate_matrix{$prev_trint_DA}->{'A'} + $rate_matrix{$prev_trint_DA}->{'C'} + 
											$rate_matrix{$prev_trint_DA}->{'G'} + $rate_matrix{$prev_trint_DA}->{'T'};
				$mut_rate_overall += $rate_row_sum_prev_trint_DA;
			}
			
			$trint_counts{$trint_DA}++;
			my $rate_row_sum_DA = $rate_matrix{$trint_DA}->{'A'} + $rate_matrix{$trint_DA}->{'C'} + 
										$rate_matrix{$trint_DA}->{'G'} + $rate_matrix{$trint_DA}->{'T'};
			$mut_rate_overall += $rate_row_sum_DA; # for example, if CpG created, may be higher
			
			if($next_trint_DA ne '') {
			#if ($mutation_site_index < (length($curr_sequence) - 2)) { # if not LAST trinucleotide
				$trint_counts{$next_trint_DA}++;
				my $rate_row_sum_next_trint_DA = $rate_matrix{$next_trint_DA}->{'A'} + $rate_matrix{$next_trint_DA}->{'C'} + 
											$rate_matrix{$next_trint_DA}->{'G'} + $rate_matrix{$next_trint_DA}->{'T'};
				$mut_rate_overall += $rate_row_sum_next_trint_DA;
			}
			
			if($verbose) { print "; new mutation rate: $mut_rate_overall\n" }
			
#			if ($track_mutations || $tracked_motif) { 
#				$generational_histories{$generations_elapsed}->{nodes} .= "$node_id\,";
#			}
			if ($track_mutations) { 
				$generational_histories{$nodes}->{$generations_elapsed}->{rate} = $mut_rate_overall;
				$generational_histories{$nodes}->{$generations_elapsed}->{count} = $mutation_count;
			}
			if ($tracked_motif) { # COMEBACK that regex work?
				my @motifs_overlapping = ($curr_sequence =~ /(?=$tracked_motif)/g); # ?= means overlapping matches
				$generational_histories{$nodes}->{$generations_elapsed}->{$tracked_motif} = scalar(@motifs_overlapping);
			}
			
			# New waiting time
			$rand_number = rand();
			$waiting_time = -(1 / $mut_rate_overall) * log($rand_number);
			#print "waiting_time=$waiting_time\n";
		} # finish a new mutation
		
		# REMOVE $curr_sequence because it is no longer needed!
		#undef($curr_sequence); # will evaporate in a second anyway
		
		# Add the remaining time in which mutation DIDN'T occur.
		$generations_elapsed += $generations;
		
		if($branch_to_tip_length == 0) {
			$branch_to_tip_length = $generations_elapsed;
		} elsif (int($branch_to_tip_length) != int($generations_elapsed)) {
			die "\n### WARNING: conflicting branch-to-tip length measures for taxon $taxon\: $branch_to_tip_length vs. $generations_elapsed\. TERMINATED.\n\n";
		}
		
		if($verbose) { print "evolution of terminal taxon complete; generations_elapsed: " . sprintf("%.3f", $generations_elapsed) . "\n" }
#		print "tree: $tree\n";
		
		# STORE TAXON HISTORY
#		print "mutation_history:\n";
		foreach my $mutated_site (sort {$a <=> $b} keys %mutation_history) {
			my $this_mutation_history = $mutation_history{$mutated_site};
			chop($this_mutation_history); # get ride of last comma (,)
			$taxa_histories{$taxon}->{$mutated_site} = $this_mutation_history;
			$mutated_sites{$mutated_site} = 1;
#			print "$mutated_site: $this_mutation_history\n";
		}
	} # END BASE CASE: a terminal taxon
	
} # END SUBROUTINE


##########################################################################################
##########################################################################################
### SUBROUTINE: TELL THE USER WHAT TO DO
##########################################################################################
##########################################################################################
sub print_usage_message {
	my ($specific_warning) = @_;
	
	print "\n################################################################################".
		"\n##                                                                            ##".
		"\n##              trevolver: Evolution On Tree Using Custom Rates!              ##".
		"\n##                                                                            ##".
		"\n################################################################################\n";
	
	print "\ntrevolver was TERMINATED because of an error in the input. Specifically:\n";
	print "\n$specific_warning\n";
	
	print "\n################################################################################\n";
	
	print "\nCALL trevolver.pl USING THE FOLLOWING OPTIONS:\n";
	print "\t--tree (REQUIRED): file containing a bifurcating evolutionary tree in newick\n" . 
			"\t\tformat with branch lengths. NO NODE NAMES OR SUPPORT VALUES AT THIS TIME.\n" .
			"\t\tOnly the first encountered tree is used.\n";
	print "\t--seed_sequence (REQUIRED): FASTA file containing starting (seed) sequence at\n" . 
			"\t\ttree root, to be evolved. Only the first sequence encountered is used.\n";
	print "\t--rate_matrix (REQUIRED): file containing 64 x 4 tab-delimited trinucleotide\n" . 
			"\t\trate matrix in alphabetical order. First row values correspond to:\n" . 
			"\t\tAAA>AAA   AAA>ACA   AAA>AGA   AAA>ATA\n";
	print "\t--branch_unit (REQUIRED): branch lengths will be multiplied by this value and\n" . 
			"\t\trounded up to the nearest integer to determine number of generations.\n";
	print "\t--random_seed (OPTIONAL): integer with which to seed the random number\n" . 
			"\t\tgenerator. If not supplied, one will be chosen and reported.\n";
	print "\t--tracked_motif (OPTIONAL): a motif to track after each mutation. For example,\n" . 
			"\t\tto report the number of CpG sites over the course of a run, specify CG.\n";
	print "\t--track_mutations (OPTIONAL): reports the mutation rate and count over time.\n";
	print "\t--excluded_taxa_list (OPTIONAL) [NOT YET SUPPORTED!]: path of file containing a\n" . 
			"\t\tlist of taxa (comma-separated) to exclude from the VCF file and the consensus\n" . 
			"\t\tsequence. This might be desirable if a small number of taxa represent outgroups,\n" . 
			"\t\tto which polymorphism in an ingroup is being compared.\n";
	print "\t--vcf_output (OPTIONAL): name of a VCF format output file to be placed in the working\n" . 
			"\t\tdirectory unless a full path name is given.\n";
	print "\t--print_consensus (OPTIONAL) [NOT YET SUPPORTED!]: prints the consensus sequence\n" . 
			"\t\t(containing the REF allele, here defined as the major allele, at each site) in a\n" . 
			"\t\tseparate output file.\n";
	print "\t--suppress_ancestral_seq (OPTIONAL): suppress printing the ancestral (seed) sequence\n" . 
			"\t\tin the output. This might be desirable if the seed sequence is very large and its\n" . 
			"\t\tinclusion in the output consumes too much disk space.\n";
	print "\t--verbose (OPTIONAL): tell trevolver to tell you EVERYTHING that happens.\n" . 
			"\t\tNot recommended except for development and debugging purposes.\n";
		
	
	print "\n################################################################################\n";
	print "### EXAMPLES:\n";
	print "################################################################################\n";
	
	print "\n### FORMAT:\n";
	
	print "\n\ttrevolver.pl --tree=<newick>.txt --seed_sequence=<seed>.fa --rate_matrix=<64x4>.txt \\\n" . 
			"\t--branch_unit=<#> --track_mutations --tracked_motif=<ACGT> --verbose > output.txt\n";
	
	print "\n### EXAMPLE USING ALL OPTIONS:\n";
	
	print "\n\ttrevolver.pl --tree=my_tree.txt --seed_sequence=my_ancestor.fa --rate_matrix=my_mutations.txt \\\n" . 
			"\t--branch_unit=144740 --random_seed=123456789 --tracked_motif=CG --track_mutations \\\n" . 
			"\t--vcf_output=SNP_report.vcf --suppress_ancestral_seq --verbose > my_output.txt\n";
	
	print "\n### EXAMPLE OF TYPICAL USAGE (program decides random seed; not verbose):\n";
	
	print "\n\ttrevolver.pl --tree=my_tree.txt --seed_sequence=my_ancestor.fa --rate_matrix=my_mutations.txt \\\n" . 
			"\t--branch_unit=144740 --tracked_motif=CG --track_mutations --vcf_output=SNP_report.vcf > my_output.txt\n";
	
	print "\n### EXAMPLE USING MINIMUM OPTIONS WITH OUTPUT TO SCREEN:\n";
	
	print "\n\ttrevolver.pl --tree=my_tree.txt --seed_sequence=my_ancestor.fa --rate_matrix=my_mutations.txt \\\n" .
			"\t--branch_unit=144740\n";
	
	print "\n################################################################################\n";
	print "################################################################################\n\n";
	
	exit;
	
}


exit;



