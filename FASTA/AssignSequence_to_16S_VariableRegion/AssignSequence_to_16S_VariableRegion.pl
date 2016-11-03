#!/usr/bin/env perl

###############################################################################

use strict;
use Getopt::Std;
use File::Basename;
use vars qw($opt_f $opt_o $opt_n $opt_t);

getopts("f:o:n:t:");

my $usage = "usage: 
$0 
	-f <fasta file>
	-o <output report filename root>
	[-n <max sequences to analyze>]
	[-t <reuse reference file name>]

	This script will take the input fasta file and determine where
	all the sequences fall on a reference 16S sequence in order to 
	determine the median coverage of the reads.

	This script will call on blast to perform the alignment
	and the results will be parsed and analyzed.

	The reference sequence is hard coded because the variable
	regions are mapped on to the reference sequence.

	if -n option is specified grab the top N sequences, else all.

	use the -t option, to specify a temporary that will can be
	reused between runs.

	Example:
		$0 \\
			-f input.fasta \\
			-o output \\
			-n 3000 \\
			-t /tmp/ref
";


if(!defined($opt_f) || !defined($opt_o)){
	die "$usage\n";
}

my $input_fasta=$opt_f;
my $output_root=$opt_o;

my $num_seq=-1;
if(defined($opt_n)){
	$num_seq=$opt_n;
}

my $tmp_ref_fname=undef;
if(defined($opt_t)){
	$tmp_ref_fname=$opt_t;
}

print STDERR "Input FASTA: $input_fasta\n";
print STDERR "Ouput Root: $output_root\n";
print STDERR "Max Sequences to Analyze: $num_seq\n";
print STDERR "Reuse Reference: $tmp_ref_fname\n";

###############################################################################

my $reference_16s="
AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGA
AGCTTGCTTCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGAGGGATAACTACTGGA
AACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATG
GGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGA
ACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGC
CGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATT
GACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAA
TTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGA
TACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACC
GGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGG
TAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACC
GCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAAT
TCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACAGAACTTTCCAGAGATGGATTGGTGCCTTCGGGAACTGT
GAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCC
TTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTC
ATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGA
CCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCA
GAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGT
AGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAAC
CTGCGGTTGGATCACCTCCTTA
";

my $reference_defline=">ecoli.O157_H7";

# Have temp reference sequence file based on input file name
my $ref_tmp="$output_root.ref.fasta";

# If user specified fixed reference file name, use it
if(defined($tmp_ref_fname)){
	$ref_tmp=$tmp_ref_fname;
}

# If reference fasta doesn't exist, create and populate it 
if(!(-e $ref_tmp)){
	print STDERR "Creating $ref_tmp...\n";
	open(REF_FH, ">$ref_tmp") || die "Could not open $ref_tmp\n";
	print REF_FH "$reference_defline";
	print REF_FH "$reference_16s";
	close(REF_FH);

	print `formatdb -p F -i $ref_tmp`;
}else{
	print STDERR "Reusing exiting $ref_tmp...\n";
}

###############################################################################

# Perform subsampling of top/first sequences if requested
my $tmp_subsmp_fasta="$output_root.subsmp.fasta";
if($num_seq>0){

	print STDERR "Grabbing top $num_seq for blast.\n";
	open(TMP_FH, ">$tmp_subsmp_fasta") || die "Could not open $tmp_subsmp_fasta\n";
	open(F_FH, "<$input_fasta") || die "Could not open $input_fasta\n";

	my $num_recs=0;
	while(<F_FH>){
		if($_=~/^>/){
			if($num_recs==$num_seq){
				last;
			}
			$num_recs++;
		}
		print TMP_FH $_;
	}

	close(TMP_FH);
	close(F_FH);

	print STDERR "Num Records: $num_recs / $num_seq\n";

	$input_fasta=$tmp_subsmp_fasta;
}

###############################################################################

# Set up and perform blast
my $bl_out_tmp="$output_root.blast_out";
my $blast_cmd="blastall -p blastn -i $input_fasta -d $ref_tmp -o $bl_out_tmp -F F -e 1e-5 -G 1 -a 2 -W 7 -q -2 -m 8";

print STDERR "Running blast...\n";
print STDERR "$blast_cmd\n";
print `$blast_cmd`;

##############################################################################

# Clean up temp reference fasta, if necessary
if(!defined($tmp_ref_fname)){
	`rm $ref_tmp`;
	`rm $ref_tmp.nhr`;
	`rm $ref_tmp.nin`;
	`rm $ref_tmp.nsq`;
}

# Remove subsample file if necessary
if($num_seq>0){
	`rm $tmp_subsmp_fasta`;
}

###############################################################################

my ($output_filename)=fileparse($output_root);

my $Rcode;

$Rcode=qq`
options(echo=F);

###############################################################################
# Read blast output
data=read.delim("$bl_out_tmp", sep="\t", header=F);
#print(dim(data));
colnames(data)=c("qry", "sbj", "perc_id", "aln_len", "mism", "gap",
               "qry_st", "qry_end", "sbj_st", "sbj_end", "eval", "bit");
#print(data);


pdf("$output_root.pdf", height=11, width=8.5);

par(oma=c(0,0,2,0));
par(mfrow=c(3,1));
par(mar=c(5.1,4.1,6.1,2.1));

ss=data[,"sbj_st"];
se=data[,"sbj_end"];
qs=data[,"qry_st"];
qe=data[,"qry_end"];
s_al=abs(ss-se);
q_al=abs(qs-qe);

max_qry=max(c(qs, qe));

###############################################################################
# Plot subject (reference) alignment
dss=density(ss);
dse=density(se);
plot(0, type="n", xlim=c(0,1542), ylim=c(0,max(c(dss\$y, dse\$y))), 
	main="Reference Alignment Positions",
	xlab="Position on Reference 16S (bp)",
	ylab="Density",
	xaxt="n"
);

vregion_pos=c(8,96,306,487,746,885,1029,1180,1372,1468);
label_pos=diff(vregion_pos)/2 + vregion_pos[1:(length(vregion_pos)-1)];
abline(v=vregion_pos, col="grey", lty=5);
axis(side=3, at=label_pos, labels=paste("V", 1:9, sep=""));
bottom_axis=seq(0, 1500, 100);
axis(side=1, at=bottom_axis, labels=bottom_axis);

points(dss, col="blue", type="l");
points(dse, col="green", type="l");

###############################################################################
# Plot query (read) alignment
dqs=density(qs);
dqe=density(qe);
max_len=max(c(qs,qe));
plot(0, type="n", xlim=c(0,max(c(qs, qe))), ylim=c(0, max(c(dqs\$y, dqe\$y))), 
	main="Subject Alignment Positions",
	xlab="Subject Position (bp)",
	ylab="Density",
	xaxt="n");
points(dqs, col="blue", type="l");
points(dqe, col="green", type="l");
bottom_axis=c(seq(0, max_len, 50), max_len);
axis(side=1, at=bottom_axis, labels=bottom_axis);

###############################################################################
# Plot histogram of alignment lengths
hist(s_al, main="Reference Alignment Lengths", xlab="Lengths (bp)");

mtext("$output_filename", side=3, line=0, outer=T, font=2, cex=1.4);

dev.off();

`;

###############################################################################

# Prep R script file
my $tmp_r_code_file="$output_root.tmp.r";
open(Rcode_FH, ">$tmp_r_code_file") || die "Could not open $tmp_r_code_file.\n";
print Rcode_FH "$Rcode";
close(Rcode_FH);

# Execute R
system("R --no-save --quiet < $tmp_r_code_file");

# Remove blast out and temporary R code
`rm $tmp_r_code_file`;
`rm $bl_out_tmp`;

print STDERR "Done.\n";

