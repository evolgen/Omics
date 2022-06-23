# Purpose: Calculation of Ka/Ks valus  of orthologous sequences
# Input: file with coding sequences
# Output: tab-separated file with Ka/Ks values
# Important: coding sequences have to be in the right frame, should not produce protein sequences with Stop codon
#	Muscle to align protein sequences of orthologos
# Paml for Ka/Ks analysis

use strict;
use warnings;
$|=1;

use Bio::Seq;
use Bio::SeqIO;
use Bio::Tools::Run::Phylo::PAML::Codeml;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Tools::Run::Alignment::MAFFT;
use Bio::Align::Utilities qw(aa_to_dna_aln); # for projecting alignments from protein to R/DNA space
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::TreeIO;

### 1 ### Variables for creating tab-delimited file with coding sequences for DB1 and DB2
# Variables for first database e.g. human coding sequences
my $DB1_infile;	# SeqIO object holding the input file in fasta format
my $DB1_seq;	# Seq object holding each sequence of the input file
my @DB1_header;	# array holding all parts of the header; were split by "|" before
my $DB1_headerPart; # variable holding each header part

# Variables for second database e.g. chimp coding sequences
my $DB2_infile; # SeqIO object holding the input file in fasta format
my $DB2_seq;    # Seq object holding each sequence of the input file
my @DB2_header; # array holding all parts of the header; were split by "|" before
my $DB2_headerPart; # variable holding each header part

### 1 ### Variables for creating tab-delimited file with coding sequences for DB1 and DB2
# Variables for first database e.g. human coding sequences
my $DB3_infile;	# SeqIO object holding the input file in fasta format
my $DB3_seq;	# Seq object holding each sequence of the input file
my @DB3_header;	# array holding all parts of the header; were split by "|" before
my $DB3_headerPart; # variable holding each header part

### 1 ### Variables for creating tab-delimited file with coding sequences for DB1 and DB2
# Variables for first database e.g. human coding sequences
my $DB4_infile;	# SeqIO object holding the input file in fasta format
my $DB4_seq;	# Seq object holding each sequence of the input file
my @DB4_header;	# array holding all parts of the header; were split by "|" before
my $DB4_headerPart; # variable holding each header part

### 3 ### Variables for truncating Stop codons and creating fasta file with OrthoPair sequences
my $OrthoPairs;		# when file with all the ortholgs get read in, every line will be an OrthoPair
my @info;       	# to hold each part of OrthoPairs: info[0]=DB1 gene name, info[1]=DB1 sequence, info[2]=DB2 gene name, info[3]=DB2 sequence
my $trunc_DB1_seq;	# will be DB1 sequence but without stop codon at the end
my $trunc_DB2_seq;	# will be DB2 sequence but without stop codon at the end
my $trunc_DB3_seq;
my $trunc_DB4_seq;
my $DB1_seq_obj; 	# sequence object of the human gene / first gene of the orthologous pair
my $DB2_seq_obj; 	# sequence object of the chimp gene / second gene of the orthologous pair
my $DB3_seq_obj;
my $DB4_seq_obj;
my $outfilename;	# to set name of outfile together out of gene name;
my $Pair_seqio_obj;	# SeqIO object for output file that will have the truncated sequences of the OrthoPair in fasta format

### 4 ### Variables for Muscle to Paml
# Variables for DNA sequence data and translated amino acid sequence data
my $DNA_seqio;
my $DNA_seq_obj;
my %DNA_seqHash;
my $AA_seq_obj;
my $AA_seq;
my @all_AA_seqs;
# Variables for the alignment
my $AA_alnio_obj;
my $aln_factory;
my $AA_aln_obj;
my $DNA_aln_obj;
my @each_DNA_aln;
my $kaks_factory;
# Variables for PAML results
my $KaKs_result;

### 3 ### This part reads in tab-delimited file with orthologous genes; names are in odd columns, sequences are in even columns
### 3 ### Puts sequences in Seq objects with the name as display ID
### 3 ### Truncates Stop codons at the end of each sequence
### 3 ### Creates a multi fasta file with the sequences that can then go through Muscle

open(KAKS_OUT, ">D://artigos/artigo_tese/scripts/Pantro3/kaks/kaks_output_pantro3.txt") ||  die("cannot open KaKs_outputfile for writing");
print KAKS_OUT "gene_ids \t H_dn \t H_dS \t H_w \t C_dn \t C_dS \t C_w \t HC_dN \t HC_ds \t HC_w \t O_dN \t O_ds \t O_W \t M_dn \t M_ds \t M_W\n";

#paml_results
my @H_Ka_Ks_W;
my @C_Ka_Ks_W;
my @H_C_Ka_Ks_W;
my @O_Ka_Ks_W;
my @M_Ka_Ks_W;
my $line;

open ORTHOLOGS, "<D://artigos/artigo_tese/scripts/Pantro3/kaks/ortho_seq_pantro3.txt" || die "Can not open Orthologs file";

#print ("What alignment program do you wish to use: 1: muscle; 2: MAFFT\n");
#my $prog_align = <>;
my $prog_align =1;
chomp $prog_align;

#print " what model for codons? 0: one, 1:b, 2:2 or more dN/dS ratios for branches \n";
#my $model=<>;
my $model=1;
chomp $model;

# goes through whole Orthologs.txt file and puts every pair in sequence objects
while ($OrthoPairs = <ORTHOLOGS>)
	{
	chomp $OrthoPairs;
	# seperates gene names and sequences
	@info = split ("\t", $OrthoPairs);
	# put sequences into variables trunc_ ... and truncates them afterwards by their stop codon; "+$" used to match only end of the sequence
	$trunc_DB1_seq = $info[1];
	$trunc_DB2_seq = $info[3];
  $trunc_DB3_seq = $info[5];
  $trunc_DB4_seq = $info[7];

	if ($trunc_DB1_seq =~ s/TAA+$//)
        	{
        	$trunc_DB1_seq =~ s/TAA+$//;
        	}
	elsif ($trunc_DB1_seq =~ s/TAG+$//)
        	{
        	$trunc_DB1_seq =~ s/TAG+$//;
        	}
	elsif ($trunc_DB1_seq =~ s/TGA+$//)
        	{
        	$trunc_DB1_seq =~ s/TGA+$//;
        	}
	# the same for DB2
	if ($trunc_DB2_seq =~ s/TAA+$//)
                {
                $trunc_DB2_seq =~ s/TAA+$//;
                }
        elsif ($trunc_DB2_seq =~ s/TAG+$//)
                {
                $trunc_DB2_seq =~ s/TAG+$//;
                }
        elsif ($trunc_DB2_seq =~ s/TGA+$//)
                {
                $trunc_DB2_seq =~ s/TGA+$//;
                }

   if ($trunc_DB3_seq =~ s/TAA+$//)
        	{
        	$trunc_DB3_seq =~ s/TAA+$//;
        	}
	elsif ($trunc_DB3_seq =~ s/TAG+$//)
        	{
        	$trunc_DB3_seq =~ s/TAG+$//;
        	}
	elsif ($trunc_DB3_seq =~ s/TGA+$//)
        	{
        	$trunc_DB3_seq =~ s/TGA+$//;
        	}

         if ($trunc_DB4_seq =~ s/TAA+$//)
        	{
        	$trunc_DB4_seq =~ s/TAA+$//;
        	}
	elsif ($trunc_DB4_seq =~ s/TAG+$//)
        	{
        	$trunc_DB4_seq =~ s/TAG+$//;
        	}
	elsif ($trunc_DB4_seq =~ s/TGA+$//)
        	{
        	$trunc_DB4_seq =~ s/TGA+$//;
        	}

	# builds sequence objects of the gene name coming from Orthologs and the truncated sequences
	$DB1_seq_obj = Bio::Seq -> new (-id => 'human', -seq => $trunc_DB1_seq);
	$DB2_seq_obj = Bio::Seq -> new (-id => 'chimp', -seq => $trunc_DB2_seq);
  $DB3_seq_obj = Bio::Seq -> new (-id => 'orangutan', -seq => $trunc_DB3_seq);
  $DB4_seq_obj = Bio::Seq -> new (-id => 'macaque', -seq => $trunc_DB4_seq);
	# generates outfilename consisting of gene name
	$outfilename = $info[0]."fa"; #."_".$info[2]."_".$info[4]."_". $info[6].".fa";
	# creates SeqIO object and writes file with fasta sequences of the OrthoPair in it
	$Pair_seqio_obj = Bio::SeqIO -> new ('-file' => ">$outfilename", '-format' => 'Fasta');
	$Pair_seqio_obj -> write_seq ($DB1_seq_obj, $DB2_seq_obj, $DB3_seq_obj, $DB4_seq_obj);

	### 4 ###

	$DNA_seqio = new Bio::SeqIO(-file   => "$outfilename",
                            -format => 'fasta');
	$outfilename =~ s/.fa+$//;
	# File that contains translated protein sequences of outfilename.fa
	open TMP_AA_File, ">$outfilename.tmp.fas" || die "Can't open tmp protein sequence fasta file";

	# process each sequence:
	while ($DNA_seq_obj = $DNA_seqio->next_seq )
		{
		# puts sequences in a hash to later on translate AA sequences back to DNA sequence
		$DNA_seqHash{$DNA_seq_obj->display_id} = $DNA_seq_obj;
		# translate them into protein and creates AA sequence object
		$AA_seq_obj = $DNA_seq_obj->translate();
		# gets AA sequence of the AA sequence object
		$AA_seq = $AA_seq_obj->seq();
		# checks AA sequence for Stop codons
		if( $AA_seq =~ /\*/ && $AA_seq !~ /\*$/ )
			{
              print $AA_seq, "\n";
          		warn("provided a CDS sequence with a stop codon, PAML will choke!");
          		#exit(0);
   	        	}
		# Tcoffee can't handle '*' even if it is trailing
		$AA_seq =~ s/\*//g;
    #print $AA_seq, "\n";
		$AA_seq_obj->seq($AA_seq);
		push @all_AA_seqs, $AA_seq_obj;
		print TMP_AA_File (">", $DNA_seq_obj->display_id, "\n");
 	  print TMP_AA_File ($AA_seq, "\n");
		}
	# checks if at least two AA sequences got provided
	if( @all_AA_seqs < 2 )
		{
        	warn("Need at least 2 CDS sequences to proceed");
        	exit(0);
		}

	# Runs MUSCLE with a multi fasta file with translated protein sequences; will create a multi fasta file with all sequences aligned and same length
# if ($prog_align == 1){

  print ("I'm running Muscle now\n");
	system "D://muscle3.6/muscle -in $outfilename.tmp.fas -out $outfilename.muscle.out.fas";
  $outfilename="$outfilename.muscle.out.fas";
# }
# elsif ($prog_align == 2){
#  print ("I'm running MAFFT now \n");
#  system "C://Programas/mafft/usr/local/bin/mafft.bat --auto $outfilename.tmp.fas > $outfilename.mafft.out.fas";
#  $outfilename="$outfilename.mafft.out.fas";
# }
# else {
#  print "please choose a valid alignment program.";
#  exit ;
# }
 
	# to convert ailgnment into AlignI object
	$AA_alnio_obj  = Bio::AlignIO -> new (-file => $outfilename,
      	                               '-format' => 'fasta');
	# Puts every protein alignment in an AlignI object and processes it through Paml
	while ( my $AA_aln_obj = $AA_alnio_obj->next_aln() )
      	{
		print ("This is the protein alignment of your genes: \n");
        	print $AA_aln_obj->length, "\n";
        	print $AA_aln_obj->percentage_identity, "\n";

		# project the protein alignment object back to CDS coordinates and creates DNA alignment object
		$DNA_aln_obj = aa_to_dna_aln($AA_aln_obj, \%DNA_seqHash);
    my $DNA_aln;
		@each_DNA_aln = $DNA_aln_obj->each_seq();

   open (ALIGN, ">C://Perl64/bin/align_paml.ph");
   print ALIGN $DNA_aln_obj->no_sequences, "\t",$DNA_aln_obj->length, "\n";

  foreach $DNA_aln(@each_DNA_aln){

    print ALIGN $DNA_aln->display_id()."  ".$DNA_aln->seq(),"\n";
  }
  close ALIGN;
}
    #Paml version 4.3
    system "C://Perl64/bin/codeml";
    print "Finished paml\n";
    
    open (paml_results, "<C://Perl64/bin/mlc");
    
@H_Ka_Ks_W=();
@C_Ka_Ks_W=();
@H_C_Ka_Ks_W=();
@O_Ka_Ks_W=();
@M_Ka_Ks_W=();
$line="";

#parses paml output
while($line=<paml_results>){

  #dn_ds: ((hg18: 0.000015, panTro2: 0.012065): 0.011358, ponAbe2: 0.019717, rheMac2: 0.106161);

  if ($line=~m /\(\(human: ([0-9\.]+), chimp: ([0-9\.]+)\): ([0-9\.]+), orangutan: ([0-9\.]+), macaque: ([0-9\.]+)\);\n/)
				{
						push(@H_Ka_Ks_W,$1);
						push(@C_Ka_Ks_W, $2);
						push(@H_C_Ka_Ks_W, $3);
						push(@O_Ka_Ks_W, $4);
            push(@M_Ka_Ks_W, $5)
		}
  #W: ((human '#0.6693' , chimp '#1.7558' ) '#0.9156' , orangutan '#0.3646' , macaque '#0.3012' );

  if ($line=~m /\(\(human '#([0-9\.]+)' , chimp '#([0-9\.]+)' \) '#([0-9\.]+)' , orangutan '#([0-9\.]+)' , macaque '#([0-9\.]+)' \);\n/)
				{
						push(@H_Ka_Ks_W,$1);
						push(@C_Ka_Ks_W, $2);
						push(@H_C_Ka_Ks_W, $3);
						push(@O_Ka_Ks_W, $4);
            push(@M_Ka_Ks_W, $5)
		}
}
  print KAKS_OUT join ("\t",  $info[0], $H_Ka_Ks_W[2], $H_Ka_Ks_W[1], $H_Ka_Ks_W[3],
                              $C_Ka_Ks_W[2], $C_Ka_Ks_W[1], $C_Ka_Ks_W[3],
                              $H_C_Ka_Ks_W[2], $H_C_Ka_Ks_W[1], $H_C_Ka_Ks_W[3],
                              $O_Ka_Ks_W[2], $O_Ka_Ks_W[1], $O_Ka_Ks_W[3],
                              $M_Ka_Ks_W[2], $M_Ka_Ks_W[1],  $M_Ka_Ks_W[3]), "\n";



   	#### 4 ### end of Paml part

}
close KAKS_OUT;

print "the End\n";
