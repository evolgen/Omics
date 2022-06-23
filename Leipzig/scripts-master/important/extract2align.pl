#!/usr/bin/env perl

#perl extract2align.pl --infile [path_to_contigs] --idlist [path_to_list] --outfile [path_to_output]

use warnings;
use strict;

use Bio::SeqIO;
use Bio::DB::Fasta;
use Getopt::Long;

#This script is designed to take in a large multisequence fasta file and return
#a file containing a subset of those sequences
#Optionally a user can get the resulting sequences aligned using either muscle or clustalw


#Notes:
#The ids in the returned file are truncated (as requested by user) 
#If you have duplicate sequences (duplicate ids), clustalw will fail. 

##################USER SETTINGS###########################

my $viewProg = "seaview";  # set to "my $viewProg = 0;" if you don't want to view alignment
my $cutIdsShort = 1;  #ids in return file are truncated at . as requested by user
                      #note that this should have happened within Bioperl function but I couldn't
		      #get it to work! See make_my_id function below
#######################################################		      
		      

my $inSeqFile = '';   
my $listOfIds = '';  
my $outSeqFile = '';
my $help = '';

my $alignMethod = '';
my $outAlignFile = '';



    
     
GetOptions ('infile|i=s' => \$inSeqFile, 'idlist|l=s' => \$listOfIds, 'outfile|o=s' => \$outSeqFile, 
            'help|h' => \$help, 'align|a=s' => \$alignMethod, 'alignfile|af=s' => \$outAlignFile);

if ($help) 
{
print <<"END"
	\n\nThis program extracts a set of sequences from a multi-fasta sequence file and will optionally align  them.\n 
           The user must have a text list of sequences to extract.\n 
	   usage:  extract2file.pl --infie = multiseqs.fasta --idlist mylist.txt --outfile = outfile.fasta\n\n
	   the options available are:
	   
	   --infile   or   --i   infilename.fasta
	   --idlist   or   --l   listfilename
	   --outfile  or   --o   outfilename
	   
	   --align    or   --a   alignment method   (Choices are clustalw or muscle)
	   --alignfile or  --af  alignment_out_filename 
	   
	   If any of the first three options are not on the command line, the user is prompted for 
	   the information. If no --align or --a is given on the command line, it assumed that you 
	   do not want to align the sequences.
	   
	   A setting in the script causes any alignment produced to be opened up in seaview when 
	   the alignment finishes. This can be edited in the script.
	   
	   Currently there is not a lot of error checking in this script, so any input files entered that
	   do not exist will cause the program to keel over.\n\n
END
;
exit 0 ;
}


($inSeqFile, $listOfIds, $outSeqFile, $alignMethod, $outAlignFile) = checkInfo($inSeqFile, $listOfIds, $outSeqFile,$alignMethod, $outAlignFile);

#read in sequence ids desired into array

open(IDLIST, $listOfIds) or die "I can't open the file $listOfIds: $!\n";
my @idList = <IDLIST>;
close(IDLIST);
#my $seqInx = Bio::Index::Fasta->new (-filename => $inxFile, -write_flag => 1);

my $db = Bio::DB::Fasta->new("$inSeqFile");#, -makeid => \&make_my_id );
#technically make_my_id should work! I know the regexp in the subroutine is working
my @ids     = $db->ids;
#print "@ids \n";



my %resultsHash = ();
my $sequence;
my $fail = 0;
my $success = 0;

 foreach my $id (@idList)
 {
 chomp($id);
 my $seqobj = $db->get_Seq_by_id($id);
#if a sequence id is found that is not found in the main file, report this to user
 eval 
 {
 	$sequence = $seqobj->seq;
 };
 if ($@) 
 {
	print "\n$id not included subsequence file because $id not found in main file";
	$fail = 1;
	
 }
  #bioperl doing something strange, putting > from next sequence onto end of this sequence
 #take off 

if ($fail == 0)
{ 
 $sequence =~ s/(\.*)>/$1/;
 $resultsHash{$id} = $sequence;
 $success = $success + 1;
}
else { $fail = 0; } 
 
 }
 
#my %readySeqs = simplify_ids(%resultsHash);

print_for_align(\%resultsHash, $outSeqFile, $success);

if ($alignMethod) 
{
  if ($success > 0 ) 
  {
	if ( $alignMethod eq "muscle") { require Bio::Tools::Run::Alignment::Muscle; }
	elsif ($alignMethod eq "clustalw") { require Bio::Tools::Run::Alignment::Clustalw; }

	my $runViewer = align_seqs($outSeqFile, $outAlignFile, $alignMethod);
	if (($viewProg) and ($runViewer))  
	{ 
		my $isInstalled = isProgInstalled($viewProg);
		if ($isInstalled) { system ("$viewProg $outAlignFile"); }
		
	}
  }

else { print "\nSorry, no sequences were extracted to be aligned. Are the ids in your list correct?\n\n"; }		
	exit 1;

}


sub checkInfo
{    
	my ($inSeqFile, $listOfIds, $outSeqFile) = @_;

	unless ($inSeqFile) 
	{	
		 my $askString = "\nPlease enter name of sequence file: "; 
		$inSeqFile = promptUser($askString);	
	}

	unless ($listOfIds) 
	{
		 my $askString = "\nPlease enter name of file with desired sequence ids: "; 
		$listOfIds = promptUser($askString);
	}	

	unless ($outSeqFile) 
	{
	 	my $askString = "\nPlease enter name of file to write sequences to: "; 
		$outSeqFile = promptUser($askString);
	}
	
	if ($alignMethod)
	{
	
		unless ($outAlignFile)
		{
	 		my $askString = "\nPlease enter name of file to write alignment to: "; 
			$outAlignFile = promptUser($askString);
		}
	
		if ($alignMethod =~ /^m/) 
			{ 
			my $isInstalled = isProgInstalled("muscle");
			
			if ($isInstalled == 1)  {$alignMethod = "muscle";}
			else {$alignMethod = "none";} 	
			}
		elsif ($alignMethod =~ /^c/) 
		{ 
			my $isInstalled = isProgInstalled("clustalw");
			if ($isInstalled == 1)  {$alignMethod = "clustalw";}
			else {$alignMethod = "none";} 		
		}

		unless (($alignMethod =~ /clustal/ ) || ($alignMethod =~ /muscle/ ) || ($alignMethod =~ /no/))
		{
			do
			{
				my $askString = "\nAlign using clustalw, muscle or no alignment [c|m|n]: ";
				$alignMethod = promptUser($askString);
			}
			while ($alignMethod ne /^[mcn]/);
		
		}
		
		if ( $alignMethod =~ /^n/ ) { print "\nNo alignment will be carried out\n\n"; }	 
		
		
	} 
	
	return ($inSeqFile, $listOfIds, $outSeqFile, $alignMethod, $outAlignFile);
}


#################Functions##############################

sub promptUser
{
#note: no error checking here!
	my $askString = shift;
	print $askString;
	my $answer = <STDIN>;
	chomp $answer;
	return $answer;
}


  
 sub make_my_id 
 {
           my $description_line = shift;
	   $description_line =~ /^>(.*)\..*/;	   
	   return $description_line;
  }

#because I can't get makeid to work in creating my database,
 #take suffixes off in new fasta file
 
 sub print_for_align 
 {
   my $inHashRef = shift;
   my $outSeqFile = shift;
   my $success = shift;
   my %newHash = ();
   my @oldKeys = keys(%$inHashRef);
   
   open(OUTFILE, ">>$outSeqFile") or warn "I can't print to $outSeqFile:$! \n";
   
   foreach (@oldKeys)
   {
   
	if ($cutIdsShort)
	{   
	my @newId = split /\./;
	print OUTFILE ">$newId[0]\n";
	}
	else { print OUTFILE ">$_\n"; }
	
	#print OUTFILE "${%$inHashRef}{$_}\n";
	print OUTFILE $$inHashRef{$_} . "\n";
   }
    close (OUTFILE);
    print "\n\nFinished adding $success sequences to $outSeqFile\n\n";
   return ;
 
 }
 
 sub isProgInstalled {
 
	my $program = shift;
 	my @PATH=split(":","$ENV{'PATH'}");     #Get users Path
	my $found = 0;
	foreach my $path (@PATH) 
	{
   	  if (-x "$path/$program") 
   	  {
   		  $found= 1;
   		  last;
    	  }
	}
	
	if (($program =~ /^clustalw/) or ($program =~ /^muscle/))
	{
        	unless ($found) {
                	print "\nThe program $program cannot be found on your system\n";
			print "Have you installed it?\n\n";
			print "\nPlease install $program or try another program\n";
			print "for alignment next time\n\n\n";
			exit;
                
		}        
 	}	
     return $found;
 
 }
 
 
 sub align_seqs 
 {
   my $seqFile = shift;
   my $outAlignFile = shift;
   my $alignMethod = shift;
   
   my $factory=();
   my $runViewer = 0; 
 
 # Build an alignment factory of the right type
	if ($alignMethod eq "muscle") 
	{
		$factory = new Bio::Tools::Run::Alignment::Muscle (-outfile_name => $outAlignFile);
	}
	elsif ($alignMethod eq "clustalw")
	{
		$factory = new Bio::Tools::Run::Alignment::Clustalw ('OUTFILE' => $outAlignFile);
	}
  # $aln is a SimpleAlign object.
  eval
  {
  my $aln = $factory->align($seqFile);
  };
  if ($@)  #well, this lurvely error catching works for muscle, but not clustalw. Grrrr.
#hence double error catching by looking for number of successes = 0 before getting into this 
#function in the first place.
  {
     print "\nThere was a problem aligning your new sequence file. \n";
     print "\nPlease check that the list of sequence ids you provided contains some sequences\n";
     print "\npresent in your original sequence file\n";
   }
  else { $runViewer = 1; }
return $runViewer;
}
