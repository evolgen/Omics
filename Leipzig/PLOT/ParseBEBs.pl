use strict;
use warnings;

#Purpose: Script to parse out the positions of positively selected sites on a full-length reference sequence (here human) 
#Call: perl ~/TFs_Katja/POTION/POTION_hack/ParseBEBs.pl Group1219 model2 /scr/k61san/henrike/TFCarlos/BranchSiteModel/Species_filter0/intermediate_files
#Outputs: File1: list of positions on reference and on trimmed seq + BEB posterior probabilities
#         File2: reference sequence with positively selected positions in lower case (at the moment standard out)

my $group = $ARGV[0];
my $info = $ARGV[1];
my $path = $ARGV[2];
#my $group = "Group1219";
#my $info = "model2";
#my $path = "/scr/k61san/henrike/TFCarlos/BranchSiteModel/Species_filter0/intermediate_files";

`cd $path/$group/`;
my $results = "$group.cluster.aa.fa.aln.nt.phy.trim.paml.${info}";
my $alignment = "$group.cluster.aa.fa.aln";
my $output = "${group}_BEB.out";

open(RESULTS, "<", $results) or die "Could not parse $results. $!\n";
open(OUT, ">", $output) or die "Could not open $output. $!\n";

my $switch = "0";
my @Input;
while(<RESULTS>){
if ($_ =~ /Bayes Empirical Bayes/){
$switch = "1";
}
if ($switch eq "1" && $_ =~ /\*/ && $_ =~ /\+/){
    chomp $_;
    #print "$_\n";
    my @array = split '\s+', $_;
    push @Input, "$array[1],$array[2],$array[3]";
    #print "$array[1]\n";
    }
}
#print "@Input\n";
#print "$Input[1]\n";
#`sed 's/,/\n/g' Group1219.trimal_cols.aa | grep '[0-9]'| head -n 40| tail -n2 `;
my @PosAlign;

foreach my $position(@Input){
    my @ALL = split ',', $position;
    my $Pos = $ALL[0];
    my $posalign = `sed -e 's/,/\\n/g' *.trimal_cols.aa | grep '[0-9]'| head -n $Pos | tail -n1| sed 's/\\s//g'`;
    chomp $posalign;
    push @PosAlign, $posalign;
}
#print "@PosAlign\n";

my $ref = `head -n2 $alignment| tail -n1`;
my $refname = `head -n1 $alignment| sed 's/>//g'`;
chomp $refname;
my $refspec = `grep -w '$refname' id_names`;
chomp $refspec;
my @refseq = split '', $ref;
#print "@ref\n";
print OUT "#$refname\t$refspec\n";
print OUT "#PosOnTrimmedSeq\tAA\tBEBpp\tAA\tPosOnReference\n";

print ">$refname\n";
my $counter = 0;
my $tracker = 0;
#print "$counter\n";
foreach my $i(0 .. $#refseq) {
    if ($refseq[$i] eq "-") {
        $counter++;
    };
    if (grep{$_ eq $i}@PosAlign) {    
        my @array2 = split ',', $Input[$tracker];
        #print "@array2\n";
        #print "$i\n";
        print OUT "$array2[0]\t$array2[1]\t$array2[2]\t$refseq[$i]\t";
        print OUT $i - $counter . "\n";
        $tracker ++;
        print lc($refseq[$i]);
    }else{print "$refseq[$i]"};

};

