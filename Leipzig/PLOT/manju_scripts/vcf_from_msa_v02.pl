#!/usr/bin/perl
use List::Util qw[min max];

$file = $ARGV[0];
open(Input,"<","$file");
open (OUT,">","$file.vcf_v02.out");
#open(Input,"<","ch");
@array=<Input>;
foreach $y(@array)
{
chomp $y;
}
print "The file is read and now processing...\n";
#print "##H_chr\t H_strd\tH_base\tC_chr\tC_strd\tC_base\tb_chr\tb_strd\tb_base\tgo_chr\tgo_strd\tgo_base\tgi_chr\tgi_strd\tgi_base\tor_chr\tor_strd\tor_base\tm_chr\tm_strd\tm_base\tr_chr\tr_strd\tr_base\n";
$k=0;
#print OUT  "score\th_name\th_sign\thst\thu\tch_name\tch_sign\tchst\tch\tbo_name\tbo_sign\tbost\tbo\tgo_name\tgo_sign\tgost\tgo\tor_name\tor_sign\torst\tor\tgi_name\tgi_sign\tgist\tgi\trh_name\trh_sign\trhst\trh\tmr_name\tmr_sign\tmrst\tmr\n";
print OUT  "score\th_name\thst\thu\tch_name\tchst\tch\tbo_name\tbost\tbo\tgo_name\tgost\tgo\tor_name\torst\tor\tgi_name\tgist\tgi\trh_name\trhst\trh\tmr_name\tmrst\tmr\n";
foreach $line(@array)
{
@H=();@bo=();@ch=();@go=();@gi=();@or=();@mr=();@rh=();
if ($array[$k] =~ m/a score/)
{
$c=0;
	$score=$array[$k];
	@sc=split /=/,$score;
		$tmp=$k;
		while($array[$tmp])
		{
		    if ($array[$tmp] =~ m/^s/ )
			{ 	$c=$c+1; 
				@V=split /\s+/,$array[$tmp];
				@name=split /\./,$V[1];
				if ($name[0] =~ m/hg19/)  { $human_start=$V[2]; @H=split //,$V[6];@hx=split /chr/,$V[1];$h_name=$hx[1];$h_sign=$V[4];} 
				if ($name[0] =~ m/bonobo/)  { $bonobo_start=$V[2]; @bo=split //,$V[6];@box=split /chr/,$V[1];$bo_name=$box[1]; $bo_sign=$V[4];} 
                                if ($name[0] =~ m/gorGor3/) { $gorilla_start=$V[2]; @go=split //,$V[6];@gox=split /chr/,$V[1];$go_name=$gox[1]; $go_sign=$V[4];} 
                                if ($name[0] =~ m/ponAbe2/) { $orangutan_start=$V[2]; @or=split //,$V[6];@orx=split /chr/,$V[1]; $or_name=$orx[1]; $or_sign=$V[4];} 
                                if ($name[0] =~ m/calJac3/) { $marmoset_start=$V[2]; @mr=split //,$V[6]; @mrx=split /chr/,$V[1];$mr_name=$mrx[1]; $mr_sign=$V[4];} 
                                if ($name[0] =~ m/rheMac3/) { $rhesus_start=$V[2]; @rh=split //,$V[6]; @rhx=split /chr/,$V[1];$rh_name=$rhx[1]; $rh_sign=$V[4];} 
                                if ($name[0] =~ m/nomLeu1/) { $gibbon_start=$V[2]; @gi=split //,$V[6]; @gix=split /\./,$V[1];$gi_name=$gix[1]; $gi_sign=$V[4];} 
                                if ($name[0] =~ m/panTro/)  { $chimp_start=$V[2]; @ch=split //,$V[6]; @chx=split /chr/,$V[1]; $ch_name=$chx[1]; $ch_sign=$V[4];} 
			}			
			$tmp++;
		}
$a=$b=$c=$d=$e=$f=$g=$h=0;
$h_first=-1;$ch_first=-1;$go_first=-1; $h_end=-1;$ch_end=-1;$go_end=-1; $bo_first=-1;$bo_end=-1;
$mr_first=-1;$mr_end=-1;$rh_first=-1;$rh_end=-1;$gi_first=-1;$gi_end=-1;$or_first=-1;$or_end=-1;

$hst=$human_start-1; $chst=$chimp_start-1;$gost=$gorilla_start-1;$bost=$bonobo_start-1;
$orst=$orangutan_start-1;$mrst=$marmoset_start-1;$rhst=$rhesus_start-1;$gist= $gibbon_start-1;
$length_array=$#H+($#H/2); 
$flag=0;
foreach $i(0..$length_array)
	{
	if ($H[$a] =~ m/[-]/ or $ch[$b]  =~ m/[-]/ or $go[$c] =~ m/[-]/ or $bo[$d] =~ m/[-]/ or $gi[$e] =~ m/[-]/ or $mr[$f] =~ m/[-]/ or $rh[$g] =~ m/[-]/ or $or[$h] =~ m/[-]/)
		{
			if ($H[$a] =~ m/[-]/)  { if ($h_first == -1)  { $h_first =$a; } else { $h_end= $a;  } if($a == $#H) { $h_end= $a; $flag=1;goto ELSE;  } } else { $hst++;}  
		        if ($ch[$b] =~ m/[-]/) { if ($ch_first == -1) { $ch_first =$b;} else { $ch_end= $b; } if($b == $#H) { $ch_end= $b;$flag=1; goto ELSE; } } else {$chst++;}
			if ($go[$c] =~ m/[-]/) { if ($go_first == -1) { $go_first =$c;} else { $go_end= $c; } if($c == $#H) { $go_end= $c;$flag=1; goto ELSE; } } else {$gost++;}
		        if ($bo[$d] =~ m/[-]/) { if ($bo_first == -1) { $bo_first =$d;} else { $bo_end= $d; } if($d == $#H) { $bo_end= $d;$flag=1; goto ELSE; } } else {$bost++;}
		        if ($gi[$e] =~ m/[-]/) { if ($gi_first == -1) { $gi_first =$e;} else { $gi_end= $e; } if($e == $#H) { $gi_end= $e;$flag=1; goto ELSE; } } else {$gist++;}
		        if ($mr[$f] =~ m/[-]/) { if ($mr_first == -1) { $mr_first =$f;} else { $mr_end= $f; } if($f == $#H) { $mr_end= $f; $flag=1;goto ELSE; } } else {$mrst++;}
		        if ($rh[$g] =~ m/[-]/) { if ($rh_first == -1) { $rh_first =$g;} else { $rh_end= $g; } if($g == $#H) { $rh_end= $g; $flag=1;goto ELSE; } } else {$rhst++;}
		        if ($or[$h] =~ m/[-]/) { if ($or_first == -1) { $or_first =$h;} else { $or_end= $h; } if($h == $#H) { $or_end= $h; $flag=1; goto ELSE; } } else {$orst++;}

		}
	else { 
		ELSE:
		if($h_first!=-1 or $ch_first!=-1 or $go_first!=-1 or $bo_first!=-1 or $gi_first!=-1 or $mr_first!=-1 or $rh_first!=-1 or $or_first!=-1 )
		{
			$h_len=$h_end-$h_start; $ch_len=$ch_end-$ch_start; $go_len=$go_end-$go_first; $gi_len=$gi_end-$gi_first;
			$bo_len=$bo_end-$bo_start; $mr_len=$mr_end-$mr_start; $rh_len=$rh_end-$rh_first; $or_len=$or_end-$or_first;
			$max_len = max($h_len,$ch_len,$go_len,$bo_len,$mr_len,$rh_len,$or_len,$gi_len);
			$end = max($h_end,$ch_end,$go_end,$bo_end,$mr_end,$rh_end,$gi_end,$or_end);
			$start=10000000000000000000;
			@names= ($h_first,$ch_first,$go_first,$bo_first,$mr_first,$rh_first,$gi_first,$or_first) ;
			foreach $y (@names){
			if( $y != -1){
				if($y < $start)
				{ $start=$y ; }
				} }
			print OUT  "$sc[1]\t";
#human
                        #if (@H){ print OUT "$h_name\t$h_sign"; $hst_tmp=$hst-1;
                        if (@H){ print OUT "$h_name"; $hst_tmp=$hst-1;
                                if($h_first != -1) {
                                                if($h_end == -1 && $end == -1 ) {
                                                         print OUT "\t$hst\t$H[$h_first-1]$H[$h_first]"; }
                                                else {
                                                        if ($h_end>$end) { print OUT "\t$hst\t";
                                                                for ($w=$h_first-1;$w<=$h_end;$w++) { print OUT "\t$H[$w]"; } }
                                                        else { $new_pos= $hst-(($end-$start)-($h_end-$h_first)); print OUT "\t$new_pos\t";
                                                                for ($w=$start-1;$w<=$end;$w++) { print OUT "$H[$w]"; } }
                                                     }
                                                   }
                                else {
                                                 if($end == -1) {
                                                         print OUT "\t$hst_tmp\t$H[$start-1]$H[$start]"; }
                                                else { $pos= $hst-($end-$start); $pos_tmp= $pos-1; print OUT "\t$pos_tmp\t";
                                                for ($w=$start-1;$w<=$end;$w++) { print OUT "$H[$w]"; }
                                                        }
                                        }
                                $a--; $hst--; } 
#chimp
                        #if(@ch) { print OUT "\t$ch_name\t$ch_sign"; $chst_tmp=$chst-1;
                        if(@ch) { print OUT "\t$ch_name"; $chst_tmp=$chst-1;
                                if($ch_first != -1) {
                                                if($ch_end == -1 && $end == -1 ) {
                                                         print OUT "\t$chst\t$ch[$ch_first-1]$ch[$ch_first]"; }
                                                else {
                                                        if ($ch_end>$end) { print OUT "\t$chst-$a-$i\t";
                                                                for ($w=$ch_first-1;$w<=$ch_end;$w++) { print OUT "$ch[$w]"; } }
                                                        else {  $new_pos= $chst-(($end-$start)-($ch_end-$ch_first)); print OUT "\t$new_pos\t";
                                                        for ($qw=$start-1;$qw<=$end;$qw++) { print OUT "$ch[$qw]"; } }
                                                        }
                                     		    }
                                else {
                                                 if($end == -1) {
                                                         print OUT "\t$chst_tmp\t$ch[$start-1]$ch[$start]"; }
                                                else { $pos= $chst-($end-$start); $pos_tmp= $pos-1; print OUT "\t$pos_tmp\t";
                                                for ($qw=$start-1;$qw<=$end;$qw++) { print OUT "$ch[$qw]"; }
                                                        }
                                     }
                                $b--; $chst--;}
#bonobo
#                        if(@bo) { print OUT "\t$bo_name\t$bo_sign"; $bost_tmp=$bost-1;
                        if(@bo) { print OUT "\t$bo_name"; $bost_tmp=$bost-1;
                                if($bo_first != -1) {
                                                if($bo_end == -1 && $end == -1 ) {
                                                         print OUT "\t$bost\t$bo[$bo_first-1]$bo[$bo_first]"; }
                                                else {
                                                        if ($bo_end>$end) { print OUT "\t$bost\t";
                                                                for ($w=$bo_first-1;$w<=$bo_end;$w++) { print OUT "$bo[$w]"; } }
                                                        else {  $new_pos= $bost-(($end-$start)-($bo_end-$bo_first)); print OUT "\t$new_pos\t";
                                                        for ($qw=$start-1;$qw<=$end;$qw++) { print OUT "$bo[$qw]"; } }
                                                        }
                                        }
                                else {
                                                 if($end == -1) {
                                                         print OUT "\t$bost_tmp\t$bo[$start-1]$bo[$start]"; }
                                                else { $pos= $bost-($end-$start); $pos_tmp= $pos-1; print OUT "\t$pos_tmp\t";
                                                for ($qw=$start-1 ;$qw<=$end;$qw++) { print OUT "$bo[$qw]"; }
                                                        }
                                        }
                                    $c--;$bost--; } 
#gorilla
#                                if(@go) {print OUT "\t$go_name\t$go_sign"; $gost_tmp=$gost-1;
                                if(@go) {print OUT "\t$go_name"; $gost_tmp=$gost-1;
                                if($go_first != -1) {
                                                if($go_end == -1 && $end == -1 ) {
                                                         print OUT "\t$gost\t$go[$go_first-1]$go[$go_first]"; }
                                                else {
                                                                if ($go_end>$end){ print OUT "\t$gost\t";
                                                                        for ($w=$go_first-1;$w<=$go_end;$w++) { print OUT "$go[$w]"; } }
                                                                else { $new_pos= $gost-(($end-$start)-($go_end-$go_first)); print OUT "\t$new_pos\t";
                                                                        for ($qw=$start-1;$qw<=$end;$qw++) { print OUT "$go[$qw]"; } }
                                                        }
                                        }
                                else {
                                                 if($end == -1) {
                                                         print OUT "\t$gost_tmp\t$go[$start-1]$go[$start]"; }
                                                else { $pos= $gost-($end-$start); $pos_tmp= $pos-1; print OUT "\t$pos_tmp\t";
                                                for ($qw=$start-1;$qw<=$end;$qw++) { print OUT "$go[$qw]"; }
                                                        }
                                        }
                                       $d--;$gost--; } 
#orangutan
#                                if(@or) { print OUT "\t$or_name\t$or_sign"; $orst_tmp=$orst-1;
                                if(@or) { print OUT "\t$or_name"; $orst_tmp=$orst-1;
                                if($or_first != -1) {
                                                if($or_end == -1 && $end == -1 ) {
                                                         print OUT "\t$orst\t$or[$or_first-1]$or[$or_first]"; }
                                                else {
                                                        if ($ch_end>$end) { print OUT "\t$orst\t";
                                                                for ($w=$or_first-1;$w<=$or_end;$w++) { print OUT "$or[$w]"; } }
                                                        else {  $new_pos= $orst-(($end-$start)-($or_end-$or_first)); print OUT "\t$new_pos\t";
                                                        for ($qw=$start-1;$qw<=$end;$qw++) { print OUT "$or[$qw]"; } }
                                                        }
                                        }
                                else {
                                                 if($end == -1) {
                                                         print OUT "\t$orst_tmp\t$or[$start-1]$or[$start]"; }
                                                else { $pos= $orst-($end-$start); $pos_tmp= $pos-1; print OUT "\t$pos_tmp\t";
                                                for ($qw=$start-1;$qw<=$end;$qw++) { print OUT "$or[$qw]"; }
                                                        }
                                        }
                                          $e--;$orst--; } 
#gibbon
#				if(@gi) { print OUT "\t$gi_name\t$gi_sign"; $gist_tmp=$gist-1;
				if(@gi) { print OUT "\t$gi_name"; $gist_tmp=$gist-1;
                                if($gi_first != -1) {
                                                if($gi_end == -1 && $end == -1 ) {
                                                         print OUT "\t$gist\t$gi[$gi_first-1]$gi[$gi_first]"; }
                                                else {
                                                        if ($gi_end>$end) { print OUT "\t$gist\t";
                                                                for ($w=$gi_first-1;$w<=$gi_end;$w++) { print OUT "$gi[$w]"; } }
                                                        else {  $new_pos= $gist-(($end-$start)-($gi_end-$gi_first)); print OUT "\t$new_pos\t";
                                                        for ($qw=$start-1;$qw<=$end;$qw++) { print OUT "$gi[$qw]"; } }
                                                        }
                                        }
                                else {
                                                 if($end == -1) {
                                                         print OUT "\t$gist_tmp\t$gi[$start-1]$gi[$start]"; }
                                                else { $pos= $gist-($end-$start); $pos_tmp= $pos-1; print OUT "\t$pos_tmp\t";
                                                for ($qw=$start-1;$qw<=$end;$qw++) { print OUT "$gi[$qw]"; }
                                                        }
                                        }
					 $f--; $gist--; } 
#rhesus
#                                if(@rh) { print OUT "\t$rh_name\t$rh_sign"; $rhst_tmp=$rhst-1;
                                if(@rh) { print OUT "\t$rh_name"; $rhst_tmp=$rhst-1;
                                if($rh_first != -1) {
                                                if($rh_end == -1 && $end == -1 ) {
                                                         print OUT "\t$rhst\t$rh[$rh_first-1]$rh[$rh_first]"; }
                                                else {
                                                        if ($rh_end>$end) { print OUT "\t$rhst\t";
                                                                for ($w=$rh_first-1;$w<=$rh_end;$w++) { print OUT "$rh[$w]"; } }
                                                        else {  $new_pos= $rhst-(($end-$start)-($rh_end-$rh_first)); print OUT "\t$new_pos\t";
                                                        for ($qw=$start-1;$qw<=$end;$qw++) { print OUT "$rh[$qw]"; } }
                                                        }
                                        }
                                else {
                                                 if($end == -1) {
                                                         print OUT "\t$rhst_tmp\t$rh[$start-1]$rh[$start]"; }
                                                else { $pos= $rhst-($end-$start);  $pos_tmp= $pos-1; print OUT "\t$pos_tmp\t";
                                                for ($qw=$start-1;$qw<=$end;$qw++) { print OUT "$rh[$qw]"; }
                                                        }
                                        }
                                         $g--;$rhst--; } 

#marmoset
#				if(@mr) { print OUT "\t$mr_name\t$mr_sign";$mrst_tmp=$mrst-1;
				if(@mr) { print OUT "\t$mr_name";$mrst_tmp=$mrst-1;
                                if($mr_first != -1) {
                                                if($mr_end == -1 && $end == -1 ) {
                                                         print OUT "\t$mrst\t$mr[$mr_first-1]$mr[$mr_first]"; }
                                                else {
                                                        if ($mr_end>$end) { print OUT "\t$mrst\t";
                                                                for ($w=$mr_first-1;$w<=$mr_end;$w++) { print OUT "$mr[$w]"; } }
                                                        else {  $new_pos= $mrst-(($end-$start)-($mr_end-$mr_first)); print OUT "\t$new_pos\t";
                                                        for ($qw=$start-1;$qw<=$end;$qw++) { print OUT "$mr[$qw]"; } }
                                                        }
                                        }
                                else {
                                                 if($end == -1) {
                                                         print OUT "\t$mrst_tmp\t$mr[$start-1]$mr[$start]"; }
                                                else { $pos= $mrst-($end-$start); $pos_tmp= $pos-1; print OUT "\t$pos\t";
                                                for ($qw=$start-1;$qw<=$end;$qw++) { print OUT "$mr[$qw]"; }
                                                        }
                                        }
					 $h--;$mrst--; } 

		print OUT "\n";
		$h_first=-1;$ch_first=-1;$h_end=-1;$ch_end=-1;$go_first=-1; $go_end=-1; $bo_first=-1;$bo_end=-1;
		$mr_first=-1;$mr_end=-1;$rh_first=-1;$rh_end=-1;$gi_first=-1;$gi_end=-1;$or_first=-1;$or_end=-1;
		if($flag==1) { goto ENDD;}
		} # closing of if after the tag ELSE:
          	else 
		{
		$h_h=$hst+1; $c_c=$chst+1; $b_b=$bost+1; $o_o=$orst+1; $g_g=$gost+1;$gi_gi=$gist+1; $r_r=$rhst+1; $m_m=$mrst+1;
#		print OUT "$sc[1]\t$h_name\t$h_sign\t$h_h\t$H[$a]\t$ch_name\t$ch_sign\t$c_c\t$ch[$b]\t$bo_name\t$bo_sign\t$b_b\t$bo[$c]\t$go_name\t$go_sign\t$g_g\t$go[$d]\t$or_name\t$or_sign\t$o_o\t$or[$e]\t$gi_name\t$gi_sign\t$gi_gi\t$gi[$f]\t$rh_name\t$rh_sign\t$r_r\t$rh[$g]\t$mr_name\t$mr_sign\t$m_m\t$mr[$h]\t \n";
	print OUT "$sc[1]\t$h_name\t$h_h\t$H[$a]\t$ch_name\t$c_c\t$ch[$b]\t$bo_name\t$b_b\t$bo[$c]\t$go_name\t$g_g\t$go[$d]\t$or_name\t$o_o\t$or[$e]\t$gi_name\t$gi_gi\t$gi[$f]\t$rh_name\t$r_r\t$rh[$g]\t$mr_name\t$m_m\t$mr[$h]\t \n";
		}
		$hst++; $chst++;$gost++;$gist++;$orst++;$mrst++;$rhst++;$bost++;
#		print "$a with $H[$a+1] \n";
		if(!$H[$a+1]) { goto ENDD;}
	     } # closing of else after the ELSE:
	$a++; $b++; $c++;$d++;$e++;$f++;$g++;$h++; 
	} # closing of the for block
}# closing of the main if block
ENDD:
$k=$k+1;
} # closing of the for loop


