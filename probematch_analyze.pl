#!/usr/bin/perl
if (!$ARGV[0]) {
	print STDERR "\n USAGE: probematch_analyze.pl mode prombematch_file db_file outfile\n\n";
	print STDERR "\n\tmodes:\n";
	print STDERR "\n\tampli  : prombematch_file db_file outfile\n";
	print STDERR "\n\tprobe : prombematch_file db_file outfile\n";
	print STDERR "\n\ttrflp : outfile_of_this_script restriction_enzyme|probe\n";
	exit;
}
$mode = $ARGV[0];
$pm_infile = $ARGV[1];
$db_infile = $ARGV[2];
$outfile = $ARGV[3];
goto TRFLP if ($mode eq 'trflp');
open(OUT1,">$outfile.summary");
open(OUT2,">$outfile.amplicons") if ($mode eq 'ampli');
open(IN2,$pm_infile) or die "ERROR: Cannot open $pm_infile\n"; #probematch out
open(IN1,$db_infile) or die "ERROR: Cannot open $db_infile\n"; #db
#goto HERE if ($mode eq 'ampli');
$cnt_db = 0;
while($line = <IN1>) {
	next if ($line !~ />/);
	chomp $line;
	$line =~ s/.+?\t//;
#	$debug = 1;
#	$debug = 1 if ($line =~ /Fuchsiella/);
	$line =~ s/\"//g;
	$cnt_db++;
	@tmp = split (/;/,$line);
#	print STDERR join "+", @tmp,"\n";
	for($l = 3;$l<=$#tmp;$l+=2) {
#		print STDERR "$tmp[$l] -> $tmp[$l-1]\n" if ($debug);
#		<STDIN> if ($debug);
		if ($tmp[$l] eq 'genus') {
			$genus_all{$tmp[$l-1]}++;
		}
		if ($tmp[$l] eq 'family') {
			$family_all{$tmp[$l-1]}++;
		}
		if ($tmp[$l] eq 'order') {
			$order_all{$tmp[$l-1]}++;
		}
		if ($tmp[$l] eq 'class') {
			$class_all{$tmp[$l-1]}++;
		}
		if ($tmp[$l] eq 'phylum') {
			$phylum_all{$tmp[$l-1]}++;
		}
		if ($tmp[$l] eq 'domain') {
			$domain_all{$tmp[$l-1]}++;
		}
	}
	print STDERR "Loading DB line $cnt_db\r" if ($cnt_db % 1000 == 0); 
}
close IN1;
print STDERR "Loading DB line $cnt_db\n\n"; 
print STDERR "DOMAINS : ",scalar keys %domain_all,"\n";
print STDERR "PHYLA   : ",scalar keys %phylum_all,"\n";
print STDERR "CLASSES : ",scalar keys %class_all,"\n";
print STDERR "ORDERS  : ",scalar keys %order_all,"\n";
print STDERR "FAMILIES: ",scalar keys %family_all,"\n";
print STDERR "GENERA  : ",scalar keys %genus_all,"\n";
#exit;

SCAN:
#seqname	desc	primer_index	primer_name	position	mismatches	seq_primer_region
$cnt_pm = 0;
%cnt = ();
while($line = <IN2>) {
	chomp $line;
	$line =~ s/\"//g;
	@tmp = split (/\t/,$line);
	$pos = $tmp[$#tmp-2];
	$cnt_pm++;
	next if ($cnt_pm == 1);
#	print STDERR "$tmp[0] -> $pos -> $tmp[2]\n";
#	<STDIN>;
	if($cnt{$tmp[0]} || $mode eq 'probe') {
		$lineage{$tmp[0]} = $tmp[2];
		@{$ampli{$tmp[0]}} = ($pos_old,$len_old,$pos,length($tmp[$#tmp])); #i.e. start of fwd, start of rev, length of rev 
		print STDERR scalar keys %ampli, " amplicons / ",$cnt_pm," probe hits\r" if ($mode eq 'ampli');
		print STDERR $cnt_pm," probe hits\r" if ($mode eq 'probe');
#		<STDIN>;
	}
	$cnt{$tmp[0]} = 1;
	$pos_old = $pos;
	$len_old = length($tmp[$#tmp]);
}
close IN2;

foreach $s (keys %lineage) {
	@tmp = split (/;/,$lineage{$s});
	for($l = 1;$l<=$#tmp;$l+=2) {
#		print STDERR "$tmp[$l] -> $tmp[$l-1]\n";
#		<STDIN>;
		if ($tmp[$l] eq 'genus') {
			$genus_cnt{$tmp[$l-1]}++;
		}
		if ($tmp[$l] eq 'family') {
			$family_cnt{$tmp[$l-1]}++;
		}
		if ($tmp[$l] eq 'order') {
			$order_cnt{$tmp[$l-1]}++;
		}
		if ($tmp[$l] eq 'class') {
			$class_cnt{$tmp[$l-1]}++;
		}
		if ($tmp[$l] eq 'phylum') {
			$phylum_cnt{$tmp[$l-1]}++;
		}
		if ($tmp[$l] eq 'domain') {
			$domain_cnt{$tmp[$l-1]}++;
		}
	}
}

print STDERR "\n";
print STDERR "DOMAINS : ",scalar keys %domain_cnt,"\n";
print STDERR "PHYLA   : ",scalar keys %phylum_cnt,"\n";
print STDERR "CLASSES : ",scalar keys %class_cnt,"\n";
print STDERR "ORDERS  : ",scalar keys %order_cnt,"\n";
print STDERR "FAMILIES: ",scalar keys %family_cnt,"\n";
print STDERR "GENERA  : ",scalar keys %genus_cnt,"\n";

print STDERR "Writing summary to $outfile.summary\n";

print OUT1 "#AMPLICONS\t",scalar keys %ampli, "\t",$cnt_db,"\n" if ($mode eq 'ampli');
print OUT1 "#MATCHES\t",scalar keys %ampli, "\t",$cnt_db,"\n" if ($mode eq 'probe');
print OUT1 "#DOMAINS\t",scalar keys %domain_cnt,"\t",scalar keys %domain_all,"\n";
print OUT1 "#PHYLA\t",scalar keys %phylum_cnt,"\t",scalar keys %phylum_all,"\n";
print OUT1 "#CLASSES\t",scalar keys %class_cnt,"\t",scalar keys %class_all,"\n";
print OUT1 "#ORDERS\t",scalar keys %order_cnt,"\t",scalar keys %order_all,"\n";
print OUT1 "#FAMILIES\t",scalar keys %family_cnt,"\t",scalar keys %family_all,"\n";
print OUT1 "#GENERA\t",scalar keys %genus_cnt,"\t",scalar keys %genus_all,"\n";

foreach $r (sort keys %domain_all) {
	print OUT1 "DOMAIN\t",$r,"\t",$domain_cnt{$r},"\t",$domain_all{$r},"\t",sprintf("%.2f",$domain_cnt{$r}/$domain_all{$r}*100),"\n" if ($domain_cnt{$r});
	print OUT1 "DOMAIN*\t$r\t0\t",$domain_all{$r},"0\n" if (!$domain_cnt{$r});
}
foreach $r (sort keys %phylum_all) {
	print OUT1 "PHYLUM\t",$r,"\t",$phylum_cnt{$r},"\t",$phylum_all{$r},"\t",sprintf("%.2f",$phylum_cnt{$r}/$phylum_all{$r}*100),"\n" if ($phylum_cnt{$r});
	print OUT1 "PHYLUM*\t$r\t0\t",$phylum_all{$r},"0\n" if (!$phylum_cnt{$r});
}
foreach $r (sort keys %class_all) {
	print OUT1 "CLASS\t",$r,"\t",$class_cnt{$r},"\t",$class_all{$r},"\t",sprintf("%.2f",$class_cnt{$r}/$class_all{$r}*100),"\n" if ($class_cnt{$r});
	print OUT1 "CLASS*\t$r\t0\t",$class_all{$r},"0\n" if (!$class_cnt{$r});
}
foreach $r (sort keys %order_all) {
	print OUT1 "ORDER\t",$r,"\t",$order_cnt{$r},"\t",$order_all{$r},"\t",sprintf("%.2f",$order_cnt{$r}/$order_all{$r}*100),"\n" if ($order_cnt{$r});
	print OUT1 "ORDER*\t$r\t0\t",$order_all{$r},"0\n" if (!$order_cnt{$r});
}
foreach $r (sort keys %family_all) {
	print OUT1 "FAMILY\t",$r,"\t",$family_cnt{$r},"\t",$family_all{$r},"\t",sprintf("%.2f",$family_cnt{$r}/$family_all{$r}*100),"\n" if ($family_cnt{$r});
	print OUT1 "FAMILY*\t$r\t0\t",$family_all{$r},"0\n" if (!$family_cnt{$r});
}
foreach $r (sort keys %genus_all) {
	print OUT1 "GENUS\t",$r,"\t",$genus_cnt{$r},"\t",$genus_all{$r},"\t",sprintf("%.2f",$genus_cnt{$r}/$genus_all{$r}*100),"\n" if ($genus_cnt{$r});
	print OUT1 "GENUS*\t$r\t0\t",$genus_all{$r},"0\n" if (!$genus_cnt{$r});
}
close OUT1;
HERE:
if ($mode eq 'ampli') {
	print STDERR "Writing amplicons to $outfile.amplicons\n";
#	@{$ampli{"S000494589"}} = (21,20,100,20);
	open(IN,$db_infile) or die "ERROR: Cannot open $db_infile\n";
	while($line = <IN>) {
		chomp $line;
		$line =~ s/\"//g;
		if ($line =~ />/) {
			if ($name && $ampli{$name}) {
				$st = $ampli{$name}[0]-$ampli{$name}[1];
				$en = $ampli{$name}[2];
				$len = $en-$st;
				$ampli = substr($seq,$st,$len); #e.g. 350 780 25 => substr from 350 length 780-350+25 = 445
				print OUT2 "$preline|$st:$en:$len\n",$ampli,"\n"; #preline is the previous line, i.e. that form which the name has been extracted
			}
			$preline = $line;
			$name = $line;
			$name =~ s/ .+//;
			$name =~ s/>//;
			$seq = '';
		} else {
			$seq .= $line;
		}
	}
	close IN;
	close OUT2;
}
exit;

TRFLP:
&sites;
open(IN,"$pm_infile") or die "No amplicon file with this name\n";
if ($site{$db_infile}) {
	$site = uc($site{$db_infile}) 
} else {
	print STDERR "Unknown enzyme\n\n";
	exit;
}
#$site = "GGG";
while ($site =~ /(.)/g) {
	$site_ok .= $deg{$1};
}
$site = $site_ok;
print "Cutting with $outfile ($site)\n";
open(OUT,">$db_infile.trflp");
while ($line = <IN>) {
	chomp $line;
	if ($line =~ />/) {
		if ($seq) {
			if ($seq =~ /$site/) {
				@parts = split (/$site/,$seq);
				$ls = $parts[0];
				$rs = $parts[$#parts];
				$lsize_cnt{length($ls)}++;
				$rsize_cnt{length($rs)}++;
				#print $preline,":",length($ls),"-",length($rs),"\n$seq\n",$ls,"\n",$rs,"\n";
				#<STDIN>;
			} else {
#				print $preline,"\t0\t0\n";
			}
		}
		$seq = '';
		$preline = $line;
	} else {
		$seq .= uc($line);
	}		
}
foreach $s (sort{$a <=> $b} keys %lsize_cnt) {
	print OUT "L\t",$s,"\t",$lsize_cnt{$s},"\n";
}
foreach $s (sort{$a <=> $b} keys %rsize_cnt) {
	print OUT "R\t",$s,"\t",$rsize_cnt{$s},"\n";
}
close OUT;

open(OUT, ">script.R");
print OUT "d<-read.table(\"$db_infile.trflp\",sep=\"\\t\")\n";
print OUT "dl<-d[d[,1]==\"L\",]\n";
print OUT "dr<-d[d[,1]==\"R\",]\n";
print OUT "pdf(\"$db_infile.trflp.pdf\",width=20,height=10)\n";
print OUT "par(mfrow=c(2,1))\n";
print OUT "plot(dl[,2],dl[,3],main=\"t-RFLP Left side\",type=\"l\",ylab=\"Intensity\",xlab = \"Position\")\n";
print OUT "plot(dr[,2],dr[,3],main=\"t-RFLP Right side\",type=\"l\",ylab=\"Intensity\",xlab = \"Position\")\n";
print OUT "dev.off()\n";
close OUT;
`R --lave --vanilla < script.R`;
unlink "script.R";
sub sites {
$site{'AclI'} = 'AACGTT';
$site{'HindIII'} = 'AAGCTT';
$site{'HindIII-HF'} = 'AAGCTT';
$site{'SspI'} = 'AATATT';
$site{'SspI-HF'} = 'AATATT';
$site{'MluCI'} = 'AATT';
$site{'PciI'} = 'ACATGT';
$site{'AgeI'} = 'ACCGGT';
$site{'AgeI-HF'} = 'ACCGGT';
$site{'BspMI'} = 'ACCTGC';
$site{'BfuAI'} = 'ACCTGC';
$site{'SexAI'} = 'ACCWGGT';
$site{'MluI'} = 'ACGCGT';
$site{'MluI-HF'} = 'ACGCGT';
$site{'BceAI'} = 'ACGGC';
$site{'HpyCH4IV'} = 'ACGT';
$site{'HpyCH4III'} = 'ACNGT';
$site{'BaeI'} = 'ACNNNNGTAYC';
$site{'BsaXI'} = 'ACNNNNNCTCC';
$site{'AflIII'} = 'ACRYGT';
$site{'SpeI-HF'} = 'ACTAGT';
$site{'SpeI'} = 'ACTAGT';
$site{'BsrI'} = 'ACTGG';
$site{'BmrI'} = 'ACTGGG';
$site{'BglII'} = 'AGATCT';
$site{'AfeI'} = 'AGCGCT';
$site{'AluI'} = 'AGCT';
$site{'StuI'} = 'AGGCCT';
$site{'ScaI-HF'} = 'AGTACT';
$site{'BspDI'} = 'ATCGAT';
$site{'ClaI'} = 'ATCGAT';
$site{'PI-SceI'} = 'ATCTATGTCGGGTGCGGAGAAAGAGGTAAT';
$site{'NsiI-HF'} = 'ATGCAT';
$site{'NsiI'} = 'ATGCAT';
$site{'AseI'} = 'ATTAAT';
$site{'SwaI'} = 'ATTTAAAT';
$site{'CspCI'} = 'CAANNNNNGTGG';
$site{'MfeI'} = 'CAATTG';
$site{'MfeI-HF'} = 'CAATTG';
$site{'Nb.BssSI'} = 'CACGAG';
$site{'BssSαI'} = 'CACGAG';
$site{'BmgBI'} = 'CACGTC';
$site{'PmlI'} = 'CACGTG';
$site{'DraIII-HF'} = 'CACNNNGTG';
$site{'AleI'} = 'CACNNNNGTG';
$site{'EcoP15I'} = 'CAGCAG';
$site{'PvuII-HF'} = 'CAGCTG';
$site{'PvuII'} = 'CAGCTG';
$site{'AlwNI'} = 'CAGNNNCTG';
$site{'BtsIMutI'} = 'CAGTG';
$site{'NdeI'} = 'CATATG';
$site{'FatI'} = 'CATG';
$site{'CviAII'} = 'CATG';
$site{'NlaIII'} = 'CATG';
$site{'MslI'} = 'CAYNNNNRTG';
$site{'FspEI'} = 'CC';
$site{'XcmI'} = 'CCANNNNNNNNNTGG';
$site{'BstXI'} = 'CCANNNNNNTGG';
$site{'PflMI'} = 'CCANNNNNTGG';
$site{'BccI'} = 'CCATC';
$site{'NcoI'} = 'CCATGG';
$site{'NcoI-HF'} = 'CCATGG';
$site{'BseYI'} = 'CCCAGC';
$site{'FauI'} = 'CCCGC';
$site{'SmaI'} = 'CCCGGG';
$site{'TspMI'} = 'CCCGGG';
$site{'XmaI'} = 'CCCGGG';
$site{'Nt.CviPII'} = 'CCD';
$site{'LpnPI'} = 'CCDG';
$site{'AciI'} = 'CCGC';
$site{'SacII'} = 'CCGCGG';
$site{'BsrBI'} = 'CCGCTC';
$site{'MspI'} = 'CCGG';
$site{'HpaII'} = 'CCGG';
$site{'ScrFI'} = 'CCNGG';
$site{'StyD4I'} = 'CCNGG';
$site{'BsaJI'} = 'CCNNGG';
$site{'BslI'} = 'CCNNNNNNNGG';
$site{'BtgI'} = 'CCRYGG';
$site{'NciI'} = 'CCSGG';
$site{'AvrII'} = 'CCTAGG';
$site{'MnlI'} = 'CCTC';
$site{'Nt.BbvCI'} = 'CCTCAGC';
$site{'Nb.BbvCI'} = 'CCTCAGC';
$site{'BbvCI'} = 'CCTCAGC';
$site{'SbfI'} = 'CCTGCAGG';
$site{'SbfI-HF'} = 'CCTGCAGG';
$site{'Bpu10I'} = 'CCTNAGC';
$site{'Bsu36I'} = 'CCTNAGG';
$site{'EcoNI'} = 'CCTNNNNNAGG';
$site{'HpyAV'} = 'CCTTC';
$site{'BstNI'} = 'CCWGG';
$site{'PspGI'} = 'CCWGG';
$site{'StyI'} = 'CCWWGG';
$site{'StyI-HF'} = 'CCWWGG';
$site{'BcgI'} = 'CGANNNNNNTGC';
$site{'PvuI'} = 'CGATCG';
$site{'PvuI-HF'} = 'CGATCG';
$site{'BstUI'} = 'CGCG';
$site{'EagI-HF'} = 'CGGCCG';
$site{'EagI'} = 'CGGCCG';
$site{'RsrII'} = 'CGGWCCG';
$site{'BsiEI'} = 'CGRYCG';
$site{'BsiWI'} = 'CGTACG';
$site{'BsiWI-HF'} = 'CGTACG';
$site{'BsmBI'} = 'CGTCTC';
$site{'Hpy99I'} = 'CGWCG';
$site{'MspA1I'} = 'CMGCKG';
$site{'MspJI'} = 'CNNR';
$site{'SgrAI'} = 'CRCCGGYG';
$site{'BfaI'} = 'CTAG';
$site{'BspCNI'} = 'CTCAG';
$site{'XhoI'} = 'CTCGAG';
$site{'PaeR7I'} = 'CTCGAG';
$site{'EarI'} = 'CTCTTC';
$site{'AcuI'} = 'CTGAAG';
$site{'PstI'} = 'CTGCAG';
$site{'PstI-HF'} = 'CTGCAG';
$site{'BpmI'} = 'CTGGAG';
$site{'DdeI'} = 'CTNAG';
$site{'SfcI'} = 'CTRYAG';
$site{'AflII'} = 'CTTAAG';
$site{'BpuEI'} = 'CTTGAG';
$site{'SmlI'} = 'CTYRAG';
$site{'AvaI'} = 'CYCGRG';
$site{'BsoBI'} = 'CYCGRG';
$site{'MboII'} = 'GAAGA';
$site{'BbsI'} = 'GAAGAC';
$site{'BbsI-HF'} = 'GAAGAC';
$site{'XmnI'} = 'GAANNNNTTC';
$site{'BsmI'} = 'GAATGC';
$site{'Nb.BsmI'} = 'GAATGC';
$site{'EcoRI'} = 'GAATTC';
$site{'EcoRI-HF'} = 'GAATTC';
$site{'HgaI'} = 'GACGC';
$site{'ZraI'} = 'GACGTC';
$site{'AatII'} = 'GACGTC';
$site{'PflFI'} = 'GACNNNGTC';
$site{'Tth111I'} = 'GACNNNGTC';
$site{'PshAI'} = 'GACNNNNGTC';
$site{'AhdI'} = 'GACNNNNNGTC';
$site{'DrdI'} = 'GACNNNNNNGTC';
$site{'Eco53kI'} = 'GAGCTC';
$site{'SacI'} = 'GAGCTC';
$site{'SacI-HF'} = 'GAGCTC';
$site{'BseRI'} = 'GAGGAG';
$site{'Nt.BstNBI'} = 'GAGTC';
$site{'PleI'} = 'GAGTC';
$site{'MlyI'} = 'GAGTC';
$site{'HinfI'} = 'GANTC';
$site{'EcoRV'} = 'GATATC';
$site{'EcoRV-HF'} = 'GATATC';
$site{'Sau3AI'} = 'GATC';
$site{'MboI'} = 'GATC';
$site{'DpnII'} = 'GATC';
$site{'DpnI'} = 'GATC';
$site{'BsaBI'} = 'GATNNNNATC';
$site{'TfiI'} = 'GAWTC';
$site{'BsrDI'} = 'GCAATG';
$site{'Nb.BsrDI'} = 'GCAATG';
$site{'BbvI'} = 'GCAGC';
$site{'BtsαI'} = 'GCAGTG';
$site{'Nb.BtsI'} = 'GCAGTG';
$site{'BstAPI'} = 'GCANNNNNTGC';
$site{'SfaNI'} = 'GCATC';
$site{'SphI'} = 'GCATGC';
$site{'SphI-HF'} = 'GCATGC';
$site{'SrfI'} = 'GCCCGGGC';
$site{'NmeAIII'} = 'GCCGAG';
$site{'NgoMIV'} = 'GCCGGC';
$site{'NaeI'} = 'GCCGGC';
$site{'BglI'} = 'GCCNNNNNGGC';
$site{'AsiSI'} = 'GCGATCGC';
$site{'BtgZI'} = 'GCGATG';
$site{'HhaI'} = 'GCGC';
$site{'HinP1I'} = 'GCGC';
$site{'BssHII'} = 'GCGCGC';
$site{'NotI'} = 'GCGGCCGC';
$site{'NotI-HF'} = 'GCGGCCGC';
$site{'Fnu4HI'} = 'GCNGC';
$site{'Cac8I'} = 'GCNNGC';
$site{'MwoI'} = 'GCNNNNNNNGC';
$site{'BmtI'} = 'GCTAGC';
$site{'BmtI-HF'} = 'GCTAGC';
$site{'NheI'} = 'GCTAGC';
$site{'NheI-HF'} = 'GCTAGC';
$site{'Nt.BspQI'} = 'GCTCTTC';
$site{'BspQI'} = 'GCTCTTC';
$site{'SapI'} = 'GCTCTTC';
$site{'BlpI'} = 'GCTNAGC';
$site{'TseI'} = 'GCWGC';
$site{'ApeKI'} = 'GCWGC';
$site{'Bsp1286I'} = 'GDGCHC';
$site{'AlwI'} = 'GGATC';
$site{'Nt.AlwI'} = 'GGATC';
$site{'BamHI-HF'} = 'GGATCC';
$site{'BamHI'} = 'GGATCC';
$site{'FokI'} = 'GGATG';
$site{'BtsCI'} = 'GGATG';
$site{'HaeIII'} = 'GGCC';
$site{'FseI'} = 'GGCCGGCC';
$site{'SfiI'} = 'GGCCNNNNNGGCC';
$site{'SfoI'} = 'GGCGCC';
$site{'KasI'} = 'GGCGCC';
$site{'PluTI'} = 'GGCGCC';
$site{'NarI'} = 'GGCGCC';
$site{'AscI'} = 'GGCGCGCC';
$site{'EciI'} = 'GGCGGA';
$site{'BsmFI'} = 'GGGAC';
$site{'ApaI'} = 'GGGCCC';
$site{'PspOMI'} = 'GGGCCC';
$site{'Sau96I'} = 'GGNCC';
$site{'NlaIV'} = 'GGNNCC';
$site{'KpnI'} = 'GGTACC';
$site{'KpnI-HF'} = 'GGTACC';
$site{'Acc65I'} = 'GGTACC';
$site{'BsaI-HF'} = 'GGTCTC';
$site{'BsaI'} = 'GGTCTC';
$site{'BsaI-HFv2'} = 'GGTCTC';
$site{'HphI'} = 'GGTGA';
$site{'BstEII-HF'} = 'GGTNACC';
$site{'BstEII'} = 'GGTNACC';
$site{'AvaII'} = 'GGWCC';
$site{'BanI'} = 'GGYRCC';
$site{'BaeGI'} = 'GKGCMC';
$site{'BsaHI'} = 'GRCGYC';
$site{'BanII'} = 'GRGCYC';
$site{'RsaI'} = 'GTAC';
$site{'CviQI'} = 'GTAC';
$site{'BstZ17I-HF'} = 'GTATAC';
$site{'BciVI'} = 'GTATCC';
$site{'SalI'} = 'GTCGAC';
$site{'SalI-HF'} = 'GTCGAC';
$site{'Nt.BsmAI'} = 'GTCTC';
$site{'BcoDI BsmAI'} = 'GTCTC';
$site{'ApaLI'} = 'GTGCAC';
$site{'BsgI'} = 'GTGCAG';
$site{'AccI'} = 'GTMKAC';
$site{'Hpy166II'} = 'GTNNAC';
$site{'Tsp45I'} = 'GTSAC';
$site{'HpaI'} = 'GTTAAC';
$site{'PmeI'} = 'GTTTAAAC';
$site{'HincII'} = 'GTYRAC';
$site{'BsiHKAI'} = 'GWGCWC';
$site{'TspRI'} = 'NNCASTGNN';
$site{'ApoI'} = 'RAATTY';
$site{'ApoI-HF'} = 'RAATTY';
$site{'NspI'} = 'RCATGY';
$site{'BsrFαI'} = 'RCCGGY';
$site{'BstYI'} = 'RGATCY';
$site{'HaeII'} = 'RGCGCY';
$site{'CviKI-1'} = 'RGCY';
$site{'EcoO109I'} = 'RGGNCCY';
$site{'PpuMI'} = 'RGGWCCY';
$site{'I-CeuI'} = 'TAACTATAACGGTCCTAAGGTAGCGAA';
$site{'SnaBI'} = 'TACGTA';
$site{'I-SceI'} = 'TAGGGATAACAGGGTAAT';
$site{'BspHI'} = 'TCATGA';
$site{'BspEI'} = 'TCCGGA';
$site{'MmeI'} = 'TCCRAC';
$site{'TaqαI'} = 'TCGA';
$site{'NruI'} = 'TCGCGA';
$site{'NruI-HF'} = 'TCGCGA';
$site{'Hpy188I'} = 'TCNGA';
$site{'Hpy188III'} = 'TCNNGA';
$site{'XbaI'} = 'TCTAGA';
$site{'BclI'} = 'TGATCA';
$site{'BclI-HF'} = 'TGATCA';
$site{'HpyCH4V'} = 'TGCA';
$site{'FspI'} = 'TGCGCA';
$site{'PI-PspI'} = 'TGGCAAACAGCTATTATGGGTATTATGGGT';
$site{'MscI'} = 'TGGCCA';
$site{'BsrGI'} = 'TGTACA';
$site{'BsrGI-HF'} = 'TGTACA';
$site{'MseI'} = 'TTAA';
$site{'PacI'} = 'TTAATTAA';
$site{'PsiI'} = 'TTATAA';
$site{'BstBI'} = 'TTCGAA';
$site{'DraI'} = 'TTTAAA';
$site{'PspXI'} = 'VCTCGAGB';
$site{'BsaWI'} = 'WCCGGW';
$site{'BsaAI'} = 'YACGTR';
$site{'EaeI'} = 'YGGCCR';

$deg{"A"} = "A";
$deg{"T"} = "T";
$deg{"C"} = "C";
$deg{"G"} = "G";
$deg{"R"} = "[AG]";
$deg{"R"} = "[AG]";
$deg{"Y"} = "[CT]";
$deg{"S"} = "[GC]";
$deg{"W"} = "[AT]";
$deg{"K"} = "[GT]";
$deg{"M"} = "[AC]";
$deg{"B"} = "[CGT]";
$deg{"D"} = "[AGT]";
$deg{"H"} = "[ACT]";
$deg{"V"} = "[ACG]";
$deg{"N"} = "[ATGC]";

}
