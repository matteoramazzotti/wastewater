#!/usr/bin/perl
&degen;
&er_site_list;

$mode = $ARGV[0];
$infile = $ARGV[1]; #a 16S collection in fasta format 
$fwd = $ARGV[2]; #written from 5' to 3'
$rev = $ARGV[3]; #written from 5' to 3'
$cut = $ARGV[4]; #optinal enzyme to simulate TRFLP

$p1 = '';
while ($fwd =~ /(.)/g) {
	$p1 .= $deg{$1};
}

$p2 = '';
$rev = reverse($rev);
$rev =~ tr/ATCGRTSWKMBDHVN/TAGCYRSWMKVHDBN/;
while ($rev =~ /(.)/g) {
	$p2 .= $deg{$1};
}

$p3 = '';
while ($cut =~ /(.)/g) {
	$p3 .= $deg{$1};
}

#print STDERR "$p1 -> $p2\n\n";
open(IN,$infile) or die "no infile $infile";
$cnt = 0;
$seq = '';
while ($line = <IN>) {
	chomp $line;
	if ($line =~ />/) {
		$cnt++;
		$name = $line;
		$name =~ s/\t.+//;
		$name =~ s/>//;
		&analyze if ($seq);
		print STDERR "\r$cnt sequences analyzed" if ($cnt%1000 == 0);
		$seq = '';
	} else {
		$seq .= uc($line);
	}
}
close IN;

sub analyze {
	$seqr = $seq;
	$seqr =~ tr/ATGC/TACG/;
	$seqr = reverse($seqr);
	$ok = 0;
	$ok = "F" if ($seq =~ /($p1.+?$p2)/);
	$ok = "R" if ($seqr =~ /($p1.+?$p2)/);
#	print STDERR "$seq\n$seqr\n$p1 -> $p2 -> $ok\n";
	if ($ok) {
		$amp_cnt++;
		$ampli = $1;
		if ($mode eq 'pcr'){
			print "$cnt\t$amp_cnt\t$name\t$ok\t",length($ampli),"\t$ampli\n";
		}
		if ($mode eq 'trflp'){
			if ($ampli =~ /$p3/g) {
				$cut_cnt++;
				$lf = $-[0];
				$rf = $+[0];
				$la = substr($ampli,0,$lf);
				$ra = substr($ampli,$rf);
				print STDERR "$cnt\t$amp_cnt\t",length($ampli),"\t$cut_cnt\t$p1\t$p2\t$p3\t$name\n$ampli\n$la\n$ra\n";
				$lsize_cnt{$lf}++;
				$rsize_cnt{$rf}++;
	 		}
		}
	}
}
if ($mode eq 'trflp'){
	foreach $s (keys %lsize) {
		print "L:\t",$s,"\t",$lsize_cnt{$s},"\n";
	}
	foreach $s (keys %lsize) {
		print "R:\t",$s,"\t",$lsize_cnt{$s},"\n";
	}
}

sub degen {
(%deg) = (
"A" => "A",
"T" => "T",
"C" => "C",
"G" => "G",
"R" => "[AG]",
"R" => "[AG]",
"Y" => "[CT]",
"S" => "[CG]",
"W" => "[AT]",
"K" => "[GT]",
"M" => "[AC]",
"B" => "[CGT]",
"D" => "[AGT]",
"H" => "[ACT]",
"V" => "[ACG]",
"N" => "[ACGT]"
);
}

sub er_site_list {
$site{'AclI'} = 'AACGTT';
$site{'HindIII HindIII-HF'} = 'AAGCTT';
$site{'SspI SspI-HF'} = 'AATATT';
$site{'MluCI'} = 'AATT';
$site{'PciI'} = 'ACATGT';
$site{'AgeI AgeI-HF'} = 'ACCGGT';
$site{'BspMI BfuAI'} = 'ACCTGC';
$site{'SexAI'} = 'ACCWGGT';
$site{'MluI MluI-HF'} = 'ACGCGT';
$site{'BceAI'} = 'ACGGC';
$site{'HpyCH4IV'} = 'ACGT';
$site{'HpyCH4III'} = 'ACNGT';
$site{'BaeI'} = 'ACNNNNGTAYC';
$site{'BsaXI'} = 'ACNNNNNCTCC';
$site{'AflIII'} = 'ACRYGT';
$site{'SpeI-HF SpeI'} = 'ACTAGT';
$site{'BsrI'} = 'ACTGG';
$site{'BmrI'} = 'ACTGGG';
$site{'BglII'} = 'AGATCT';
$site{'AfeI'} = 'AGCGCT';
$site{'AluI'} = 'AGCT';
$site{'StuI'} = 'AGGCCT';
$site{'ScaI-HF'} = 'AGTACT';
$site{'ClaI BspDI'} = 'ATCGAT';
$site{'PI-SceI'} = 'ATCTATGTCGGGTGCGGAGAAAGAGGTAAT';
$site{'NsiI NsiI-HF'} = 'ATGCAT';
$site{'AseI'} = 'ATTAAT';
$site{'SwaI'} = 'ATTTAAAT';
$site{'CspCI'} = 'CAANNNNNGTGG';
$site{'MfeI MfeI-HF'} = 'CAATTG';
$site{'Nb.BssSI'} = 'CACGAG';
$site{'BssSαI'} = 'CACGAG';
$site{'BmgBI'} = 'CACGTC';
$site{'PmlI'} = 'CACGTG';
$site{'DraIII-HF'} = 'CACNNNGTG';
$site{'AleI'} = 'CACNNNNGTG';
$site{'EcoP15I'} = 'CAGCAG';
$site{'PvuII PvuII-HF'} = 'CAGCTG';
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
$site{'NcoI NcoI-HF'} = 'CCATGG';
$site{'BseYI'} = 'CCCAGC';
$site{'FauI'} = 'CCCGC';
$site{'SmaI'} = 'CCCGGG';
$site{'XmaI TspMI'} = 'CCCGGG';
$site{'Nt.CviPII'} = 'CCD';
$site{'LpnPI'} = 'CCDG';
$site{'AciI'} = 'CCGC';
$site{'SacII'} = 'CCGCGG';
$site{'BsrBI'} = 'CCGCTC';
$site{'HpaII MspI'} = 'CCGG';
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
$site{'SbfI SbfI-HF'} = 'CCTGCAGG';
$site{'Bpu10I'} = 'CCTNAGC';
$site{'Bsu36I'} = 'CCTNAGG';
$site{'EcoNI'} = 'CCTNNNNNAGG';
$site{'HpyAV'} = 'CCTTC';
$site{'BstNI'} = 'CCWGG';
$site{'PspGI'} = 'CCWGG';
$site{'StyI StyI-HF'} = 'CCWWGG';
$site{'BcgI'} = 'CGANNNNNNTGC';
$site{'PvuI PvuI-HF'} = 'CGATCG';
$site{'BstUI'} = 'CGCG';
$site{'EagI-HF EagI'} = 'CGGCCG';
$site{'RsrII'} = 'CGGWCCG';
$site{'BsiEI'} = 'CGRYCG';
$site{'BsiWI BsiWI-HF'} = 'CGTACG';
$site{'BsmBI'} = 'CGTCTC';
$site{'Hpy99I'} = 'CGWCG';
$site{'MspA1I'} = 'CMGCKG';
$site{'MspJI'} = 'CNNR';
$site{'SgrAI'} = 'CRCCGGYG';
$site{'BfaI'} = 'CTAG';
$site{'BspCNI'} = 'CTCAG';
$site{'XhoI PaeR7I'} = 'CTCGAG';
$site{'EarI'} = 'CTCTTC';
$site{'AcuI'} = 'CTGAAG';
$site{'PstI PstI-HF'} = 'CTGCAG';
$site{'BpmI'} = 'CTGGAG';
$site{'DdeI'} = 'CTNAG';
$site{'SfcI'} = 'CTRYAG';
$site{'AflII'} = 'CTTAAG';
$site{'BpuEI'} = 'CTTGAG';
$site{'SmlI'} = 'CTYRAG';
$site{'AvaI BsoBI'} = 'CYCGRG';
$site{'MboII'} = 'GAAGA';
$site{'BbsI BbsI-HF'} = 'GAAGAC';
$site{'XmnI'} = 'GAANNNNTTC';
$site{'BsmI'} = 'GAATGC';
$site{'Nb.BsmI'} = 'GAATGC';
$site{'EcoRI EcoRI-HF'} = 'GAATTC';
$site{'HgaI'} = 'GACGC';
$site{'ZraI'} = 'GACGTC';
$site{'AatII'} = 'GACGTC';
$site{'Tth111I PflFI'} = 'GACNNNGTC';
$site{'PshAI'} = 'GACNNNNGTC';
$site{'AhdI'} = 'GACNNNNNGTC';
$site{'DrdI'} = 'GACNNNNNNGTC';
$site{'Eco53kI'} = 'GAGCTC';
$site{'SacI SacI-HF'} = 'GAGCTC';
$site{'BseRI'} = 'GAGGAG';
$site{'Nt.BstNBI'} = 'GAGTC';
$site{'PleI'} = 'GAGTC';
$site{'MlyI'} = 'GAGTC';
$site{'HinfI'} = 'GANTC';
$site{'EcoRV EcoRV-HF'} = 'GATATC';
$site{'Sau3AI MboI DpnII'} = 'GATC';
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
$site{'SphI SphI-HF'} = 'GCATGC';
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
$site{'NotI NotI-HF'} = 'GCGGCCGC';
$site{'Fnu4HI'} = 'GCNGC';
$site{'Cac8I'} = 'GCNNGC';
$site{'MwoI'} = 'GCNNNNNNNGC';
$site{'BmtI BmtI-HF'} = 'GCTAGC';
$site{'NheI NheI-HF'} = 'GCTAGC';
$site{'Nt.BspQI'} = 'GCTCTTC';
$site{'BspQI SapI'} = 'GCTCTTC';
$site{'BlpI'} = 'GCTNAGC';
$site{'TseI ApeKI'} = 'GCWGC';
$site{'Bsp1286I'} = 'GDGCHC';
$site{'AlwI'} = 'GGATC';
$site{'Nt.AlwI'} = 'GGATC';
$site{'BamHI-HF BamHI'} = 'GGATCC';
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
$site{'KpnI KpnI-HF'} = 'GGTACC';
$site{'Acc65I'} = 'GGTACC';
$site{'BsaI-HF BsaI BsaI-HFv2'} = 'GGTCTC';
$site{'HphI'} = 'GGTGA';
$site{'BstEII-HF BstEII'} = 'GGTNACC';
$site{'AvaII'} = 'GGWCC';
$site{'BanI'} = 'GGYRCC';
$site{'BaeGI'} = 'GKGCMC';
$site{'BsaHI'} = 'GRCGYC';
$site{'BanII'} = 'GRGCYC';
$site{'RsaI'} = 'GTAC';
$site{'CviQI'} = 'GTAC';
$site{'BstZ17I-HF'} = 'GTATAC';
$site{'BciVI'} = 'GTATCC';
$site{'SalI SalI-HF'} = 'GTCGAC';
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
$site{'ApoI ApoI-HF'} = 'RAATTY';
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
$site{'NruI NruI-HF'} = 'TCGCGA';
$site{'Hpy188I'} = 'TCNGA';
$site{'Hpy188III'} = 'TCNNGA';
$site{'XbaI'} = 'TCTAGA';
$site{'BclI BclI-HF'} = 'TGATCA';
$site{'HpyCH4V'} = 'TGCA';
$site{'FspI'} = 'TGCGCA';
$site{'PI-PspI'} = 'TGGCAAACAGCTATTATGGGTATTATGGGT';
$site{'MscI'} = 'TGGCCA';
$site{'BsrGI BsrGI-HF'} = 'TGTACA';
$site{'MseI'} = 'TTAA';
$site{'PacI'} = 'TTAATTAA';
$site{'PsiI'} = 'TTATAA';
$site{'BstBI'} = 'TTCGAA';
$site{'DraI'} = 'TTTAAA';
$site{'PspXI'} = 'VCTCGAGB';
$site{'BsaWI'} = 'WCCGGW';
$site{'BsaAI'} = 'YACGTR';
$site{'EaeI'} = 'YGGCCR';

}
