#!/usr/bin/perl
if (!$ARGV[0]) {
print STDERR "\nUSAGE: mgnify.pl prep|down|collect [MGnifile] folder\n\n";
print STDERR "       prep: reads the tabular file obtained form the MGnify \"Download results as CSV\" button.\n";
print STDERR "             and creates the mgnify.links file in the specified folder\n";
print STDERR "       down: reads the mgnify.links file form the specified folder and starts downloading\n";
print STDERR "             taxonomic abundance tables that are accumulated in the mgnify.out file in the specified folder\n";
print STDERR "    collect: reads the mgnify.out file in the specified folder and prepares the various mgnify.rank files\n";
print STDERR "       grep: reads the mgnify.out filtering matching lines\n";
print STDERR "     grepan: reads the mgnify.out file and displays tax tables of specific analyses\n\n";
exit;
}

$mode = $ARGV[0];

if ($mode eq "prep") { #usage: MGnify.pl prep infile folder
	open(IN,$ARGV[1]); #the input file from MGnify 
	#the input file derives form the default MGnify "Download results as CSV" button pressed on the alyses tab. 
	#it forms lines that looks like this: "MGYA00087407","3.0","ERS1462490","MGYS00001338","amplicon",,"ERR1744699",
	mkdir $ARGV[2] if (!-d $ARGV[2]); #a destination folder
	open(OUT,">$ARGV[2]/mgnify.links");
	#the output should be like this: MGYA00067182	MGYS00001166	https://www.ebi.ac.uk/metagenomics/api/v1/analyses/MGYA00067182/file/SRR2001789_MERGED_FASTQ_otu.tsv
	$cnt = 0;
	while ($line = <IN>) {
		$cnt++;
		next if ($cnt == 1);
		chomp $line;
		$line =~ s/\"//g;
		($analysis,$version,$sample,$study,$exptype,undef,$run) = split (/,/,$line);
		print OUT $analysis,"\t",$study,"\t","https://www.ebi.ac.uk/metagenomics/api/v1/analyses/",$analysis,"/file/",$run,"_MERGED_FASTQ_\n";
	}
	close IN;
	close OUT;
}

if ($mode eq "down") { #usage: Mgnify.pl down folder
	open(IN,"$ARGV[1]/mgnify.out.log");
	while ($line = <IN>) {
		chomp $line;
		@parts = split (/\t/,$line);
		$used{$parts[0]} = 1 if ($parts[$#parts] ne "no");
	}
	close IN;
	open(IN,"$ARGV[1]/mgnify.links");
	open(OUT,">>$ARGV[1]/mgnify.out");
	open(LOG,">>$ARGV[1]/mgnify.out.log");
	@case = qw /_MERGED_FASTA_SSU_OTU_TABLE_JSON.biom _MERGED_FASTQ_SSU_OTU_TABLE_JSON.biom _FASTA_SSU_OTU_TABLE_JSON.biom _FASTQ_SSU_OTU_TABLE_JSON.biom _MERGED_FASTQ_otu.tsv _FASTQ_otu.tsv _MERGED_FASTA_otu.tsv _FASTQ_otu.tsv/;
	while ($line = <IN>) {
		next if ($line =~ /^#/);
		chomp $line;
#		$mode = 1;
		#a: analysis, s: study, l: link (partial)
		($a,$s,$l) = split (/\t/,$line);
		$tmp = '';
		$run = $l;
		$run =~ s/.+\///;
#		print STDERR "$s @ $l:";
		print STDERR "$a @ $s: $run ";
		print STDERR " already done.\n" if ($used{$a});
		next if ($used{$a});
		$end = 'no';
		foreach $mode (0..$#case) {
			$link = $l.$case[$mode];
			#print STDERR "\n TESTING $link...";
			$tmp = `wget -q -O- $link`;
			#print STDERR "$mode." if $tmp =~ /\w/ && $tmp !~ /DOCTYPE html/;
			if ($tmp =~ /\w/ && $tmp !~ /DOCTYPE html/) {
				$end = $case[$mode];
				last;
			}
#			<STDIN>;
		}
		#a parser for the biom file
		if ($tmp =~ /\w/ && $link =~ /JSON/) {
			@tax = ();%tax = ();
			$tmp =~ /\"data\": \[\[(.+?)\]\]/;
			$data = ",[$1";
			while ($tmp =~ /{\"id\": \"(\d+)\", \"metadata\": \{\"taxonomy\": \[(.+?)\]/g) {
				push(@tax,$1);
				$tax{$1} = $2;
			}
			$cnt = 0;
			$tmp = '';
			while($data =~ /\[\d+,\d+,(.+?)\]/g) {
				$abund[$cnt] = $1;
				$tmp .= $tax[$cnt]."\t".$abund[$cnt]."\t".$tax{$tax[$cnt]}."\n";
				$cnt++;
			}
		}
		print STDERR "$end\n";		
		print LOG "$a\t$s\t$run\t$end\n";
		print OUT $a,"\t",$s,"\n",$tmp,"\n";
#		<STDIN>;
#		sleep 2;
		print STDERR "Compressing mgnify.out...";
		`"gzip $ARGV[1]/mgnify.out`;	
		print STDERR "done.\n";		
	}
	close IN;
	close OUT;
	close LOG;
}
if ($mode eq "collect") {
	open(IN,"zcat $ARGV[1]/mgnify.out.gz|");
	$cnt = 0;
	while ($line = <IN>) {
		chomp $line;
		next if ($line =~ /^#/);
		if ($line =~ /^MGYA/) {
			print STDERR $cnt{$exp},"\n" if (%exp); 
			$exp = $line;
			$exp{$line} = 1;
			push @exp, $line;
			$cnt++;
			print STDERR "$cnt. Parsing $line: ";
		} else {
			$line =~ /k__(.*?); p__(.*?); c__(.*?); o__(.*?); f__(.*?); g__(.*?); s__(.*?)/;
#			print $6,"\n";
			$k{$1} = 1;$p{$2} = 1;$c{$3} = 1;$o{$4} = 1;$f{$5} = 1;$g{$6} = 1;$s{$6." ".$7} = 1;
			$kx{$1."@".$exp} = 1;$px{$2."@".$exp} = 1;$cx{$3."@".$exp} = 1;$ox{$4."@".$exp} = 1;$fx{$5."@".$exp} = 1;$gx{$6."@".$exp} = 1;$sx{$6." ".$7."@".$exp} = 1;
			$cnt{$exp}++;
		}
#		last if $cnt > 100;
	}
	close IN;
	print STDERR scalar keys %cnt," / ",scalar keys %exp," experiments with info\n";
	print STDERR "Kingdom\t",scalar keys %k,"\n";
	print STDERR "Phylum\t",scalar keys %p,"\n";
	print STDERR "Class\t",scalar keys %c,"\n";
	print STDERR "Order\t",scalar keys %o,"\n";
	print STDERR "Family\t",scalar keys %f,"\n";
	print STDERR "Genus\t",scalar keys %g,"\n";
	print STDERR "Species\t",scalar keys %s,"\n";

  	open(OUT,">$ARGV[1]/mgnify.kingdom");
	print OUT "Name\t",(join "\t",sort keys %exp),"\t#exp\n";
	foreach $v (sort keys %k) {
		print OUT $v;
		$tot = 0;
		foreach $e (sort keys %exp) {
			print OUT "\t1" if ($kx{$v."@".$e});	
			$tot++ if ($kx{$v."@".$e});	
			print OUT "\t0" if (!$kx{$v."@".$e});	
		}
		print OUT "\t$tot\n";
	}
	close OUT;

  	open(OUT,">$ARGV[1]/mgnify.phyla");
	print OUT "Name\t",(join "\t",sort keys %exp),"\t#exp\n";
	foreach $v (sort keys %p) {
		print OUT $v;
		$tot = 0;
		foreach $e (sort keys %exp) {
			print OUT "\t1" if ($px{$v."@".$e});	
			$tot++ if ($px{$v."@".$e});	
			print OUT "\t0" if (!$px{$v."@".$e});	
		}
		print OUT "\t$tot\n";
	}
	close OUT;

 	open(OUT,">$ARGV[1]/mgnify.class");
	print OUT "Name\tnexp\t",(join "\t",sort keys %exp),"\t#exp\n";
	foreach $v (sort keys %c) {
		print OUT $v;
		$tot = 0;
		foreach $e (sort keys %exp) {
			print OUT "\t1" if ($cx{$v."@".$e});	
			$tot++ if ($cx{$v."@".$e});	
			print OUT "\t0" if (!$cx{$v."@".$e});	
		}
		print OUT "\t$tot\n";
	}
	close OUT;

 	open(OUT,">$ARGV[1]/mgnify.order");
	print OUT "Name\tnexp\t",(join "\t",sort keys %exp),"\t#exp\n";
	foreach $v (sort keys %o) {
		print OUT $v;
		$tot = 0;
		foreach $e (sort keys %exp) {
			print OUT "\t1" if ($ox{$v."@".$e});	
			$tot++ if ($ox{$v."@".$e});	
			print OUT "\t0" if (!$ox{$v."@".$e});	
		}
		print OUT "\t$tot\n";
	}
	close OUT;

 	open(OUT,">$ARGV[1]/mgnify.family");
	print OUT "Name\tnexp\t",(join "\t",sort keys %exp),"\t#exp\n";
	foreach $v (sort keys %f) {
		print OUT $v;
		$tot = 0;
		foreach $e (sort keys %exp) {
			print OUT "\t1" if ($fx{$v."@".$e});	
			$tot++ if ($fx{$v."@".$e});	
			print OUT "\t0" if (!$fx{$v."@".$e});	
		}
		print OUT "\t$tot\n";
	}
	close OUT;

 	open(OUT,">$ARGV[1]/mgnify.genus");
	print OUT "Name\t",(join "\t",sort keys %exp),"\t#exp\n";
	foreach $v (sort keys %g) {
		print OUT $v;
		$tot = 0;
		foreach $e (sort keys %exp) {
			print OUT "\t1" if ($gx{$v."@".$e});	
			$tot++ if ($gx{$v."@".$e});	
			print OUT "\t0" if (!$gx{$v."@".$e});	
		}
		print OUT "\t$tot\n";
	}
	close OUT;

 	open(OUT,">$ARGV[1]/mgnify.species");
	print OUT "Name\t",(join "\t",sort keys %exp),"\t#exp\n";
	foreach $v (sort keys %s) {
		print OUT $v;
		$tot = 0;
		foreach $e (sort keys %exp) {
			print OUT "\t1" if ($sx{$v."@".$e});	
			$tot++ if ($sx{$v."@".$e});	
			print OUT "\t0" if (!$sx{$v."@".$e});	
		}
		print OUT "\t$tot\n";
	}
	close OUT;
}

if ($mode eq "grep") { #this is a filter for the output of collect
	open(IN,"$ARGV[1]/mgnify.links");
	while ($line = <IN>) {
		next if ($line =~/^#/);
		chomp $line;
		($s,$p,undef) = split (/\t/,$line);
		$proj{$s} = $p;
	}
	close IN;
	open(IN,"zcat $ARGV[1]/mgnify.out.gz|");
	while ($line = <IN>) {
		chomp $line;
		next if ($line =~ /^#/);
		$current = $line if ($line =~ /^MGYA/);
		print $current,"\t",$proj{$current},"\t",$line,"\n" if ($line =~ /$ARGV[2]/i);
	}
	close IN;
}

if ($mode eq "grepan") { #this is a filter for the output of collect, specifically writes the tax table of a specific analysis
	open(IN,"zcat $ARGV[1]/mgnify.out.gz|");
	$print = 0;
	while ($line = <IN>) {
		next if ($line =~ /^#/);
		$print = 1 if ($line =~ /^MGYA/ && $line =~ /$ARGV[2]/i);
		$print = 0 if ($line =~ /^MGYA/ && $line !~ /$ARGV[2]/i);
		print $line if ($print);
	}
	close IN;
}

if ($mode eq "count") { #this is a filter for the output of collect
	open(IN,"$ARGV[1]/mgnify.links");
	while ($line = <IN>) {
		next if ($line =~/^#/);
		chomp $line;
		($s,$p,undef) = split (/\t/,$line);
		$proj{$s} = $p;
	}
	close IN;
	open(IN,"zcat $ARGV[1]/mgnify.out.gz|");
	while ($line = <IN>) {
		chomp $line;
		next if ($line =~ /^#/);
		next if ($line !~ /\w/);
		$current = $line if ($line =~ /^MGYA/);
		next if ($line =~ /^MGYA/);
		$cnt{$current}++;
	}
	close IN;
	foreach $s (sort {$cnt{$b} <=> $cnt{$a}} keys %cnt) {
		print $s,"\t",$proj{$s},"\t",$cnt{$s},"\n";
	}
}
