$mode = $ARGV[0];
if ($mode eq "grep") {
	open(IN,"mgnify.links");
	while ($line = <IN>) {
		next if ($line =~/^#/);
		chomp $line;
		($s,$p,undef) = split (/\t/,$line);
		$proj{$s} = $p;
	}
	close IN;
	open(IN,"mgnify.out");
	while ($line = <IN>) {
		chomp $line;
		next if ($line =~ /^#/);
		$current = $line if ($line =~ /^MGYA/);
		print $current,"\t",$proj{$current},"\t",$line,"\n" if ($line =~ /$ARGV[1]/i);
	}
	close IN;
}
if ($mode eq "collect") {
	open(IN,"mgnify.out");
	$cnt = 0;
	while ($line = <IN>) {
		chomp $line;
		next if ($line =~ /^#/);
		if ($line =~ /^MGYA/) {
			$exp = $line;
			$exp{$line} = 1;
			push @exp, $line;
			$cnt++;
			print STDERR "$cnt. Parsing $line\n";
		} else {
			$line =~ /k__(.*?); p__(.*?); c__(.*?); o__(.*?); f__(.*?); g__(.*?); s__(.*?)/;
			$k{$1} = 1;$p{$2} = 1;$c{$3} = 1;$o{$4} = 1;$f{$5} = 1;$g{$6} = 1;$s{$5." ".$6} = 1;
			$kx{$1."@".$exp} = 1;$px{$2."@".$exp} = 1;$cx{$3."@".$exp} = 1;$ox{$4."@".$exp} = 1;$fx{$5."@".$exp} = 1;$gx{$6."@".$exp} = 1;$sx{$5." ".$6."@".$exp} = 1;
		}
#		last if $cnt > 100;
	}
	close IN;
	print STDERR "Kingdom\t",scalar keys %k,"\n";
	print STDERR "Phylum\t",scalar keys %p,"\n";
	print STDERR "Class\t",scalar keys %c,"\n";
	print STDERR "Order\t",scalar keys %o,"\n";
	print STDERR "Family\t",scalar keys %f,"\n";
	print STDERR "Genus\t",scalar keys %g,"\n";
	print STDERR "Species\t",scalar keys %s,"\n";

  	open(OUT,">mgnify.kinkdom");
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

  	open(OUT,">mgnify.phyla");
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

 	open(OUT,">mgnify.class");
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

 	open(OUT,">mgnify.order");
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

 	open(OUT,">mgnify.family");
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

 	open(OUT,">mgnify.genus");
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

 	open(OUT,">mgnify.species");
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

if ($mode eq "down") {
	open(IN,"mgnify.out");
	while ($line = <IN>) {
		chomp $line;
		$used{$line}=1 if ($line =~ /^MGYA/);
	}
	close IN;
	open(IN,"mgnify.links");
	open(OUT,">>mgnify.out");
	while ($line = <IN>) {
		next if ($line =~ /^#/);
		chomp $line;
		$mode = 1;
		#s: sample, l:link,undef is the study, that is not of interest here...
		($s,undef,$l) = split (/\t/,$line);
		print STDERR "$s...";
		print STDERR "done.\n" if ($used{$s});
		next if ($used{$s});
		$tmp = `wget -q -O- $line` if ($mode == 1);
		print STDERR "1...";
		DOWN:
		$line =~ s/MERGED_// if ($mode == 2);
		$tmp = `wget -q -O- $line` if ($mode == 2);
		if ($tmp !~ /\w/ && $mode == 1) {
			$mode = 2;
			print STDERR "2...";
			goto DOWN;
		}
		print OUT $s,"\n",$tmp,"\n";
		print STDERR "\n";
		sleep 2;
	}
	close IN;
}
