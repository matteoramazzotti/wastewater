probematch="/media/matteo/500gb/installers/RDPTools/ProbeMatch.jar"

#release 11.5
wget http://rdp.cme.msu.edu/download/current_Bacteria_unaligned.fa.gz
wget http://rdp.cme.msu.edu/download/current_Archaea_unaligned.fa.gz

zcat current_Archaea_unaligned.fa.gz && grep -c ">" current_Archaea_unaligned.fa
#160767
zcat current_Bacteria_unaligned.fa.gz && grep -c ">" current_Bacteria_unaligned.fa
#3196041

#to create revcom using iupacs (n.b. probematch needs both primers to be written as forward)
#xclip -o | perl -ne 'tr/ATGCRYKMBV/TACGYRMKVB/;print' | rev

#IMPORTANT: matches are intended as in the forward strand, so if e.g. a fish probe has tobe searched, its revcom match must be searched
#=> single probe check: revcom probe
#=> pcr primer check: dir ok, revcom of rev   

###trying primers from BMR (Takahashi 2014)
java -jar $probematch BMR/primers_BMR.probematch -n 0 -o BMR/archea.pair current_Archaea_unaligned.fa
perl probematch_analyze.pl ampli BMR/archea.pair current_Archaea_unaligned.fa BMR/archea.pair
java -jar $probematch BMR/primers_BMR.probematch -n 0 -o BMR/bacteria.pair current_Bacteria_unaligned.fa
perl probematch_analyze.pl ampli BMR/bacteria.pair current_Bacteria_unaligned.fa BMR/bacteria.pair

#hmmsearch -E 0.00001 --domtblout BMR/probematch_bacteria_BMR.pair.16Shmm --noali --cpu 2 -o /dev/null ../16S_bact_for3.hmm BMR/probematch_bacteria_BMR.pair.amplicons

###trying primers from BMR but with the fwd modified for increasing Brocadiales coverage
java -jar $probematch broca_opt/primers_broca_opt.probematch -n 0 -o broca_opt/archea.pair current_Archaea_unaligned.fa
perl probematch_analyze.pl ampli broca_opt/archea.pair current_Archaea_unaligned.fa broca_opt/archea.pair
java -jar $probematch broca_opt/primers_broca_opt.probematch -n 0 -o broca_opt/bacteria.pair current_Bacteria_unaligned.fa
perl probematch_analyze.pl ampli broca_opt/bacteria.pair current_Bacteria_unaligned.fa broca_opt/bacteria.pair

#hmmsearch -E 0.00001 --domtblout broca_opt/probematch_bacteria_broca_opt.pair.16Shmm --noali --cpu 2 -o /dev/null ../16S_bact_for3.hmm broca_opt/probematch_bacteria_broca_opt.pair.amplicons

###trying primers from Albertsen et al. 2015
java -jar $probematch albertsen_v1v3/albertsen_v1v3.probematch -n 0 -o albertsen_v1v3/archea.pair current_Archaea_unaligned.fa
perl probematch_analyze.pl ampli albertsen_v1v3/archea.pair current_Archaea_unaligned.fa albertsen_v1v3/archea.pair
java -jar $probematch albertsen_v1v3/albertsen_v1v3.probematch -n 0 -o albertsen_v1v3/bacteria.pair current_Bacteria_unaligned.fa
perl probematch_analyze.pl ampli albertsen_v1v3/bacteria.pair current_Bacteria_unaligned.fa albertsen_v1v3/bacteria.pair

java -jar $probematch albertsen_v3v4/albertsen_v3v4.probematch -n 0 -o albertsen_v3v4/archea.pair current_Archaea_unaligned.fa
perl probematch_analyze.pl ampli albertsen_v3v4/archea.pair current_Archaea_unaligned.fa albertsen_v3v4/archea.pair
java -jar $probematch albertsen_v3v4/albertsen_v3v4.probematch -n 0 -o albertsen_v3v4/bacteria.pair current_Bacteria_unaligned.fa
perl probematch_analyze.pl ampli albertsen_v3v4/bacteria.pair current_Bacteria_unaligned.fa albertsen_v3v4/bacteria.pair

java -jar $probematch albertsen_v4/albertsen_v4.probematch -n 0 -o albertsen_v4/archea.pair current_Archaea_unaligned.fa
perl probematch_analyze.pl ampli albertsen_v4/archea.pair current_Archaea_unaligned.fa albertsen_v4/archea.pair
java -jar $probematch albertsen_v4/albertsen_v4.probematch -n 0 -o albertsen_v4/bacteria.pair current_Bacteria_unaligned.fa
perl probematch_analyze.pl ampli albertsen_v4/bacteria.pair current_Bacteria_unaligned.fa albertsen_v4/bacteria.pair

