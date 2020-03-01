# wastewater
scripts develoepd and used for searching and evaluating improved primers for 16S surveys in wastewater


<b>mgnify.pl</b>
parser for EBI MGnigy taxonomy tables. It reads input from the file produced by the "Download" button of the MGnify interface and produces a link file (e.g. mgnify.links) that ius further used to downoad, collect, and reorganize results.

<b>mgnify.links</b>
input for mgnify.pl (see above)

<b>analysis.sh</b>
bash script used to simulate PCRs results on RDP bacteria and archaea databases

<b>probematch_analyze.pl</b>
analysis and reorganization of the output of probeMatch

<b>degen_pcr.pl</b>
scipt used to simulate a PCR with degenerate primers. Tested on RDP sequences
