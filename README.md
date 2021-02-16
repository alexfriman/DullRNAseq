# DullRNAseq
The aligner was written with a solely purpose to understand how RNA sequencing is working. It does not account Q scores, it does not distinguish mapped and not-mapped reads when computing RPKM.
Both gene sequnces and reads are stored as strings.<br>
The author learned how NOT to write an aligner.<br>
Requires MPI compiler.<br>
Building sources:<br>
$ ./configure<br>
$ make<br>
Binary file DullRNAseq appears in ./<br>
Usage:<br>
$ mpirun -np /threads_number/ DullRNAseq genome.fasta reads.fastq [reads2.fastq ...]
