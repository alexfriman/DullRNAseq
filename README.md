# DullRNAseq
The aligner was written with a solely purpose to understand how RNA sequencing is working. It does not account Q scores, it does not distinguish mapped and not-mapped reads when computing RPKM.
Both gene sequnces and reads are stored as strings.
Requires MPI compiler.
Usage:
mpirun -np /threads_number/ DullRNAseq genome.fasta reads.fastq [reads2.fastq ...]
