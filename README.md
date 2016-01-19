#MeGAMerge
MeGAMerge (A tool to merge assembled contigs, long reads from metagenomic sequencing runs)

##Description
MeGAMerge is a perl based wrapper/tool that can accept any number of sequence (FASTA) files containing assembled contigs of any length in Multi-FASTA format to produce an improved contig set based on OLC based assembly.  All overlap parameters (Minimum Overlap Length, Identity, etc) are user-declarable at runtime. It is written to run on Linux.

##Requirements:
You will need to have the following tools installed and in $PATH, or added to $binpath in the tool:

- Newbler (specifically runAssembly)
- Minimus2 (part of AMOS, also requires MUMmer)

##Usage:

MeGAMerge-1.0.pl [options] output_directory <list of fastas>

##Options:

-overlap=NN            Parameter for minimum overlap length in minimus2/Newbler (default = 80)

-minID=NN              Minimum % identity for overlap in minimus2/Newlber (default 98)

-conserr=NN            Maximum conservation error for minimus2 (default 0.06)

-cpu=NN                Number of CPU for Newbler (default 4)

-bindir=directory      Directory containing MUMmer executables and AMOS executables

-newblerdir=direcoty   Directory for newbler executable (runAssembly)

-o=outfile             Name of final file to output in output_directory (default MergedContigs.fasta)

-minLen=NN             Minimum length to include in newbler assemblies (default 150)

-minIncludeLen=NN      Minimum length to include in minimus assembly (default, 200)

-d                     Turns on debug information

-force                 Force overwrite of previous runs.  

-single_genome=1       Runs assuming single genome, reducing auto-options
                       (one newbler run, exclude fewer contigs, overrides -minLen and minIncludeLen)


##Citation
Please cite:

Scholz, M., Lo, C.-C., & Chain, P. S. G. (2014). Improved Assemblies Using a Source-Agnostic Pipeline for MetaGenomic Assembly by Merging (MeGAMerge) of Contigs. Scientific Reports, 4, 6480. Retrieved from http://dx.doi.org/10.1038/srep06480

if you use this software for your publications
