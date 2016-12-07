#MeGAMerge
MeGAMerge (A tool to merge assembled contigs, long reads from metagenomic sequencing runs)

##Description
MeGAMerge is a perl based wrapper/tool that can accept any number of sequence (FASTA) files containing assembled contigs of any length in Multi-FASTA format to produce an improved contig set based on OLC based assembly.  All overlap parameters (Minimum Overlap Length, Identity, etc) are user-declarable at runtime. It is written to run on Linux.

##Requirements:
You will need to have the following tools installed and in $PATH, or added to $binpath in the tool:

- [Newbler (specifically runAssembly)] (http://www.454.com/products/analysis-software/)
- [Minimus2 (part of AMOS, also requires MUMmer)] (http://amos.sourceforge.net/wiki/index.php/AMOS)

And a Perl module for Parallel nucmer script we developed for speeding up minimus2.
- [Parallel::ForkManager module from CPAN] (http://search.cpan.org/~yanick/Parallel-ForkManager-1.19/lib/Parallel/ForkManager.pm)

###Installation notes:
MUMmer:

For larger genome projects, the MUMmer package must be compiled in 64 bit mode.  This can be accomplished using:
make all CPPFLAGS="-O3 -DSIXTYFOURBITS"

AMOS:

For installation of AMOS, AMOS tools must be able to find nucmer, delta-filter and show-coords as compiled above, either by adding it  to the path before running ./configure
Or by specifying variables:
NUCMER,DELTAFILTER, and SHOWCOORDS when running ./configure in the amos directory. 
Example:
./configure NUCMER=/usr/local/bin/nucmer/bin/nucmer --prefix /usr/local/amos

Minimus2 (Required changes):

To take advantabe of the mutiple threads for the nucmer alignments step of Minimus2, we provide a custom minimus2 script, [minimus2_nucmer_multithreads] (https://github.com/LANL-Bioinformatics/MeGAMerge/blob/dev/minimus2_nucmer_multithreads), which needs modify the AMOS installation path before using it.
The following lines are path to above tools, users need to change the path to corresponding installed path. Usually, the tool path can be found by `which` command if they are in your environment $PATH variable. (ex: `which nucmer`)

    line 1: #!/PATH/to/AMOS-3.1.0/bin/runAmos -C

    line 46: BINDIR=/PATH/to/AMOS-3.1.0/bin
    line 47: NUCMER=/PATH/to/nucmer
    line 48: DELTAFILTER=/PATH/to/delta-filter
    line 49: SHOWCOORDS=/PATH/to/show-coords
 

##Usage:

MeGAMerge-1.2.pl [options] output_directory <list of fastas>

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
