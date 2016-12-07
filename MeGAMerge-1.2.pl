#!/usr/bin/perl -w

#This program was prepared by Los Alamos National Security, LLC at Los Alamos National Laboratory (LANL)
#under contract No. DE-AC52-06NA25396 with the U.S. Department of Energy (DOE). All rights in the program
#are reserved by the DOE and Los Alamos National Security, LLC.  Permission is granted to the public to copy
#and use this software without charge, provided that this Notice and any statement of authorship are reproduced on all copies.
#Neither the U.S. Government nor LANS makes any warranty, express or implied, or assumes any liability or responsibility for the use of this software.

#License restrictions:
#Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
#•	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#•	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#•	Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
#THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


# Requires Newbler (runAssembly) and Amos (minimus2 and toAmos) to be in $bindir
# Algorithm:
#   sequences < 1800 bp and > $minLen (default 100 bp) run by Newbler
#   perform minimus2 on sequences >= 1800bp plus the above newbler assembly and singletons >= 200bp
#
#Author: mscholz@lanl.gov
#Modified by chienchi@lanl.gov  @ 20110502
#Modified by mscholz@msu.edu @20140210

use strict;
use File::Basename;
use Getopt::Long;
use Cwd;
use FindBin qw($RealBin);

umask 0002;

#$ENV{PATH} ="/usr/local/bin:/usr/bin:/bin:/gorgonzola2/bin:/opt/apps/AMOS/bin:/opt/apps/bin";
my $debug = 0;
my $force;

#Declare variables and default value
my ($dir, $iteration, $combine1, $combine2,
    $minimusHold, $singleGenome,
    $bnk, @remove, $cmd);
my $cpu=4;
my $overlap=80;
my $minID=98;
my $conserr=0.06;
my $minLen=150;
my $minInclude=200;
my $covCutoff = 3;
my $outfile="MergedContigs.fasta";
my $newblerdir = "/opt/apps/bin/";
my $AMOSdir = "/opt/apps/AMOS-short/bin";
my $bindir = "/opt/apps/AMOS/bin";
my ($minimusCount, $newblerCount, $smallCount, $filecount, $iter,$smallCutoff,$newblerCutoff) = (0,0,0,0,1,300,1800);

#get options
#defaults will work for our system, until we update directories for AMOS and Newbler.
GetOptions(
 "overlap=s"       => \$overlap,
 "minID=s"         => \$minID,
 "cpu=i"           => \$cpu,
 "conserr=s"       => \$conserr,
 "bindir=s"        => \$bindir,
 "newblerdir=s"    => \$newblerdir,
 "o=s"             => \$outfile,
 "minLen=s"        => \$minLen,
 "minIncludeLen=s" => \$minInclude,
 "d"               => \$debug,
 "single_genome=s" => \$singleGenome,
 "covCutoff=s"     => \$covCutoff,
 "force"				 => \$force,
 "help|?"          => sub{Usage()},
 "bnk"             => \$bnk
);

#my $minimus2_options = "-D CONSERR=$conserr -D MINID=$minID -D OVERLAP=$overlap";
my $minimus2 = "$RealBin/minimus2_nucmer_multithreds ";
my $minimus2_options = "-D CONSERR=$conserr -D MINID=$minID -D OVERLAP=$overlap -D THREADS=$cpu";

$ENV{PATH} .= ":$newblerdir" if (defined $newblerdir);
$ENV{PATH} .= ":$bindir" if (defined $bindir);

if (!@ARGV) { &Usage;}

#declare variables for location and minimum lengths for include and for consideration
if($minID > 100 || $minID < 0){print "problem with minID:$minID\n"; Usage()}
if($conserr > 1 || $conserr < 0){print "problem with conserr:$conserr\n";Usage();}
if($newblerdir){if($newblerdir !~ /\/$/){$newblerdir .= "/";} }
if($singleGenome){$minLen=0; $smallCutoff = 0}
#my $toAmos = "$bindir/toAmos";
#my $minimus2 = "$bindir/minimus2";
my $currentdir=getcwd;

#number of sequences to generate new file for minimus, or to use -large flag for newblerMax
my $newblerMax = 10000000;

#####Begin code using inputs#########
#first line in command must be directory
$dir = shift @ARGV;

if ( -e $dir && ! -d $dir){
    print "The $dir is an existed file, not a directory."; Usage();
}

$outfile= "$dir/". basename($outfile);

print qx/ps -o args $$/ if ($debug);
print "The Merged FASTA will be stored in $outfile\n" if ($debug);

unless (-e $dir){
	mkdir $dir;
}

#declare output files for fastas
my $minimusFile = "$dir/minimus.fasta";
my $largeFile   = "$dir/largefile.fasta";
my $NewblerFile = "$dir/newblerIn.fasta";

if($force){
	if (-e $minimusFile) {unlink $minimusFile};
	if (-e $largeFile)  {unlink $largeFile};
	if (-e $NewblerFile) {unlink $NewblerFile};
}
else{
	if (-e $minimusFile) {die "file $minimusFile exists in target directory, please remove, or use -force"};
	if (-e $largeFile)  {die "file $largeFile exists in target directory, please remove, or use -force"};
	if (-e $NewblerFile) {die "file $NewblerFile exists in target directory, please remove, or use -force"};
	if (-e $outfile){die "file $outfile exists in target directory, please remove, rename target output or use -force"};
}

push @remove, ($NewblerFile,$largeFile);

# all input fasta files
my @temp = @ARGV;

foreach(@temp){
	$filecount++;
	if ($debug){print "Reading $_\n";}
 	my $offset = &sortFiles($NewblerFile, $largeFile, $newblerCutoff, $_, $newblerCount);
 	$newblerCount += $offset;
}


#assemble them using newbler
print "Running Newbler assembly with $newblerCount sequences\n";
my $newblerOut = &newblerRun($NewblerFile,"$dir/newbler", $newblerCount);
push @remove, ("$dir/newbler");

#dynamically determine length of sequence to exclude from minimus2 combine
#this is necessary because of the large number of small contigs generated by SOAPdenovo
unless($minInclude){
#	$minInclude = &histogram($newblerOut);
	$minInclude = 200;
}

## run minimus2
my $seq_num_for_minimus2 = &openPrintFastas(0, $minimusFile, $largeFile, $newblerOut);
print "Running Minimus2 with $seq_num_for_minimus2 sequences\n";
&runMinimus($dir, $minimusFile);

#open the last file to print out to the renamed directory
my $final_contigs_num = &openPrintFastas($minInclude, $outfile, "$minimusFile.out");

print "\n Final Contigs number: $final_contigs_num \n";
if ($bnk){
	rename ("$dir/minimus.bnk", "$outfile.bnk");
#	my $cmd = "mv $dir/minimus.bnk $outfile.bnk";
#	system($cmd);
}

#remove old directories
&cleanup ("$dir/newbler",@remove);

&print_run_time;

########################## SUBS
sub Usage
{
	print <<END;
Usage: perl $0 [options] output_directory <list of fastas>

Options:
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
-single_genome=1       Runs assuming single genome, reducing auto-options
                       (one newbler run, exclude fewer contigs, overrides -minLen and minIncludeLen)
-force						Force $0 to run in non-empty directories (overwrite old run)

END

exit;
}

###
sub coverageTest{
	my ($test,$cutoff) = @_;
	if($test =~ /cvg/){
		my @temp =split /_/, $test;
		if($temp[1] >= $cutoff){return 0;}
		else {return 1;}
	}
	elsif($test =~ /numreads/){
		return 0;
	}
	else{
		return 0;
	}
}

###
sub sortFiles{
	my ($small, $large, $newblerCutoff, $filename,$offset) = @_;
	my @histogram;
	open NEWBLER, ">>$small" or die "trying to open $small $!";
	open MINIMUS, ">>$large" or die $!;
	open FILE, $filename or die "File $filename not found $!";
	my $count = 1;
	$/ = ">";
	my $iter1=0;
	while (<FILE>){
		my ($keys, $seq, $seq1) = ();
		my @seq = ();
		$_ =~ s/\>//g;
		unless ($_){next;}
		($keys, @seq) = split /\n/, $_;
		my $test = &coverageTest($keys,$covCutoff);
		if ($test > 0){next;}
		$keys =~ s/ /_/g;
		$seq = join "", @seq;
		undef @seq;
		my $len = length $seq;
		$count++;
		unless($len>= $minLen){next;}

		#if($seq =~ /N+/){
		#	@seq = split /N+/, $seq;
		#}
		#else{
			push @seq, $seq;
		#}

 		foreach $seq1 (@seq){
                   next if (! $seq1);
 		   my $len2 = length ($seq1);
 		   if( $len2 >= $newblerCutoff){
			  print MINIMUS ">Contig\_$minimusCount\n$seq1\n";
			  $histogram[$len2]++;
			  $minimusCount++;
		   }else{
	 		  print NEWBLER ">Newbler_contig$offset\n$seq1\n";
	 		  $offset++;
	 	   }

		}
	}
	close FILE;
	$/ = "\n";
	return $offset;
}


###
sub runMinimus{
	my ($dir,$next) = @_;
	die "problem finding $next" unless (-e $next);
	my $cmd = "toAmos -s $next -o $next.afg ";
	$cmd .= "1>/dev/null 2>/dev/null" unless ($debug);
	#chdir $dir;
	if(system ($cmd)){print "$cmd failed:\n$!"; die;}

	$cmd = "cd $currentdir; $minimus2 $minimus2_options $next ";
	$cmd .= "1>/dev/null 2>/dev/null" unless ($debug);
	print "$cmd\n";
	if(system ($cmd)){print "$cmd failed:\n$!"; die "$!";}

	my $offset = &openPrintFastas(0,"$next.out", "$next.fasta", "$next.singletons.seq");
	return $next;
}

###
sub delete_die{
	print "$cmd failed: $!";
	print "exiting nicely";
#  	&cleanup("$dir/newbler","$dir/small",@remove);
	exit;
}

sub newblerRun{
	my ($newblers, $outdir, $numberSeqs) = @_;
	#check to see if this point has been reached
	if(-e "$outdir/All.fasta" && ! $force){print "Newbler finished\n"; return "$outdir/All.fasta";}

	#declare local variables
	my $files;
	#my $large = "";
	#if ($numberSeqs > $newblerMax && !$singleGenome){
	#	print "$numberSeqs sequences, using -large flag\n\n";
	#	$large = "-large";
	#}

	my $cmd =  "runAssembly -force -large -rip -mi $minID -ml $overlap -pairt -cpu $cpu -a $minInclude -o $outdir $newblers";
	unless ($debug){ $cmd .= " 1>/dev/null 2>/dev/null";}
	print "$cmd\n";
	system($cmd);
	if (! -e "$outdir/454AllContigs.fna"){
		my $rerun = 1;
		until ($rerun >= 5){
			print "newbler failed $rerun times\n";
			system($cmd);
			if(! -e "$outdir/454AllContigs.fna"){$rerun ++;}
			else {$rerun = 10;}

		}
	}

	$files = "$outdir/454AllContigs.fna";
	#if (-e "$outdir/All.fasta"){unlink "$outdir/All.fasta";}
	open OUT, ">>$outdir/All.fasta";
	my @seq = &retrieveNewblerSingletons($outdir, $newblers);
	foreach (@seq){print OUT "$_";}
	open IN, "$files" or die $!;
	while (<IN>){print OUT "$_";}
	close IN;
	close OUT;

	return ("$outdir/All.fasta");
}


###
sub retrieveNewblerSingletons{
	my ($readdir,$filename) = @_;
	my ($sequence,$seqcount);
	my @return;
	my %getSeqs;

	open FILE, "$readdir/454ReadStatus.txt" or &delete_die($!) ;
	while (<FILE>){
		chomp $_;
		my @all = split /\t/, $_;
		if ($all[1] eq 'Singleton' || $all[1] eq 'Outlier'){
			unless($all[0] =~ /split/){
				$getSeqs{$all[0]} = 1;
			}
		}

	}
	close FILE;

	if($debug){print scalar(keys %getSeqs)." newbler singletons sequences\n";}
	if($debug){print "Checking 454PairAlign.txt\n";}
	open FILE, "$readdir/454PairAlign.txt" or  &delete_die($!);
	while (<FILE>){
	     chomp;
	     my @fields = split /\t/,$_;
	     next if ($_ !~ /\d/);
	     #$fields[10]=~ s/-//g;
	     #$fields[11]=~ s/-//g;
             # identical sequence
	     #if ($fields[3] == $fields[8] && $fields[7] == $fields[8])
	     #{
	     #    $getSeqs{$fields[4]} = 0;
	     #}
             # identical sub sequence
	     if ($fields[3] >= $fields[7] && $fields[7] == $fields[8])
	     {
	         $getSeqs{$fields[4]} = 0;
	     }
             # identical sub sequence
	     if ($fields[7] >= $fields[3] && $fields[3] == $fields[8])
	     {
	         $getSeqs{$fields[1]} = 0;
	     }
	}
	close FILE;

	my $key;
	open FILE, $filename or &delete_die($!) ;
	while (<FILE>){
		chomp $_;
		if($_ =~ /\>/){
			$_ =~ s/>//g;
			if ($sequence){
                if ($getSeqs{$key} && length($sequence) >= $minInclude ){
			       $seqcount++;
				   push @return, ">NewblerSingleton_$seqcount\n$sequence\n";
                }
				$sequence = "";
			}
			$key=$_;
		}
		else{
			$sequence .= $_;
		}
	}

	if($sequence && $getSeqs{$key} && length($sequence) >= $minInclude){
		$seqcount++;
		push @return, ">NewblerSingleton_$seqcount\n$sequence";
	}
	close FILE;

	if($debug){print "Unique Newbler Singletons (>=$minInclude bp): $seqcount\n";}
	return @return;
}

###

sub openPrintFastas{
	my $minIncluded = shift @_;
	my $file = shift @_;
	my @arry;
	my $return = 0;
	my $filenum = 0;
	my $count = 1;
	if($debug){print "opening $file for writing\n";}
	open OUT, ">$file" or die "can't open $file because: ".$! ;
	while(@_){
		$filenum++;
		my $filename = shift @_;
		if($debug){print "opening $filename for reading\n";}
		open FILE, $filename or next;
		$/ = ">";
		while (<FILE>){
			$_ =~ s/\>//g;
			unless ($_){next;}
			my ($name, @seq) = split /\n/, $_;
			if($seq[0]){
				my $seq = join "", @seq;
                                $seq =~ s/\?//g;
                                $seq = uc $seq;
				if(length $seq < $minIncluded){next;}
				print OUT ">Contig\_$return\n";
				$return++;
				undef @seq;
				if(length $seq > 100){$seq =~ s/(.{100})/$1\n/g;}
				print OUT "$seq";
			unless ($seq =~ /\n$/){print OUT "\n";}
			}
			else{next;}

		}
		$/ = "\n";
		close FILE;
	}
	close OUT;
	return $return;
}

sub cleanup{
	my $newblerdir = shift @_;
	foreach (@_){
		if(-e $_){unlink $_;}
	}
	
	$cmd = "rm -rf $newblerdir";
	system $cmd;
	$cmd = "rm -rf $dir/minimus*";
	system $cmd;
	return 0;
}

sub print_run_time {
  # Print runtime #
  my $run_time = time() - $^T;
  printf("\n Total running time: %02d:%02d:%02d\n\n", int($run_time / 3600), int(($run_time % 3600) / 60),
  int($run_time % 60));
}
