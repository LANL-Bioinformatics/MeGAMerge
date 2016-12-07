#! /usr/bin/perl
# 20110824   chienchi at lanl dot gov

use Getopt::Long;
use Parallel::ForkManager;
use File::Basename;
use Cwd qw /abs_path/;
use strict;
use warnings;

use FindBin qw($Bin);

$ENV{PATH}="$Bin/../bin:$Bin:$ENV{PATH}";

#my $debug=1;
my ($prefix, $minimusOverlap, 
    $thread ,$NUCMER,
    $query, $ref,
    $breaklen);
my $nosimplify="";
### options default values ###
$thread=1;
$minimusOverlap=60;
$breaklen=250;
$prefix="Output";
#$NUCMER="/opt/apps/bin/nucmer";

#my $subset_seq_size=15000000;#15M
my $subset_seq_size=100000000; #100M
######################
GetOptions(   
            'overlap=i' => \$minimusOverlap,
            'thread=i'  => \$thread,
            'breaklen=i'=> \$breaklen,
            'nosimplify'=> \$nosimplify,
            'query=s'   => \$query,
            'ref=s'     => \$ref,
            'prefix=s'  => \$prefix,
	    'subset=f'  => \$subset_seq_size,
            'help|?' =>  sub{Usage()}
);
$subset_seq_size +=0;

&Usage if (!$query and !$ref);
my $max_thread=`grep -c ^processor /proc/cpuinfo`;
if ( $thread < 1 || $thread >$max_thread)  {
  die("-thread value must be between than 1 and $max_thread.\n");
}
my $self_nucmer=&check_self_nucmer($ref,$query);

my ($file_name, $file_path, $suffix)=fileparse("$prefix", qr/\.[^.]*/);
$ref=abs_path($ref);
$query=abs_path($query);

my $dir= "${file_path}minimusTemp$$";

if (-e "$dir") {system ("rm -rf $dir");}
mkdir ("$dir");

$nosimplify = "--nosimplify " if ($nosimplify);

my (@query_files,@ref_files);
@query_files=&splitFasta_by_size($dir,$query,"query",$subset_seq_size);
@ref_files=&splitFasta_by_size($dir,$ref,"ref",$subset_seq_size);  

  my $pm = new Parallel::ForkManager($thread);
  foreach my $ref_file_i(0..$#ref_files)
  { 
    foreach my $query_file_i(0..$#query_files)
    {
       next if ($ref_file_i > $query_file_i && $self_nucmer);
       $pm->start($ref_file_i) and next;
       my $cmd1="nucmer --maxmatch $nosimplify -b $breaklen -p $dir/$ref_file_i.$query_file_i -c $minimusOverlap $ref_files[$ref_file_i] $query_files[$query_file_i] ";
       if (system($cmd1)){print STDERR "error on $cmd1\n"}  
       $pm->finish;
    }
  }
    $pm->wait_all_children;

## concate delta files
open (OUT,">$prefix.delta") or die $!;
print OUT "$ref $query\nNUCMER\n";
foreach my $ref_file_i(0..$#ref_files)
{ 
   foreach my $query_file_i(0..$#query_files)
   {
        next if ($ref_file_i > $query_file_i && $self_nucmer);
        open (FILE, "< $dir/$ref_file_i.$query_file_i.delta") or die $!;
        my $tmp=<FILE>; 
        my $tmp2=<FILE>;
        while(<FILE>){print OUT $_;}
        close FILE;
        unlink "$dir/$ref_file_i.$query_file_i.delta";
   }
}
close OUT;
# clean up
system ("rm -rf $dir");
exit;

sub Usage
{
	print <<END;
Usage: $0 [options] -ref <fasta> -query <fasta>
This script is for minimus2 on large data. It implemented with --maxmatch, -c, -b, --nosimplify and -p of nucmer and
It splits the input by size 100M bp and run numcer on each splitted dataset. 
If your input is < 100M bp, it always run with 1 thread only.
Options:
-overlap        Parameter for minimum overlap length in minimus2 (default = 60)
-breaklen       the distance an alignment extension will attempt to
                extend poor scoring regions before giving up
-prefix         Name of final file to output in output_directory (default Output)
-thread         Threads usage (default 1)
-nosimplify     Turn this option off if aligning a sequence to itself to look for repeats
-subset		subset size. (default 100e6)
-help           Show this usage

END
exit;
} 

sub splitFasta_by_size
{
     my ($dir,$file,$type,$seq_size_in_file)=@_;
     my $file_count=0;
     my $seq_count=0;
     my @file_list;
     my $len=0;
     my $seq;
     my $head;
     open FILE, ">> $dir/$type.$file_count" or die $!;
     push @file_list,"$dir/$type.$file_count";
     open (IN,"$file");
     while(<IN>)
     {
           chomp;
           if (/>/)
           {
              if ($seq){
                 $seq_count++;
                 $len = $len + length ($seq);
                 if ($len < $seq_size_in_file)
                 {
                   print FILE "$head\n$seq\n"; 
                 }
                 else
                 {
                   $file_count++;
                   $len=0;
                   close FILE;
                   open FILE, ">>$dir/$type.$file_count" or die $!;
                   push @file_list, "$dir/$type.$file_count";
                   print FILE "$head\n$seq\n";
                 }
              }
              $head=$_;
              $seq="";
           }
           else
           {
              $_ =~ s/[\n\t\f\r_0-9\s]//g;
              $seq.=$_;
           }
     }
           if ($seq){
                 $seq_count++;
                 $len = $len + length ($seq);
                 if ($len < $seq_size_in_file)
                 {
                   print FILE "$head\n$seq\n"; 
                 }
                 else
                 {
                   $file_count++;
                   $len=0;
                   close FILE;
                   open FILE, ">>$dir/$type.$file_count" or die $!;
                   push @file_list,"$dir/$type.$file_count";
                   print FILE "$head\n$seq\n";
                 }
              }
     close FILE;         
     close IN;
     
    
     return (@file_list);
     
}

sub check_self_nucmer 
{
    my $ref=shift;
    my $query=shift;
    my $self_nucmer=1;
    open (REF,$ref) or die "$!\n";
    open (QUERY,$query) or die "$!\n";
    
    my $check_line_number = 1000;
    for my $i (1..10000)
    {
       my $ref_line=<REF>;
       my $query_line=<QUERY>;
       next if (!defined $ref_line);
       if ($ref_line ne $query_line)
       {
           $self_nucmer=0;
           last;
       }
    }
    
    return ($self_nucmer);
}
