#!/usr/bin/perl -w


use strict;
use warnings;
no warnings 'uninitialized';	#### Added to original to run without error
use Getopt::Long qw(:config pass_through);
# &Getopt::Long::Configure qw(pass_through);
use IO::File;
use FindBin qw($Bin);

if (!@ARGV){
   print "Usage: perl $0 -i <the raw reads file>  -m <Barcodes mapping file>  -s <indexed genome name>  -d <length_disrupt_percent (max=1)> [-operon -c <operon_probability_cutoff (max=1)>] [-arrayed]\n";
   print "Required argument: ";
   print "-i gives the input raw reads file, -m gives the mapping file with the barcode sequence and name for each sample in the format as <barcode>\t<Sample name>, -s gives the name of the indexed genome that the reads should map to";
   print "Optional argument: \n";
   print "-d gives the region of the gene in which insertions are expected to disrupt gene function. The default is 1, which means when insertion falls anywhere in the gene (100%), the gene function will be disrupted. Setting the -d argument to 0.9, for example, would exclude insertions in the distal 10% of the gene when calculating the total number of reads/insertions for that gene. -operon (no argument) specifies that putative downstream (polar) insertions should be calculated based on a user-provided operon probability file. -c is the cutoff for operon probability, default is 1, which means only when the probability of two genes being in an operon is equal to 1 (100%), the polar effect will be considered for the downstream gene. Setting the -c argument to 0.8, for example, will calculate a polar effect for the downstream gene if the probability of the genes being in an operon is at least 0.8 (80%). -arrayed (no argument) is the option for the arrayed library. Refer to README for more information\n";
   exit;
}

my $datasource;   #### The raw reads file which will be analyzed
my $mismatch1=1;  #### Mismatches allowed in finding the transposon sequence
my $mismatch2=1;  #### Mismatches allowed when mapping the sequences to reference using bowtie
my $mapping;      #### Mapping file name
my $poolname="INSEQ_experiment";    #### The prefix of all the analyzed files
my $index;                          #### The name of the index fold in the indexes
my $operon=0;                       #### The option if analyze the sequence using operon information
my $cutoff=100;                     #### The cutoff for operon probability, default is 100
my $disruption=100;                 #### The disruptable percentage of genes. Only when transposon was inserted into this proximal region of gene, the gene is considered interrupted since the insertion in the distal region of a gene may not affect the function of gene at all. The default is 100.
my $arrayed=0;     #### The option if the data is from an arrayed library
my $path=$Bin;     #### Get the path for this analysis package


GetOptions (
  "i|input=s"=>\$datasource,
  "s|species=s"=>\$index,
  "m|mappingfile=s"=> \$mapping,
  "operon"=>\$operon,
   "d|disruption=s"=>\$disruption,
   "c|cutoff=i"=>\$cutoff,
   "arrayed"=>\$arrayed,
);

my $bowtie_path="";

open BT,"$path/config.txt" || die "Error, please check configuration file as config.txt in the package";
while (<BT>){
    chomp;
    if ($_=~m/\=\"(.*)\"/){
         $bowtie_path=$1;
         if ($1 eq ""){
            die "please specify where the bowtie directory is\n";
         }
    }
}


my $outputdir=`pwd`;
chomp $outputdir;
#print "$outputdir \n";


my $input=$poolname.".scarf";
if (-e $input){
 `rm $input`;
}

`ln -s $datasource $input`;

mkdir "bcsortedseqs";     
mkdir "results";

# Step 1
# Use the mapping file to assign each read to a barcode and store the output file as inputfile_assigned.txt", and store the statistics in the log file

my $barcode_assigned=$input."_assigned.txt";
my $logfile=$input.".log";

open CODES,$mapping || die "Error: can't open the mapping file, check the README for more details about mapping file\n";
my @codes;
my %codes_hash; #$codes_hash{$code} = sampleID
my %codes_number; #$codes_number{$code} = sample number (within a given sequencing lane)
my $sample_count=0;
my %has_barcode; #$has_barcode{$code} = number of reads with barcode $code 
my $total=0;
my $total_mapped=0;
my $code_length;

while (my $line=<CODES>){
  chomp $line;
  my @codes_array = split (/\s+/, $line);
  $code_length=length ($codes_array[0]);   
  $codes_hash{$codes_array[0]}=$codes_array[1];
  $codes_number{$codes_array[0]}=$sample_count;
  $sample_count++;
}
close CODES;

open OUT,">$barcode_assigned";
open INPUT,$input;
while (my $line=<INPUT>) {
  chomp $line;
  $total++;
  my @temp_array = split (/:/, $line);	####splits line into array of strings based on separator (:)
  my $seq_barcode = substr ($temp_array[9], 0, $code_length);   ####retreives barcode after 9th colon
  if (exists $codes_hash{$seq_barcode}){ #if first bases of read are one of the barcodes
    $total_mapped++;
    unless (exists $has_barcode{$seq_barcode}) {  ####what's this unless for???
      $has_barcode{$seq_barcode}=1;
    } else {
      $has_barcode{$seq_barcode}++;
    }
    my $line2 = <INPUT>;
    my $seq = substr ($line2,0,33);  			####retreives sequence
    print OUT ">".$codes_hash{$seq_barcode}.":".$seq_barcode."\n".$seq."\n";
  }
}
close OUT;

my %filehandlehash;
foreach my $a (keys %codes_hash) {
        if ($has_barcode{$a}){
           my $fh = IO::File->new(">bcsortedseqs\/$input\_$codes_hash{$a}\.fas");
           $filehandlehash{$a} = $fh;
        }
}

open LOG,">$logfile";
print LOG "Number of total reads in $input is $total\n";
print LOG "Total mapped $total_mapped\t".100*$total_mapped/$total." %\n";
print LOG "Barcode\tSample\tReads\tPercent\n";
foreach my $key (sort {$codes_number{$a} <=> $codes_number{$b}} keys %codes_number){
  if ($has_barcode{$key}){
    my $pct = 100*$has_barcode{$key}/$total;
    print LOG "$key\t$codes_hash{$key}\t$has_barcode{$key}\t$pct\n";
  }
}

#Step 2: Trimmed reads to remove transposon, append 16bp reads with a 5'N 

open IN,$barcode_assigned;

my @tn_array=split (//,"ACAGGTTG");
$total=0;
my $count=0;
my $exact_match_count=0;
my $mismatch1_count=0;
my $mismatch2_count=0;
my $header; 
my %tn_match; #$tn_match{sampleID} = number of reads that have TN sequence at allowable #mismatches
my @header; 
my $bc;
my $number;
my $length;

while (my $line =<IN>){ #go through each line
  chomp $line;
  if ($line =~m/^>/){ #if a header row
    $header = $line;
    @header = split (/:/,$line);
    $bc=$header[1];
    unless (exists $tn_match{$header[0]}){
      $tn_match{$header[0]}=0;
    }
  } else { #if a sequence row
    my $seq = $line;
    my $pos1_match=0;
    my $pos2_match=0;
    #does read contain transposon at bp 20 or 21?
    my $tn_pos1=substr($seq, 16, 8);
    my $tn_pos2=substr($seq, 17, 8);
    if ($tn_pos1 eq 'ACAGGTTG') {
      $pos1_match=1;
      $exact_match_count++;
      $count++;
    } else {
      if ($tn_pos2 eq 'ACAGGTTG') {
        $pos2_match=1;
        $exact_match_count++;
        $count++;
      }
    } 
    if ($mismatch1 == 1) { #if 1bp tn mismatches allowed
      my $pos1_score=0;
      my $pos2_score=0;
      my @pos1_array = split (//,$tn_pos1);
      my @pos2_array = split (//,$tn_pos2);
      for (my $n=0; $n<8;$n++){
        if ($pos1_array[$n] eq $tn_array[$n]) {
          $pos1_score++;
        }
        if ($pos2_array[$n] eq $tn_array[$n]) {
          $pos2_score++;
        }
      }
      if ($pos1_score == 7) {
        $pos1_match=1;
        $mismatch1_count++;
        $count++;
      } 
      if ($pos2_score == 7) {
        $pos2_match=1;
        $mismatch1_count++;
        $count++;
      }
    }
    if ($mismatch1 == 2) { #if 2bp tn mismatches allowed
      my $pos1_score=0;
      my $pos2_score=0;
      my @pos1_array = split (//,$tn_pos1);
      my @pos2_array = split (//,$tn_pos2);
      for (my $n=0; $n<8;$n++){
        if ($pos1_array[$n] eq $tn_array[$n]) {
          $pos1_score++;
        }
        if ($pos2_array[$n] eq $tn_array[$n]) {
          $pos2_score++;
        }
      }
      if ($pos1_score == 7) {
        $pos1_match=1;
        $mismatch1_count++;
        $count++;
      } 
      if ($pos2_score == 7) {
        $pos2_match=1;
        $mismatch1_count++;
        $count++;
      }
      if ($pos1_score == 6) {
        $pos1_match=1;
        $mismatch2_count++;
        $count++;
      } 
      if ($pos2_score == 6) {
        $pos2_match=1;
        $mismatch2_count++;
        $count++;
      }
    }
    
     #if match found
    if ($pos1_match==1){
      #trim transposon sequence  
      my $trimmed_seq = substr ($seq, 0, 16);
      $tn_match{$header[0]}++;
      $number->{$bc}->{$trimmed_seq}++;
      $length->{$bc}->{$trimmed_seq}=16;
    }
    if ($pos2_match==1){
      #trim transposon sequence
      my $trimmed_seq = substr ($seq, 0, 17);
      $tn_match{$header[0]}++;
      $number->{$bc}->{$trimmed_seq}++;
      $length->{$bc}->{$trimmed_seq}=17;
    }
  }
}
close OUT;

foreach my $sample (sort keys %$number){   
       my $fh=$filehandlehash{$sample};
       foreach my $read (sort keys %{$number->{$sample}}){
          $header=$sample.":".$length->{$sample}->{$read}.":".$number->{$sample}->{$read};
          print $fh ">$header\n";
          print $fh "$read\n";
      }
}

foreach (keys %codes_hash){
   close ($filehandlehash{$_});
}

print LOG "In the trimming process, Here is the statistics for each sample\n";
print LOG "Sample\tTrimmed\tPercentage\n";
foreach my $key (sort {$codes_number{$a} <=> $codes_number{$b}} keys %codes_number){
  if ($has_barcode{$key}){
    my $trimkey=">".$codes_hash{$key};
    my $percentage=sprintf("%.1f",$tn_match{$trimkey}/$has_barcode{$key}*100);
    print LOG "$key\t$tn_match{$trimkey}\t$percentage\n";
  }
}

close LOG;

chdir "$path/indexes/$index";
my @ptt=`ls *.ptt`;

chdir $outputdir;

my $clean_up="clean_up.sh";
open CL,">$clean_up";
print CL "rm $barcode_assigned\n",
         "rm mappingjobs_*\n",
         "rm -r bcsortedseqs\n";


#Create job files for each sample   
####Change to split jobs into two files, allowing MATLAB filtering in between? Or generate additional "_job2" for remapping after MATLAB filter?
foreach (keys %codes_hash) {
        if ($has_barcode{$_}){
          open (TASKSFORBC, ">$outputdir\/mappingjobs\_$input\_$_\.job") || die "Error: Can't create $outputdir\/mappingjobs\_$input\_$_\.job\n\n";
          print TASKSFORBC  "$bowtie_path/bowtie -m 1 --best --strata -a --fullref -n $mismatch2 -l 17  $path\/indexes\/$index\/$index -f bcsortedseqs\/$input\_$codes_hash{$_}\.fas results\/$input\_$codes_hash{$_}\.bowtiemap\n";
          print TASKSFORBC  "perl $path/process_bowtie_output.pl results\/$input\_$codes_hash{$_}\.bowtiemap \n";
          foreach my $ptt (@ptt){
               chomp $ptt;
               if ($ptt =~m/(.*)\.ptt/){
               print TASKSFORBC  "perl $path/normalize_processed_filter.pl results\/$input\_$codes_hash{$_}\.bowtiemap_processed.txt_$1\n";
               if (!$arrayed){

                  if ($operon==0){
                      print TASKSFORBC  "perl $path/map_genes_v2.pl $path/indexes\/$index\/$ptt results\/$input\_$codes_hash{$_}\.bowtiemap_processed.txt_$1_filter_cpm.txt $disruption\n";
                  }elsif ($operon==1){
                      print TASKSFORBC "perl $path/map_genes.pl $path/indexes\/$index\/$ptt $path/indexes\/$index\/$1.operons  $disruption $cutoff results\/$input\_$codes_hash{$_}\.bowtiemap_processed.txt_$1_filter_cpm.txt \n";
                  }
               }else{
               }  
              }
         }
            close TASKSFORBC;
       }
#         system("qsub -l h_vmem=4G mappingjobs\_$input\_$_\.job");
}


#### Generate additional "_job2" for remapping after MATLAB filtering  (PNM, Dec.2018)
	#### omit first half of job file; start with "normalize_process_filter

foreach (keys %codes_hash) {
        if ($has_barcode{$_}){
          open (TASKSFORBC, ">$outputdir\/mappingjobs\_$input\_$_\.job2") || die "Error: Can't create $outputdir\/mappingjobs\_$input\_$_\.job2\n\n";
          ####print TASKSFORBC  "$bowtie_path/bowtie -m 1 --best --strata -a --fullref -n $mismatch2 -l 17  $path\/indexes\/$index\/$index -f bcsortedseqs\/$input\_$codes_hash{$_}\.fas results\/$input\_$codes_hash{$_}\.bowtiemap\n";
          ####print TASKSFORBC  "perl $path/process_bowtie_output.pl results\/$input\_$codes_hash{$_}\.bowtiemap \n";
          foreach my $ptt (@ptt){
               chomp $ptt;
               if ($ptt =~m/(.*)\.ptt/){
               print TASKSFORBC  "perl $path/normalize_processed_filter.pl results\/$input\_$codes_hash{$_}\.bowtiemap_processed.txt_$1\n";
               if (!$arrayed){

                  if ($operon==0){
                      print TASKSFORBC  "perl $path/map_genes_v2.pl $path/indexes\/$index\/$ptt results\/$input\_$codes_hash{$_}\.bowtiemap_processed.txt_$1_filter_cpm.txt $disruption\n";
                  }elsif ($operon==1){
                      print TASKSFORBC "perl $path/map_genes.pl $path/indexes\/$index\/$ptt $path/indexes\/$index\/$1.operons  $disruption $cutoff results\/$input\_$codes_hash{$_}\.bowtiemap_processed.txt_$1_filter_cpm.txt \n";
                  }
               }else{
               }  
              }
         }
            close TASKSFORBC;
       }
#         system("qsub -l h_vmem=4G mappingjobs\_$input\_$_\.job2");  ####what's this, why is it hashed out? Also above for '.job'
}



