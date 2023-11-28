#!/usr/bin/env perl
# Author: M.R. Weigand [yrh8@cdc.gov]

use strict;
use Getopt::Long;
use File::Basename;

&GetOptions(
	'dir=s' => \my$dir,	
	'samplesheet=s' => \my$samplesheet,
	'out=s' => \my$out);
($dir and $samplesheet and $out) or &HELP_MESSAGE;

#read list of sample names
my%hash=();
open IN, "$samplesheet";
while(my$i = <IN>){
    chomp$i;
    unless( $i =~ /fastq_1/ ){
        chomp$i;
        my@si = split(",",$i);

        #initiate data structure
        $hash{$si[0]}{'totreads'}=0;
        $hash{$si[0]}{'hscount'}=0;
        $hash{$si[0]}{'hsunmap'}=0;
        $hash{$si[0]}{'hsfract'}=0;
        $hash{$si[0]}{'bpcount'}=0;
        $hash{$si[0]}{'bpfract'}=0;
        $hash{$si[0]}{'uncount'}=0;
        $hash{$si[0]}{'unfract'}=0;
        $hash{$si[0]}{'depth'}=0;
        $hash{$si[0]}{'avgread'}=0;
        $hash{$si[0]}{'readlen'}=0;
        $hash{$si[0]}{'readset'}=0;
        @{$hash{$si[0]}{'breadth'}}=(0,0,0,0);
        @{$hash{$si[0]}{'mlst'}}=("NA","NA");
        @{$hash{$si[0]}{'asm'}}=("NA","NA","NA","NA","NA");

    }
}
close IN;

#now process results files
my@files = glob( $dir . '/*' );

foreach my$f (@files){

    my$id=basename($f);

    if($f =~ /hs-unmapped.flagstat/){ #flagstat results from *.hs-unmapped.bam
        $id =~ s/.hs-unmapped.flagstat//;
        my$hscount=0;
        $hscount = qx( grep -m 1 "primary" $f | awk '{print \$1}' );
        chomp$hscount;

        $hash{$id}{'hsunmap'}=$hscount;
       # print "$id\ths-uncount:$hscount\n";

    }elsif($f =~ /bp-mapped.flagstat/){ #flagstat results from *.bp-mapped.bam
        $id =~ s/.bp-mapped.flagstat//;
        my$bpcount=0;
        $bpcount = qx( grep -m 1 "primary" $f | awk '{print \$1}' );
        chomp$bpcount;

        $hash{$id}{'bpcount'}=$bpcount;

    }elsif($f =~ /bedcov.tsv/){ #mlst bedcov, only if allele refs provided with --mlst
        $id =~ s/.mlst-bedcov.tsv//;
        my@mlst=(0,0);

        open MLST, "$f";
        while(my$m = <MLST>){
            chomp$m;
            my@sm=split("\t",$m);
            if($sm[-1] >= 20){ $mlst[0]++ };
            if($sm[-1] >= 50){ $mlst[1]++ };
        }
        close MLST;
        @{$hash{$id}{'mlst'}}=@mlst;        
    
    }elsif($f =~ /depth.tsv/){ #samtools depth results from *.bp-mapped.bam
        $id =~ s/.bp-depth.tsv//;
        my@breadth=(0,0,0,0);

        open DEPTH, "$f";
        my$len=0;
        my$totcov=0;
        while(my$d = <DEPTH>){
            chomp$d;
            $len++;
            my@sd=split("\t",$d);
            $totcov+=$sd[-1];
            if($sd[-1] >= 20){ $breadth[0]++ };
            if($sd[-1] >= 50){ $breadth[2]++ };
        }
        close DEPTH;

        my$depth = sprintf("%.2f", $totcov/$len);
        $breadth[1] = sprintf("%.4f", $breadth[0]/$len);
        $breadth[3] = sprintf("%.4f", $breadth[2]/$len);

        $hash{$id}{'depth'}=$depth;
        @{$hash{$id}{'breadth'}}=@breadth;

    }elsif($f =~ /transposed_report.tsv/){ #combined QUAST assembly stats
        open QUAST, "$f";
        while(my$q = <QUAST>){
            chomp$q;
            unless($q =~ /^Assembly/){
                my@sq = split("\t",$q);
                $sq[0] =~ s/.skesa$//;

                @{$hash{$sq[0]}{'asm'}}=($sq[13],$sq[15],$sq[19],$sq[14],$sq[38]); #contigs,tot-len,N50,longest,ref-frac
            }
        }
        close QUAST;
    }elsif($f =~ /multiqc_general_stats.txt/){ #basic stats from multiqc
        open STATS, "$f";
        while(my$stat = <STATS>){
            chomp$stat;
            unless($stat =~ /^Sample/){
                my@ss = split("\t",$stat);
                $ss[0] =~ s/_[12]$//;
                $hash{$ss[0]}{'totreads'}+=$ss[-1];
                $hash{$ss[0]}{'readlen'}+=$ss[3];
                $hash{$ss[0]}{'readset'}++;

            }
        }
        close STATS;
    }

}
#now print output file
my@header = qw(Sample
    Reads-total
    Reads-Hs
    Frac-Hs 
    Reads-Bp
    Frac-Bp 
    Reads-Unmap
    Frac-Unmap
    Depth
    Breadth-20x-len
    Breadth-20x-frac
    Breadth-50x-len
    Breadth-50x-frac
    MLST-20x
    MLST-50x
    Read-len
    Asm-contigs
    Asm-length
    Asm-N50
    Asm-longest
    Asm-Ref-frac);

#print join(",",@header)."\n";
open OUT, ">$out";
print OUT join("\t",@header)."\n";
foreach my$sample (sort keys %hash){
    
    $hash{$sample}{'hscount'}=$hash{$sample}{'totreads'}-$hash{$sample}{'hsunmap'};
    $hash{$sample}{'uncount'}=$hash{$sample}{'hsunmap'}-$hash{$sample}{'bpcount'};
    $hash{$sample}{'hsfract'}=sprintf("%.4f",$hash{$sample}{'hscount'}/$hash{$sample}{'totreads'});
    $hash{$sample}{'bpfract'}=sprintf("%.4f",$hash{$sample}{'bpcount'}/$hash{$sample}{'totreads'});
    $hash{$sample}{'unfract'}=sprintf("%.4f",$hash{$sample}{'uncount'}/$hash{$sample}{'totreads'});
    $hash{$sample}{'avgread'}=sprintf("%.2f",$hash{$sample}{'readlen'}/$hash{$sample}{'readset'});
    
    print OUT join("\t", (
        $sample,
        $hash{$sample}{'totreads'},
        $hash{$sample}{'hscount'},
        $hash{$sample}{'hsfract'},
        $hash{$sample}{'bpcount'},
        $hash{$sample}{'bpfract'},
        $hash{$sample}{'uncount'},
        $hash{$sample}{'unfract'},
        $hash{$sample}{'depth'},
        @{$hash{$sample}{'breadth'}},
        @{$hash{$sample}{'mlst'}},
        $hash{$sample}{'avgread'},
        @{$hash{$sample}{'asm'}}) )."\n";

}

#########################
sub HELP_MESSAGE { die "
.Description:
   Combine and summarize CIWGS workflow outputs in tsv file.

.Usage: $0 -dir -out -samplesheet

   [mandatory]
   -dir	<dir>           Directory of collected results files.
   -samplesheet <file>  Workflow samplesheet.csv
   -out	<name>          Filename of output summary.
   
" }
