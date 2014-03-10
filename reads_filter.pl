##!/usr/bin/perl -w
##########################################################
#-this progrome is used for high quality reads filtering-#
##########################################################
my $mode1="";
my $l=0;
my $f=2;
my $Q1=20;
my $Q2=10;
my $N=0.3; 
my $S=0.3;
parse_command_line();
open FIG,    "<$input";
open OUT,    ">$output";
open out,    ">$output1";
$cnt=0;
$t1="";$t2="";$t3="";$t4="";
$all=0;$high_qual=0;
for($i=1;$i<=$S*$l;$i++){$t1="$t1"."A";$t2="$t2"."T";$t3="$t3"."G";$t4="$t4"."C";} 		#$1由参数指定
foreach $word(<FIG>){
         chomp($word=$word);
         @file=split(/\s+/,$word);
         $cnt=$cnt+1;
         if($cnt%4==1){$name=substr($word,1);}
         elsif($cnt%4==2){$read=$word;}
         elsif($cnt%4==0){$qual=$word;
              $len=length $read;
              $qual_sum=0;
              $N_sum   =0;
              if($f==1){;}
              if($f==2){
                   $new_word="";
                   for($i=0;$i<$len;$i++){
                        $temp=substr($read,$len-1-$i,1); 
                        if(   $temp eq "A"){$new_word="$new_word"."T";}  
                        elsif($temp eq "T"){$new_word="$new_word"."A";}   
                        elsif($temp eq "G"){$new_word="$new_word"."C";}  
                        elsif($temp eq "C"){$new_word="$new_word"."G";} 
                        else {$new_word="$new_word"."$temp";} 
                   }
              }
              for($i=0;$i<$len;$i++){
                            $tt=substr($qual,0,$i);
                            if($tt eq "N"){$N_sum++;}
                            $tt=ord($tt)-33;  //ord 命令   将括号内的数转化成二进制
                            if($tt<$Q1){$qual_sum++;}               
              }
              if($read =~/$t1/||$read =~/$t2/||$read =~/$t3/||$read =~/$t4/|| $qual_sum>=$Q2||$N_sum>=$l*$N){;}
              else{ 
                         $tag=0;
                        
                         if($read=~/$base/){
                             while($read=~/$base/g ){ 
                                  $end=pos($read);
                                  $ss=substr($read,$end-$l,$l);
                                  $tag++;
                             }
                             print OUT ">$name\n$ss\n";
                         }
                         $read1=$new_word;
                         if($tag==0 && $f==2 && $read1=~/$base/){
                             while($read1=~/$base/g ){ 
                                  $end=pos($read1);
                                  $ss=substr($read1,$end-$l,$l);
                                  $tag++;
                             }
                             print OUT ">$name\n$ss\n";
                         }
                         if($tag==0){print out ">$name\n$ss\n";}
                                 
              }                     
         }
}
           
sub parse_command_line {  //下面参数的来源
    while (@ARGV) {
	$_ = shift @ARGV;
	if       ($_ =~ /^-i$/)  { $input   = shift  @ARGV; }
        elsif    ($_ =~ /^-o$/)  { $output  = shift  @ARGV; }
        elsif    ($_ =~ /^-o1$/)  { $output1  = shift  @ARGV; }
        elsif    ($_ =~ /^-b$/)  { $base    = shift  @ARGV; }
        elsif    ($_ =~ /^-l$/)  { $l       = shift  @ARGV; }
        elsif    ($_ =~ /^-f$/)  { $f       = shift  @ARGV; }
        elsif    ($_ =~ /^-Q1$/) { $Q1      = shift  @ARGV; }
        elsif    ($_ =~ /^-Q2$/) { $Q2      = shift  @ARGV; }
        elsif    ($_ =~ /^-S$/)  { $S       = shift  @ARGV; }
        else {
	       print STDERR "Unknown command line option: '$_'\n";
	       usage();
	}
        
    }
}

sub usage {
    print STDERR <<EOQ; 
    perl reads_filter.pl -i file path  -o file path -o1 trash -b -f  -T  -Q1  -Q2 -N -[-h]
    i    :input file 
    o    :outputfile
    o1   :trash file
    b    :target restriction site for BsaXI:[ATGC]{9}AC[ATGC]{5}CTCC[ATGC]{7}
    l    :length of read
    f    :f=1 palindromic structure;f=2 no palindromic structure [2] //palindromic structure 翻译
    Q1   :threshold for low quality score [20].
    Q2   :maximum no.of low-quality bases[10].  
    S    :discard reads with homopolymers >(Sxread length)[0.3].
    h    :display the help information.
EOQ
exit(0);
}
