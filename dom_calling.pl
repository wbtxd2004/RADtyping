##!/usr/bin/perl -w

#####################################################################
#-this progrome is used to genotyping dominant and codomiant marker
#####################################################################
my $dep1=4;
my $dep2=4; 
my $p   =0.8;
opendir (DIR, "reads_mapping");
@dire = readdir DIR;
@dire = sort @dire;
my $out    ="genotype/all_dom";
my $out1   ="genotype/poly_dom";
parse_command_line();
do_ref();
Do_genotype();
filter();
system("rm genotype/TEMP");
sub do_ref(){
        open FIG,"<ref/HQ_ref_dom";
        open OUT,">$out";
        foreach $word(<FIG>){
               chomp($word=$word);
               print OUT "$word\n";
        }
        close FIG;
        close OUT;
}

sub filter(){
    for($C=1;$C<=200;$C++){
             $record[$C][0]=exp(-$C);
             $here=0;$tag=0;$bd[$C]=0;
             for($j=1;$j<=200;$j++){
                             $record[$C][$j]=$record[$C][$j-1]*$C/($j);
                             $here=$here+$record[$C][$j];
                             if($here<=0.05){$bd[$C]=$j;}   
             }
    }
    open FIG,"<$out";
    open OUT,">$out1";
    $do_m=0;
    foreach $word(<FIG>){
        chomp($word=$word);
        @file=split(/\s+/,$word);
        $none=0;$have=0;$sum=0;$num=0;
        $len=@file;
        for($i=4;$i<$len;$i++){
               if(@file[$i]>0){$num++;$sum=$sum+@file[$i];}
        } 
        if($num>0 && @file[2]*@file[3]==0){
           $hh=int($sum/$num);
           if($hh>$bd[$dep2/2] && $hh>=4){    #.........change 5/9
               $ave=int($sum/$num)+1;$str="@file[0] @file[1] @file[2] @file[3]";
               for($i=4;$i<$len;$i++){
                         if($bd[$ave]<-1){
                               if(@file[$i]>$bd[$ave]){$str=$str."  @file[$i]";}
                               else{$str=$str."  --";}
                         }
                         else{
                               if(@file[$i]>$bd[$ave]){$str=$str."  @file[$i]";}
                               elsif(@file[$i]>0){$str=$str."  --";}
                               else{$str=$str."  0";}
                         }
               }
               if($d1>=$pro_num*$p){print OUT "$str\n";} 
          } 
       }      
    }
}
  
sub Do_genotype(){
    for($kk=0;$kk<@dire;$kk++){
                   print "@dire[$kk]\n";
                   if(@dire[$kk] eq "."|| @dire[$kk] eq ".." ||  @dire[$kk] eq "P1" ||  @dire[$kk] eq "P2"){;}
                   else{
                           %hash={};
                           open FIG,"<soap/@dire[$kk]";
                           foreach $word(<FIG>){
                                   chomp($word=$word);
                                   @file=split(/\s+/,$word);
                                   $hash{@file[7]}++;
                           }
                           close FIG;
                           open FIG,"<$out";
                           open OUT,">genotype/TEMP";
                           foreach $word(<FIG>){
                                   chomp($word=$word);
                                   @file=split(/\s+/,$word);
                                   $str=substr(@file[0],1);
                                   $a=0;
                                   if(exists $hash{$str}){$a=$hash{$str};}
                                   print OUT "$word $a\n";
                           }
                          close FIG;
                          close OUT;
                          open FIG,"<genotype/TEMP";
                          open OUT,">$out";
                          foreach $word(<FIG>){
                                  chomp($word=$word);
                                  print OUT "$word\n";
                          }
                   }
                   close FIG;
                   close OUT;
    }
}

sub hashValueDescendingNum {$hello{$b} <=> $hello{$a};}


sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
        if    ($_ =~ /^-p$/)   { $p        = shift  @ARGV; }
        elsif    ($_ =~ /^-o1$/)  { $out     = shift  @ARGV; }
        elsif    ($_ =~ /^-o2$/)  { $out1     = shift  @ARGV; }
        else     {
                  print STDERR "Unknown command line option: '$_'\n";
	          usage();
        }
    }
}

sub usage {
    print STDERR <<EOQ; 
    perl dom_calling.pl -p  -o1  -o2  [-h]
    p  :least percentage of genotyped progenies[0.8]
    o1  :output file of all dominant genotypes[genotype/all_dom]
    o2  :output file of polymorphic dominant markers[genotype/poly_dom]
    h  :display the help information.
EOQ
exit(0);
}

