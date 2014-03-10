##!/usr/bin/perl -w

###################################################################
#-this progrome is used to genotyping codomiant marker
###################################################################
###################################################################
my $len         ; 
my $pro_num=0   ; 
my $PP     =3.84;
my $p      =0.8 ;
mkdir("genotype");
my $out    ="genotype/all_codom";
my $out1   ="genotype/poly_codom";
my $out2   ="genotype/LOG";
parse_command_line();
if($a==0.05){$PP=3.84;}
else        {$PP=5.9;}
opendir (DIR, "reads_mapping");
@dire = readdir DIR;
@dire = sort @dire;
print "@dire\n";
here();
SNP_genotype("P1");
SNP_genotype("P2");
print $out2 "the number of progeny from left to right\n";
for($kk=0;$kk<@dire;$kk++){
    print  "@dire[$kk]\n";
    if(@dire[$kk] eq "P1"||@dire[$kk] eq "P2"|| @dire[$kk] eq "."|| @dire[$kk] eq ".."){;}
    else{SNP_genotype(@dire[$kk]);}
}
filter();
sub filter(){
    open FIG,"<$out";
    open OUT,">$out1";
    foreach $word(<FIG>){
            chomp($word=$word);
            @file=split(/\s+/,$word);
            $len=@file;
            $child_num=$len-6;
            $female=@file[4];$male=@file[5];
            $a=substr($female,0,1);$b=substr($female,1,1);
            $c=substr($male,0,1)  ;$d=substr($male,1,1)  ;
            if(($a ne $b) && ($c ne $d)&& ($female eq $male)){
                             $AA=0;$AT=0;$TT=0;
                             for($i=6;$i<$len;$i++){
                                 if(@file[$i] eq $a."$a"){$AA++;}
                                 elsif(@file[$i] eq $a."$b"){$AT++;}
                                 elsif(@file[$i] eq $b."$b"){$TT++;}
                             }
            
                            $sum=($AA+$AT+$TT)/4;
                            if($sum>=$child_num*$p/4 && ($AA>2 || $AT>2||$TT>2)){print OUT "$word\n";}
           }
           if(($a eq $b) && ($c ne $d)&&($a eq $c || $a eq $d)){
                              
                              $AA=0;$AT=0;$TT=0;
                              for($i=6;$i<$len;$i++){
                                 if(@file[$i] eq $female){$AA++;}
                                 elsif(@file[$i] eq $male)  {$AT++;} 
                              }
                              $sum=($AA+$AT)/2;
                              if($sum>$child_num*$p/2 && ($AA>2 || $AT>2)){print OUT "$word\n";}

          }
          if(($a ne $b) && ($c eq $d)&&($a eq $c || $b eq $c)){
                            $AA=0;$AT=0;$TT=0;
                             for($i=6;$i<$len;$i++){
                                 if(@file[$i] eq $female){$AT++;}
                                 elsif(@file[$i] eq $male)  {$AA++;}
                             }
                             $sum=($AA+$AT)/2;
                             if($sum>$child_num*$p/2 && ($AA>2 || $AT>2)){print OUT "$word\n";}

         } 
   }
}


sub here(){
    open FIG,"<ref/HQ_ref_codom";
    open OUT,">$out";
    $cnt=0;
    foreach $word(<FIG>){
            chomp($word=$word);
            $cnt++;
            @file=split(/\s+/,$word);
            $len=length @file[1];
            $ref[$cnt]=@file[1];
            $female{@file[2]}++;
            $male{@file[3]}++;
            @dep_1[$cnt]=@file[2];
            @dep_2[$cnt]=@file[3];
            for($i=1;$i<=$len;$i++){
                   $a=substr(@file[1],$i-1,1);
                   print OUT "@file[0] @file[1]  $i $a \n";
            }
    }
    close FIG;
}


    
#---------------------SNP-------------------------------------------
sub SNP_genotype(@_){
    my $name=@_[0];   my %depth={}; my %SNP={}; my $tol=0;$go=0; 
    open FIG,"<reads_mapping/$name";
    for(my $i=1;$i<=$len;$i++){
                       $dir[$i]->{A}=0;
                       $dir[$i]->{T}=0;
                       $dir[$i]->{G}=0;
                       $dir[$i]->{C}=0;
                       $dir[$i]->{N}=0;
    }
 
    foreach $word(<FIG>){
        chomp($word=$word);
        @file=split(/\s+/,$word);
        if($word=~/ref/){
                for($i=1;$i<=$len;$i++){
                       $A=$dir[$i]->{A};
                       $T=$dir[$i]->{T};
                       $G=$dir[$i]->{G};
                       $C=$dir[$i]->{C};
                       $n=$A+$T+$G+$C;
                       %hello=(A=>$A,T=>$T,G=>$G,C=>$C);
                       $tag=0;$genotype="--";$cnt=0;
                       foreach $key (sort hashValueDescendingNum (keys(%hello))) { @s[$cnt]=$key;$cnt++;}
                       $n1=$hello{@s[0]};$n2=$hello{@s[1]};$n3=$hello{@s[2]};$n4=$hello{@s[3]};
                       $genotype="--";
                       $t1=500;
                       if(($n1+$n2>=4 && $n1+$n2<$t1)){
                              $e=0.01;
                              $LRT1=(1-$e*3/4)**$n1*($e/4)**$n2/((0.5-$e/4)**($n1+$n2));
                              $LRT1=2*log(abs($LRT1)+0.000001);
                              if($LRT1<0 && -$LRT1>$PP){
                                          if(@s[0] le @s[1]){$genotype="@s[0]@s[1]";}
                                          else{$genotype="@s[1]@s[0]";}
                                          $tol++;
                              }
                              elsif($LRT1>$PP){$genotype="@s[0]@s[0]";}
                       }
                       $SNP{@file[0]}->[$i]=$genotype; 
                }
                for($i=1;$i<=$len;$i++){
                       $dir[$i]->{A}=0;
                       $dir[$i]->{T}=0;
                       $dir[$i]->{G}=0;
                       $dir[$i]->{C}=0;
                       $dir[$i]->{N}=0;
                }

        }
        else{
                for($i=1;$i<=$len;$i++){$dir[$i]->{substr(@file[0],$i-1,1)}++;}
        }
    } 
    my $cnt=0;
    open FIG,"<ref/HQ_ref_codom";
    foreach $word(<FIG>){
            chomp($word=$word);
            $cnt++;
            @file=split(/\s+/,$word);
            $ref[$cnt]=@file[1];
    }
    close FIG;
    open FIG,"<$out";
    open OUT,">genotype/TEMP";
    $cnt=0;$new=0;$het=0;
    foreach $word(<FIG>){
            chomp($word=$word);
            $a=int($cnt/$len);
            $str=$ref[$a+1];
            if(exists $SNP{$str}){
                    $SS=$SNP{$str}->[$cnt-$len*$a+1];
                    print  OUT"$word $SS\n";
                    if($SS ne "--"){$new++;}
                    if(substr($SS,0,1) ne substr($SS,1,1)){$het++;}
            }
            else{ print  OUT "$word --\n";}
            $cnt++;
    }
    close FIG;
    close OUT;
    open FIG,"<genotype/TEMP";
    open OUT,">$out";
    foreach $word(<FIG>){
            chomp($word=$word);
            print OUT "$word\n";
    }
    close FIG;
    close OUT;
}


sub hashValueDescendingNum {$hello{$b} <=> $hello{$a};}
 

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if       ($_ =~ /^-a$/)  { $a    = shift  @ARGV; }
        elsif    ($_ =~ /^-p$/)  { $p    = shift  @ARGV; }
        elsif    ($_ =~ /^-o1$/) { $out  = shift  @ARGV; }
        elsif    ($_ =~ /^-o2$/) { $out1 = shift  @ARGV; }
        else                     { 
                  print STDERR "Unknown command line option: '$_'\n";
	          usage();
        }
    }
}

sub usage {
    print STDERR <<EOQ; 
    perl codom_calling.pl -a -p -o1 -o2[-h]
    a  :significance level of LRT test               [0.05]
    p  :least percentage of genotyped progenies      [0.8]
    o1 :output file of all codominant genotypes      [genotype/all_codom]
    o2 :output file of polymorphic codominant markers[genotype/poly_codom]
    h :display the help information.
EOQ
exit(0);
}
