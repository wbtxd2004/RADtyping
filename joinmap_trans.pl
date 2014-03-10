##!/usr/bin/perl -w
my $r=0.8;
my $d_dir="genotype/poly_dom";
my $d_out="genotype/dom_JM";
my $c_dir="genotype/poly_codom";
my $c_out="genotype/codom_JM";
parse_command_line();
open FIG, "<$d_dir";
open OUT, ">$d_out";
foreach $word(<FIG>){
            chomp($word=$word);
            @file=split(/\s+/,$word);
            $len=@file;
            $str="@file[0] @file[1]";
            if(@file[2]==0){
                      $have=0;$none=0;
                      $str="$str <nnxnp>";
                      for($i=4;$i<$len;$i++){
                              if(@file[$i] eq "--"){$str="$str --";}    #change 5/9
                              elsif(@file[$i]==0){$str="$str nn";$none++;}
                              elsif(@file[$i]>0){$str="$str np";$have++;}
                      }       
            } 
            else{
                      $have=0;$none=0;
                      $str="$str <lmxll>";
                      for($i=4;$i<$len;$i++){
                              if(@file[$i] eq "--"){$str="$str --";}
                              elsif(@file[$i]==0){$str="$str ll";$none++;}
                              elsif(@file[$i]>0){$str="$str lm";$have++;}
                      }   
            }
            $sum=($none+$have)/2;
            $ka=(($none-$sum)**2+($have-$sum)**2)/$sum;
            print OUT  "0 $ka $str\n";
}  
close FIG;
close OUT;
open FIG, "<$c_dir";
open OUT, ">$c_out";
foreach $word(<FIG>){
            chomp($word=$word);
            @file=split(/\s+/,$word);
            $len=@file;
            $child_num=$len-6;
            $female=@file[4];$male=@file[5];
            $a=substr($female,0,1);$b=substr($female,1,1);
            $c=substr($male,0,1)  ;$d=substr($male,1,1)  ;
            $str=""; 
            if(($a ne $b) && ($c ne $d)&& ($female eq $male)){
                             $AA=0;$AT=0;$TT=0;$str="<hkxhk>";
                             for($i=6;$i<$len;$i++){
                                 if(@file[$i] eq $a."$a"){$AA++;$str="$str hh";}
                                 elsif(@file[$i] eq $a."$b"){$AT++;$str="$str hk";}
                                 elsif(@file[$i] eq $b."$b"){$TT++;$str="$str kk";}
                                 else{$str="$str --";}
                             }
            
                            $sum=($AA+$AT+$TT)/4;
                            if($sum>=$child_num*$r/4){
                               $ka=($AA-$sum)**2/$sum+($AT-$sum*2)**2/(2*$sum)+($TT-$sum)**2/$sum;
                              
                               if($ka<1000){
                                     print OUT  "1 $ka @file[0] @file[1] @file[2] @file[3] $str\n";
                               }
                            }
           }
           elsif(($a eq $b) && ($c ne $d)&&($a eq $c || $a eq $d)){
                              $AA=0;$AT=0;$TT=0;$str="<nnxnp>";
                              for($i=6;$i<$len;$i++){
                                 if(@file[$i] eq $female){$AA++;$str="$str nn";}
                                 elsif(@file[$i] eq $male)  {$AT++;$str="$str np";} 
                                 else{$str="$str --";}
                              }
                              $sum=($AA+$AT)/2;
                              if($sum>$child_num*$r/2){
                                 $ka=(($AA-$sum)**2+($AT-$sum)**2)/$sum;
                                 if($ka<=1000){ 
                                           print OUT  "0 $ka @file[0] @file[1] @file[2] @file[3] $str\n";
                                  }
                              }
          }
          elsif(($a ne $b) && ($c eq $d)&&($a eq $c || $b eq $c)){
                           $AA=0;$AT=0;$TT=0;$str="<lmxll>";
                           for($i=6;$i<$len;$i++){
                                 if(@file[$i] eq $female){$AT++;$str="$str lm";}
                                 elsif(@file[$i] eq $male)  {$AA++;$str="$str ll";}
                                 else{$str="$str --";}
                           }
                           $sum=($AA+$AT)/2;
                           if($sum>$child_num*$r/2){
                              $ka=(($AA-$sum)**2+($AT-$sum)**2)/$sum;
                              if($ka<=1000){ 
                                           print OUT  "0 $ka @file[0] @file[1] @file[2] @file[3] $str\n";
                              }
                           }
         }
         else{ print "$word\n";} 
   
        
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if    ($_ =~ /^-c$/)  { $c_dir  = shift  @ARGV; }
        elsif ($_ =~ /^-d$/)  { $d_dir  = shift  @ARGV; }
        elsif ($_ =~ /^-c1$/)  { $c_out  = shift  @ARGV; }
        elsif ($_ =~ /^-d1$/)  { $d_out  = shift  @ARGV; }
        else  {
	       print STDERR "Unknown command line option: '$_'\n";
	       usage();
	}
    }
}

sub usage {
    print STDERR <<EOQ; 
    perl joinmap_trans.pl -c file path  -d file path -o1 file path -o2 file path  -r  [-h]
    c     :input file path of codominant [genotype/poly_codom].
    d     :input file path of dominant   [genotype/poly_dom].
    c1    :output file path of codominant[genotype/codom_JM].
    d1    :output file path of dominant  [genotype/dom_JM].
    h     :display the help information.
EOQ
exit(0);
}

