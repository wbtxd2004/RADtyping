##!/usr/bin/perl -w

###################################################################
#-this progrome is used to construct high quality reference-------#
###################################################################
system("mkdir soap");
parse_command_line();
ref_filter();
sub ref_filter(){
    open FIG,  "<ref/ref_codom";
    $cnt=0;
    %hash={};
    foreach $word(<FIG>){
            chomp($word=$word);
            @file=split(/\s+/,$word);     #........change 5..8...................
            $cnt++;
            $len=length @file[1];
            $ref[$cnt]=@file[1];
            $female{@file[2]}++;
            $male{@file[3]}++;
            @dep_1[$cnt]=@file[2];
            @dep_2[$cnt]=@file[3]; 
            $hash{@file[2]}++;
            $hash{@file[3]}++;
    }
    close FIG;
    #--------Mixed normal/poisson mode to determine threshold of composite cluster--
    @dep_1[0]=$d1;
    @pare1=norm_fit_define(@dep_1);
    $upN_1=(2+@pare1[2]/@pare1[3])*@pare1[1];
    @dep_2[0]=$d2;
    @pare2=norm_fit_define(@dep_2);
    $upN_2=(2+@pare2[2]/@pare2[3])*@pare2[1];
    print "norm $upN_1  $upN_2\n";




    $mean1=0;$mean2=0;$std1=0;$std2=0;
    for($i=1;$i<=@dep_1;$i++){
            $mean1=$mean1+@dep_1[$i];
            $mean2=$mean2+@dep_2[$i];
    }
    $mean1=$mean1/@dep_1;
    $mean2=$mean2/@dep_2;
    for($i=1;$i<=@dep_1;$i++){
            $std1=$std1+(@dep_1[$i]-$mean1)**2;
            $std2=$std2+(@dep_2[$i]-$mean2)**2;
    }
    $std1=sqrt($std1/@dep_1);
    $std2=sqrt($std2/@dep_2);
    $upT_1=$mean1+2*$std1;
    $upT_2=$mean2+2*$std2;
    
    if($m eq "P" || $m eq "N"){$up_1=$upN_1;$up_2=$upN_2;}
    elsif($m eq "T"){$up_1=$upT_1;$up_2=$upT_2;}
    else{
                   @left=sort{$hash{$b}<=>$hash{$a}} keys %hash;
                   $up_1=@left[0];$up_2=@left[0];
    }
    close FIG; 
    
    open FIG,  "<ref/ref_codom";
    open OUT,  ">ref/HQ_ref_codom";
    open OUT1, ">soap/ref";
    $cnt=0;
    foreach $word(<FIG>){
            chomp($word=$word);
            $cnt++;
            @file=split(/\s+/,$word);
            $len=length @file[1];
            $all=@file[2]+@file[3];
            if($all>=8 && $all<=$up_1+$up_2){ 
                          print OUT "$word\n";
                          print OUT1"@file[0]\n";
                          print OUT1 "@file[1]"."AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
            }
    }
    close FIG;
    close OUT;
    close OUT1;
    open FIG,  "<ref/ref_dom";
    open OUT,  ">ref/HQ_ref_dom";
    open OUT1, ">>soap/ref";
    $cnt=0;
    foreach $word(<FIG>){
            chomp($word=$word);
            $cnt++;
            @file=split(/\s+/,$word);
            $len=length @file[1];
            $all=@file[2]+@file[3];
            if($all>=4){ 
                          print OUT "$word\n";
                          print OUT1"@file[0]\n";
                          print OUT1 "@file[1]"."AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
            }
    }
    close FIG;
    
}

sub norm_fit_define(@_){
    $constant=@_[0]; #########13/5/8
    @num=@_;
    print "@pare\n";
    $cnt=scalar @num;
    $a1=0.8;$u1=20;$sigma1=10;$a2=1-$a1;$u2=2*$u1;$sigma2=10; 
    for($time=1;$time<=1000;$time++){
            $sum1=0;$sum2=0;$sum3=0;
            for($i=1;$i<=$cnt;$i++){
                    @temp[0]=@num[$i];$temp[1]=$u1;$temp[2]=$sigma1;
                    if($sigma1<10**(-50)){$one=0;}
                    else{
                               $one=my_norm(@temp);
                    }

                    @temp[0]=@num[$i];$temp[1]=2*$u1;$temp[2]=$sigma2;
                    if($sigma2<10**(-50)){$two=0;}
                    else{
                               $two=my_norm(@temp);
                    }
                    $dd=$a1*$one+(1-$a1)*$two;

                    if($dd<10**(-50)){$w1[$i]=0;$w2[$i]=0;}
                    else{
                           $w1[$i]=$a1*$one/($a1*$one+(1-$a1)*$two);
                           $w2[$i]=1-$w1[$i];
                    }
                    
            }
            $sum1=0;$sum2=0;$sum3=0;$summ1=0;$summ2=0;$summ3=0;
            for($i=1;$i<=$cnt;$i++){
                    $sum1=$sum1+$w1[$i];
                    $sum2=$sum2+$w1[$i]*@num[$i];
                    $sum3=$sum3+$w1[$i]*(@num[$i]-$u1)**2;
                    $summ1=$summ1+$w2[$i];
                    $summ2=$summ2+$w2[$i]*@num[$i];
                    $summ3=$summ3+$w2[$i]*(@num[$i]-2*$u1)**2;
            }
            $a1_new=$sum1/$cnt;
            $u1_new=$sum2/$sum1;
            $sigma1_new=$sum3/$sum1;
                
            $a2_new=$summ1/$cnt;
            $u2_new=$summ2/$summ1;
            $sigma2_new=$summ3/$summ1;

            if(abs($a1-$a1_new)>10**(-3)){
                   $a1=$a1_new;$u1=$constant;$sigma1=$sigma1_new;  
                   $a2=$a2_new;$u2=$constant;$sigma2=$sigma2_new; 
            }
            else{
                   @pare[0]=$a1;@pare[1]=$constant;@pare[2]=sqrt($sigma1);@pare[3]=sqrt($sigma2);
                   return @pare;
            }
            
    }
}

sub norm_fit(@_){
    @num=@_;
    print "@pare\n";
    $cnt=scalar @num;
    $a1=0.8;$u1=20;$sigma1=10;$a2=1-$a1;$u2=2*$u1;$sigma2=10; 
    for($time=1;$time<=1000;$time++){
            $sum1=0;$sum2=0;$sum3=0;
            for($i=1;$i<=$cnt;$i++){
                    @temp[0]=@num[$i];$temp[1]=$u1;$temp[2]=$sigma1;
                    if($sigma1<10**(-50)){$one=0;}
                    else{
                               $one=my_norm(@temp);
                    }

                    @temp[0]=@num[$i];$temp[1]=2*$u1;$temp[2]=$sigma2;
                    if($sigma2<10**(-50)){$two=0;}
                    else{
                               $two=my_norm(@temp);
                    }
                    $dd=$a1*$one+(1-$a1)*$two;

                    if($dd<10**(-50)){$w1[$i]=0;$w2[$i]=0;}
                    else{
                           $w1[$i]=$a1*$one/($a1*$one+(1-$a1)*$two);
                           $w2[$i]=1-$w1[$i];
                    }
                    
            }
            $sum1=0;$sum2=0;$sum3=0;$summ1=0;$summ2=0;$summ3=0;
            for($i=1;$i<=$cnt;$i++){
                    $sum1=$sum1+$w1[$i];
                    $sum2=$sum2+$w1[$i]*@num[$i];
                    $sum3=$sum3+$w1[$i]*(@num[$i]-$u1)**2;
                    $summ1=$summ1+$w2[$i];
                    $summ2=$summ2+$w2[$i]*@num[$i];
                    $summ3=$summ3+$w2[$i]*(@num[$i]-2*$u1)**2;
            }
            $a1_new=$sum1/$cnt;
            $u1_new=$sum2/$sum1;
            $sigma1_new=$sum3/$sum1;
                
            $a2_new=$summ1/$cnt;
            $u2_new=$summ2/$summ1;
            $sigma2_new=$summ3/$summ1;

            if(abs($a1-$a1_new)>10**(-3)){
                   $a1=$a1_new;$u1=$u1_new;$sigma1=$sigma1_new;  
                   $a2=$a2_new;$u2=$u2_new;$sigma2=$sigma2_new; 
            }
            else{
                   @pare[0]=$a1;@pare[1]=$u1;@pare[2]=sqrt($sigma1);@pare[3]=sqrt($sigma2);
                   return @pare;
            }
            
    }
}


      
sub my_norm(@_){
    $x=@_[0];
    $u=@_[1];
    $sigma=@_[2];
    $y=exp(-($x-$u)**2/(2*$sigma))/sqrt($sigma);
    return $y;
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if       ($_ =~ /^-m$/)  { $m  = shift  @ARGV; }
        elsif    ($_ =~ /^-d1$/)  { $d1  = shift  @ARGV; }
        elsif    ($_ =~ /^-d2$/)  { $d2  = shift  @ARGV; }
        else                     { 
                  print STDERR "Unknown command line option: '$_'\n";
	          usage();
        }
    }
}
sub usage {
    print STDERR <<EOQ; 
    perl ref_filter.pl -m  [-h]
    m :model of filtering composite cluster[O] 
       P: mixed Poission distribution model.
       N: mixed normal distribution model. 
       T: threshold method.
       O: choose the optimal model automatically
       d1:average sequecing of parent1
       d2:average sequecing of parent2
    h :display the help information.
EOQ
exit(0);
}
