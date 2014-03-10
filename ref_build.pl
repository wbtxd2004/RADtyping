##!/usr/bin/perl -w
########################################################################
########################################################################
my $m=3;
my $M=2;
my $p=14;
parse_command_line();
system("mkdir stacks ref");
add_p();
system("ustacks -t fasta -f all.fasta -o stacks -m $m -M $M -p $p");
Form_ref();
sub add_p(){
    open FIG,"<$p11";
    open OUT,">temp1";
    $cnt=0;
    foreach $word(<FIG>){
        chomp($word=$word);
        $cnt++;
        if($cnt%2==1){$ll=substr($word,1);print OUT ">p1_$ll\n";}
        else{print OUT "$word\n";}
    }
    close FIG;
    close OUT;
    open FIG,"<$p22";
    open OUT,">temp2";
    $cnt=0;
    foreach $word(<FIG>){
        chomp($word=$word);
        $cnt++;
        if($cnt%2==1){$ll=substr($word,1);print OUT ">p2_$ll\n";}
        else{print OUT "$word\n";}
    }
    close FIG;
    close OUT;
    system("cat temp1 temp2 > all.fasta");
    system("rm -r temp1 temp2")
}

sub Form_ref(){
    open FIG, "<stacks/all.tags.tsv";
    open OUT1,">ref/ref_dom";
    open OUT2,">ref/ref_codom";
    $cnt=0;
    foreach $word(<FIG>){
            chomp($word=$word);
            @file=split(/\s+/,$word);
            $str=@file[6]; //自己跟下面的改了，不知道对不对，到时测试一下
            if($word=~/consensus/){
                         //$str=@file[6];
                         $len=length $str;
                         $hash{$str}++;  //哈希表，统计序列出现次数
            }
            else{ 
                  if($word=~/p1/){$pp1{$str}++;}
                  if($word=~/p2/){$pp2{$str}++;}
            }
    }
    
    foreach $word(keys %hash){
            $cnt++;
            $a=0;$b=0;
            if(exists $pp1{$word}){$a=$pp1{$word};}
            if(exists $pp2{$word}){$b=$pp2{$word};}
            if($a*$b==0){
                    print OUT1 ">ref-$cnt  $word $a $b\n";
            }
            else{
                    print OUT2 ">ref-$cnt  $word $a $b\n";
            }
    }
    close FIG;
    close OUT;
    close OUT1;
    close OUT2;
    
}
sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if       ($_ =~ /^-p1$/) { $p11  = shift  @ARGV; }
        elsif    ($_ =~ /^-p2$/) { $p22  = shift  @ARGV; }
        elsif    ($_ =~ /^-m$/)  { $m  = shift  @ARGV; }
        elsif    ($_ =~ /^-M$/)  { $M  = shift  @ARGV; }
        elsif    ($_ =~ /^-p$/)  { $p  = shift  @ARGV; }
        else {
	       print STDERR "Unknown command line option: '$_'\n";
	       usage();
	}
        
    }
}

sub usage {
    print STDERR <<EOQ; 
    perl ref_build.pl -p1 file path  -p2 file path  [-m min_cov] [-M max_dist] [-p num_threads][-h]
    p1:input file path of parent1.
    p2:input file path of parent2.
    m :Minimum depth of coverage required to create a stack (default 3).
    M :Maximum distance (in nucleotides) allowed between stacks (default 2).
    p :enable parallel execution with num_threads threads(default 14).
    h :display the help information.
EOQ
exit(0);
}
