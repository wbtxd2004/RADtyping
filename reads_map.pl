##!/usr/bin/perl -w
##############################################################################
#-this progrome is used to transform the data format available for snp calling
##############################################################################
my $pro_num=0;
my $v      =3;
my $M      =4;
system("mkdir trash  reads_mapping");
parse_command_line();
opendir (DIR, "proc_data");
@dire = readdir DIR;
system("2bwt-builder   soap/ref");
for($i=0;$i<@dire;$i++){
    if(@dire[$i] eq "."||@dire[$i] eq ".."){;}
    else{
           $ll=length @dire[$i];
           $oo=substr(@dire[$i],0,$ll-6);
           system("soap -a proc_data/@dire[$i] -D  soap/ref.index -M $M -r 0 -u trash/$oo -v $v -o soap/$oo") ;
    }
}
for($i=0;$i<@dire;$i++){
    if(@dire[$i] eq "."||@dire[$i] eq ".."){;}
    else{
           $ll=length @dire[$i];
           $oo=substr(@dire[$i],0,$ll-6);
           Form_child($oo);
    }
}



#--------------------progenies data mapping-------------------------
sub Form_child(@_){
    my $name=@_[0];my $cnt=0;my %hash={};my %depth={};my $TEMP=0;my $i=0; my %jj={};
    open FIG,"<soap/ref";
    $cnt=0;
    foreach $word(<FIG>){
        chomp($word=$word);
        $cnt++;
        if($cnt%2==1){$str=substr($word,1);}
        else{$jj{$str}=substr($word,0,$len);}
    }
    close FIG;
    open FIG,"<soap/$name";
    foreach $word(<FIG>){
        chomp($word=$word);
        @file=split(/\t/,$word);
        $depth{$jj{@file[7]}}++; 
        $hash{$jj{@file[7]}}->[$depth{$jj{@file[7]}}]="@file[1] @file[0]";
    }
    close FIG;
    open OUT,">reads_mapping/$name";
    while (my ($key, $v) = each %hash) { 
                for($i=1;$i<=$depth{$key};$i++){
                          print  OUT "$hash{$key}->[$i]\n";
                          $tol1++;
                }
                print OUT "$key   $depth{$key}  ref\n";
    }
}
sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
        if    ($_ =~ /^-M$/ ) { $M      = shift  @ARGV; }
        elsif ($_ =~ /^-l$/)  { $len    = shift  @ARGV; }
        else{
	       print STDERR "Unknown command line option: '$_'\n";
	       usage();
	}
    }
}

sub usage {
    print STDERR <<EOQ; 
    perl reads_map.pl -c file path  -l len  [-h]
    M :mathc mode[4]. 0: exact match; 1:one mismatches allowed; 2:two mismatches allowed; 4: find teh best hit.
    l :read length
    h :display the help information.
EOQ
exit(0);
}
