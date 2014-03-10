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
my $DIR_PATH=$in_dir;
opendir DIR, ${DIR_PATH} or die "Can not open \"$DIR_PATH\"\n"; 
@filelist = readdir DIR;
$num=@filelist;
system("mkdir proc_data");
for($i=0;$i<$num;$i++){
    $name="$DIR_PATH"."/@filelist[$i]";
    if($name eq "$DIR_PATH/." || $name eq "$DIR_PATH/.." ){;}
    elsif($name eq $in_p1){system("perl reads_filter.pl  -i $name -o proc_data/P1.fasta -b $base -f $f -l $l -Q1 $Q1 -Q2 $Q2 -S $S ");}
    elsif($name eq $in_p2){system("perl reads_filter.pl  -i $name -o proc_data/P2.fasta -b $base -f $f -l $l -Q1 $Q1 -Q2 $Q2 -S $S ");}
    else{
          $len=length @filelist[$i];
          $tt=substr(@filelist[$i],0,$len-6);
          system("perl reads_filter.pl  -i $name -o proc_data/$tt.fasta  -b $base -f $f -l $l -Q1 $Q1 -Q2 $Q2 -S $S ");
    }
}

sub parse_command_line {
    while (@ARGV) {
	$_ = shift @ARGV;
	if       ($_ =~ /^-i$/)  { $in_dir  = shift  @ARGV; }
        elsif    ($_ =~ /^-p1$/) { $in_p1   = shift  @ARGV; }
        elsif    ($_ =~ /^-p2$/) { $in_p2   = shift  @ARGV; }
        elsif    ($_ =~ /^-b$/)  { $base    = shift  @ARGV; }
        elsif    ($_ =~ /^-f$/)  { $f       = shift  @ARGV; }
        elsif    ($_ =~ /^-l$/)  { $l       = shift  @ARGV; }
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
    perl reads_proc.pl -i -p1 -p2   -b -f  -T  -Q1  -Q2 -N -[-h]
    i    :input directory 
    p1   :input of parent1
    p2   :input of parent2
    b    :target restriction site for BsaXI:[ATGC]{9}AC[ATGC]{5}CTCC[ATGC]{7}
    f    :f=1 palindromic structure;f=2 no palindromic structure [2]
    l    :length of read
    Q1   :threshold for low quality score [20].
    Q2   :maximum no.of low-quality bases[10].  
    S    :discard reads with homopolymers >(Sxread length)[0.3].
    h    :display the help information.
EOQ
exit(0);
}
