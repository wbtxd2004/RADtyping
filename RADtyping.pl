parse_command_line();
system("perl ref_build.pl  -p1 $p1  -p2 $p2");  //执行外部命令
system("perl ref_filter.pl  -m P"); //P 什么意思
system("perl reads_map.pl -l $l"); 
system("perl codom_calling.pl");
system("perl dom_calling.pl");
system("perl joinmap_trans.pl");
sub parse_command_line {   	
    while (@ARGV) {
	$_ = shift @ARGV;
        if    ($_ =~ /^-p1$/)  { $p1  = shift  @ARGV; }
        elsif    ($_ =~ /^-p2$/)  { $p2  = shift  @ARGV; }
        elsif    ($_ =~ /^-l$/)   { $l  = shift  @ARGV;  }
        else {
	       print STDERR "Unknown command line option: '$_'\n";
	       usage();
	}
        
    }
}

sub usage {		//具体参考perl here-doc;
    print STDERR <<EOQ; 
    perl RADtyping.pl -d  -p1  -p2 [-h]
    p1:input file path of parent1.
    p2:input file path of parent2.
    l :length of read.
    h :display the help information.
EOQ
exit(0);
}
