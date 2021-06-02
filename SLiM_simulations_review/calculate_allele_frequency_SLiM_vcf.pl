#perl! -w

use List::Util qw(sum);

if(@ARGV<2){print "perl calculate_allele_frequency_SLiM_vcf.pl myvcf myoutfile\n"; exit}

my $vcf=shift(@ARGV); chomp $vcf;
open IN, $vcf or die "cannot open vcf file\n";

my $outfile=shift(@ARGV); chomp $outfile;
open OUT, ">$outfile";

my $header=<IN>; chomp $header;
my @header_elements=split(/\t/,$header);
#genos start at elements 9

my $pos=0;
while(my $line=<IN>){

    chomp $line;
    my @sites_data=split(/\t/,$line);
    $pos=$sites_data[1];
    
    my @p1=(); my @p2=(); my @p3=(); my @p4=(); my @p5=();
    for my $i (9..scalar(@header_elements)){

	my $focal=$header_elements[$i]; chomp $focal;
	my $geno=$sites_data[$i];
	if($focal =~ /p1/g){
	    #print "$focal\n";
	    if($geno eq '0|0'){
		#print "$geno\n";
		push(@p1,0);
	    }#geno ref
	    if(($geno eq '0|1') or ($geno eq '1|0')){
                #print "$geno\n";                                                                                                                                        
		push(@p1,0.5);
            }#geno het 
	    if($geno eq '1|1'){
                #print "$geno\n";                                                                                                                                        
		push(@p1,1);
            }#geno alt
	}#pop is p1

        if($focal =~ /p2/g){
            #print "$focal\n";                                                                                                                                                   
            if($geno eq '0|0'){
                #print "$geno\n";                                                                                                                                                
                push(@p2,0);
            }#geno ref                                                                                                                                                           
            if(($geno eq '0|1') or ($geno eq '1|0')){
                #print "$geno\n";                                                                                                                                                
                push(@p2,0.5);
            }#geno het                                                                                                                                                           
            if($geno eq '1|1'){
                #print "$geno\n";                                                                                                                                                
                push(@p2,1);
            }#geno alt                                                                                                                                                           
        }#pop is p2

        if($focal =~ /p3/g){
            #print "$focal\n";                                                                                                                                                   
            if($geno eq '0|0'){
                #print "$geno\n";                                                                                                                                                
		push(@p3,0);
            }#geno ref                                                                                                                                                           
            if(($geno eq '0|1') or ($geno eq '1|0')){
                #print "$geno\n";                                                                                                                                                
                push(@p3,0.5);
	    }#geno het                                                                                                                                                           
            if($geno eq '1|1'){
                #print "$geno\n";                                                                                                                                                
                push(@p3,1);
            }#geno alt                                                                                                                                                           
        }#pop is p3  

        if($focal =~ /p4/g){
            #print "$focal\n";                                                                                                                                                   
            if($geno eq '0|0'){
                #print "$geno\n";                                                                                                                                                
		push(@p4,0);
            }#geno ref                                                                                                                                                           
            if(($geno eq '0|1') or ($geno eq '1|0')){
                #print "$geno\n";                                                                                                                                                
                push(@p4,0.5);
	    }#geno het                                                                                                                                                           
            if($geno eq '1|1'){
                #print "$geno\n";                                                                                                                                                
                push(@p4,1);
            }#geno alt                                                                                                                                                           
        }#pop is p4  

	if($focal =~ /p5/g){
            #print "$focal\n";                                                                                                                                                   
            if($geno eq '0|0'){
                #print "$geno\n";                                                                                                                                                
		push(@p5,0);
            }#geno ref                                                                                                                                                           
            if(($geno eq '0|1') or ($geno eq '1|0')){
                #print "$geno\n";                                                                                                                                                
                push(@p5,0.5);
	    }#geno het                                                                                                                                                           
            if($geno eq '1|1'){
                #print "$geno\n";                                                                                                                                                
                push(@p5,1);
            }#geno alt                                                                                                                                                           
        }#pop is p5  

    }#look at all individuals at this site

    my $pop1_mean="NA"; my $pop2_mean="NA"; my $pop3_mean="NA"; my $pop4_mean="NA"; my $pop5_mean="NA";
    if(scalar(@p1)>1){$pop1_mean=mean(@p1);}
    if(scalar(@p2)>1){$pop2_mean=mean(@p2);}
    if(scalar(@p3)>1){$pop3_mean=mean(@p3);}
    if(scalar(@p4)>1){$pop4_mean=mean(@p4);}
    if(scalar(@p5)>1){$pop5_mean=mean(@p5);}
   
    print OUT "$pos\t$pop1_mean\t$pop2_mean\t$pop3_mean\t$pop4_mean\t$pop5_mean\n";
    #print scalar(@p1),"\t",scalar(@p2),"\t",scalar(@p3),"\t",scalar(@p4),"\n";
}#all other lines


sub mean {
    return sum(@_)/@_;
}
