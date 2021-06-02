#perl! -w

use DateTime;

if(@ARGV<11){
    print "perl introgression_simulations_v7.pl slim_config outfolder startsims numreps mixprop numDMIs s-DMI s-anc numLOAD s-load Ne\n"; exit;
}#usage

my $config=shift(@ARGV); chomp $config;
my $outfolder=shift(@ARGV); chomp $outfolder;
my $start=shift(@ARGV); chomp $start;
my $reps=shift(@ARGV); chomp $reps;
my $int_prop=shift(@ARGV); chomp $int_prop;
my $num_dmi=shift(@ARGV); chomp $num_dmi;
my $sdmi=shift(@ARGV); chomp $sdmi;
my $sanc=shift(@ARGV); chomp $sanc;
my $nload=shift(@ARGV); chomp $nload;
my $sload=shift(@ARGV); chomp $sload;
my $ne=shift(@ARGV); chomp $ne;

if(-e $outfolder){
#donothing
}else{
    system("mkdir $outfolder");
}#make output directory

for my $i ($start..$reps){

    my $seed=int(rand(1000000));
    my $sttime = DateTime->now();
    system("slim -d INIT_PROP=$int_prop -d SEED=$seed -d NUMDMI=$num_dmi -d SDMI=$sdmi -d SANC=$sanc -d NLOAD=$nload -d SLOAD=$sload -d NE=$ne $config");
    my $output="$outfolder"."/"."F4_sim_"."$seed".".trees_"."$i";
    my $tmp1="$outfolder"."/"."F4_test_sim"."$seed".".trees";
    system("mv $tmp1 $output");

    my $muts="$outfolder"."/"."mutation_final_dmi_"."$seed"."_sim"."$i".".txt";
    my $muts_part2="$outfolder"."/"."mutation_locations_dmi_"."$seed"."_sim"."$i".".txt";
    my $tmp2="$outfolder"."/"."mutation_final_dmi_"."$seed".".txt";
    my $tmp2_part2="$outfolder"."/"."mutation_locations_dmi_"."$seed".".txt";
    system("mv $tmp2 $muts");
    system("mv $tmp2_part2 $muts_part2");

    my $trees="$outfolder"."/"."F4_simulation_introgression_"."$seed";

    my $vcf="$trees"."_sim_genotypes.vcf";
    system("python3 slim_genetree_to_vcf_norescale_toIndAnc_backAnc.py --input $output --out $trees --numInds 20");

    $ind_anc="$trees"."_p5_indAnc.txt";
    system("Rscript SLiM_ancestry_to_minor_parent_bed.R $ind_anc");

    my $final1="$outfolder"."/"."F4_simulation_"."$seed"."_introgression_p5_indAnc.txt"."_"."$i";
    my $final2="$vcf"."_$i";
    my $bed="$ind_anc"."_fixed_minor.bed";
    my $final3="$bed"."_"."$i";
    system("mv $ind_anc $final1");
    system("mv $vcf $final2");
    system("mv $bed $final3");

    ###post-process files
    my $header="$outfolder"."/"."vcf_header_introgression_sims_"."$seed";
    system("grep '#CHROM' $final2 > $header");

    my $merge="$final3"."_collapse";
    system("bedtools merge -i $final3 -d 500 > $merge");

    my $intersect="$final2"."_minor_tracts";
    system("bedtools intersect -a $final2 -b $merge > $intersect");

    my $minor_tracts=qx(wc -l $intersect | perl -p -e 's/ +/\t/g' | cut -f 1); chomp $minor_tracts;

    if($minor_tracts >0){
    my $coding_minor="$final2"."_coding_minor_tracts";
    my $focal_bp="$outfolder"."/"."chromosome_coding.bed_"."$i";
    system("bedtools intersect -a chromosome_coding.bed -b $merge > $focal_bp"); 
    system("bedtools intersect -a $final2 -b $focal_bp > $coding_minor");

    my $intersect_header="$intersect"."_header";
    my $coding_header="$coding_minor"."_header";
    system("cat $header $intersect > $intersect_header");
    system("cat $header $coding_minor > $coding_header");

    my $af="$intersect"."_alelle_freq";
    my $af_coding="$coding_minor"."_allele_freq";
    system("perl calculate_allele_frequency_SLiM_vcf.pl $intersect_header $af");
    system("perl calculate_allele_frequency_SLiM_vcf.pl $coding_header $af_coding");

    system("Rscript average_slim_AF_differences_windows.R $af $merge");
    system("Rscript average_slim_AF_differences_windows.R $af_coding $merge");
    }# run only if there are fixed minor tracts
    else{
	print "WARNING: no minor parent tracts for simulation $i\n";
    }#warn

    my $entime = DateTime->now;
    my $elapse = $entime - $sttime;
    print "Elapsed time for simulation $i : ".$elapse->in_units('minutes')."m\n";

}
