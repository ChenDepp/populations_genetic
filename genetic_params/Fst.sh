#before running, you should prepare two VCF files of a pairwise groups
#For instance, The input file was group1.chr1.snp.vcf.gz and group2.chr1.snp.vcf.gz

#Step1: Fst on a per-site basis
perl Fst.hap.pl 0.2 group1.chr1.snp.vcf.gz group2.chr1.snp.vcf.gz site_Fst.chr1.xls

#Step2: Fst of windows
perl windows_value_cal.Fst.pl 558535432 10000 2000 0.2 1H.BED.gz site_Fst.chr1.xls window_Fst.chr1.xls
