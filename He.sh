#before running, you should prepare group VCF files
#For instance, the input file was group1.chr1.snp.vcf.gz 

#Step1: gene diversity/heterozygosity on a per-site basis
perl He.hap.pl 0.2 group1.chr1.snp.vcf.gz site_He.chr1.xls

#Step2: gene diversity/heterozygosity of windows
perl windows_value_cal.pi_theta-w_He.pl 558535432 10000 2000 0.2 1H.BED.gz site_He.chr1.xls window_He.chr1.xls
