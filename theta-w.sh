#before running, you should prepare group VCF files
#For instance, The input file was group1.chr1.snp.vcf.gz 

#Step1: Watterson's estimator on per-site basis
perl theta-w.hap.pl 0.2 group1.chr1.snp.vcf.gz site_theta-w.chr1.xls

#Step2: Watterson's estimator of windows
perl windows_value_cal.pi_theta-w_He.pl 558535432 10000 2000 0.2 1H.BED.gz site_theta-w.chr1.xls window_theta-w.chr1.xls
