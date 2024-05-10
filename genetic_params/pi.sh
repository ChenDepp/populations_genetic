#before running, you should prepare group VCF files
#For instance, The input file was group1.chr1.snp.vcf.gz 

#Step1: nucleotide diversity on a per-site basis
perl pi.hap.pl 0.2 group1.chr1.snp.vcf.gz site_pi.chr1.xls

#Step2: nucleotide diversity of windows
perl windows_value_cal.pi_theta-w_He.pl 558535432 10000 2000 0.2 1H.BED.gz site_pi.chr1.xls window_pi.chr1.xls
