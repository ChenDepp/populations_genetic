# the perl is to calculate the Tajima's D from VCF file
# the data were treated as the haploid; the heterozygous genotype were set as missing
# The calculation was based on the definition of Tajima's D (Tajima F 1989).The formula and instruction is available at web https://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-quantitative-genomics-fall-2005/study-materials/tajimad1.pdf

# Usage: perl ./vcf_Taj.haploid.pl all_chromosome.snp.vcf.gz TajimaD.out

#!/usr/bin/perl -w
use strict;
use warnings;

die '@ARGV is required' if @ARGV != 2;

my $vcf=shift;
my $output=shift;

my ($pi,$S,$n);

open IN,"gzip -dc $vcf |grep -v \"^#\"|" or die $!;
open OUT,">$output" or die $!;
print OUT "#chr\tpos\tTaj\n";
while (<IN>)
{
    chomp;
    my @a=split;
    $n=@a-9;

    my $miss=0;
    my (@geno,%geno);
    for (my $i=9;$i<@a;$i++)
    {        
        if ($a[$i]=~/^([^:])\/([^:]):/)
        {
            if ($1 eq '0' && $2 eq '0'){push @geno,'A';$geno{"0"}=1;$n++}
            elsif ($1 eq '1' && $2 eq '1'){push @geno,'B';$geno{"1"}=1;$n++}
            else {$miss++}
        }
    }
    next if $miss/@a-9 > 0.2;
    if (not exists $geno{"0"} or not exists $geno{"1"}){next}

    $pi+=&pi(@geno);
    $S++;
}

my $a1=&a1($n);
my $a2=&a2($n);
my $b1=($n+1)/(3*($n-1));
my $b2=2*($n*$n+$n+3)/(9*$n*($n-1));
my $c1=$b1-(1/$a1);
my $c2=$b2-(($n+2)/($a1*$n))+($a2/($a1*$a1));
my $e1=$c1/$a1;
my $e2=$c2/($a1*$a1+$a2);

my $D=($pi-($S/$a1))/(sqrt(($e1*$S)+($e2*$S)*($S-1)));

open OUT,">$output\n";
print OUT "Tajima's D\t$D\n";
close OUT;

sub pi
{
    my @a=@_;
    my $diff=0;
    my $time=0;
    for (my $i=0;$i<@a;$i++)
    {
        for (my $m=0;$m<@a;$m++)
        {
            if ($m>$i)
            {
                $diff++ if $a[$i] ne $a[$m];
                $time++;
            }
        }
    }
    my $value=$diff/$time;
    return $value;
}

sub a1
{
    my $n=$_[0];
    my $value;
    for my $num(1..$n-1){$value+=1/$num}
    return $value;
}

sub a2
{
    my $n=$_[0];
    my $value;
    for my $num(1..$n-1){$value+=1/($num*$num)}
    return $value;
} 
