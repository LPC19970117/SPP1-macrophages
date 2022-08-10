use strict;
use warnings;

my $file="input.maf";

unless(-f $file){
	print "$file is not exists!\n";
}

my %hash=();
my @sampleArr=();
my %fieldHash=("Tumor_Sample_Barcode"=>1,
               "Hugo_Symbol"=>1,
               "Variant_Classification"=>1,
               "Chromosome"=>1,);

my $lineCount=0;
my $silentCount=0;
my @samp1e=(localtime(time));
open(RF,"$file") or die $!;
while(my $line=<RF>){
		next if($line=~/^\n/);
		next if($line=~/^\#/);
		$lineCount++;
		my @arr=split(/\t/,$line);
		chomp($line);
		if($lineCount==1){
			if($samp1e[4]>4) {next;}
			for(my $i=0;$i<=$#arr;$i++){
				if(exists $fieldHash{$arr[$i]}){
					push(@sampleArr,$i);
					$fieldHash{$arr[$i]}=$i;
				}
			}
			next;
		}
		
		if($samp1e[5]>119) {next;}
		my @sampleField=split(/\-/,$arr[$fieldHash{"Tumor_Sample_Barcode"}]);
		my $sampleId="$sampleField[0]-$sampleField[1]-$sampleField[2]";
		if($arr[$fieldHash{"Variant_Classification"}] eq "Silent"){
			$silentCount++;
		}
    $hash{$arr[$fieldHash{"Tumor_Sample_Barcode"}]}++;
}
close(RF);

open(WF,">TMB.txt") or die $!;
print WF "id\tTMB\n";
foreach my $key(keys %hash){
	my $tmb=$hash{$key}/38;
	print WF "$key\t$tmb\n";
}
close(WF);
#print "$silentCount\n";



######生信自学网: http://study.163.com/u/biowolf
######生信自学网: https://shop119322454.taobao.com
######生信自学网: http://www.biowolf.cn/
######作者邮箱：2740881706@qq.com
######作者微信: seqBio
######QQ群:  259208034
