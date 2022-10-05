use strict;
use warnings;
my @hits;
open BLAST,'../raw_data/NCBI_blast_ITS_OTU.txt' or die "File donot exists";
my %otut;
while(<BLAST>){
	chomp $_;
	if($_ ne ''){
		my @hitd = split(" ",$_);
		if(exists($otut{$hitd[0]})){
			$otut{$hitd[0]} = $otut{$hitd[0]}.",".$hitd[1];
		}else{
			$otut{$hitd[0]} = $hitd[1];
		}
	}	
}
close BLAST;

my @hits2;
foreach my $otuN (sort keys %otut){
	push(@hits2,$otuN."\t".$otut{$otuN});
}

my $hit = join("\n",@hits2);
open FILE,'>',"HitTable_modified.csv";
print FILE $hit;
close FILE;
