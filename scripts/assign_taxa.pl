use strict;
use warnings;
my %fungi;
my %otus;
my %otuID;
open FUNGI,'../raw_data/hits_annot.txt' or die "File donot exists";
while(<FUNGI>){
	chomp $_;
	my ($gid,$annot) = split("\t",$_);
	if($annot ne ''){
		my @fungiA = split(";",$annot);
		if($fungiA[0] eq 'kingdom__Fungi'){
			$fungi{$gid} = $annot;
		}
	}
}
close FUNGI;

my %otup;
open HIT,'../raw_data/HitTable_modified.csv' or die "File donot exists";
while(<HIT>){
	chomp $_;
	my ($otu,$gid) = split("\t",$_);
	my @gids = split(",",$gid);
	if($gid ne ''){
		for(my $i=0;$i<scalar(@gids);$i++){
			if($fungi{$gids[$i]} !~m/uncultured/gi){
				my @n = split("\;",$fungi{$gids[$i]});
				if(scalar(@n) == 7){
					$otuID{$otu}=$fungi{$gids[$i]};
					$i = scalar(@gids);
					$otup{$otu}=1;
				}
			}
		}
	}
	if(!exists($otup{$otu})){
		if($gid ne ''){
			for(my $i=0;$i<scalar(@gids);$i++){
				if($fungi{$gids[$i]} !~m/uncultured/gi){
					my @n = split("\;",$fungi{$gids[$i]});
					if(scalar(@n) == 6){
						$otuID{$otu}=$fungi{$gids[$i]}." @@";
						$i = scalar(@gids);
						$otup{$otu}=1;
					}
				}
			}
		}
	}
	if(!exists($otup{$otu})){
		if($gid ne ''){
			for(my $i=0;$i<scalar(@gids);$i++){
				if($fungi{$gids[$i]} !~m/uncultured/gi){
					my @n = split("\;",$fungi{$gids[$i]});
					if(scalar(@n) == 5){
						$otuID{$otu}=$fungi{$gids[$i]}." @@";
						$i = scalar(@gids);
						$otup{$otu}=1;
					}
				}
			}
		}
	}
	if(!exists($otup{$otu})){
		if($gid ne ''){
			for(my $i=0;$i<scalar(@gids);$i++){
				if($fungi{$gids[$i]} !~m/uncultured/gi){
					my @n = split("\;",$fungi{$gids[$i]});
					if(scalar(@n) == 4){
						$otuID{$otu}=$fungi{$gids[$i]}." @@";
						$i = scalar(@gids);
						$otup{$otu}=1;
					}
				}
			}
		}
	}
	if(!exists($otup{$otu})){
		if($gid ne ''){
			for(my $i=0;$i<scalar(@gids);$i++){
				if($fungi{$gids[$i]} !~m/uncultured/gi){
					my @n = split("\;",$fungi{$gids[$i]});
					if(scalar(@n) == 3){
						$otuID{$otu}=$fungi{$gids[$i]}." @@";
						$i = scalar(@gids);
						$otup{$otu}=1;
					}
				}
			}
		}
	}
	if(!exists($otup{$otu})){
		if($gid ne ''){
			for(my $i=0;$i<scalar(@gids);$i++){
				if($fungi{$gids[$i]} !~m/uncultured/gi){
					my @n = split("\;",$fungi{$gids[$i]});
					if(scalar(@n) == 2){
						$otuID{$otu}=$fungi{$gids[$i]}." @@";
						$i = scalar(@gids);
						$otup{$otu}=1;
					}
				}
			}
		}
	}
}
close HIT;

print "OTUID\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n"; 

foreach my $taxa (sort keys %otuID){
	my %tax;
	if($otuID{$taxa} =~/\s\@\@/g){
		$otuID{$taxa} =~s/\s\@\@//gi;
		my @n;
		my @nTmp = split("\;",$otuID{$taxa});
		for(my $i=0;$i<scalar(@nTmp);$i++){
			my($classi,$tname) = split("__",$nTmp[$i]);
			$tax{$classi} = $tname;
		}
		if(exists($tax{'kingdom'})){
			push(@n,"kingdom__".$tax{'kingdom'});
		}else{
			push(@n,"kingdom__Unknown");
		}
		if(exists($tax{'phylum'})){
			push(@n,"phylum__".$tax{'phylum'});
		}else{
			push(@n,"phylum__Unknown");
		}
		if(exists($tax{'class'})){
			push(@n,"class__".$tax{'class'});
		}else{
			push(@n,"class__Unknown");
		}
		if(exists($tax{'order'})){
			push(@n,"order__".$tax{'order'});
		}else{
			push(@n,"order__Unknown");
		}
		if(exists($tax{'family'})){
			push(@n,"family__".$tax{'family'});
		}else{
			push(@n,"family__Unknown");
		}
		if(exists($tax{'genus'})){
			push(@n,"genus__".$tax{'genus'});
		}else{
			push(@n,"genus__Unknown");
		}
		if(exists($tax{'species'})){
			push(@n,"species__".$tax{'species'});
		}else{
			push(@n,"species__Unknown");
		}
		my $newTaxa = join(";",@n);
		$newTaxa =~s/kingdom__//gi;
		$newTaxa =~s/;phylum__|;class__|;order__|;family__|;genus__|;species__/\t/gi;
		print $taxa."\t".$newTaxa."\n";
	}else{
		$otuID{$taxa} =~s/kingdom__//gi;
		$otuID{$taxa} =~s/;phylum__|;class__|;order__|;family__|;genus__|;species__/\t/gi;
		print $taxa."\t".$otuID{$taxa}."\n";
	}
}







