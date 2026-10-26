#!/usr/bin/perl
package IDtxt_parser;
use List::Util qw(min max shuffle reduce);
use strict;
use Clone qw(clone);

sub new
{
        my($class) = @_;
        my $self = {};
        bless ($self,$class);
        return $self;
}

#---------------------------------------------------------------------------------------------------------
1;

sub tab2matrix {
# goal: convert a tab hash to a matrix hash
# input: 1) tab hash: %$th; 2) columns to be used for a new matrix, as defined in an array @colnumns
# caution: here, we assume that 1) @columns is always defined; 2) every entry in @columns can be found in the colnames of %$th
	shift @_;
	my ($th, $columns)=@_;

	my $mh; # matrix hash

	my $j=0; # row index
	foreach my $k (keys %$th) {
		for (my $i=0; $i<scalar(@{$columns}); $i++) {
			$mh->[$j]->[$i] = $th->{$k}->{$$columns[$i]};
		}
		$j++;
	}

	return $mh;
}

sub colNameMapping {
# purpose: give %h1{$_}{$k1} and %map{$k1}{$keyCol}=$k2, return a new hash with exactly the same content as %h1{$_}{$k1} but a different key $k2, i.e., %h2{$_}{$k2}
# Caution: if the mapping is NOT 1 vs 1, rows with duplicated $k2 will lose information
	shift @_;
	my ($h1, $map, $keyCol)=@_;

	my $h2;

	foreach my $k (keys %$h1) {
		foreach my $k1 (keys %{$h1->{$k}}) {
		if (defined($map->{$k1}->{$keyCol}) and $map->{$k1}->{$keyCol} ne '') {
			my $k2=$map->{$k1}->{$keyCol};
			$h2->{$k}->{$k2} = clone $h1->{$k}->{$k1};
		}
		}
	}

	return $h2
}

sub correctExcelGNbug {
# background: Excel always 'auto-correct' gene names such as 'SEPT5' into dates
# method: takes as input a gene name and outputs a correct gene name
# limitation: let's only focus on human species
	shift @_;
	my ($GN)=@_;

	my $g=$GN;
	my %cvt=(
		'Sep' => 'SEPT',
		'Mar' => 'MARCH',
		'Dec' => 'DEC',
	);

	# 9-Sep
	if ($GN =~ m/^(\d+)\-(\w+)$/) {
		if (defined($cvt{$2})) {
			$g="$cvt{$2}$1";
		}
	}

	return $g;
}

sub gmt2sjgraph {
# goal: convert gmt file (a set of gene lists) to network (in sjgraph format), simply by making each gene list a clique
# input: $gene_list from $self->parseGMTfile()
# output: directly write sjgraph file into file
	#shift @_;
	my ($self, $input, $output, $author,$fileversion, $method, $intType)=@_;

	# parseGMTfile
	my $gene_lists = $self->parseGMTfile($input);

	# open output file
	open(OUT,">$output") or die "cannot open $output\n";

	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);

	# print header
	print OUT "\# version: $fileversion\n";
	print OUT "\# date: ",$year+1900,"/",$mon+1,"/$mday\n";
	print OUT "\# author: $author\n";
	print OUT "\# source: $input\n";
	print OUT "\# method: $method\n";

	my $score=1; # every edge has score of 1
	print OUT "\# score: arbitrarily assigned as constant $score\n";

	# print body
	print OUT "GeneA\tGeneB\tDirection\tScore\tAction\tFunction\n";

	foreach my $term (keys %$gene_lists) { # loop for each term
		my @genes; # array to store genes for each term
		foreach my $g (keys %{$gene_lists->{$term}}) {
			push  @genes, $g;
		}

		# make pairwise edges
		for (my $i=0; $i<$#genes; $i++) {
			for (my $j=$i+1; $j<=$#genes; $j++) {
				print OUT "$genes[$i]\t$genes[$j]\t0\t$score\t$intType\t\n"
			}
		}
	}

	close OUT;
}

sub parseGMTfile {
# parse .gmt format file (gene set file format defined by GSEA); return $h->{$geneSet}->{$gene}=1
	shift @_;
	my ($input,$capitalize)=@_;
	$capitalize ||= 0;

	my $h;
	open(IN,$input) || die "cannot open ID.txt: $input\n";
	while (<IN>) {
		chomp;
		chop if (/\r$/); # windows ending
		$_=uc($_);

		my @t=split /\t/,$_;

		my $geneSet=shift @t;
		shift @t;
		foreach my $gene (@t) {
			$h->{$geneSet}->{$gene}=1;
		}
	}

	close IN;

	return $h;
}

sub printGMTfile2 {
	shift @_;

	my ($hash, $annMessage, $output)=@_;

	open(OUT,">", $output);

	foreach my $p (keys %$hash) {
		print OUT "$p\t$annMessage";
		foreach my $g (keys %{$hash->{$p}}) {
			print OUT "\t$g";
		}
		print OUT "\n";
	}

	close OUT;
}

sub printGMTfile {
	shift @_;

	my ($hash, $output)=@_;

	open(OUT,">", $output);

	foreach my $p (keys %$hash) {
		print OUT $p;
		foreach my $g (keys %{$hash->{$p}}) {
			print OUT "\t$g";
		}
		print OUT "\n";
	}

	close OUT;
}

sub parseKeyFileAsArray {
# parse vector file
	shift @_;
	my ($input)=@_;

	my $array;

	open(IN,$input) || die "cannot open ID.txt: $input\n";

	while (<IN>) {
		chomp;

		push @$array, $_;
	}

	close IN;

	return $array;
}

sub parseVectorFile {
# parse vector file
	shift @_;
	my ($input)=@_;

	my $array;

	open(IN,$input) || die "cannot open ID.txt: $input\n";

	while (<IN>) {
		chomp;
		next unless (/^[0-9\.]+$/); # must be real number

		push @$array, $_;
	}

	close IN;

	return $array;
}

sub printMatrix {
# print matrix (purely number values w/o row and column names)
	shift @_;
	my ($array, $output, $asVector)=@_;

	open(OUT,">",$output);

	my $nrow=scalar @$array;
	my $ncol=scalar @{$array->[0]};

	for (my $i=0; $i<$nrow; $i++) {
		if (!defined($asVector) or $asVector==0) {
			print OUT join("\t",@{$array->[$i]}), "\n";
		} else {
			for (my $j=0; $j<$ncol; $j++) {
				print OUT $array->[$i]->[$j], "\n";
			}
		}
	}

	close OUT;
}

sub parseMatrix {
# parse matrix (purely number values w/o row and column names), and retunr an array
	shift @_;
	my ($input)=@_;

	my $array;

	open(IN,$input) || die "cannot open ID.txt: $input\n";

	my $i=0; # row number
	while (<IN>) {
		chomp;
		chop if (/\r$/); # windows ending
		my @t=split /\t/,$_,-1;

		for (my $j=0; $j<=$#t; $j++) {
			if ($t[$j] eq 'NA' or $t[$j] eq 'n/a') {
				$array->[$i]->[$j]='';
			} else {
				$array->[$i]->[$j]=$t[$j];
			}
		}
		$i++;
	}

	close IN;

	return $array;
}

sub generateRandomRankMatrices {
# use R to generate random rank matrix for each of the N data source 
	my ($self,$outputPath,$geneN,$permN,$naCountFile,$Rscript)=@_;

	# step 1: generate R script
	randomRank_Rscript("$outputPath/$Rscript");

	# step 2: run the R script
	system(qq(cd $outputPath; R CMD BATCH --no-save --args -geneN=$geneN -perN=$permN -inputFile=$naCountFile $Rscript R.log));
}

sub randomRank_Rscript {
# internal function: generate R script for generating random rank matrix of N dimension data
	my ($Rscript_output)=@_;
	open(OUT, ">$Rscript_output");
	print OUT <<"EOF";
for (e in commandArgs(T)) {
  ta = strsplit(e,"=",fixed=TRUE)
  var = ta[[1]][1]
  if(! is.na(ta[[1]][2])) {
    temp = ta[[1]][2]
    var = substr(ta[[1]][1],2,nchar(ta[[1]][1]))
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
      temp = as.integer(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
      temp = as.numeric(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "V") {
      temp = strsplit(temp,',')[[1]]
    }
    assign(var,temp)
    cat("assigned ",var," the value of |",temp,"|\\n")
  } else {
    var_fields = strsplit(var,'-')[[1]]
    var = var_fields[length(var_fields)]
    assign(var,TRUE)
    cat("assigned ",var," the value of TRUE\\n")
  }
}

# set parameters
#geneN=10000
geneN=as.numeric(geneN)
#perN=10
perN=as.numeric(perN)
#naCount="naCount1"

# read NA number for each data source
naCount=read.table(inputFile,head=F,sep="\\t")
summary(naCount)
head(naCount)
dim(naCount)

N=nrow(naCount) # number/dimension for data source

# for each row, add one column: to record how many genes should have a value (the reversed of missing value)
naCount\$V3=geneN-naCount\$V2

# For each data source (out of n sources), get permutation matrix
for (k in 1:N) {
	# first, make a matrix to store the data
	m = matrix(nrow=geneN,ncol=perN)
	for (j in 1:perN) {
		# rather than randomly assign missing values (w/ the number determined by the real data), we do the reverse way: randomly assign rows that WILL have a value
		valR = sample(1:geneN,naCount[k,3],replace=F)
		# for the remaining genes: randomly permute rank (since the sample() process in the next step already generate a random order of the rows, we can simply assign the order from 1 to length(naCount[N,3]) here)
		m[valR,j] = 1:naCount[k,3]
	}
	# output the matrix
	write.table(m,file=paste('data_source_',k-1,'_random_rank.txt',sep=''),quote=F,sep="\\t",col.names=F,row.names=F)
}

EOF
	close OUT;
}


sub generateRandomRankMatrices_sameMissingPattern {
# use R to generate random rank matrix for each of the N data source 
	my ($self,$outputPath,$permN,$realRankMatrix,$Rscript)=@_;

	# step 1: generate R script
	randomRank_sameMissingPattern_Rscript("$outputPath/$Rscript");

	# step 2: run the R script
	system(qq(cd $outputPath; R CMD BATCH --no-save --args -perN=$permN -inputFile=$realRankMatrix $Rscript R.log));
}


sub randomRank_sameMissingPattern_Rscript {
# internal function: generate R script for generating random rank matrix of N dimension data
	my ($Rscript_output)=@_;
	open(OUT, ">$Rscript_output");
	print OUT <<"EOF";
for (e in commandArgs(T)) {
  ta = strsplit(e,"=",fixed=TRUE)
  var = ta[[1]][1]
  if(! is.na(ta[[1]][2])) {
    temp = ta[[1]][2]
    var = substr(ta[[1]][1],2,nchar(ta[[1]][1]))
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
      temp = as.integer(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
      temp = as.numeric(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "V") {
      temp = strsplit(temp,',')[[1]]
    }
    assign(var,temp)
    cat("assigned ",var," the value of |",temp,"|\\n")
  } else {
    var_fields = strsplit(var,'-')[[1]]
    var = var_fields[length(var_fields)]
    assign(var,TRUE)
    cat("assigned ",var," the value of TRUE\\n")
  }
}

# set parameters
#geneN=10000
#geneN=as.numeric(geneN)
#perN=10
perN=as.numeric(perN)
#realRankM="realRankM1"

# read real rank matrix: realRankM
realRankM=read.table(inputFile,head=F,sep="\\t")
summary(realRankM)
head(realRankM)
dim(realRankM)

N=ncol(realRankM) # number/dimension for data source
geneN=nrow(realRankM) # number of genes

# for each source, get coordinates of non-missing values: nonMissPos (list)
nonMissPos=list()
for (k in 1:N) {
	nonMissPos[[k]] = which(!is.na(realRankM[,k]))
}

# For each data source (out of n sources), get permutation matrix
for (k in 1:N) {
	# first, make a matrix to store the data
	m = matrix(nrow=geneN,ncol=perN)
	for (j in 1:perN) {
		#m[nonMissPos[[k]],j] = 1:length(nonMissPos[[k]])
		m[nonMissPos[[k]],j] = sample( 1:length(nonMissPos[[k]]), length(nonMissPos[[k]]), replace=F)
	}
	# output the matrix
	write.table(m,file=paste('data_source_',k-1,'_random_rank.txt',sep=''),quote=F,sep="\\t",col.names=F,row.names=F)
}

EOF
	close OUT;
}


sub combinedProbability {
# input: probability array of individual evidence (S[i])
# output: combined probability = 1 - multiply(1 - S[i])
# assumption: the evidence is independent of each other
	shift @_;
	my @s=@_;

	# check range
	foreach (@s) {
		die "ERROR: unexpected value for combinedProbability(): $_\n Probability range of [0, 1] is expected.\n" unless (0<=$_ and $_<=1);
	}

	my $p=1; # variable for the final result

	foreach (@s) {
		$p*=1-$_;
	}

	return 1-$p;
}

sub reciprocal {
# add a new field (reciprocal) to an existing hash according to a specified field
	shift @_;
	my ($hash,$field,$newField)=@_;

	foreach my $k (keys %$hash) {
		next unless (defined($hash->{$k}->{$field}) and $hash->{$k}->{$field} ne 'NA' and $hash->{$k}->{$field} ne '');
		if ($hash->{$k}->{$field} != 0) {
			$hash->{$k}->{$newField} = 1/$hash->{$k}->{$field};
		} else {
			die "ERROR: 0 ranking: cannot genereate reciprocal for $k!\n(program exit)\n\n";
		}
	}
}

sub rank {
# add a new field (ranking) to an existing hash according to a specified field
# if there is a tie, the top most rank will be assigned to all members with that same value
	shift @_;
	my ($hash,$field,$de_in,$newField,$zeroBased, $abs)=@_;

	# build unique value hash
	# for details, see /home/yli4/development/JUMPrank/v1.0.0/tie_ranks/rank_test2.pl
	my @keys = grep {$hash->{$_}->{$field} ne ''} keys %$hash; # non blank values
	my %uniqValue; # %uniqValue{$value}=$key;

	# take absolute values before ranking?
	if (defined($abs) and $abs==1) {
		push @{$uniqValue{abs($hash->{$_}->{$field})}}, $_ for @keys;
	} else {
		push @{$uniqValue{$hash->{$_}->{$field}}}, $_ for @keys;
	}

	# get sorted values
	my @sortedValues;

	if ($de_in eq 'increase') {
		@sortedValues= sort {$a<=>$b} keys %uniqValue;
	} elsif ($de_in eq 'decrease') {
		@sortedValues= sort {$b<=>$a} keys %uniqValue;
	}

	# assign rank in hash
	my $k=1; # the rank
	$k=0 if (defined($zeroBased) and $zeroBased==1); # zero based ranking
	for (my $i=0; $i<scalar(@sortedValues); $i++) {
		my @genes=@{$uniqValue{$sortedValues[$i]}};  # gene name array

		foreach (@genes) {
			$hash->{$_}->{$newField} = $k;
		}

		$k+=scalar(@genes);
	}

	# for missing values, assign ''
	foreach (keys %$hash) {
		$hash->{$_}->{$newField} = '' if (!defined($hash->{$_}->{$newField}));
	}
}

sub orderStatistic_wrapper {
# input: rank ratio hash
# method: for each gene/row, calculate Q and assign as a new field
	my ($self,$rrhash)=@_;

	foreach my $g (keys %$rrhash) {
		my @r; # vector of rank ratio for this gene

		#for ($i=0; $i<$N; $i++) {
		foreach my $f (keys %{$rrhash->{$g}}) {
			#push @r $rrhash->{$}->{"rank$i"} if ($rrhash->{$}->{"rank$i"} ne '');
			push(@r,$rrhash->{$g}->{$f}) if ($f =~ /^rank(\d+)$/ and $rrhash->{$g}->{$f} ne '');
		}

		$rrhash->{$g}->{Q} = $self->orderStatistic_Q(@r);
	}
}

sub orderStatistic_Q {
# calculate order statistic Q based on Aerts et al., Nat. Biotech. 2006
# input: vector of n rank ratios
# output: Q
	shift @_;
	my @r=sort {$a <=> $b} @_; # vector of n rank ratios, sorted ascending

	my $n=scalar(@r); # dimension of data source
	return '' if ($n==0);

	# ensure 0<=rank ratio<=1
	foreach (@r) {
		die "ERROR: rank ratio out of bound: $_\n" unless (0<=$_ and $_<=1);
	}

	# calculate Vn
	my @V; 
	$V[0]=1;
	for (my $k=1; $k<=$n; $k++) {
		my $sign=1;
		my $ifact=1;
		$V[$k]=0; # Vk to be calculated; here initialized as 0
		for (my $i=1; $i<=$k; $i++) {
			$ifact=$ifact*$i;
			my $t=$sign * $r[$n-$k]**$i * $V[$k-$i] / $ifact;    # one term within the sum for Vk; @r range from 0 to $n-1
			$V[$k]+=$t; # sum for Vk
			$sign=$sign * (-1);
		}
	}

	# return Q=n!Vn
	return factorial($n)*$V[$n];
}

sub rankratio {
# based on %rankhash, return the rank ratio hash: %rrhash
	my ($self, $rankhash,$denominatorOption)=@_;

	my $rrhash;
	my $geneN=scalar(keys %$rankhash); # total gene number

	# get numbers of missing values for each source
	my ($naCount,$completeN)=$self->countNA($rankhash);

	# update rank based on the original rank (the original rank could be not accurate due to redundant IDs)
	foreach my $f (keys %$naCount) {
		# zero based ranking: this method puts too much weight on genes as long as it's top ranked in one source 
		#$self->rank($rankhash,$f,'increase',"new$f",1); # how about missing values? missing values will be still blanck

		# 1-based ranking
		$self->rank($rankhash,$f,'increase',"new$f"); # how about missing values? missing values will be still blanck
	}

	# calculate rank ratio
	foreach my $g (keys %$rankhash) {
		#foreach my $f (keys %{$rankhash->{$g}}) {
		foreach my $f (keys %$naCount) {

			if ($rankhash->{$g}->{$f} eq '') { # if n/a, assign blank
				$rrhash->{$g}->{$f}='';
				next;
			}

			if ($denominatorOption==1) { # option 1: varies across resource
				$rrhash->{$g}->{$f}=$rankhash->{$g}->{"new$f"}/($geneN-$naCount->{$f});
			} elsif ($denominatorOption==2) { # option 2: number of genes w/o missing values (Aerts et al., Nat. Biotech. 2006)
				#$rrhash->{$g}->{$f}=$rankhash->{$g}->{$f}/$completeN; # in this case, rank ratio > 1; will figure out later
			}
		}
	}

	return $rrhash;
}

sub countNA {
# 1) calculate number of rows w/ missing values for EACH column. Output: $countsCol->{$field}= # of blank in that column
# 2) count # rows w/o any n/a. Output: $completeCases = # of rows w/o any blank
# 3) for each ROW: count how many blanks. Output: $countsRow->{$key}= # of blank for that row
	shift @_;
	my ($hash)=@_;

	my ($countsCol,$completeCases,$countsRow);

	# initialize %countsCol: to store n/a # for each column
	foreach my $k (keys %$hash) {
		foreach my $f (keys %{$hash->{$k}}) {
			$countsCol->{$f}=0;
		}
		last;
	}

	$completeCases=0; # # rows w/o any n/a

	# start to count
	foreach my $k (keys %$hash) {
		my $k_na=0; # n/a # for this row
		$countsRow->{$k}=0;

		foreach my $f (keys %{$hash->{$k}}) {
			if ($hash->{$k}->{$f} eq '') {
				$countsCol->{$f}++;
				$k_na++;
				$countsRow->{$k}++;
			}
		}
		$completeCases++ unless ($k_na);
	}

	return $countsCol,$completeCases,$countsRow;
}

sub isGeneName {
# goal: check whether a string is gene name
	shift @_;
	my ($s)=@_;

	$s =~ s/\-//g; # we tolerate '-'

	return 0 if ($s eq '');

	my $ck=($s =~ /\W/)?0:1; # if match any non-"word" (word: letters, numbers, _) characters, it does not look like a normal gene name
	$ck=0 if ($s =~ /^\d+$/); # if only match to number, it  does not look like a normal gene name

	return $ck;
}

sub parseIdmappingDat {
# goal: parse HUMAN_9606_idmapping.dat, generate %hash{$uniprot}{$DB}{$id}=1
	shift @_;
	my ($inputFile)=@_;

	my $hash;

	open(IN,$inputFile) or die "cannot open $inputFile\n";

	while (<IN>) {
		chomp;
		my ($k,$db,$id)=split /\t/;
		$hash->{$k}->{$db}->{$id}=1;
	}

	close IN;

	return $hash;
}

=head
sub printDepth2Hash {
# goal: print out hash of IdmappingDat, which has depth of 2
	shift @_;
	my ($hash, $separator, $output)=@_;

	open(OUT,">",$output) or die "cannot open $output\n";

	# print header

	# print body
	foreach my 

	close OUT;
}
=cut

sub sjgraphIDconvert {
# goal: id convert for sjgraph format
# input: sjgraph file, %id1_to_id2
# output: sjgraph file with new IDs
	shift @_;
	my ($inputFile,$id1_to_id2,$outputFile,$field)=@_;

	open(IN,$inputFile) or die "cannot open $inputFile\n";
	open(OUT,">",$outputFile) or die "cannot open $outputFile\n";

	while (<IN>) {
		chomp;
		if (/^#/ or /^GeneA/ and /Function$/) { 
			print OUT "$_\n"; 
			next; 
		}

		my ($item_id_a,$item_id_b,$Direction,$score,$mode,$action)=split /\t/;

		if (defined($field)) {
			$item_id_a=(defined($id1_to_id2->{$item_id_a}))?$id1_to_id2->{$item_id_a}->{$field}:'n/a';
			$item_id_b=(defined($id1_to_id2->{$item_id_b}))?$id1_to_id2->{$item_id_b}->{$field}:'n/a';
		} else {
			$item_id_a=(defined($id1_to_id2->{$item_id_a}))?$id1_to_id2->{$item_id_a}:'n/a';
			$item_id_b=(defined($id1_to_id2->{$item_id_b}))?$id1_to_id2->{$item_id_b}:'n/a';
		}

		print OUT "$item_id_a\t$item_id_b\t$Direction\t$score\t$mode\t$action\n";
	}

	close IN;
	close OUT;
}

sub subhash2newhash {
# goal: we have a hash: %h1{$k}{$f1}{$id1/$id2/etc.} and also %h1{$k}{$f2}=$v
#       output: a new hash: %h2{$id1/$id2/etc.}=$v
	shift @_;
	my ($h1,$f1,$f2)=@_;

	my $h2; # a new hash

	foreach my $k (keys %$h1) {
		foreach my $id (keys %{$h1->{$k}->{$f1}}) {
			if (!defined($h2->{$id}) or $h2->{$id} eq 'n/a') {
				$h2->{$id}=$h1->{$k}->{$f2};
			}
		}
	}

	return $h2;
}

sub publicPPI2sjgraph {
# goal: convert public PPI DB format to sjgraph format
	shift @_;
	my ($inputFormat,$inputFile,$author,$fileversion,$outputFile)=@_;

	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);

	open(IN,$inputFile) or die "cannot open $inputFile\n";
	open(OUT,">",$outputFile) or die "cannot open $outputFile\n";

	print OUT "\# version: $fileversion\n";
	print OUT "\# date: ",$year+1900,"/",$mon+1,"/$mday\n";
	print OUT "\# author: $author\n";
	print OUT "\# source: $inputFile\n";

	if ($inputFormat eq 'bioplex') {
		print OUT "\# method: Immunopurification of the tagged protein and detection of associated proteins by mass spectrometry\n";
		print OUT "\# score: based on 'pInt' column (range 0.75-1)\n";
	}

	print OUT "GeneA\tGeneB\tDirection\tScore\tAction\tFunction\n";

	if ($inputFormat eq 'bioplex') {
		my $line=<IN>; # header
		while (<IN>) {
			chomp;
			chop if (/\r$/); # windows ending
			my ($GeneA,$GeneB,$UniprotA,$UniprotB,$SymbolA,$SymbolB,$pW,$pNI,$pInt)=split /\t/;
			print OUT "$SymbolA\t$SymbolB\t0\t$pInt\tbinding\t\n";
			#print OUT "$UniprotA\t$UniprotB\t0\t$pInt\tbinding\t\n";
		}
	} elsif ($inputFormat eq 'string') {
		my $line=<IN>; # header
		while (<IN>) {
			chomp;
			my ($item_id_a,$item_id_b,$mode,$action,$a_is_acting,$score)=split /\t/;
			my $Direction=0;
			$Direction=1 if ($mode eq 'activation' || $mode eq 'inhibition' || $mode eq 'ptmod' || $mode eq 'catalysis');
			#$item_id_a = (split /\./,$item_id_a)[1] if ($item_id_a =~ /\./);
			#$item_id_b = (split /\./,$item_id_b)[1] if ($item_id_b =~ /\./);
			print OUT "$item_id_a\t$item_id_b\t$Direction\t$score\t$mode\t$action\n";
		}
	} elsif ($inputFormat eq 'in_web') {
		while (<IN>) {
			chomp;
			my @t=split /\t/;

			# extract gene name
			my $geneA=(split /\|/,$t[4])[0];
			$geneA =~ /uniprotkb:([\w\-]+)\(gene name\)/;
			$geneA = $1;

			my $geneB=(split /\|/,$t[5])[0];
			$geneB =~ /uniprotkb:([\w\-]+)\(gene name\)/;
			#$geneB =~ /uniprotkb:(\w+)\(gene name\)/;
			$geneB = $1;

			# extact edge score
			my $score=(split /\|/,$t[14])[0];

			print OUT "$geneA\t$geneB\t0\t$score\tphysical_interaction\t\n";
		}
	} else {
		die "$inputFormat format is not supported!!\n";
	}

	close IN;
	close OUT;
}

sub parse_params
{
	shift @_;
        my ($par,$parahash)=@_;
        open(IN,$par);
        while (<IN>)
        {
                s/^\s+//;
                next if (/^#/);
                chomp;
                next unless (/ = /);

                s/\s*([;\#].*)$//; # delete comments
                my ($key,$value) = split(' = ',$_);

                $$parahash{$key}=$value;
        }
        close IN;
}

sub moduleVsGeneset {
# input: 1) Coexpressed_modules_annotated.txt (from -n); 2) gene list (e.g., an interesting pathway that we want to focus)
# purpose: check for each module (file 1), what genes are in the pathway (file 2)
# output:
#   add two columns in Coexpressed_modules_annotated.txt: 1) # of genes in the list; 2) genes in the list
#   genes in the list (file 2) but not in Coexpressed_modules_annotated.txt (file 1)
#   both information will be print out to Standard output
	my ($self,$moduleFile, $genelistFile)=@_;

        # parse inputs
        my @columnOrder;
        my $geneset=$self->parseKeyFile($genelistFile,1,-1,0);
        my $modules=$self->parseTab($moduleFile,0,\@columnOrder);

        # split gene list into subhash, with initial value of 0
        $self->splitValue2Hash($modules,'proteins','genes');

	# foreach module: overlap with %geneset
	# if match: set value as 1 for terms in %modules
	my $modulegenes;
        foreach my $m (keys %$modules) {
                $self->ifMatch($modules->{$m}->{genes},$geneset,'',1,0);
                foreach (keys %{$modules->{$m}->{genes}}) {
                        $modulegenes->{$_}=$modules->{$m}->{genes}->{$_};
                }
        }

        $self->ifMatch($geneset,$modulegenes,'',1,0);
        my $WB_proteins_remain=join(",",grep {$geneset->{$_}==0} keys %$geneset);

        # make a new field to store WB proteins for each module
        foreach my $m (keys %$modules) {
                my @t=grep {$modules->{$m}->{genes}->{$_}>0} keys %{$modules->{$m}->{genes}};
                $modules->{$m}->{'genes_in_list_n'}=scalar @t;
                $modules->{$m}->{'genes_in_list'}=join(",",@t);
        }
        push @columnOrder,'genes_in_list_n','genes_in_list';

        # print out %modules
        $self->printTab_order($modules,0,\@columnOrder,0);

        # print out WB proteins not in modules:
        print "genes not in modules\tn=",scalar(grep {$geneset->{$_}==0} keys %$geneset),"\t$WB_proteins_remain\n";

}

sub splitValue2Hash {
# purpose: suppose the %hash has a term ($hash->{$key}->{$columnToBeSplit}) that is a gene list, with gene names seperated by comma ($separator).
#          Here, we want to split this gene list and create a new subhash to store each gene name, with the following format:
#          $hash->{$key}->{$newKeyName}->{$geneName}=$initialValue
        shift @_;
        my ($hash,$columnToBeSplit,$newKeyName,$separator,$initialValue)=@_;

        $separator ||= ',';
        $initialValue ||= 0;

        foreach my $k (keys %$hash) {
                my @genes=split $separator,$hash->{$k}->{$columnToBeSplit};
                foreach my $g (@genes) {
                        $hash->{$k}->{$newKeyName}->{$g}=$initialValue;
                }
        }
}


sub hashKey2string {
# goal: take a hash as input, and concatenate all keys by $separator
# return: a string
	shift @_;
	my ($hash,$separator)=@_;
	$separator ||= ',';
	my @s;
	foreach my $k (keys %$hash) {
		push @s, $k;
	}
	return join($separator,@s);
	#return join(",",@s);
}

sub kappa {
# calculate  Cohen's Kappa according to https://en.wikipedia.org/wiki/Cohen%27s_kappa
#
# inputs: %h{$keys}=''
# 1) gene list 1
# 2) gene list 2
# 3) background
# output: kappa score
#
# procedures:
# 1) only focus on genes in the background (sub filterHash): %g1 -> %g1f, %g2 -> %g2f
# 2) calculate 4 numbers: 
#   $a=overlapped
#   $b=specific to %g1f
#   $c=specific to %g2f
#   $d=not in both lists
# 3) marginal probability:
#   $ma=(($a+$b)*($a+$c))/($a+$b+$c+$d)
#   $mb=(($c+$d)*($b+$d))/($a+$b+$c+$d)
# 4) he overall probability of random agreement is $pe = ($ma+$mb)/($a+$b+$c+$d)
# 5) The observed proportionate agreement is p0=($a+$d)/($a+$b+$c+$d)
# 6) applying our formula for Cohen's Kappa: $kp=($p0-$pe)/(1-$pe)
	my ($self, $g1, $g2, $bg)=@_;
	my ($g1f, $g2f, $kp, $a, $b, $c, $d, $p0, $pe);

	#1) only focus on genes in the background (sub filterHash): %g1 -> %g1f, %g2 -> %g2f
	$g1f=$self->filterHash($g1, $bg);
	$g2f=$self->filterHash($g2, $bg);

	# 2) calculate 4 numbers:
	my $ovlp=$self->filterHash($g1f, $g2f);
	$a=scalar(keys %$ovlp); #correct?
	$b=scalar(keys %$g1f)-$a;
	$c=scalar(keys %$g2f)-$a;
	$d=scalar(keys %$bg)-$a-$c-$b;

	# 3) marginal probability:
	my $ma=(($a+$b)*($a+$c))/($a+$b+$c+$d);
	my $mb=(($c+$d)*($b+$d))/($a+$b+$c+$d);

	# 4) he overall probability of random agreement
	$pe = ($ma+$mb)/($a+$b+$c+$d);

	# 5) The observed proportionate agreement
	$p0=($a+$d)/($a+$b+$c+$d);

	# 6) applying our formula for Cohen's Kappa
	$kp=($p0-$pe)/(1-$pe);

	return $kp;
}

sub rmRow_byPattern {
# filter row if a specific column match (or not match) a specific pattern
	shift @_;
	my ($h,$valueCol,$pattern,$match,$exact)=@_;

	foreach my $k (keys %$h) {
		my $m=0; # match or not
		if ($exact==1 and $h->{$k}->{$valueCol} eq $pattern 
		or $exact==0 and $h->{$k}->{$valueCol} =~ m/$pattern/) {
			$m=1;
		}

		if ($match and $m or !$match and !$m) {
			delete $h->{$k};
		}
	}
}

sub ifMatch {
# Purpose: 
# default usage: if $h2->{$k} is defined (note: $k is a key in %h1), then $h1->{$k}->{$newCol}=$trueValue; else $h1->{$k}->{$newCol}=$falseValue
# alernative usage: define $altKeyCol; if $h2->{$h1->{$k}->{$altKeyCol}} is defined (note: $k is a key in %h1), then $h1->{$k}->{$newCol}=$trueValue; else $h1->{$k}->{$newCol}=$falseValue
	shift @_;
	my ($h1,$h2,$newCol,$trueValue,$falseValue,$altKeyCol)=@_;

	foreach my $k (keys %$h1) {
		if (!defined($altKeyCol)) {
			if (defined($h2->{$k})) {
				$h1->{$k}->{$newCol}=$trueValue;
			}else{
				$h1->{$k}->{$newCol}=$falseValue;
			}
		} else {
			my $kk=$h1->{$k}->{$altKeyCol};
			if (defined($h2->{$kk})) {
				$h1->{$k}->{$newCol}=$trueValue;
			}else{
				$h1->{$k}->{$newCol}=$falseValue;
			}
		}
	}

	return $h1;
}

sub updateValues {
# Purpose: update values in target file ($tg) by those in source file ($sr), for specific columns (defined by $colIndex_s and $colIndex_t)
# Note: $colOrder_s and $colOrder_t are arrays of @column_order[$i]=$column_name, as defined by parseTab()
	shift @_;
	my ($sr,$colOrder_s,$colIndex_s,$tg,$colOrder_t,$colIndex_t)=@_;

	if (scalar(@{$colIndex_t}) != scalar(@{$colIndex_s})) {
		print "The target file and source file column index have different length!\n updateValues() is skipped!\n";
		return;
	}

	foreach my $id (keys %$tg) {
		# check whether this $id also found in $sr
		next unless (defined($sr->{$id}));

		# update values
		# Here, it's assumned that $colOrder_s and $colOrder_t are correct (i.e., match to $sr and $tg); but may not be always the case
		# May need more robust method
		for (my $i=0; $i<scalar(@{$colIndex_t}); $i++) {
			$tg->{$id}->{$colOrder_t->[$$colIndex_t[$i]]}=$sr->{$id}->{$colOrder_s->[$$colIndex_s[$i]]};
		}
		
	}

	return $tg;
}

sub capitalizeKey {
# take a hash input, and return a new hash with all keys capitalize
	shift @_;
	my ($hash)=@_;

	my $nwhash;

	foreach my $k (keys %$hash) {
		$nwhash->{uc($k)}=$hash->{$k};
	}

	return $nwhash;
}

sub filterInput {
# Purpose: filter a table (i.e., $input here) with a 3rd party file (i.e., $keyhash here)
# Function: read each row of $input: if the key term belong to (keys %keyhash) -> print this row 
# extra function: if $eliminate is defined and equals 1, it will print rows that does NOT match any of %keyhash
	shift @_;
	my ($input,$headerlines,$keyCol,$keyhash,$output,$cap,$eliminate)=@_;

	open(IN,$input);
	open(OUT,">",$output);

	for (my $i=0; $i<$headerlines; $i++) {
		my $line=<IN>;
		print OUT "$line";
	}

	while (<IN>) {
		chomp;
		chop if (/\r$/);
		my @t=split /\t/,$_,-1;
		if (defined($cap) and $cap==1) {
			$t[$keyCol]=uc($t[$keyCol]);
		}

		if (defined($keyhash->{$t[$keyCol]}) and (!defined($eliminate) or $eliminate==0)
		or !defined($keyhash->{$t[$keyCol]}) and defined($eliminate) and $eliminate==1 ) {
			print OUT "$_\n";
		}
	}

	close IN;
	close OUT;
}

sub condenseHash {
# old function: similar as twoLayer2single():
# 1) use a new key, e.g. $proID -> $GN
# 2) if select multiple isoforms: use $parameterKey to judge
	shift @_;
	my ($hash,$newkey,$parameterKey,$largerBetter,$oldkeyName)=@_;

	my $newhash;

	foreach my $k (keys %$hash) {
		my $newK=$hash->{$k}->{$newkey};
		if (!defined($newhash->{$newK}) 
		or $largerBetter and $newhash->{$newK}->{$parameterKey}<$hash->{$k}->{$parameterKey}
		or $largerBetter==0 and $newhash->{$newK}->{$parameterKey}>$hash->{$k}->{$parameterKey}) {
			$newhash->{$newK}->{$oldkeyName}=$k;
			foreach my $f (keys %{$hash->{$k}}){
				$newhash->{$newK}->{$f}=$hash->{$k}->{$f};
			}
		} 
	}

	return $newhash;
}

sub mergeHash {
# purpose: simply concatanate "rows" from two tables (hashes)
	shift @_;
	my ($hash1,$hash2)=@_;

	my $hash3;

	$hash3=clone $hash1;
	foreach (keys %$hash2) {
		$hash3->{$_}=clone $hash2->{$_};
	}

	return $hash3;
}

sub splitHash {
	shift @_;

	my ($largehash,$fghash)=@_;

	my ($fg,$bg);

	foreach (keys %$largehash) {
		if (defined($fghash->{$_})) {
			$fg->{$_}=clone $largehash->{$_};
		} else {
			$bg->{$_}=clone $largehash->{$_};
		}
	}

	return ($fg,$bg);
}

sub pValue {
	shift @_;
	my ($obs,@null)=@_;

	my $p;

	my $k=0; # number of times a null value is larger than the oberved
	foreach (@null) {
		if ($_>$obs) {
			$k++;
		}
	}
	$p=$k/scalar(@null);

	return $p;
}

sub array2hash {
	shift @_;
	my ($array)=@_;

	my $hash;

	foreach (@$array) {
		$hash->{$_}='';
	}

	return $hash;
}

sub array2hash_v2 {
	shift @_;
	my (@array)=@_;

	my $hash;

	foreach (@array) {
		$hash->{$_}='';
	}

	return $hash;
}

sub sample {
	shift @_;
	my ($n,@array)=@_;

	my @ra;

	# Shuffled list of indexes into @deck
	my @shuffled_indexes = shuffle(0..$#array);
	
	# Get just N of them.
	my @pick_indexes = @shuffled_indexes[ 0 .. $n - 1 ];  
	
	# Pick cards from @deck
	@ra = @array[ @pick_indexes ];

	return @ra;
}

sub filterCol {
# column filter: based on hash1, only keep columns with column names specified in hash2
	shift @_;
	my ($hash1,$hash2)=@_;

	my $hash3;
	foreach my $k (keys %$hash1) {
		foreach my $f (keys %$hash2) {
			if (defined($hash1->{$k}->{$f})) {
				$hash3->{$k}->{$f}=$hash1->{$k}->{$f};
			}
		}
	}

	return $hash3;
}

sub preorderedRows {
# goal: 
# methods: make a wrapper to achieve the goal
# step 1: parse the key file as $keyHash, and build a small hash using $smallHash=$self->filterHash($bigHash, $keyHash)
# step 2: parse the key file as $keyArray, and add a new index to $smallHash->{preorder}
# output: the resulting $smallHash, which can be printed out using $self->printTab_order by sorting $smallHash->{preorder}
	my ($self, $keyFile, $bigHash)=@_;

	my $keyHash=$self->parseKeyFile2($keyFile);
	my $smallHash=$self->filterHash($bigHash, $keyHash);

	my $keyArray=$self->parseKeyFileAsArray($keyFile);

	for (my $i=0; $i<scalar(@$keyArray); $i++) {
		my $k=$keyArray->[$i];
		$smallHash->{$k}->{preorder}=$i if (defined($smallHash->{$k}));
	}

	return $smallHash;
}

sub filterHash {
# row filter: based on hash1, only keep items with key stores in hash2
	shift @_;
	my ($hash1,$hash2,$reverse)=@_;
	#my ($hash1,$hash2)=@_;

	my ($hash3,$hash3_c); # hash3: keys in hash2; hash3_c: keys not in hash2
	$hash3=clone $hash1;
	$hash3_c=clone $hash1;
	foreach (keys %$hash3) {
		if (!defined($hash2->{$_})) {
			delete $hash3->{$_};
		} else {
			delete $hash3_c->{$_};
		}
	}

	if (defined($reverse) and $reverse) {
		return $hash3_c;
	} else {
		return $hash3;
	}

}

sub sortByScore {
	shift @_;
	my ($score)=@_;

	my $sortedGenes;

	foreach (sort {$score->{$b} <=> $score->{$a}} keys %$score) {
		push @$sortedGenes, $_;
	}

	return $sortedGenes;
}

sub enrichmentScoreFast {
	shift @_;
	# S, genes sorted by score
	my ($DEgenes,$sortedGenes)=@_;

	# step 0. configure parameters
	my $Nh=scalar(keys %$DEgenes);
	my $N=scalar(@$sortedGenes);
	my $increment=sqrt(($N-$Nh)/$Nh);
	my $decrement=sqrt($Nh/($N-$Nh));

	# step 3. starting form highest score: random walk => @KS
	my ($KS_current, @KS);
	$KS_current=0;
	push @KS, $KS_current;
	my ($inN,$deN);
	$inN=$deN=0;

	for (my $i=0;$i<scalar(@$sortedGenes);$i++) {
		if (defined($DEgenes->{$sortedGenes->[$i]})) {
			$KS_current+=$increment;
			$inN++;
		} else {
			$KS_current-=$decrement;
			$deN++;
		}
		push @KS, $KS_current;
	}

	# step 4. ES=max(abs(@KS))
	return (max(max(@KS),(-1)*min(@KS)),max(@KS),min(@KS),$Nh,$N,@KS);
}

sub enrichmentScore {
#	shift @_;
	# L, S, scores of whole data
	my ($self,$geneList,$DEgenes,$scoresAllGenes)=@_;

	# step 0. configure parameters
	# N=length(L)
	# Nh=length(S)
	# increase=sqrt((N-Nh)/Nh)
	# decrease=sqrt(Nh/(N-Nh))

	# step 1. intersect: only keep scores with gene names in $geneList 
	my $scoresL=$self->filterHash($scoresAllGenes,$geneList);
	my $DEgenesL=$self->filterHash($DEgenes,$scoresL);

	my $Nh=scalar(keys %$DEgenesL);
	my $N=scalar(keys %$scoresL);
	my $increment=sqrt(($N-$Nh)/$Nh);
	my $decrement=sqrt($Nh/($N-$Nh));

	# step 2. sort $scoresL by score
	# step 3. starting form highest score: random walk => @KS
	my ($KS_current, @KS);
	$KS_current=0;
	push @KS, $KS_current;
	my ($inN,$deN);
	$inN=$deN=0;
	foreach (sort {$scoresL->{$b} <=> $scoresL->{$a}} keys %$scoresL) {
		if (defined($DEgenes->{$_})) { # foreground
			$KS_current+=$increment;
			$inN++;
		} else { # background
			$KS_current-=$decrement;
			$deN++;
		}
		#push @KS, abs($KS_current);
		push @KS, $KS_current;
	#	print "$KS_current,";
	}
	#print "$Nh,",$N-$Nh,",$inN,$deN,",$KS_current,",",max(@KS),"\n";

	# step 4. ES=max(abs(@KS))
	return (max(max(@KS),(-1)*min(@KS)),max(@KS),min(@KS),$Nh,$N,@KS);
	#return (max(abs(@KS)),$Nh,$N,@KS);
}

sub print2layerTab {
# Goal: print out one file for each key value on the 1st layer of the hash
	#shift @_;
	my ($self,$hash,$outputSuffix)=@_;

	# print out one file for each key value on the 1st layer of the hash
	foreach my $k (keys %$hash) {
		my $output="$k$outputSuffix";
		$self->printTab($hash->{$k},$output);
	}
}

sub parse2layerTab {
	shift @_;
	my ($input,$colkey1, $colkey2)=@_;
	#$colkey1 ||= 0;
	#$colkey2 ||= 1; # if $colkey2 is set as 0, this statement will wrongly re-assign it as 1, because perl interprets '0' as 'false' hence perfer '1' (true) in the '||' operation
	$colkey1=(!defined($colkey1))?0:$colkey1;
	$colkey2=(!defined($colkey2))?1:$colkey2;

	my $hash;

	open(IN,$input) || die "cannot open file: $input\n";
	my $header=<IN>; chomp($header);
	#chop($header) if ($header =~ m/\r$/); # windows ending
	$header =~ s/\r$//g; # windows ending
	#$header =~ s/\#//g; # no '#'
	my @colnames=split /\t/,$header;
	while (<IN>) {
		chomp;
		#chop if (/\r$/); # windows ending
		s/\r$//g; # windows ending
		my @t=split /\t/,$_;

		next if ($t[$colkey1] eq $colnames[$colkey1]); # advoid repeated headers (e.g., using >1 annotation DBs)

		for (my $i=0; $i<=$#t; $i++) {
			#next if ($i==$colkey1 or $i==$colkey2);
			$hash->{$t[$colkey1]}->{$t[$colkey2]}->{$colnames[$i]}=$t[$i];
		}
	}
	close IN;

	return $hash;
}

sub parse3layerTab {
	shift @_;
	my ($input,$colkey1, $colkey2, $colkey3)=@_;

	my $hash;

	open(IN,$input) || die "cannot open file: $input\n";
	my $header=<IN>; chomp($header);
	#chop($header) if ($header =~ m/\r$/); # windows ending
	$header =~ s/\r$//g; # windows ending
	#$header =~ s/\#//g; # no '#'
	my @colnames=split /\t/,$header;
	while (<IN>) {
		chomp;
		#chop if (/\r$/); # windows ending
		s/\r$//g; # windows ending
		my @t=split /\t/,$_;

		next if ($t[$colkey1] eq $colnames[$colkey1]); # advoid repeated headers (e.g., using >1 annotation DBs)

		for (my $i=0; $i<=$#t; $i++) {
			#next if ($i==$colkey1 or $i==$colkey2);
			$hash->{$t[$colkey1]}->{$t[$colkey2]}->{$t[$colkey3]}->{$colnames[$i]}=$t[$i];
		}
	}
	close IN;

	return $hash;
}

sub twoLayer2single {
# compress a 2-layer hash to a 1-layer hash: 
# 1) will use the 1st layer key as new key, hence the 2nd layer key will be transfered to a normal column;
# 2) accordingly, only one representative item of the 2nd layer will be selected;
# 3) the selection criteria is determined by the 3 following parameters:
# 	a) $columnName: a column with values that determine the selection outcome
# 	b) $chooseLarge: 1=prefer large value; 0=prefer small value
# 	c) $abs: take an absolute value before selection

	shift @_;
	my ($h1, $columnName, $chooseLarge, $abs)=@_;
	my $h2; # 1-layer hash to be returned

	foreach my $k1 (keys %$h1) {
		my $s_k2=''; # record the selected key of 2nd layer
		my $s_value=''; # record the best value

		foreach my $k2 (keys %{$h1->{$k1}}) {
			my $value=(defined($abs) and $abs)?abs($h1->{$k1}->{$k2}->{$columnName}):$h1->{$k1}->{$k2}->{$columnName};

			if ($s_value eq '') { # first value
				$s_value=$value; $s_k2=$k2;
			} elsif ($chooseLarge==1 and $s_value<$value) { # prefer larger and is large
				$s_value=$value; $s_k2=$k2;
			} elsif ($chooseLarge==0 and $s_value>$value) { # prefer small and is small
				$s_value=$value; $s_k2=$k2;
			}
		}

		# let %h2 record the new item
		foreach my $f (keys %{$h1->{$k1}->{$s_k2}}) {
			$h2->{$k1}->{$f}=$h1->{$k1}->{$s_k2}->{$f};
		}
	}

	return $h2;
}

sub column2hash {
# convert a particular column content to a new hash
	shift @_;
	my ($hash,$col,$col2, $col2_value)=@_;
	my $colhash;
	$col2 ||= '';
	$col2_value ||= '';

	foreach my $h (values %$hash) {
		if (defined($h->{$col})) {
			if ($col2 ne '') {
				if (defined($h->{$col2}) and $h->{$col2} eq $col2_value) {
					$colhash->{$h->{$col}}=1;
				}
			} else {
				$colhash->{$h->{$col}}=1;
			}
		}
	}

	return $colhash;
}

sub attachCol {
# Default (both hashes shared keys): add %h2->{$key}->{$newKey} to %h1->{$key}, otherwise 'NA'
# Alternatively, specify which column in %h1 that correspond to keys of %h2
	shift @_;
	my ($h1,$h2,$newKey,$connectKey_h1,$altNewKey)=@_;
	$altNewKey ||= $newKey;

	foreach my $k(keys %$h1) {
		my $key=''; # the key to connect to %h2
		if (defined($connectKey_h1) and $connectKey_h1 ne '0') {
			$key=$h1->{$k}->{$connectKey_h1};
		} else {
			$key=$k;
		}

		if (defined($h2->{$key}->{$newKey})) {
			$h1->{$k}->{$altNewKey}=$h2->{$key}->{$newKey};
		} elsif (!defined($h1->{$k}->{$altNewKey})) {
			$h1->{$k}->{$altNewKey}='NA';
		}
	}
}

sub phoNorm {
	shift @_;
	my ($whl,$pho,$proKey,$tmtCols)=@_;

	#foreach my $hash (values %$pho) {
	foreach my $p (keys %$pho) {
		my $hash=$pho->{$p};
		my $pro=$hash->{$proKey};
		if (defined($whl->{$pro})) {
			foreach my $c (keys %$tmtCols) {
				if (defined($hash->{$c}) and defined($whl->{$pro}->{$c})) {
					$hash->{$c}=$hash->{$c}-$whl->{$pro}->{$c}; # norm. ratio = pho - whl
				}else {
					die "TMT chanennels $c undefined for $pro\n";
				}
			}
		} else {
			#print "$pro not detected in whl proteome\n";
			#delete $hash; # delete pho items w/o whl
			#undef $hash; # delete pho items w/o whl
			delete $pho->{$p};
		}
	}
}

sub applyExponential {
	shift @_;
	my ($hash, $numKeys, $logBase)=@_;
	$logBase ||= 2;

	foreach my $k (keys %$hash) {
		foreach my $f (keys %$numKeys) {
			$hash->{$k}->{$f}= $logBase ** $hash->{$k}->{$f};
		}
	}
}

sub applyScale {
	shift @_;
	my ($hash, $numKeys, $scale)=@_;

	foreach my $k (keys %$hash) {
		foreach my $f (keys %$numKeys) {
			$hash->{$k}->{$f}= $hash->{$k}->{$f}*$hash->{$k}->{$scale};
		}
	}
}

sub applyLogRatio {
	shift @_;
	my ($hash, $numKeys, $ctrlKey, $logBase)=@_;
	$logBase ||= 2;

	foreach my $k (keys %$hash) {
		foreach my $f (keys %{$hash->{$k}}) {
			if (defined($numKeys->{$f}) and $hash->{$k}->{$ctrlKey} !=0) {
				$hash->{$k}->{$f}=log($hash->{$k}->{$f}/$hash->{$k}->{$ctrlKey})/log($logBase);
			}
		}
	}
}

sub applyMean {
	shift @_;
	my ($hash, $numKeys, $newKey)=@_;

	foreach my $k (keys %$hash) {
		my @t;
		foreach my $f (keys %{$hash->{$k}}) {
			if (defined($numKeys->{$f})) {
				push @t,$hash->{$k}->{$f};
			}
		}
		$hash->{$k}->{$newKey}=mean(@t);
	}
}

sub applySum {
	shift @_;
	my ($hash, $numKeys, $newKey)=@_;

	foreach my $k (keys %$hash) {
		my @t;
		my $sum=0;
		foreach my $f (keys %{$hash->{$k}}) {
			if (defined($numKeys->{$f})) {
				#push @t,$hash->{$k}->{$f};
				$sum+=$hash->{$k}->{$f};;
			}
		}
		$hash->{$k}->{$newKey}=$sum;
	}
}

sub applyMax {
	shift @_;
	my ($hash, $numKeys, $newKey)=@_;

	foreach my $k (keys %$hash) {
		my @t;
		foreach my $f (keys %{$hash->{$k}}) {
			if (defined($numKeys->{$f})) {
				push @t,$hash->{$k}->{$f};
			}
		}
		$hash->{$k}->{$newKey}=max(@t);
	}
}

sub applyMin {
	shift @_;
	my ($hash, $numKeys, $newKey)=@_;

	foreach my $k (keys %$hash) {
		my @t;
		foreach my $f (keys %{$hash->{$k}}) {
			if (defined($numKeys->{$f})) {
				push @t,$hash->{$k}->{$f};
			}
		}
		$hash->{$k}->{$newKey}=min(@t);
	}
}

sub addColnamePrefix {
	shift @_;
	my ($hash,$prefix,$colOrder)=@_;

	foreach my $k (keys %$hash) {
		foreach my $f (keys %{$hash->{$k}}) {
			my $new_f="$prefix$f";
			$hash->{$k}->{$new_f}=$hash->{$k}->{$f};
			delete $hash->{$k}->{$f};
		}
	}

	# also update @colOrder if defined
	if (defined($colOrder)) {
		for (my $i=0; $i<scalar(@{$colOrder}); $i++) {
			$$colOrder[$i]="$prefix$$colOrder[$i]";
		}
	}
}

sub parseMultiTab {
# a wrapper to merge columns of multiple tables, with assumptions below:
# 1. all the keys is in column '0'
# 2. all column names are unique across ALL files; redundant columns (i.e. with same column names) will be randomly collapsed

	my ($self,$inlist)=@_;
	my @inputs=@$inlist;

	my $mergedhash;

	my $input=shift @inputs;
	my $mergedhash=$self->parseTab($input);

	foreach $input (@inputs) {
		my $h1=$mergedhash;
		my $h2=$self->parseTab($input);
		$mergedhash=$self->mergeTabHash($h1,$h2);
	}

	return $mergedhash;
}

sub parseMultiTab_order {
	my ($self,$inlist,$colOrder)=@_;
	my @inputs=@$inlist;

	my $mergedhash;

	my $input=shift @inputs;
	my $mergedhash=$self->parseTab($input,0,$colOrder);

	foreach $input (@inputs) {
		my $h1=$mergedhash;

		my @colOrder2;
		my $h2=$self->parseTab($input,0,\@colOrder2);

		$mergedhash=$self->mergeTabHash($h1,$h2);
		for (my $i=0; $i<=$#colOrder2; $i++) {
			push @$colOrder, $colOrder2[$i];
		}
	}

	return $mergedhash;
}

sub mergeTabHash {
# to merge columns of two tables
	shift @_;
	my ($hash1, $hash2)=@_;
	my $hash3;

	# clone hash1 to hash3; also get field names
	my %fields;
	foreach my $k (keys %$hash1) {
		foreach my $f (keys %{$hash1->{$k}}) {
			$hash3->{$k}->{$f}=$hash1->{$k}->{$f};
			$fields{$f}='';
		}
	}

	# add $hash2 to $hash3; for hash 1 unique rows, assigned ''
	foreach my $k (keys %$hash2) {
		# if not in hash1: initialize field from hash1
		if (!defined($hash3->{$k})) {
			foreach my $f (keys %fields) {
				$hash3->{$k}->{$f}='';
			}
		}
		# add hash2 fields
		foreach my $f (keys %{$hash2->{$k}}) {
			$hash3->{$k}->{$f}=$hash2->{$k}->{$f};
			$fields{$f}='';
		}
	}

	# add blank hash2 fields to hash1 specific rows of hash3
	foreach my $k (keys %$hash1) {
		foreach my $f (keys %fields) {
			if (!defined($hash3->{$k}->{$f})) {
				$hash3->{$k}->{$f}='';
			}
		}
	}

	return $hash3;
}

sub parseMultiTab_specificCol {
# a wrapper to merge specific columns of multiple tables, with assumptions below:
# 1. all the keys is in column '0'
# 2. every file has one column that exactly match $colName, and this column should be in column '1'
# 3. if non-unique ID exist, only one record will be retained, with the criteria defined by two paramaters: $chooseLarge,$abs (see details of these two parameters in sub twoLayer2single())

	my ($self,$inlist,$colName,$chooseLarge,$abs)=@_;
	my @inputs=@$inlist;

	my $mergedhash;

	my $input=shift @inputs;
	my $mergedhash_0 = $self->parse2layerTab($input,0,1);
	my $mergedhash = $self->twoLayer2single($mergedhash_0, $colName, $chooseLarge, $abs);

	my $index=0; # counting the table; used for naming new columns: "$colName$index"
	foreach $input (@inputs) {
		my $h1=$mergedhash;
		#my $h2=$self->parseTab($input);
		my $h2_0 = $self->parse2layerTab($input,0,1);
		my $h2 = $self->twoLayer2single($h2_0, $colName, $chooseLarge, $abs);

		$index++;
		$mergedhash=$self->mergeTabHash_specificCol($h1,$h2,$colName,$index);
	}

	return $mergedhash;
}

sub mergeTabHash_specificCol {
# to merge specific columns of two tables: either exactly match $colName or "$colName(\d+)"
	shift @_;
	my ($hash1, $hash2, $colName, $index)=@_;
	my $hash3;

	# clone hash1 to hash3; also get field names
	# note: hash 1 may be a novel hash, or it's already a combined hash from previous merging steps
	my %fields;
	foreach my $k (keys %$hash1) {
		foreach my $f (keys %{$hash1->{$k}}) {
			if ($f eq $colName) { # hash 1 is a novel hash
				my $ff=$f."0";
				$hash3->{$k}->{$ff}=$hash1->{$k}->{$f};
				$fields{$ff}='';
			} elsif ($f =~ /^$colName(\d+)$/) { # hash 1 is already a combined hash
				$hash3->{$k}->{$f}=$hash1->{$k}->{$f};
				$fields{$f}='';
			}
		}
	}

	# add $hash2 to $hash3; for hash 1 unique columns, assigned ''
	foreach my $k (keys %$hash2) {
		# if not in hash1 (hash 2 specific rows): initialize field from hash1
		if (!defined($hash3->{$k})) {
			foreach my $f (keys %fields) {
				$hash3->{$k}->{$f}='';
			}
		}
		# add hash2 fields
		foreach my $f (keys %{$hash2->{$k}}) {
			if ($f eq $colName) {
				$hash3->{$k}->{"$f$index"}=$hash2->{$k}->{$f};
				$fields{"$f$index"}='';
			}
		}
	}

	# add blank hash2 fields to hash1 columns of hash3
	foreach my $k (keys %$hash1) {
		foreach my $f (keys %fields) {
			if (!defined($hash3->{$k}->{$f})) {
				$hash3->{$k}->{$f}='';
			}
		}
	}

	return $hash3;
}

sub printIDtxt {
        shift @_;
        my ($scanhash,$output)=@_;

        open(OUT,">$output");
        foreach my $outfile (keys %$scanhash) {
                print OUT $scanhash->{$outfile}->{IDtxt_line},"\n";
        }
        close OUT;
}

sub printGeneList {
	shift @_;
	my ($hash, $output, $printN)=@_;

	open(OUT,">$output");

	foreach my $list (keys %$hash) {
		print OUT $list,"\t";
		foreach my $gene (keys %{$hash->{$list}}) {
			print OUT "$gene,";
		}
		if (defined($printN) and $printN) {
			print OUT "\t",scalar(keys %{$hash->{$list}});
		}
		print OUT "\n";
	}

	close OUT;
}

sub print2layerTab_to_matrix {
# goal: given a 2-layer hash ( $h->{$k1}->{$k2}->{$f} ) and a field ($valueCol) that repreasent the content/value of the matrix, this subroutine will generate a matrix:
# nrow = number of keys of layer 1
# ncol = number of keys of layer 2 (the dimension is define by scalar(@{$colnames}))
	shift @_;
	my ($h, $f, $colnames, $output)=@_;

	open(OUT,">$output");

	# print colnames
	for (@{$colnames}) {
		print OUT "\t$_";
	}
	print OUT "\n";

	# then print each row
	foreach my $k1 (keys %$h) {
		print OUT "$k1";

		for my $k2 (@{$colnames}) {
			print OUT "\t";
			if (defined($h->{$k1}->{$k2}->{$f})) {
				print OUT $h->{$k1}->{$k2}->{$f};
			}
		}

		print OUT "\n";
	}

	close OUT;
}

sub print2layerGeneList {
# goal (just guess; may doing totally different things): designed for PPI module, with the following format: moduleName, moduleSize, pro1,pro2,...,proN
	shift @_;
	my ($hash, $output)=@_;

	open(OUT,">$output");

	# print header
	print OUT "module";
	foreach my $list (keys %$hash) {
		foreach my $cluster (sort keys %{$hash->{$list}}) {
			print OUT "\t$cluster";
		}
		foreach my $cluster (sort keys %{$hash->{$list}}) {
			print OUT "\t$cluster";
		}
		foreach my $cluster (sort keys %{$hash->{$list}}) {
			print OUT "\t$cluster";
		}
		last;
	}
	print OUT "\n";

	foreach my $list (sort keys %$hash) {
		print OUT $list;
		# protein number of each cluster
		my $n=0;
		foreach my $cluster (sort keys %{$hash->{$list}}) {
			my $wpcn=scalar(keys %{$hash->{$list}->{$cluster}});
			$n+=$wpcn;
			print OUT "\t",$wpcn;
		}
		# protein % of each cluster
		#print OUT "\t";
		foreach my $cluster (sort keys %{$hash->{$list}}) {
			printf OUT "\t%.2f",scalar(keys %{$hash->{$list}->{$cluster}})*100/$n;
		}
		# print protein names
		foreach my $cluster (sort keys %{$hash->{$list}}) {
			print OUT "\t";
			foreach my $g (keys %{$hash->{$list}->{$cluster}}) {
				print OUT "$g,";
			}
		}
		print OUT "\n";
	}

	close OUT;
}

sub geneList2KeyHash {
# convert GeneList hash (2 layered) to a single layer hash:
# key: gene names
# remove list names
	shift @_;
	my ($hash)=@_;

	my $keyhash;
	foreach my $h (values %$hash) {
		foreach my $g (keys %$h) {
			$keyhash->{$g}='';
		}
	}

	return $keyhash;
}

sub FisherExactTest_Rscript {
# generate R script for FET
	shift @_;
	my ($output)=@_;
	open(OUT, ">$output");
	print OUT <<"EOF";
for (e in commandArgs(T)) {
  ta = strsplit(e,"=",fixed=TRUE)
  var = ta[[1]][1]
  if(! is.na(ta[[1]][2])) {
    temp = ta[[1]][2]
    var = substr(ta[[1]][1],2,nchar(ta[[1]][1]))
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
      temp = as.integer(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
      temp = as.numeric(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "V") {
      temp = strsplit(temp,',')[[1]]
    }
    assign(var,temp)
    cat("assigned ",var," the value of |",temp,"|\\n")
  } else {
    var_fields = strsplit(var,'-')[[1]]
    var = var_fields[length(var_fields)]
    assign(var,TRUE)
    cat("assigned ",var," the value of TRUE\\n")
  }
}
# load input file (output from combineTables.pl)
tb=read.table(inputFile,head=T,sep="\\t")
summary(tb)
head(tb)
dim(tb)

# fold enrichment
tb\$foldEnrich=(tb\$overlapped.n/tb\$list2.n)/(tb\$list1.n/tb\$total.n)

# perform FET row by row
tb\$p=rep(1,nrow(tb))
for (i in 1:nrow(tb)) {
        m=matrix(nrow=2,ncol=2)
        m[1,1]=tb\$overlapped.n[i]
        m[1,2]=tb\$list2.n[i]-m[1,1]
        m[2,1]=tb\$list1.n[i]-m[1,1]
        m[2,2]=tb\$total.n[i]-m[1,1]-m[1,2]-m[2,1]
        tb\$p[i]=fisher.test(m)\$p.value
}

# BH FDR
tb\$fdr=p.adjust(tb\$p,method='BH')
write.table(tb,file=outputFile,quote=F,sep="\\t",row.name=F)

EOF
	close OUT;
}

sub FisherInputTable {
# prepare Fisher's Exact Test input matrix: list1 name, list2 name, list1 n, list2 n, overlapped n, total n
	my ($self, $genelist1, $genelist2, $totalN, $output)=@_;

	open(OUT, ">$output");

	print OUT "list1 name	list2 name	list1 n	list2 n	overlapped n	total n\n";
	
	foreach my $list1 (keys %$genelist1) {
		my $n1=scalar(keys %{$genelist1->{$list1}});

		foreach my $list2 (keys %$genelist2) {
			my $n2=scalar(keys %{$genelist2->{$list2}});
			my $cmn=$self->filterHash($genelist1->{$list1}, $genelist2->{$list2});
			my $ovlpN=scalar(keys %{$cmn});

			print OUT "$list1\t$list2\t$n1\t$n2\t$ovlpN\t$totalN\n";
		}
	}

	close OUT;
}

sub parseGeneList {
# convert gene list format (i.e., input for Fisher exact test) to hash
	shift @_;
	my ($input,$keyCol, $genelistCol, $withHeader)=@_;

	$keyCol ||= 0;
	$genelistCol ||= 1;
	$withHeader ||= 0;

	my $hash;

	open(IN,$input)  or die "Cannot open $input!!!\n";

	if ($withHeader) {
		my $line=<IN>;
	}

	while (<IN>) {
		chomp;
		#my ($listname,$genelist)=split(/\s/,$_);
		my @t=split(/\s/,$_);
		my ($listname,$genelist)=($t[$keyCol],$t[$genelistCol]);
		my @genes=split /,/,$genelist;
		foreach my $gene (@genes) {
			$hash->{$listname}->{$gene}='';
		}
	}

	close IN;

	return $hash;
}

sub filterDtas {
	shift @_;
	my ($hash,$dtas,$output)=@_;

	open(IN,$dtas)  or die "Cannot open $dtas!!!\n";
	open(OUT,">",$output);

	while(<IN>) {
		my $headerline=$_;
		my @t=split(/\s/,$_);
		my $line;
		#$line=<IN>; my @mz=split(/\s/,$line);
		#$line=<IN>; my @int=split(/\s/,$line);
		my $mzline=<IN>;
		my $intline=<IN>;

		# t20.10000.1.2.dta
		print "  parsing $t[0]\r";
		my $outfile=$t[0];
		$outfile =~ s/\.dta//;

		if (defined($hash->{$t[0]}) || defined($hash->{"$outfile.out"}) || defined($hash->{"$outfile.spout"}) || defined($hash->{$outfile}) ) {
			#print OUT "$t[1] $t[2]\n";
			#for (my $i=0; $i<scalar(@int); $i++) { print OUT "$mz[$i] $int[$i]\n"; }			
			print OUT "$headerline$mzline$intline";
		}
	}
	print "Finish parsing .dtas\n";

	close IN;
	close OUT;
}

sub parseKeyFile2 {
	shift @_;
	my ($input,$header)=@_;
	$header ||= 0; # exist header?
	my $hash;

	open(IN,$input) || die "Cannot open $input\n";
	if ($header) {
		my $line=<IN>;
	}
	while (<IN>) {
		chomp;
		$hash->{$_}='';
	}
	close IN;

	return $hash;
}

sub parseKeyFile {
	shift @_;
	my ($input,$header,$valueCol)=@_;
	$header ||= 0; # exist header?
	$valueCol ||= -1; # exist value column? which one?
	my $hash;

	open(IN,$input) || die "Cannot open $input\n";
	if ($header) {
		my $line=<IN>;
	}
	while (<IN>) {
		chomp;
		my @t=split(/\s/,$_);
		if ($valueCol>=0) {
			$hash->{$t[0]}=$t[$valueCol];
		} else {
			$hash->{$t[0]}='';
		}
	}
	close IN;

	return $hash;
}

sub printSimpleHash {
# print hash, w/ structure: $h->{$key}=$value
	shift @_;
	my ($hash,$output,$keyOnly)=@_;
	open(OUT,">",$output);

	foreach my $k (keys %$hash) {
		if (defined($keyOnly) and $keyOnly) {
			print OUT "$k\n";
		} else {
			print OUT "$k\t$hash->{$k}\n";
		}
	}

	close OUT;
}

sub printTab {
# 3 printing mode of sorting/printing columns:
# 1) default: sort by column name ($pattn and $column_order are blank)
# 2) only print columns with column name that match to $pattn
# Future development: 2) and 3) can be combined, but not available yet
	shift @_;
	my ($hash,$output,$pattn)=@_;
	open(OUT,">",$output);

	# print header
	#my $sampleK=shift (keys %$hash);
	my $sampleK=(keys %$hash)[0];
	print OUT "key";
	#foreach my $colname (keys %{$hash->{$sampleK}}) {
	foreach my $colname (sort {$a cmp $b} keys %{$hash->{$sampleK}}) {
		if (!defined($pattn) or $pattn eq '' or $colname =~ m/$pattn/) {
			print OUT "\t",$colname;
		}
	}
	print OUT "\n";

	# print table content
	foreach my $k (keys %$hash) {
		print OUT $k;
#		foreach my $value (values %{$hash->{$k}}) {
#			print OUT "\t",$value;
#		}
		foreach my $colname (sort {$a cmp $b} keys %{$hash->{$sampleK}}) {
			if (!defined($pattn) or $pattn or $colname =~ m/$pattn/) {
				print OUT "\t",$hash->{$k}->{$colname};
			}
		}
		print OUT "\n";
	}
	close OUT;
}

sub printTab_order {
# 3) keep original column orders (if $column_order is defined)
	shift @_;
	my ($hash,$output,$column_order,$noHeader,$continueWrite,$printKey, $sortCol, $increasing, $NA_label, $nonNumericSortCol)=@_;
	$NA_label ||= '';

	if (defined($continueWrite) and $continueWrite==1) {
		open(OUT,">>",$output);
	} else {
		open(OUT,">",$output);
	}

	# print header
	unless (defined($noHeader) and $noHeader==1) {
		# print key column
		if (defined($printKey) and $printKey==1) {
			print OUT "key\t";
		}

		print OUT $$column_order[0];
		for (my $i=1; $i<scalar(@{$column_order}); $i++) {
			print OUT "\t",$$column_order[$i];
		}
		print OUT "\n";
	}

	# print table content
	# decide how to sort the order of rows (if specified in the parameters)
	my @keys = keys %$hash;
	if (defined($sortCol)) {
		if (defined($increasing) and $increasing) {
			if (defined($nonNumericSortCol) and $nonNumericSortCol) {
				@keys = sort {$hash->{$a}->{$sortCol} cmp $hash->{$b}->{$sortCol}} @keys;
			} else {
				@keys = sort {$hash->{$a}->{$sortCol} <=> $hash->{$b}->{$sortCol}} @keys;
			}
		} else {
			if (defined($nonNumericSortCol) and $nonNumericSortCol) {
				@keys = sort {$hash->{$b}->{$sortCol} cmp $hash->{$a}->{$sortCol}} @keys;
			} else {
				@keys = sort {$hash->{$b}->{$sortCol} <=> $hash->{$a}->{$sortCol}} @keys;
			}
		}
	}

	# then print each row (with row orders defined above)
	foreach my $k (@keys) {

		# then print other columns
		my @tmp;
		my $naCount=0;
		#print OUT $hash->{$k}->{$$column_order[0]};
		for (my $i=0; $i<scalar(@{$column_order}); $i++) {
			push @tmp, $hash->{$k}->{$$column_order[$i]};

			$naCount++ if ($tmp[$#tmp] eq '');

			# how to print blank/missing data?
			if ($NA_label ne '' and  $tmp[$#tmp] eq '') {
				$tmp[$#tmp]=$NA_label;
			}
		}

		if ($naCount<scalar(@{$column_order})) { # if not a totally empty row, let's print it
			# print key column if required
			if (defined($printKey) and $printKey==1) {
				print OUT "$k\t";
			}

			# print other required fields
			print OUT join("\t",@tmp);
			print OUT "\n";
		}
	}

	close OUT;

}

sub get_colnames {
# extract the colnames of a hash (that correponds to a dataframe)
	shift @_;
	my ($h) = @_;

	my @colnames;
	my $k=(keys %$h)[0]; # 1st element
	for (keys %{$h->{$k}}) {
		push @colnames, $_;
	}

	return sort(@colnames);
}

sub printTab_hashParam {
# 3) keep original column orders (if $column_order is defined)
	shift @_;
	my ($params) = @_;

	my $hash = $params->{hash};
	my $output = $params->{output};
	my $column_order = $params->{column_order};
	my $noHeader = defined($params->{noHeader})?$params->{noHeader}:0;
	my $continueWrite = defined($params->{continueWrite})?$params->{continueWrite}:0;
	my $printKey = defined($params->{printKey})?$params->{printKey}:0;
	my $sortCol = defined($params->{sortCol})?$params->{sortCol}:'';
	my $increasing = defined($params->{increasing})?$params->{increasing}:0;
	my $NA_label = defined($params->{NA_label})?$params->{NA_label}:'';
	my $nonNumericSortCol = defined($params->{nonNumericSortCol})?$params->{nonNumericSortCol}:0;

	#my ($hash,$output,$column_order,$noHeader,$continueWrite,$printKey, $sortCol, $increasing, $NA_label, $nonNumericSortCol)=@_;
	#$NA_label ||= '';

	if (defined($continueWrite) and $continueWrite==1) {
		open(OUT,">>",$output);
	} else {
		open(OUT,">",$output);
	}
	#print "writing output to $output ...\n";

	# print header
	unless (defined($noHeader) and $noHeader==1) {
		# print key column
		if (defined($printKey) and $printKey==1) {
			print OUT "key\t";
		}

		print OUT $$column_order[0];
		for (my $i=1; $i<scalar(@{$column_order}); $i++) {
			print OUT "\t",$$column_order[$i];
		}
		print OUT "\n";
	}

	# print table content
	# decide how to sort the order of rows (if specified in the parameters)
	my @keys = keys %$hash;
	if (defined($sortCol) and $sortCol ne '') {
		if (defined($increasing) and $increasing) {
			if (defined($nonNumericSortCol) and $nonNumericSortCol) {
				@keys = sort {$hash->{$a}->{$sortCol} cmp $hash->{$b}->{$sortCol}} @keys;
			} else {
				@keys = sort {$hash->{$a}->{$sortCol} <=> $hash->{$b}->{$sortCol}} @keys;
			}
		} else {
			if (defined($nonNumericSortCol) and $nonNumericSortCol) {
				@keys = sort {$hash->{$b}->{$sortCol} cmp $hash->{$a}->{$sortCol}} @keys;
			} else {
				@keys = sort {$hash->{$b}->{$sortCol} <=> $hash->{$a}->{$sortCol}} @keys;
			}
		}
	}

	# then print each row (with row orders defined above)
	foreach my $k (@keys) {

		# then print other columns
		my @tmp;
		my $naCount=0;
		#print OUT $hash->{$k}->{$$column_order[0]};
		for (my $i=0; $i<scalar(@{$column_order}); $i++) {
			push @tmp, $hash->{$k}->{$$column_order[$i]};

			$naCount++ if ($tmp[$#tmp] eq '');

			# how to print blank/missing data?
			if ($NA_label ne '' and  $tmp[$#tmp] eq '') {
				$tmp[$#tmp]=$NA_label;
			}
		}

		if ($naCount<scalar(@{$column_order})) { # if not a totally empty row, let's print it
			# print key column if required
			if (defined($printKey) and $printKey==1) {
				print OUT "$k\t";
			}

			# print other required fields
			print OUT join("\t",@tmp);
			print OUT "\n";
		}
	}

	close OUT;

}


sub parseTab
{
	shift @_;
	my ($input,$colkey,$columnOrder,$noheader, $noDuplicateKeyCol, $skipRowN)=@_;
	$colkey ||= 0;
	$noheader ||= 0;
	my $hash;

	open(IN,$input) || die "cannot open ID.txt: $input\n";
	#print "Parsing $input ...\n";

	my $header;
	if (defined($skipRowN) and $skipRowN) {
		for (my $i=0; $i<$skipRowN; $i++) {
			$header=<IN>;
		}
	}

	$header=<IN>; chomp($header);
	chop($header) if ($header =~ m/\r$/); # windows ending

	my @colnames;
	if ($noheader == 1) {
		my @t=split /\t/,$header,-1;
		for (my $i=0; $i<=$#t; $i++) {
			$colnames[$i]="V$i"; # if no header, define column names as 'V1', 'V2', etc.
			$hash->{$t[$colkey]}->{$colnames[$i]}=$t[$i]; # get the value
		}
	} else {
		@colnames=split /\t/,$header;
	}

	# return @columnOrder
	if (defined($columnOrder)) {

		for (my $i=0; $i<=$#colnames; $i++) {
			$$columnOrder[$i]=$colnames[$i];
		}

		if (defined($noDuplicateKeyCol) and $noDuplicateKeyCol) {
			shift @{$columnOrder};
		}
	}

	# parse each row and build hash
	while (<IN>) {
		chomp;
		chop if (/\r$/); # windows ending
		my @t=split /\t/,$_,-1;
		unless ($#colnames == $#t) {
			die "$input is not a matrix: $#colnames colnames, but $#t columns found for row with key $t[$colkey]!!!\n";
		}

		next if ($t[$colkey] eq ''); # ignore blank keys

		for (my $i=0; $i<=$#t; $i++) {
			next if (defined($noDuplicateKeyCol) and $noDuplicateKeyCol and $i==$colkey);
			$hash->{$t[$colkey]}->{$colnames[$i]}=$t[$i];
		}
	}
	close IN;

	return $hash;
}

sub parseArray
{
# goal: parse a file and store in an array
# output: @array->[0..$#row_num]->[0..$#col_num]=$value
	shift @_;
	my ($input,$noheader,@arrayIndex)=@_;
	$noheader ||= 0;
	my $array;

	open(IN,$input) || die "cannot open ID.txt: $input\n";

	unless ($noheader) {
		my $header=<IN>; chomp($header);
		chop($header) if ($header =~ m/\r$/); # windows ending
	}

	# parse each row and build array
	my $k=0;
	while (<IN>) {
		chomp;
		chop if (/\r$/); # windows ending
		my @tmp=split /\t/,$_,-1;

		if (scalar(@arrayIndex)>0) {
			@{$array->[$k]}=@tmp[@arrayIndex];
		} else {
			@{$array->[$k]}=@tmp;
		}

		$k++;
	}
	close IN;

	return $array;
}

sub parse_IDtxt
{
	#shift @_;
	my ($self,$idtxt,$wDecoy)=@_;
	$wDecoy ||= 0;

	my ($scanhash,$peptidehash,$proteinhash);

	open(IN,"$idtxt") || die "cannot open ID.txt: $idtxt\n";
	while(<IN>)
	{
		chomp;
		next if (/^Database/ || /^Peptide/);
		my ($Peptide,$Protein,$Outfile,$measuredMH,$calcMH,$ppm,$XCorr,$dCn,$Ions,$red,$group,$subgroup,$unique,$tryptic,$pos,$precursor_peak_intensity_percentage)=split(/\;/,$_);
		my $pep=$self->intpep($Peptide);			# pure peptide sequence
		#my $pep=intpep($Peptide);			# pure peptide sequence
		my $newscan=(defined($scanhash->{$Outfile}))?0:1; 	# new scan or duplicated scan for shared proteins?

		# scanhash
		$scanhash->{$Outfile}->{IDtxt_line}=$_;
		$scanhash->{$Outfile}->{measuredMH}=$measuredMH;
		$scanhash->{$Outfile}->{calcMH}=$calcMH;
		$scanhash->{$Outfile}->{ppm}=$ppm;
		$scanhash->{$Outfile}->{XCorr}=$XCorr;
		$scanhash->{$Outfile}->{dCn}=$dCn;
		$scanhash->{$Outfile}->{Ions}=$Ions;
		$scanhash->{$Outfile}->{peptide}=$pep;
		$scanhash->{$Outfile}->{full_peptide}=$Peptide;
		$scanhash->{$Outfile}->{proteins}->{$Protein}='';
		unless ($wDecoy) {
			$scanhash->{$Outfile}->{tryptic}=$tryptic;
			$scanhash->{$Outfile}->{precursor_peak_intensity_percentage}=$precursor_peak_intensity_percentage;
		}

		# peptidehash
		$peptidehash->{$pep}->{outfiles}->{$Outfile}='';
		$peptidehash->{$pep}->{full_peptide}=$Peptide;
		$peptidehash->{$pep}->{proteins}->{$Protein}='';
		if (!defined($peptidehash->{$pep}->{XCorr}) || $XCorr>$peptidehash->{$pep}->{XCorr})
		{
			$peptidehash->{$pep}->{top_score_outfile}=$Outfile;
			
			$peptidehash->{$pep}->{measuredMH}=$measuredMH;
			$peptidehash->{$pep}->{calcMH}=$calcMH;
			$peptidehash->{$pep}->{ppm}=$ppm;
			$peptidehash->{$pep}->{XCorr}=$XCorr;
			$peptidehash->{$pep}->{dCn}=$dCn;
			$peptidehash->{$pep}->{Ions}=$Ions;

			unless ($wDecoy) {
				$peptidehash->{$pep}->{tryptic}=$tryptic;
				$peptidehash->{$pep}->{precursor_peak_intensity_percentage}=$precursor_peak_intensity_percentage;
			}
		}

		# proteinhash
		# future development
	}
	close IN;

	return ($scanhash,$peptidehash);
}

sub intpep
{
	shift @_;
	my ($pep)=@_;

	if ($pep =~ m/\./)
	{
		chop($pep); chop($pep); $pep=reverse($pep);
		chop($pep); chop($pep); $pep=reverse($pep);
	}

	return $pep;
}

sub nomod_pep
{
	my ($self,$pep)=@_;
	my $nomod=$self->intpep($pep);
	$nomod =~ s/[\@\#\*\^\~\$\%]//g;
	return $nomod;
}

# internal functions, i.e., cannot be applied by programs/modules outside
# basic statistics functions cp from http://www.perlmonks.org/?node_id=1089988
sub mean {
    my (@data) = @_;
    my $sum;
    foreach (@data) {
        $sum += $_;
    }
    return ( $sum / @data );
}
sub median {
    my (@data) = sort { $a <=> $b } @_;
    if ( scalar(@data) % 2 ) {
        return ( $data[ @data / 2 ] );
    } else {
        my ( $upper, $lower );
        $lower = $data[ @data / 2 ];
        $upper = $data[ @data / 2 - 1 ];
        return ( mean( $lower, $upper ) );
    }
}
sub std_dev {
    my (@data) = @_;
    my ( $sq_dev_sum, $avg ) = ( 0, 0 );

    $avg = mean(@data);
    foreach my $elem (@data) {
        $sq_dev_sum += ( $avg - $elem )**2;
    }
    return ( sqrt( $sq_dev_sum / ( @data - 1 ) ) );
}
sub factorial {
# https://perlmaven.com/factorial-in-perl
	my $n=shift @_;
	my $fact = reduce { $a * $b } 1 .. $n;
	return $fact;
}
# end of cp
