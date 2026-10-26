#!/bin/env perl

use strict;
use warnings;
use Cwd;
use File::Basename;
use FindBin qw($Bin);
use lib "$Bin"; # $Bin = '/hpcf/authorized_apps/proteomics_apps/pipeline/release/version1.13.003/JUMP/bin'


my $require_file = $Bin."/batch/BatchLib.pl"; # use the file $Bin/batch/BatchLib.pl
require $require_file;

# check input
if (scalar(@ARGV)<1)
{
	die "Usage: jump -batch-id jump_batchID.params\n";
}

# current datapath (e.g. /home/zyuan1/testBatch)
my $cwd_dir = getcwd();
if (index($cwd_dir, '/home.nfs/')==0) {
	my @new_dirs = split('/',$cwd_dir);
	$new_dirs[1]='home';# '/home.nfs/' to '/home/'
	$cwd_dir = join('/',@new_dirs);
}

# get nbatch from batchID, and check if they are the same
my $batchid_param_fullfile = $cwd_dir."/".$ARGV[0];
my ($nbatch_ID,$nbatch_fj_params,$mods_inID) = get_nbatch_batchID($batchid_param_fullfile);
if ($nbatch_ID!=$nbatch_fj_params) {
	die "In ".$ARGV[0].", some batches in 'input_path_batch' do not have filter results.\n";
}

# if input memory
my $input_mem=0;
if (scalar(@ARGV)>=2 and index($ARGV[1], 'mem')!=-1) {
	my $value1 = get_value($ARGV[1]);
	if (index($value1, 'M')!=-1) {
		my ($number1,$value2) = split('M',$value1);
		$input_mem = $number1/1000;
	} elsif (index($value1, 'G')!=-1) {
		my ($number1,$value2) = split('G',$value1);
		$input_mem = $number1*1;
	} elsif (index($value1, 'T')!=-1) {
		my ($number1,$value2) = split('T',$value1);
		$input_mem = $number1*1000;
	} else {
		$input_mem = $value1/1000;
	}
}

# write_step1_fbatch
my $fbatch_sh_fullfile = $cwd_dir."/ParameterFiles/Batch/run1_fbatch.sh";
write_step1_fbatch($fbatch_sh_fullfile,$cwd_dir,$Bin,$ARGV[0]);

# first get mem and queue, and then run step1_fbatch
my $pepXMLsize = get_pepXMLsize_from_fbatch($batchid_param_fullfile);
my $nMEM1 = calculate_mem($pepXMLsize,4); # choice: 4,fbatch
if ($input_mem>0) {
	$nMEM1 = $input_mem;
}
submit_jobs($nMEM1,$fbatch_sh_fullfile);