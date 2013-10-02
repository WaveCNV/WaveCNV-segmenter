#!/usr/bin/perl

#Hack to get around FindBin error where broken Carp is preloaded
BEGIN {
    if(@ARGV == 1 && $ARGV[0] eq 'findbin'){
        eval 'require FindBin';
        print $FindBin::RealBin;
        exit;
    }

    my $Bin = `$0 findbin`;
    eval "use lib '$Bin/../inc/lib'";
    eval "use lib '$Bin/../inc/perl/lib'";
    eval "use lib '$Bin/../../perl/lib'";
    eval "use lib '$Bin/../../lib'";
    eval "use lib '$Bin/../src/inc/lib'";
    eval "use lib '$Bin/../src/inc/perl/lib'";
    eval "use lib '$Bin/../src/lib'";
    eval "use lib '$Bin/../perl/lib'";
    eval "use lib '$Bin/../lib'";
}

use strict;
use warnings;
use FindBin;
use POSIX;
use File::Copy;
use Getopt::Long;

my $usage = "
USAGE:
     $0 <pileup_file>
     $0 <pileup_file> <chr>

DESCRIPTION:
     This script takes a samtools pileup file and generates segments for CNV analysis
     using Wavelette theory.

OPTIONS:
     p_value <FLOAT> Threshold for initial segment KS merging (DEFAULT: 1e-10)

     help|?          Prints this usage statement
   
";

my ($MCC, $MCRROOT, $h5fromtxt);
eval "require WaveCNV::ConfigData";
if(!$@){
    $MCC = WaveCNV::ConfigData->config('MCC');
    $MCRROOT = WaveCNV::ConfigData->config('MCRROOT');
    $h5fromtxt = WaveCNV::ConfigData->config('h5fromtxt');
}

my $level = 11;
my $threshmethod = 'sqtwolog';
my $hdf_datafield = 'GposRd';
#my $output_file_string;
my $slice_toggle = 0;
my $start_string = '[]';
my $end_string = '[]';
my $scaling_string = 1;
my $p_value_string = 1e-10;
#my $chr_string;
my $bpoint_only_string = 1;
my $cpus = 1;

GetOptions("cpus=i"      => \$cpus,
           "level=i"     => \$level,
           "p_value=f"   => \$p_value_string,
	   "MCRROOT=s"   => \$MCRROOT,
           "help|?"      => sub{print $usage; exit(0);});

my $file = shift;
my $chr_select = shift || '';

if(! $file){
    print $usage;
    exit(0);
}

die "ERROR: File $file does not exists\n" if(! -f $file);
if($p_value_string <=0){
    warn "WARNING: P-value must be >= 0. Setting value to 1e-10\n";
    $p_value_string = 1e-10;
}
if($cpus < 1){
    warn "WARNING: cpus must be set to a positive integer. Setting value to 1\n";
    $cpus = 1;;
}
if($cpus > 1){
    warn "NOTE: cpus is set to $cpus. Make sure you have at least 25Gb of free memory for each CPU\n";
}

#now do analyis on each chromosome
my ($name) = $file =~ /([^\/]+)$/;
my $outdir = "$name.segmenter.output";
my $h5dir = "$outdir/h5";
my $segdir = "$outdir/seg";
mkdir $outdir unless(-d $outdir);
mkdir $h5dir  unless(-d $h5dir);
mkdir $segdir unless(-d $segdir);

#read in pileup
my %finished;
foreach my $f (<$h5dir/*.h5>){
    my ($id) = $f =~ /([^\/]+)\.h5$/;
    next if($chr_select && $id ne $chr_select);
    $finished{$id}++;
}

unless($finished{$chr_select}){
    my $tell = 30000000000;
    if(-f "$outdir/tell"){
	open(TELL, "<$outdir/tell");
	my $tell = <TELL>;
	close(TELL);
	chomp $tell;
    }

    my $last_ch;
    my $last_p;
    my $OUT;
    open(IN, "<$file");
    seek(IN,$tell,0) if(defined($tell));
    if($tell == 30000000000){ #temp
	<IN>; #temp
    } #temp
    while(my $line = <IN>){
	chomp $line;
	my @F = split(/\t/, $line);
	next unless(@F >= 4);
	
	my ($chr, $pos, $depth) = ($F[0], $F[1], $F[3]);
	next if($finished{$chr});
	next if($chr_select && $chr ne $chr_select);
	
	if(!$last_ch){
	    unlink("$h5dir/$chr.h5.tmp") if(-f "$h5dir/$chr.h5.tmp");
	    open($OUT, '>', "$h5dir/$chr.h5.tmp"); 

	    if(defined $tell){
		open(TELL, ">$outdir/tell");
		print TELL $tell;
		close(TELL);
	    }
	}
	elsif($chr ne $last_ch){
	    close($OUT);
	    undef $OUT;
	    File::Copy::move("$h5dir/$last_ch.h5.tmp", "$h5dir/$last_ch.h5");
	    $finished{$last_ch}++;
	    undef $last_ch;
	    undef $last_p;
	    unlink("$h5dir/$chr.h5.tmp") if(-f "$h5dir/$chr.h5.tmp");
	    if(-f "$h5dir/$chr.h5"){
		$finished{$chr}++;
		next;
	    }
	    open($OUT, '>', "$h5dir/$chr.h5.tmp");

	    if(defined $tell){
		open(TELL, ">$outdir/tell");
		print TELL $tell;
		close(TELL);
	    }
	}

	if($last_p && $last_p != $F[1]-1){
	    my $d = ($F[1]-$last_p) + 1;

	    if($d > 2000){
		#print left side
		for(my $i = $last_p+1; $i <= $last_p+1000; $i++){
		    print $OUT "$i,0\n";
		}

		#print right side
		for(my $i = $F[1]-1000-1; $i < $F[1]; $i++){
		    print $OUT "$i,0\n";
		}
	    }
	    else{
		for(my $i = $last_p+1; $i < $F[1]; $i++){
		    print $OUT "$i,0\n";
		}
	    }
	}
	
	print $OUT "$F[1],$F[3]\n";
	$last_ch = $chr;
	$last_p = $F[1];
	$tell = tell(IN);
    }
    close(IN);

    if($OUT){
	close($OUT);
	File::Copy::move("$h5dir/$last_ch.h5.tmp", "$h5dir/$last_ch.h5");
	$finished{$last_ch}++;
    }
}

#set up parameters
my $EXE;
if($MCC){ #use compiled from source app
    $ENV{DYLD_LIBRARY_PATH}  =  '.:'.join(':', grep {-d $_} <$MCRROOT/runtime/*>);
    $ENV{DYLD_LIBRARY_PATH} .=  ':'.join(':', grep {-d $_} <$MCRROOT/bin/*>);
    $ENV{DYLD_LIBRARY_PATH} .=  ':'.join(':', grep {-d $_} <$MCRROOT/sys/os/*>);
    $ENV{LD_LIBRARY_PATH}  =  '.:'.join(':', grep {-d $_} <$MCRROOT/runtime/*>);
    $ENV{LD_LIBRARY_PATH} .=  ':'.join(':', grep {-d $_} <$MCRROOT/bin/*>);
    $ENV{LD_LIBRARY_PATH} .=  ':'.join(':', grep {-d $_} <$MCRROOT/sys/os/*>);
    $ENV{XAPPLRESDIR} = "$MCRROOT/X11/app-defaults";
    ($EXE) = grep {-f $_} (<$FindBin::Bin/../arch/other/wavecnv_segmenter>,
			   <$FindBin::Bin/../arch/other/*/wavecnv_segmenter>,
			   <$FindBin::Bin/../arch/other/*/*/wavecnv_segmenter>,
			   <$FindBin::Bin/../arch/other/*/*/*/wavecnv_segmenter>);
}
else{ #use precompiled binary
    if(!$MCRROOT){
	my @libs = '/Applications/MATLAB/MATLAB_Compiler_Runtime/v82'
	    if(-d '/Applications/MATLAB/MATLAB_Compiler_Runtime/v82');
	if($ENV{DYLD_LIBRARY_PATH}){
	    push(@libs, grep {/\/sys\/os\//} split(/\:/, $ENV{DYLD_LIBRARY_PATH}));
	}
	if($ENV{LD_LIBRARY_PATH}){
	    push(@libs, grep {/\/sys\/os\//} split(/\:/, $ENV{LD_LIBRARY_PATH}));
	}
	@libs = map {s/\/sys\/os\/.*//; $_} @libs;
	$MCRROOT = shift @libs;
    }
    
    die "ERROR: Can't find MATLAB MCR.\n" if(!$MCRROOT);
    die "ERROR: MCR directory $MCRROOT is invalid.\n" if(! -d $MCRROOT);

    my ($OS, $ARC) = (POSIX::uname())[0,4];
    die "ERROR: This does not appear to be MCR version 8.2"
	if($OS eq 'Darwin' && !<$MCRROOT/runtime/*/libmwmclmcrrt.8.2.dylib>);
    die "ERROR: This does not appear to be MCR version 8.0"
	if($OS eq 'Linux' && !<$MCRROOT/runtime/*/libmwmclmcrrt.so.8.0>);

    if($OS eq 'Darwin' && $ARC eq 'x86_64'){
	$ENV{DYLD_LIBRARY_PATH}  = ".:$MCRROOT/runtime/maci64";
	$ENV{DYLD_LIBRARY_PATH} .= ":$MCRROOT/bin/maci64";
	$ENV{DYLD_LIBRARY_PATH} .= ":$MCRROOT/sys/os/maci64";
	$ENV{XAPPLRESDIR} = "$MCRROOT/X11/app-defaults";
	$EXE = "$FindBin::Bin/../arch/maci64/wavecnv_segmenter.app/Contents/MacOS/wavecnv_segmenter";
    }
    elsif($OS eq 'Linux' && $ARC eq 'x86_64'){
	$ENV{LD_LIBRARY_PATH}  = ".:$MCRROOT/runtime/glnxa64";
	$ENV{LD_LIBRARY_PATH} .= ":$MCRROOT/bin/glnxa64";
	$ENV{LD_LIBRARY_PATH} .= ":$MCRROOT/sys/os/glnxa64";
	$ENV{XAPPLRESDIR} = "$MCRROOT/X11/app-defaults";
	$EXE = "$FindBin::Bin/../arch/glnxa64/wavecnv_segmenter";
    }
}

#segmenter here
my @results;
foreach my $chr (keys %finished){
    next if($chr eq 'chrMT' || $chr eq 'chrM');
    next if(-f "$segdir/$chr.seg_short");
    my @args = ("$h5dir/$chr.h5",
		$level,
		$threshmethod,
		$hdf_datafield,
		"$segdir/$chr.seg",
		$slice_toggle,
		$start_string,
		$end_string,
		$scaling_string,
		$p_value_string,
		$chr,
		$bpoint_only_string);
    

    system($EXE, @args) && die "ERROR: Segmenter failed on $chr\n";
    push(@results, "$segdir/$chr.seg");
}
