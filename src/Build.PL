#!/usr/bin/perl
eval 'exec /usr/bin/env perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/inc/lib"; #libary just needed for configuration
use lib "$FindBin::Bin/inc/perl/lib"; #libary just needed for configuration
use lib "$FindBin::Bin/lib"; #libary just needed for configuration
use lib "$FindBin::Bin/../perl/lib"; #libary just needed for configuration
use lib "$FindBin::Bin/../lib"; #libary just needed for configuration
use INSTALL::Build;

chdir($FindBin::Bin); #alway go here for safety

#configure Build
my $build = INSTALL::Build->new(
     module_name => 'WaveCNV',
     dist_author => 'Carson Holt <carson.holt@oicr.on.ca>',
     dist_version_from  => "$FindBin::Bin/version",
     dist_abstract => ('segmenter.pl takes a pileup file and produces segments as input to the CNV caller'
		       ),
     dynamic_config => 1,
     license  => 'artistic_2',
     auto_configure_requires => 1,
     requires => {
	 'perl'               => '5.8.8',
	 #'Other::Module' => '>= 1.2, != 1.5, < 2.0',
     },
     exe_requires => {
	 #'samtools' => 'samtools',
     },
     lib_requires => {
	 #'samtools bam.h' => 'bam.h',
     },
     build_requires => {
	 #'Test::More'         => '0',
     },
     config_requires => {
	 #'File::Which'        => '0',
     },
     recommends => {
	 'Carp'               => '1.24',
	 #'Other::Module' => '>= 1.2, != 1.5, < 2.0',
     },
     install_base => "$FindBin::Bin/../",
     install_base_relpaths => {
	 arch   => [qw(perl/lib)],
	 lib    => [qw(perl/lib)],
	 bindoc => [qw(perl/man)],
	 libdoc => [qw(perl/man)],
	 data   => [qw(data)],
     },
     PL_files => {#'inc/build/Bio-SamTools-mod_of_1.33/Build.PL' => [],
     },
     create_packlist => 0,
);

print "\n*You can use precompiled binaries of the WaveCNV segmenter or compile from source code.\n".
      "Running the precompiled binaries requirs MATLAB MCR 8.0 on a 64-bit Linux machine.\n".
      "Compiling from source requires any version of MATLAB after 2010\n\n";
my $from_source = $build->y_n("Do you prefer to compile from source (requires that you have MATLAB and MCR installed):", 'N');

if($from_source){ #compile from source
    my $MCC = File::Which::which('mcc');
    
    $MCC = $build->prompt("\nPlease specify the path to the MATLAB mcc compiler (leave blank to cancel):", $MCC);
    $MCC =~ s/^\s|\s$//g if($MCC);
    while($MCC && ! -f $MCC){
	$MCC = $build->prompt("\nCannot find mcc\n.".
			      "Please specify the path to the MATLAB mcc compiler (leave blank to cancel):", '');
	$MCC =~ s/^\s|\s$//g if($MCC);
	last if(!$MCC);
    }
    die "ERROR: Failed to specify path to mcc compiler\n" if(!$MCC);
    $build->config_data(MCC => $MCC);

    my @libs;
    my ($OS, $ARC) = (POSIX::uname())[0,4];
    if($OS eq 'Darwin'){
	@libs = grep {-d $_} </Applications/MATLAB/MATLAB_Compiler_Runtime/*>;
	push(@libs, grep {/\/sys\/os\//} split(/\:/, $ENV{DYLD_LIBRARY_PATH}))
	    if($ENV{DYLD_LIBRARY_PATH});
    }
    elsif($OS eq 'Linux'){
	push(@libs, grep {/\/sys\/os\//} split(/\:/, $ENV{LD_LIBRARY_PATH}))
	    if($ENV{LD_LIBRARY_PATH});
    }
    @libs = map {s/\/sys\/os\/.*//; $_} @libs;

    my $MCRROOT = shift @libs || '';
    $MCRROOT = $build->prompt("\nPlease specify the MATLAB MCR root directory (leave blank to cancel):", $MCRROOT);
    $MCRROOT =~ s/^\s|\s$//g if($MCRROOT);

    while($MCRROOT && (! grep {-f $_} <$MCRROOT/runtime/*/libmwmclmcrrt.*>)){
	$MCRROOT = $build->prompt("\nCannot find required MCR files\n.".
				  "Please specify the MATCLAB MCR root directory (leave blank to cancel):", '');
	$MCRROOT =~ s/^\s|\s$//g if($MCRROOT);
    }
    die "ERROR: Failed to specify MCR root directory\n" if(!$MCRROOT);
    $build->config_data(MCRROOT => $MCRROOT);
}
else{ #use MCR for preconfigured binaries
    my @libs;
    my ($OS, $ARC) = (POSIX::uname())[0,4];
    if($OS eq 'Darwin'){
	print "\n*You are using Darwin on Mac OSX, MCR version 8.2 is required\n\n";
	@libs = '/Applications/MATLAB/MATLAB_Compiler_Runtime/v82';
	push(@libs, grep {/\/sys\/os\//} split(/\:/, $ENV{DYLD_LIBRARY_PATH}))
	    if($ENV{DYLD_LIBRARY_PATH});
    }
    elsif($OS eq 'Linux'){
	print "\n*You are using Linux, MCR version 8.0 is required\n\n";
	push(@libs, grep {/\/sys\/os\//} split(/\:/, $ENV{LD_LIBRARY_PATH}))
	    if($ENV{LD_LIBRARY_PATH});
    }
    @libs = map {s/\/sys\/os\/.*//; $_} @libs;

    my $MCRROOT = shift @libs || '';
    $MCRROOT = $build->prompt("\nPlease specify the MATLAB MCR root directory (leave blank to cancel):", $MCRROOT);
    $MCRROOT =~ s/^\s|\s$//g if($MCRROOT);

    my $file = ($OS eq 'Darwin') ? 'libmwmclmcrrt.8.2.dylib' : 'libmwmclmcrrt.so.8.0';
    while($MCRROOT && (! grep {-f $_} <$MCRROOT/runtime/*/$file>)){
	$MCRROOT = $build->prompt("\nCannot find required MCR files\n.".
				  "Please specify the MATCLAB MCR root directory (leave blank to cancel):", '');
	$MCRROOT =~ s/^\s|\s$//g if($MCRROOT);
    }
    die "ERROR: Failed to specify MCR root directory\n" if(!$MCRROOT);
    $build->config_data(MCRROOT => $MCRROOT);
}

$build->config_data(PATH => $ENV{PATH});
$build->config_data(build_done => 0);
$build->create_build_script;
$build->log_info("\n\nThe file 'Build' has been created for you to finish installation.\n");
$build->status;
