#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use POSIX qw(ceil floor);
use Data::Dumper;
use File::Basename;
use Cwd;

my $verbose = 1; my $label = ''; 
my ($dataset,$dbsql,$filelist,$filedir,$castor,$filesperjob,$jobs,$pretend,$args,$evjob,$triangular,$rrb,$customize,$inlinecustomize,$maxfiles,$maxevents,$skipfiles,$json,$fnal,$AAA,$addparents,$randomize);
my ($bash,$lsf,$help,$byrun,$bysize,$nomerge,$evperfilejob,$evperfile,$eosoutdir,$outdir);
my $monitor="/afs/cern.ch/user/g/gpetrucc/pl/cmsTop.pl";#"wc -l";
my $report= "/afs/cern.ch/user/g/gpetrucc/sh/report";   #"grep 'Events total'";
my $maxsyncjobs = 99;
my $subprocesses = 0;
my $firstlumiblock = 1;
my $condormem = 2000;
my %replaces;
my @job2run;
my @job2evts;

GetOptions(
    'args'=>\$args,
    'dbsql=s'=>\$dbsql,
    'dataset=s'=>\$dataset,
    'filelist=s'=>\$filelist,
    'files=s'=>\$filedir,
    'castor=s'=>\$castor,
    'files-per-job|fj=i'=>\$filesperjob,
    'jobs|n=i'=>\$jobs,
    'json|j=s'=>\$json,
    'pretend|dry-run'=>\$pretend,
    'verbose|debug+'=>\$verbose,
    'bash'=>\$bash,
    'byrun'=>\$byrun,
    'bysize'=>\$bysize,
    'monitor-script=s'=>\$monitor,
    'report-script=s'=>\$report,
    'lsf=s'=>\$lsf,
    'label=s'=>\$label,
    'triangular'=>\$triangular,
    'round-robin|rrb'=>\$rrb,
    'events-per-job|ej=i'=>\$evjob,
    'events-per-file-job|efj=i'=>\$evperfilejob,
    'events-per-file|ef=i'=>\$evperfile,
    'max-sync-jobs=i'=>\$maxsyncjobs,
    'help|h|?'=>\$help,
    'customize|c=s'=>\$customize,
    'inline-customize=s'=>\$inlinecustomize,
    'maxevents=i'=>\$maxevents,
    'maxfiles=i'=>\$maxfiles,
    'skipfiles=i'=>\$skipfiles,
    'subprocess=i'=>\$subprocesses,
    'nomerge'=>\$nomerge,
    'fnal'=>\$fnal,
    'AAA'=>\$AAA,
    'add-parents'=>\$addparents,
    'randomize'=>\$randomize,
    'first-lumi-block|flb=i'=>\$firstlumiblock,
    'eosoutdir=s'=>\$eosoutdir,
    'outdir=s'=>\$outdir,
    'condor-memory=s'=>\$condormem,
    'replace=s'=>\%replaces
);

sub usage() {
print <<EOF

usage: 
  $0 cfg.py split_options [ input ] [ options ] [ cfg.py_args]

it will create cmssw jobs   
  - cfg_job<N>.py     for all the N jobs, each producing outputs renamed as "X.root" => "X_job<N>.root"
  - cfg_merge_<X>.py  for all enabled output modules X define in the cfg.py,
                      it merges the output "X_job<N>.root" back into "X.root"
                      (enabled means in a cms.EndPath; doesn't look at the cms.Schedule yet)
  - cfg_merge_TFileService.sh to merge TFileService output (if any is there)
  - cfg_cleanup.sh    script that deletes all intermediate outputs files (but not the final output root files)

split_options:
  * split per files: use one of these two options
       --jobs:          specify the number of jobs (files in input are split evenly among jobs)
                        shortcut: -n
       --files-per-job: specify the number of files per job (number of jobs is computed automatically)
                        shortcut: --fj
  * split per events: use *both* these options
       --jobs:          specify the number of jobs 
                        shortcut: -n
       --events-per-job: specify the number of events per job
                        shortcut: --ej
      note: this will not work if you use other PoolSource parameters like 'lumisToProcess' or 'eventsToProcess'
      note 2: it won't change the random seeds unless you add the --randomize option
  * split per events from files: 
       --events-per-file-job: specify the number of events per job
                        shortcut: --efj
      note: this will be very slow, as it will have to query the number of events in each file to edmFileUtil or DAS
      note 2: some jobs will process less than this number of events 

input:
  default:    takes as input the files from the cfg.py
  --dbsql:    executes DBS query, and takes as input the filenames in the output
  --dataset:  takes as input the files in the given dataset and available at cern
  --filelist: takes as input the root files contained in this file. 
              it works also if the file contains more columns than the file names (but no more than one filename per line)
              it does not work if the file contains duplicates (not yet, at least)
  --files:    takes as input the contents of a local directory, plus an optional glob for the file
              (e.g. /data/gpetrucc/<dir>/  or /data/gpetrucc/<dir>/<pattern>)
  --json:     takes as input a JSON file to apply

options:
  --verbose: print out debug information; can be used more than once for extra verbosity
  --pretend: just print out most of what this script should do, and execute nothing.
  --byrun:   split files by run (works faster for prompt reco only, where the run can be inferred by the file name).
  --bash:    creates also a bash script cfg_local.sh that spawns the <N> jobs locally, waits for them and then merges the outputs.
  --lsf:     creates also a bash script that will submit jobs to the specified LSF queue.
             warning: this uses my private cmsRun LSF wrapper, it might not work for you.
  --args:    handle cfg.py files that take command line arguments 
  --label:   rename all intermediate and output files, inserting an extra "_<label>" before "_job<N>", "_local" & so on...
   --triangular: instead of making jobs of uniform size, scale them linearly from zero to twice the average size
                 it will give you longer latency for some jobs, but also a fast feedback on the first ones  
  --monitor-script: to be used with "--bash", changes the command used to monitor logfiles (default is "wc -l")
  --report-script:  to be used with "--bash", changes the command used to make a final report of logfiles (default is "grep 'Events total'")
  --customize: append specified python fragment to the cfg
  --replace A=B: replace A with B (can be specified multiple times)
  --subprocess N: assume the subprocess of order N is the one doing the real output. Only N=0 (no subprocess) or N=1 are supported now.
  --fnal: use fnal xrootd redirector
  --AAA:  use cern global xrootd redirector
  --nomerge:  don't merge the outputs
  --eosoutdir: EOS directory to store the jobs output files
EOF
}


if ($help) { usage(); exit(0); }
if (not(@ARGV)) { print " *** You have to specify an input file *** \n"; usage(); exit(0); }

my $filename = shift(@ARGV);
my $basename = $filename; $basename =~ s/\.py$//;
if ($label) { $label = "_$label"; }
my $runme    = "python -i $filename " . join(' ', map("'$_'", @ARGV)) . " 2> /dev/null";
if ($ENV{'CMSSW_VERSION'} =~ /1[2-9]_\d+_.*/) {
    $runme =~ s/^python /python3 /;
}
if ($customize) {
    if (!(-f $customize)) {  print " *** File $customize not found *** \n"; exit(1); }
    $runme = "cat $customize - | $runme";
}


my @cleanup = ();


my $THEPROCESS = "process";
if ($subprocesses == 1) {
    $THEPROCESS = 'getattr(process.subProcess, "_SubProcess__process")';
}

#===============================================================
# query inputs from python file 
use File::Temp qw/ :POSIX /;
my $py_out_file = tmpnam();
my $queryPythonFile = <<EOF;
cmsSplit_output_file = open("$py_out_file","w")

## INPUT
for x in process.source.fileNames: 
    cmsSplit_output_file.write("IN\\t\%s\\t-\\n" % x)

## OUTPUT
class EPVisitor:
    def enter(self,visitee):
        global cmsSplit_output_file
        #print "Visiting type ",type(visitee),": ",visitee
        if type(visitee) == cms.OutputModule:
            cmsSplit_output_file.write("OUT\\t\%s\\t\%s\\t\%s\\n" % (visitee.label_(), visitee.fileName.value(), visitee.type_()));
        if ("Dump" in visitee.label_()) and (type(visitee) == cms.EDAnalyzer) and hasattr(visitee, "outName"):
            cmsSplit_output_file.write("OUT\\t\%s\\t\%s\\t\%s\\n" % (visitee.label_(), visitee.outName.value(), "Dump"));
    def leave(self,visitee): 
        pass

epv = EPVisitor()

if getattr(process,"schedule",None):
    for p in process.schedule:
        if type(p) == cms.EndPath:
            p.visit(epv)
else:
    for ep in [V for N,V in $THEPROCESS.endpaths_().items()]:
       ep.visit(epv)

## TFileService
if hasattr(process,"TFileService") and type($THEPROCESS.TFileService) == cms.Service: 
    cmsSplit_output_file.write("TFS\t"+$THEPROCESS.TFileService.fileName.value()+"\\tX\\tX\\n")

cmsSplit_output_file.close()
EOF
if ($inlinecustomize) {
    $queryPythonFile = "## Inline customize begin\n$inlinecustomize\n## Inline customize end\n" . $queryPythonFile;
}

my @pythonCrap = qx{ echo '$queryPythonFile'  | $runme 2>&1 };
open PYTHONFILEINFO, "$py_out_file" or die "Python inspection didn't produce output.\nIt shouted ".join('',@pythonCrap)."\n";
my @pythonFileInfo = <PYTHONFILEINFO>;
print @pythonCrap if $verbose > 1;
print @pythonFileInfo if $verbose > 1;
close PYTHONFILEINFO;

my @files = ();
if (defined($dbsql)) {
    print "Using input files from DBS Query $dbsql\n" if $verbose;
    $maxfiles = 999999 unless defined($maxfiles);
    @files = grep(m/\/store.*.root/, qx(dasgoclient --query='$dbsql' --limit $maxfiles));
} elsif (defined($dataset)) {
    print "Using input files from DBS Query for dataset $dataset \n" if $verbose;
    $maxfiles = 999999 unless defined($maxfiles);
    @files = grep(m/\/store.*.root/, qx(dasgoclient --limit $maxfiles --query='file dataset=$dataset'));
} elsif (defined($filelist)) {
    print "Using input files from file $filelist\n" if $verbose;
    open FILELIST, "$filelist" or die "Can't read file list $filelist\n";
    foreach (<FILELIST>) {
        my ($f) = m/.*?(\S+\.root(\?svcClass=\w+)?).*/ or next;
        if ($f =~ m{^/castor}) { $f = 'rfio:'.$f; }
        elsif ($f !~ m{^(/store|\w+:)}) { $f = 'file:'.$f; }
        push @files, $f;
    }
}  elsif (defined($filedir)) {
    print "Using input files from $filedir\n" if $verbose;
    if ($filedir =~ m{^/eos/cms} && !(-d '/eos/cms')) {
        $filedir =~ s{^/eos/cms}{};
    }
    if ($filedir =~ m{^/store/}) {
        my $stdir  = dirname($filedir);
        my $stglob = basename($filedir);
        if ($stglob =~ /.*(\*|\?).*|.*\.root/) {
            $stglob =~ s/\./\\./; # convert a 
            $stglob =~ s/\?/./;   # regex to 
            $stglob =~ s/\*/.*/;  # a glob.
        } else {
            $stdir = $filedir;
            $stglob = ".*";
        }
        @files = ();
        print "$stdir\n";
        foreach my $line (qx(eos ls $stdir)) {
            my ($thisfile) = ($line =~ m/(\S+\.root)/) or next;
            if (basename($thisfile) =~ m{$stglob}) {
                push @files, "$stdir/$thisfile";
            }
        }
    } else {
        if (-d $filedir) { $filedir .= "/*.root"; }
        @files = map( "file:$_", glob($filedir) );
    }
} else {
    print "Using files in existing cfg file\n" if $verbose;
    foreach (@pythonFileInfo) { 
        my ($what,$arg,$arg2) = split(/\s+/,$_); 
        if ($what eq "IN") { push @files, $arg; } 
    }
}
chomp @files;
if ($AAA) {
    print "Using cms-xrd-global.cern.ch redirector\n" if $verbose;
    foreach (@files) { s{^/store/}{root://cms-xrd-global.cern.ch//store/}; }
} elsif ($fnal) {
    print "Using cmsxrootd.fnal.gov redirector\n" if $verbose;
    foreach (@files) { s{^/store/}{root://cmsxrootd.fnal.gov//store/}; }
}
if (defined($skipfiles)) {
    @files = @files[$skipfiles .. $#files];
}
if (defined($maxfiles) and (scalar(@files) > $maxfiles)) {
    @files = @files[0 .. ($maxfiles-1)];
}

if (defined($evjob)) {
    die "If you specify the number of events per job, you must specify the number of jobs and not the number of files per job" unless defined($jobs) and not defined($filesperjob);
    $filesperjob = scalar(@files);
} elsif (defined($filesperjob)) {
    die "Can't use 'files-per-job and jobs at the same time.\n" if defined($jobs);
    $jobs = ceil(scalar(@files)/$filesperjob);
} elsif (defined($evperfilejob)) {
    $filesperjob = scalar(@files);
    $jobs = 1;
    print "This will be slow\n" if $verbose and not $evperfile;
} else {
    die "Please specify the number of jobs (parameter --jobs or -n, or --files-per-job or -nj).\n" unless defined($jobs);
    if ($jobs > scalar(@files)) { $jobs = scalar(@files); }
    $filesperjob = ceil(scalar(@files)/$jobs);
}


print "Found ".scalar(@files)." files, jobs set to $jobs, files per job $filesperjob.\n" if $verbose;
print "List of files: \n\t" . join("\n\t", @files, '') . "\n" if $verbose > 1;

#===============================================================
# pool output module
my @outputModules = ();
my @nanoModules = ();
my @dumpModules = ();
my %mergeList = ();
my %mergeNano = ();
my %mergeDump = ();
foreach (@pythonFileInfo) {
    chomp;
    #print STDERR "[[$_]]\n";
    my ($type,$module,$file,$otype) = split(/\s+/) or next;
    next unless $type eq "OUT";
    #$file =~ s/\.root$/$label.root/;
    my $ofile = $file; $ofile =~ s/\.root$/$label.".root"/e;
    if ($otype eq "NanoAODOutputModule") {
        push @nanoModules, [$module,$file];
        print "Found enabled $otype (NANO) output module $module producing $file\n" if $verbose > 0;
        unless ($nomerge) {
            $mergeNano{$module} = {'outfile'=>$ofile, 'infiles'=>[]};
        }
    } elsif ($otype eq "Dump") {
        $ofile =~ s/\.dump$/$label.".dump"/e;
        push @dumpModules, [$module,$file];
        print "Found enabled $otype output module $module producing $file\n" if $verbose > 0;
        unless ($nomerge) {
            $mergeDump{$module} = {'outfile'=>$ofile, 'infiles'=>[]};
        }
    } else {
        push @outputModules, [$module,$file];
        print "Found enabled output $otype (EDM) module $module producing $file\n" if $verbose > 0;
        unless ($nomerge) {
            $mergeList{$module} = {'outfile'=>$ofile, 'infiles'=>[]};
        }
    }
}

#===============================================================
# TFileService
my $tfsFile;
foreach (@pythonFileInfo) { 
    my ($what,$arg,$arg2,$arg3) = split(/\s+/,$_); 
    if ($what eq "TFS") { 
        $tfsFile = $arg;
        print "Must handle TFileService producing $tfsFile\n" if $verbose > 0;
    }
}
my @tfsMerge = ();


#===============================================================
# handling of jobs that take command line arguments
sub fixArgs {
    my $cfgtxt = shift(@_);
    my $argsAsTxt = "[" . join(", ", map("'$_'", "cmsRun", $filename, @ARGV)) . "]";
    if ($cfgtxt =~ m/^import\s+sys\s*;?\s*$/m) {
        $cfgtxt =~ s/^import\s+sys\s*;?\s*$/import sys\nsys.argv = $argsAsTxt/m;
    } elsif ($cfgtxt =~ m/^from\s+sys\s+import\s+argv\s*;?\s*$/m) {
        $cfgtxt =~ s/^from\s+sys\s+import\s+argv\s*;?\s*$/from sys import argv\nsys.argv = $argsAsTxt/m;
    } else {
        die "I can't see where you import sys.argv in the cfg file. hints are:\n\t" . join("\t",grep(/sys|arg/, split(/^/m,$cfgtxt))) ."\n";
    }
    return $cfgtxt;
}
#===============================================================
# make jobs

sub split_even {
    my @x = @files;
    my @ret = ();
    my $min = floor(scalar(@x)/$jobs);
    my $rem = scalar(@x) - $min*$jobs;
    foreach my $i ( 1 .. $jobs ) {
        my @this = ();
        foreach ( 1 .. $min ) { push @this, shift(@x); }
        if ($i <= $rem) { push @this, shift(@x); }
        push @ret, [@this];
    }
    return [@ret];
}
sub split_rrb {
    my @x = @files;
    my @ret = ();
    foreach my $j ( 1 .. $jobs ) {
        push @ret, [];
    }
    foreach my $i ( 0 .. $#files ) {
        my $j = ($i % $jobs);
        push @{$ret[$j]}, shift(@x);
    }
    return [@ret];
}
sub split_triang {
    my @x = @files;
    my $scale = scalar(@x)/($jobs*($jobs+1)/2);
    my @ret = (); my $got = 0;
    foreach my $i ( 1 .. $jobs ) {
        my @this = ();
        my $limit = ceil($i*($i+1)/2 * $scale);
        do {
            push @this, shift(@x); 
            $got++;
        } while ($got < $limit);
        push @ret, [@this];
        last unless @x;
    }
    return [@ret];
}

open SRC, $filename; my $src = join('',<SRC>); close SRC;
if ($customize) {
    $src .= "###\n### Begin customize using $customize\n###\n"; 
    open CUST, $customize; $src .= join('',<CUST>); close CUST;
    $src .= "###\n### END customize using $customize\n###\n"; 
}

my $splits;
if ($byrun) {
    my @allfiles = @files;
    my %run2file = ();
    foreach my $f (@allfiles) {
        my @runs; 
        if ($f =~ m{/store/data/[^/]+/[^/]+/[^/]+/v\d+/000/(\d\d\d)/(\d\d\d)/[0-9A-F]+/[0-9A-F\-]+.root}) {
            push @runs, "$1$2";
        } else {
            my @edmlsout = qx{edmFileUtil -e $f};
            foreach (@edmlsout) {
                m/^\s+(\d\d\d\d\d\d)\s+0\s+0\s+\d+\s+\(Run\).*/ and push @runs, $1;
            }
            print "File: $f, runs = @runs\n";
        }
        if (($#runs > 0) and $json) { die "Can't use both a json file and a split by run with multiple runs per input file, sorry.\n"; }
        #print "File: $f, runs = @runs\n";
        foreach my $run (@runs) {
            $run2file{$run} = [] unless defined $run2file{$run};
            push @{$run2file{$run}}, $f;
        }
    }
    #print Dumper(\%run2file);
    my @alljobs = ();
    foreach my $run (sort(keys(%run2file))) {
        @files = @{$run2file{$run}};
        $jobs = ceil(scalar(@files)/$filesperjob);
        print "Run $run: Will make $jobs jobs for " . scalar(@files) . " files\n";
        foreach my $filelist ( @{ split_even() } ) {
            push @alljobs, $filelist;
            push @job2run, $run;
        }
    }
    $splits = \@alljobs;
    $jobs = scalar(@alljobs);
} elsif ($bysize) {
    die "--bysize makes no sense if making just one job\n" if ($jobs == 1);
    my %file2size = ();
    foreach my $f (@files) {
        # skip if we already know the size
        next if defined $file2size{$f};
        # it's an eos file
        if ($f =~ m{^(/eos/cms)?/store/.*.root}) {
            my $dir = dirname($f);
            $dir =~ s{^(/eos/cms)}{};
            my @eosls = qx{ls -l /eos/cms$dir};
            foreach (@eosls) {
                my (undef,undef,undef,undef,$size,undef,undef,undef,$base) = split(/\s+/) or next;
                $file2size{"$dir/$base"} = $size;
            }
        } elsif ($f =~ m{^((?:file:)?)(.+\.root)$}) {
            my $dir = dirname($2);
            if ($dir) {
                my @ls = qx{ls -l $dir};
                foreach (@ls) {
                    m/^total\s+\d+/ and next;
                    my (undef,undef,undef,undef,$size,undef,undef,undef,$base) = split(/\s+/) or print;
                    $file2size{"$1$dir/$base"} = $size;
                }
            } else {
                my @ls = qx{ls -l $2};
                foreach (@ls) {
                    my (undef,undef,undef,undef,$size,undef,undef,undef,$f2) = split(/\s+/) or next;
                    $file2size{"$1$f2"} = $size;
                }
            }
        } else {
            die "Not supported yet\n";
        }
    }
    my $tot = 0; 
    foreach my $f (@files) {
        die "Could not find size for file $f " unless $file2size{$f};
        $tot += $file2size{$f};
    }
    printf ("Total size: %.3f Mb, approximate size per job %.3f Mb \n", $tot/1024.0/1024.0, $tot/1024.0/1024.0/$jobs);
    my @fsorted = sort { $file2size{$b} cmp $file2size{$a} } @files;
    my %taken = ();
    my $cut = $tot / $jobs;
    my @alljobs = ();
    while (scalar(keys(%taken)) < scalar(@files)) {
        my $subtot = 0;
        my @jobfiles = ();
        foreach my $f (@fsorted) {
            next if $taken{$f};
            if ($subtot == 0 or $subtot + $file2size{$f} < $cut) {
                push @jobfiles, $f; $taken{$f} = 1;
                $subtot += $file2size{$f};
            }
        }
        push @alljobs, [ @jobfiles ];
        printf ("Job %d, %d files, %.3f Mb\n", scalar(@alljobs), scalar(@jobfiles), $subtot/1024.0/1024.0); 
    }
    #print Dumper(\@alljobs);
    $splits = \@alljobs;
    $jobs = scalar(@alljobs);
    #die "Fin qui tutto bene\n";
} elsif ($evperfilejob) {
    my %file2events = ();
    if ($evperfile) {
        foreach my $f (@files) {
            $file2events{$f} = $evperfile;
        }
    } else {
        my $flatfiles = join(' ', @files);
        my @epj = qx{ ~gpetrucc/py/edmlsevents.py $flatfiles };
        foreach (@epj) {
            my ($f,$e) = m/(^\S+)\s+(\d+)/ or next;
            $file2events{$f} = $e;
        }
    }
    my @fsorted = sort { $file2events{$b} cmp $file2events{$a} } @files;
    my %taken = ();
    my @alljobs = ();
    while (scalar(keys(%taken)) < scalar(@files)) {
        my $subtot = 0;
        my @jobfiles = ();
        foreach my $f (@fsorted) {
            next if $taken{$f};
            if ($subtot == 0 or $subtot + $file2events{$f} < $evperfilejob) {
                push @jobfiles, $f; $taken{$f} = 1;
                $subtot += $file2events{$f};
            }
        }
        printf ("Step 0: job %d, %d files, %d events\n", scalar(@alljobs), scalar(@jobfiles), $subtot); 
        if ($subtot < 1.5 * $evperfilejob) {
            printf ("Step 1: job %d, %d files, %d events\n", scalar(@alljobs), scalar(@jobfiles), $subtot); 
            push @alljobs, [ @jobfiles ];
            push @job2evts, [ 0, -1 ];
        } else {
            my $estart = 0;    
            while ($estart < $subtot) {
                printf ("Step 1: job %d, %d files, max %d events, skip %d\n", scalar(@alljobs), scalar(@jobfiles), $evperfilejob, $estart); 
                push @alljobs, [ @jobfiles ];
                push @job2evts, [ $estart, $evperfilejob ];
                $estart += $evperfilejob;
            }
        }
    }
    $splits = \@alljobs;
    $jobs = scalar(@alljobs);
} elsif ($triangular) {
    $splits = split_triang(); 
} elsif ($rrb) {
    $splits = split_rrb(); 
} elsif ($evjob) {  # in this case, put all files in all jobs
    $splits = [ map [@files], ( 1 .. $jobs ) ];
} else {
    $splits = split_even();
}

if ($json) {
    if ($json =~ /^https?:.*/) {
        system("wget -N $json");
        $json = basename($json);
    }
    if ($json !~ /^\/.*/) {
        $json = getcwd() . "/" . $json;
    }
    print "Will use JSON $json\n";
}
foreach my $j (1 .. $jobs) {
    my $pyfile = $basename . $label . "_job$j.py";
    print "Will create job $j, source $pyfile\n" if $verbose;
    my $postamble = ""; my @myfiles = ();
    if (defined($evjob)) {
        @myfiles = @files;
        my $inputfiles = join('', map("\t'$_',\n", @myfiles));
        $postamble .= "process.source.fileNames = [\n$inputfiles]\n" if (@files);
        $postamble .= "if process.source.type_() != 'EmptySource':\n";
        $postamble .= "    process.source.skipEvents = cms.untracked.uint32(".(($j-1)*$evjob).")\n";
        $postamble .= "else:\n";
        $postamble .= "    process.source.firstLuminosityBlock = cms.untracked.uint32(".($j +  $firstlumiblock).")\n";
        $postamble .= "process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32($evjob))\n";
    } else {
        @myfiles = @{$splits->[$j-1]};
        my $inputfiles = join('', map("\t'$_',\n", @myfiles));
        $postamble .= "process.source.fileNames = [\n$inputfiles]\n";
        if (defined($evperfilejob)) {
            $postamble .= "process.source.skipEvents = cms.untracked.uint32(".($job2evts[$j-1]->[0]).")\n";
            $postamble .= "process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(".($job2evts[$j-1]->[1])."))\n";
        } elsif (defined($maxevents)) {
            $postamble .= "process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32($maxevents))\n";
        } else {
            $postamble .= "process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))\n";
        }
    }
    if ($addparents) {
        my %parents = ();
        print " querying parents for all ".scalar(@myfiles)." files in job - this can take some time.\n" if $verbose;
        foreach my $child (@myfiles) {
            die "--add-parents works only if the files are lfns " unless $child =~ m{^/store/.*root};
            my @myparents = grep(m/\/store.*.root/, qx(das_client.py --limit 10000 --query='parent file=$child'));
            chomp @myparents;
            print STDERR "ERROR, no parents found for $child" unless @myparents;
            foreach my $parent (@myparents) {
                $parents{$parent} = 1;
            }
        }
        print " found ".scalar(keys(%parents))." parent files\n" if $verbose;
        my $parentfiles = join('', map("\t'$_',\n", keys(%parents)));
        $postamble .= "process.source.secondaryFileNames = cms.untracked.vstring()\n";
        $postamble .= "process.source.secondaryFileNames = [\n$parentfiles]\n";
    }
    foreach (@outputModules) {
        my ($n,$f) = @$_;
        $f =~ s/\.root$/$label ."_job$j.root"/e;
        if ($jobs == 1) { $f =~ s/_job1//; }
        if (defined($outdir)) { $f = $outdir . "/" . basename($f); }
        $postamble .= "$THEPROCESS.$n.fileName = '$f'\n";
        unless($nomerge or ($jobs == 1)) {
            push @{$mergeList{$n}->{'infiles'}}, $f;
            push @cleanup, $f;
        }
    }
    foreach (@nanoModules) {
        my ($n,$f) = @$_;
        $f =~ s/\.root$/$label ."_job$j.root"/e;
        if ($jobs == 1) { $f =~ s/_job1//; }
        if (defined($outdir)) { $f = $outdir . "/" . basename($f); }
        $postamble .= "$THEPROCESS.$n.fileName = '$f'\n";
        unless($nomerge or ($jobs == 1)) {
            push @{$mergeNano{$n}->{'infiles'}}, $f;
            push @cleanup, $f;
        }
    }
    foreach (@dumpModules) {
        my ($n,$f) = @$_;
        $f =~ s/\.dump$/$label ."_job$j.dump"/e;
        if ($jobs == 1) { $f =~ s/_job1//; }
        if (defined($outdir)) { $f = $outdir . "/" . basename($f); }
        $postamble .= "$THEPROCESS.$n.outName = '$f'\n";
        unless($nomerge or ($jobs == 1)) {
            push @{$mergeDump{$n}->{'infiles'}}, $f;
            push @cleanup, $f;
        }
    }
    if ($tfsFile) {
        my $f = $tfsFile; 
        $f =~ s/\.root$/$label ."_job$j.root"/e;
        if ($jobs == 1) { $f =~ s/_job1//; }
        if (defined($outdir)) { $f = $outdir . "/" . basename($f); }
        $postamble .= "$THEPROCESS.TFileService.fileName = '$f'\n";
        unless($nomerge or ($jobs == 1)) {
            push @tfsMerge, $f;
            push @cleanup, $f;
        }
    }
    if ($json) {
        $postamble .= "import FWCore.PythonUtilities.LumiList as LumiList\n";
        $postamble .= "process.source.lumisToProcess = LumiList.LumiList(filename = '$json').getVLuminosityBlockRange()\n";
    } elsif ($byrun) {
        my $run = $job2run[$j-1];
        $postamble .= "process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('$run:1-$run:9999999',)\n";
    }
    if ($randomize) {
        $postamble .= "## Scramble\n";
        $postamble .= "import random\n";
        $postamble .= "rnd = random.SystemRandom()\n";
        $postamble .= "for X in process.RandomNumberGeneratorService.parameterNames_():\n";
        $postamble .= "   if X != 'saveFileName': getattr(process.RandomNumberGeneratorService,X).initialSeed = rnd.randint(1,99999999)\n";
    }
    if ($inlinecustomize) {
        $postamble = "## Inline customize begin\n$inlinecustomize\n## Inline customize end\n" . $postamble;
    }
    print " and will append postamble\n$postamble\n" if $verbose > 1;
    my $text = $src . "\n### ADDED BY cmsSplit.pl ###\n" . $postamble;
    if ($args) {
        $text = fixArgs($text);
        print "and it did change this lines about cmd line args:\n\t" . join("\t",grep(/sys|arg/, split(/^/m,,$text))) ."\n" if $verbose > 1;
    }
    if (%replaces) {
        foreach my $s(keys(%replaces)) {
            my $r = $replaces{$s};
            $text =~ s{$s}{$r}g;
        }
    }

    next if $pretend;
    open OUT, "> $pyfile" or die "Can't write to $pyfile\n"; push @cleanup, $pyfile;
    print OUT $text;
    close OUT;
}

## make merge jobs
foreach my $m (sort(keys(%mergeList))) {
    my $pyfile = $basename . $label . "_merge_$m.py";
    print "Will create merge job $m, source $pyfile\n" if $verbose;
    print "Merge output file ",$mergeList{$m}->{'outfile'}," from:\n\t",join("\n\t",@{$mergeList{$m}->{'infiles'}},''),"\n" if $verbose > 1;

    next if $pretend;
    open OUT, "> $pyfile" or die "Can't write to $pyfile\n";  push @cleanup, $pyfile;
    my $out = $mergeList{$m}->{'outfile'};
    if (defined($outdir)) { $out = $outdir . "/" . basename($out); }
    my $in  = join("\n",map("\t'file:$_',",@{$mergeList{$m}->{'infiles'}}));
    print OUT <<EOF;
import FWCore.ParameterSet.Config as cms
process = cms.Process('cmsMerge')

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [\n$in\n]

process.out = cms.OutputModule("PoolOutputModule",fileName = cms.untracked.string('$out'))
process.end = cms.EndPath(process.out)   
EOF
    close OUT;
}

## make merge jobs (Nano)
foreach my $m (sort(keys(%mergeNano))) {
    my $pyfile = $basename . $label . "_merge_$m.sh";
    print "Will create merge job $m, source $pyfile\n" if $verbose;
    print "Merge output file ",$mergeNano{$m}->{'outfile'}," from:\n\t",join("\n\t",@{$mergeNano{$m}->{'infiles'}},''),"\n" if $verbose > 1;
    next if $pretend;
    open OUT, "> $pyfile" or die "Can't write to $pyfile\n";  push @cleanup, $pyfile;
    my $out = $mergeNano{$m}->{'outfile'};
    my $in  = join(" ",@{$mergeNano{$m}->{'infiles'}});
    print OUT <<EOF;
#!/bin/bash
if which haddnano.py > /dev/null 2>&1; then
    MERGE="python3 \$(which haddnano.py)";
elif which wget > /dev/null 2>&1 && wget https://raw.githubusercontent.com/cms-nanoAOD/nanoAOD-tools/master/scripts/haddnano.py -O \$CMSSW_BASE/bin/\$SCRAM_ARCH/haddnano.py -q; then
    chmod +x \$CMSSW_BASE/bin/\$SCRAM_ARCH/haddnano.py
    echo "Retrieved haddnano.py from cms-nanoAOD github master"
    MERGE="python3 \$(which haddnano.py)";
else
    MERGE="hadd -ff";
    echo "WARNING: haddnano.py not available.";
fi
\$MERGE $out $in
EOF
    close OUT;
}

## make merge jobs (Dump)
foreach my $m (sort(keys(%mergeDump))) {
    my $pyfile = $basename . $label . "_merge_$m.sh";
    print "Will create merge job $m, source $pyfile\n" if $verbose;
    print "Merge output file ",$mergeDump{$m}->{'outfile'}," from:\n\t",join("\n\t",@{$mergeDump{$m}->{'infiles'}},''),"\n" if $verbose > 1;
    next if $pretend;
    open OUT, "> $pyfile" or die "Can't write to $pyfile\n";  push @cleanup, $pyfile;
    my $out = $mergeDump{$m}->{'outfile'};
    my $in  = join(" ",@{$mergeDump{$m}->{'infiles'}});
    print OUT <<EOF;
#!/bin/bash
cat $in > $out
EOF
    close OUT;
}




## make TFileSerivce merge script
if ($tfsFile and not($nomerge)) {
    my $pyfile = $basename . $label . "_merge_TFileService.sh";
    print "Will create TFileService merge job, source $pyfile\n" if $verbose;
    my $tfsOut = $tfsFile; $tfsOut =~ s/\.root$/$label.root/;
    if (defined($outdir)) { $tfsOut = $outdir . "/" . basename($tfsOut); }

    if (!$pretend) {
        open OUT, "> $pyfile" or die "Can't write to $pyfile\n";  push @cleanup, $pyfile;
        print OUT "#!/bin/sh\n";
        print OUT "hadd -ff $tfsOut " . join(" ",@tfsMerge) . "\n";
        close OUT;
        chmod 0755, $pyfile;
    }
}

## make report script
my $reportfile = "$basename$label.report.txt";
my @outs = grep(/^(OUT|TFS)/, @pythonFileInfo);
if (scalar(@outs) == 1) {
    if ($tfsFile) {
        $reportfile = $tfsFile;
    } else {
        my ($w,$m,$file) = split(/\s+/,$outs[0]);
        $reportfile = $file;
    }
    $reportfile =~ s/\.root$/$label.report.txt/;
}

if (($bash or $lsf) and not ($pretend)) {
    my $pyfile = $basename . $label . "_report.sh";
    print "Will create reort script $pyfile\n" if $verbose;

    print STDERR "Report file will be called $reportfile\n";

    my $jlglob = $basename . $label . "_job[0-9]*.log*";

    open OUT, "> $pyfile" or die "Can't write to $pyfile\n";  push @cleanup, $pyfile;
    print OUT "#!/bin/bash\n";
    print OUT "$report $jlglob | tee $reportfile;\n";
    close OUT;
}

## make bash driver scripts
if ($bash and not($pretend)) {
    my $pyfile = $basename . $label . "_local.sh";
    print "Will create bash driver script  source $pyfile\n" if $verbose;

    open OUT, "> $pyfile" or die "Can't write to $pyfile\n";  push @cleanup, $pyfile;
    print OUT "#!/bin/bash\n";

    my $jgrep  = $basename . $label . "_job[0-9]\\+.py";
    my $jlglob = $basename . $label . "_job[0-9]*.log";

    foreach my $j (1 .. $jobs) {
        my $jfile = $basename . $label . "_job$j.py";
        my $lfile = $basename . $label . "_job$j.log";
        print OUT "(cmsRun $jfile > $lfile 2>&1 &)\n"; push @cleanup, $lfile;

        if ($j % $maxsyncjobs == 0) {
            print OUT <<EOF;
sleep 10;
while ps x | grep -q "cmsRun $jgrep"; do
    clear;
echo "At \$(date), jobs are still running...";
top -n 1 -bc -u \$UID | grep 'COMMAND\\|cmsRun $jgrep' | grep -v grep;
$monitor $jlglob;
sleep 5;
done;
EOF
        }
    }

    print OUT <<EOF;
sleep 10;
while ps x | grep -q "cmsRun $jgrep"; do
    clear;
    echo "At \$(date), jobs are still running...";
    top -n 1 -bc -u \$UID | grep 'COMMAND\\|cmsRun $jgrep' | grep -v grep;
    $monitor $jlglob;
    sleep 5;
    done;
echo "Jobs done.";
$report $jlglob | tee $reportfile;
EOF
    unless ($nomerge) {
        print OUT "echo 'Doing the merge'\n";
        foreach my $m (sort(keys(%mergeList))) {
            my $mfile = $basename . $label . "_merge_$m.py";
            my $mlog  = $basename . $label . "_merge_$m.log";
            print OUT "cmsRun $mfile >  $mlog 2>&1;\n"; push @cleanup, $mlog;
        }
        foreach my $m (sort(keys(%mergeNano))) {
            my $mfile = $basename . $label . "_merge_$m.sh";
            my $mlog  = $basename . $label . "_merge_$m.log";
            print OUT "bash $mfile >  $mlog 2>&1;\n"; push @cleanup, $mlog;
        }
        foreach my $m (sort(keys(%mergeDump))) {
            my $mfile = $basename . $label . "_merge_$m.sh";
            my $mlog  = $basename . $label . "_merge_$m.log";
            print OUT "bash $mfile >  $mlog 2>&1;\n"; push @cleanup, $mlog;
        }
        if ($tfsFile) {
            my $tfsfile = $basename . $label . "_merge_TFileService.sh";
            my $tfslog  = $basename . $label . "_merge_TFileService.log";
            print OUT "bash $tfsfile > $tfslog 2>&1; \n"; push @cleanup, $tfslog;
        }
        print OUT "echo 'All merge jobs done.'\n";
    }
    close OUT;
    chmod 0755, $pyfile;
}

## make lsf driver scripts
if ($lsf and not($pretend)) {
    if (!defined($eosoutdir)) { die "This only works with EOS for this cmsSplit."; }
    my $pyfile = $basename . $label . "_bsub.sh";
    if ($lsf =~ m/condor.*/) {
        print "Condor!\n";
        $pyfile = $basename . $label . "_condor.sub";
        print "Will create condor submit file $pyfile\n" if $verbose;
        die "Malformed --lsf condor option, it should be just 'condor' or 'condor-(number)[n](unit)', with unit = m|h|d|w." unless ($lsf =~ m/condor(-(\d+)n?([mhdw]))?/);
    } else {
        print "Will create bash driver script source $pyfile\n" if $verbose;
    }
    my $runner = $ENV{"CMSSW_BASE"}."/src/FastPUPPI/NtupleProducer/python/scripts/cmsRunBatchEOS";
    next if $pretend;
    open OUT, "> $pyfile" or die "Can't write to $pyfile\n"; push @cleanup, $pyfile;
    if ($lsf =~ m/condor(-(\d+)n?([mhdw]))?/) {
        my ($numtime, $timeunit) = ("8", "h");
        if ($1) { $numtime = $2; $timeunit = $3 };
        my $time = $numtime;
        if ($timeunit eq "m")    { $time *= 60; }
        elsif ($timeunit eq "h") { $time *= 60*60; }
        elsif ($timeunit eq "d") { $time *= 60*60*24; }
        elsif ($timeunit eq "w") { $time *= 60*60*24*7; }
        else { die "Not understood job time: $2$3\n"; }
        my $args = `pwd`; chomp($args);
        $args .= "  $eosoutdir";
print OUT <<EOF;
Universe = vanilla
Executable = $runner 
use_x509userproxy = \$ENV(X509_USER_PROXY)

Log        = ${basename}${label}_\$(Job).condor
Output     = ${basename}${label}_\$(Job).out
Error      = ${basename}${label}_\$(Job).error
getenv      = True
request_memory = $condormem
transfer_output_files = ""
+MaxRuntime = $time

Arguments = $args ${basename}${label}_\$(Job).py
Queue Job from (
EOF
        foreach my $j (1 .. $jobs) { print OUT "   job$j\n"; }
        print OUT ")\n"
    } else {
        print OUT "#!/bin/sh\n";
        foreach my $j (1 .. $jobs) {
            my $jfile = $basename . $label . "_job$j.py";
            my $lfile = $basename . $label . "_job$j.log";
            print OUT "$runner -$lsf $eosoutdir $jfile\n"; push @cleanup, "$lfile.[0-9]*";
        }
    }
    close OUT;
    chmod 0755, $pyfile;
}

# make cleanup scripts
if (not($pretend)) {
    my $pyfile = $basename . $label . "_cleanup.sh";
    print "Will create cleanup script $pyfile\n" if $verbose;
    open OUT, "> $pyfile" or die "Can't write to $pyfile\n";  push @cleanup, $pyfile;
    print OUT "#!/bin/sh\n";
    foreach my $f (@cleanup) {
        print OUT "rm $f\n";
    }
    close OUT;
    chmod 0755, $pyfile;
}

unlink ($py_out_file) if (-f $py_out_file);
#print "Please delete $py_out_file\n" if (-f $py_out_file);
