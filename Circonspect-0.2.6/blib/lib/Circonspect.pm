package Circonspect;

use 5.006;
use strict;
use warnings;
use Cwd;
use File::Basename;
use File::Spec;
use Math::Random::MT qw(srand rand);
use Bio::Assembly::Singlet;
use Bio::Assembly::Tools::ContigSpectrum;

use mco; # Module for creation and manipulation of Circonspect "database"

#use Parallel::MPI qw(:all);

our $VERSION = '0.2.6';

=head1 NAME

Circonspect - Control In Research on CONtig SPECTra

=head1 SYNOPSIS

Run this in a terminal: Circonspect --help

=head1 DESCRIPTION

Circonspect (Control In Research on CONtig SPECTra) was designed to create
contig spectra based on the assembly of random subsets of metagenomes. The
contig spectra can then be modeled downstream to get diversity estimates (for
example with PHACCS (http://sourceforge.net/projects/phaccs). Different kind of
contig spectra can be created: the standard contig spectrum for alpha-diversity,
the mixed contig spectrum, for gamma-diversity, and the cross-contig spectrum,
for beta diversity. Building a cross-contig spectrum requires building a mixed
contig spectrum, which requires building normal contig spectra based on several
selected metagenomes (FASTA files). Circonspect uses a indexed sequence
database to store temporary data. Any assembly program that support a minimum
overlap length and similarity is adapted; at this time, Minimo, TIGR
Assembler, LIGR Assembler, CAP3 and Newbler are supported. When multiple
repetitions are done (bootstrapping), the average contig spectrum is
calculated. The user can control the size of the sequences by trimming them to
a given length and removing sequences smaller than a threshold.


=head1 AUTHOR

Florent ANGLY, C<< <florent.angly at gmail.com> >>

=head1 BUGS

All complex software has bugs lurking in it, and this program is no exception.
If you find a bug, please report it on the SourceForge Tracker for Grinder:
L<https://sourceforge.net/tracker/?group_id=231834&atid=1084346>

Bug reports, suggestions and patches are welcome. Circonspect's is developed on
Sourceforge under Git revision control (L<https://sourceforge.net/scm/?type=git&group_id=231834>).
To get started with a patch, do:

   git clone git://circonspect.git.sourceforge.net/gitroot/circonspect/circonspect

=head1 COPYRIGHT

Copyright 2010 Florent ANGLY

=head1 LICENSE

Circonspect is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License (GPL) as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
Circonspect is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with Circonspect. If not, see <http://www.gnu.org/licenses/>.

=cut

#------------------------------------------------------------------------------#

sub Circonspect {
  # Circonspect main routine
  my ( $fastafiles, $qualfiles, $discard, $trim, $outdir, $outprefix, $outbasic,
    $outdetails, $outstats, $outace, $samplesize, $metapercent, $samplecompo,
    $repnumber, $mincoverage, $seed, $asmprog, $minidentity, $minoverlap, $mixed,
    $cross, $cplxthres, $cplxwindow ) = @_;

  # Do basic checks on sampling parameters
  ($mixed, $cross, $asmprog, $samplecompo, $outdir, $outprefix, $discard) =
    basic_param_check( $discard, $trim, $asmprog, $samplecompo, $fastafiles,
    $mixed, $cross, $outdir, $outprefix, $outstats, $outace);

  # Need some filehandles
  my $basic_fh = \*STDOUT;
  if ($outbasic) {
    my $outbasic_file = File::Spec->catfile($outdir, $outprefix.'.csp');
    open($basic_fh, ">$outbasic_file") || die ("Error: Cannot write output ".
      "file $outbasic_file\n");
  }
  my $detail_fh;
  if ($outdetails) {
    my $outdetails_file = File::Spec->catfile($outdir, $outprefix.'-details.txt');
    open($detail_fh, ">$outdetails_file") || die ("Error: Cannot write output ".
      "file $outdetails_file\n");
  }

  print $detail_fh "Circonspect $VERSION\n\n" if $detail_fh;

  # Display input parameters
  #print $detail_fh params_str($fastafiles, $qualfiles, $discard, $trim,
  #  $outdir, $outprefix, $outdetails, $outstats, $outace, $samplesize,
  #  $metapercent, $samplecompo, undef, $repnumber, $mincoverage, $seed,
  #  $asmprog, $minidentity, $minoverlap, $mixed, $cross) if $detail_fh;

  # Find QUAL file corresponding to each FASTA file
  if ($qualfiles) {
    $qualfiles = getqual($fastafiles, $qualfiles);
  }

  # Seed the random number generator
  $seed = srand($seed);

  # Create struct to hold BioPerl objects
  my $mco = mco->new();

  # Import general metagenome info and produce unique metagenome names
  my $metanames = $mco->import_metainfo($fastafiles, $qualfiles, $outdir,
    $outprefix);

  # Import original sequences (and quality scores)
  my $imported = $mco->import_sequences($fastafiles, $qualfiles, $discard,
    $cplxthres, $cplxwindow, $outdir, $outprefix);
  print $detail_fh $imported if ($detail_fh);

  # Get some stats
  my ($stats, $report) = $mco->get_size_stats($fastafiles, $trim);
  print $detail_fh $report."\n" if ($detail_fh);

  # More input validty tests and compute sampling parameters
  ($samplesize, $metapercent, $repnumber, $mincoverage, my $compolist) =
    extended_param_check($samplesize, $metapercent, $repnumber, $mincoverage,
    $samplecompo, $mixed, $cross, $fastafiles, $stats);

  # Display processed parameters
  print $detail_fh params_str($fastafiles, $qualfiles, $discard, $trim, $outdir,
    $outprefix, $outdetails, $outstats, $outace, $samplesize, $metapercent,
    $samplecompo, $compolist, $repnumber, $mincoverage, $seed, $asmprog,
    $minidentity, $minoverlap, $mixed, $cross) if $detail_fh;

  # Real computation begins here    
  # For all repetitions
  print $detail_fh "Computing...\n" if ($detail_fh);
  my %csp_collection;
  #MPI_Init();
  #my $mpi_size = MPI_Comm_size(MPI_COMM_WORLD); # number of parallel executions
  #my $mpi_rank = MPI_Comm_rank(MPI_COMM_WORLD); # rank of the execution
  for (my $rep = 1 ; $rep <= $repnumber ; $rep++) { # non-parallel loop
  #for (my $rep = $mpi_rank+1 ; $rep <= $repnumber ; $rep += $mpi_size) { # parallel loop

    if ($repnumber > 1) {
      print $detail_fh "Repetition $rep\n" if ($detail_fh);
      #print $detail_fh "Repetition $rep (MPI proces $mpi_rank)\n" if ($detail_fh);
    }

    # Assemble and calculate contig spectra
    my $csps = process( $fastafiles, $compolist, $metanames, $mco, $samplesize,
      $trim, $asmprog, $minidentity, $minoverlap, $mixed, $cross, $detail_fh,
      $outstats, $outace, $outdir, $outprefix, $rep );

    # Save the results
    for my $sample ( @$metanames ) {
      $csp_collection{$sample}{$rep}          = $$csps{$sample};
      my $dissolved = 'dissolved '.$sample;
      $csp_collection{$dissolved}{$rep}       = $$csps{$dissolved};
    }
    $csp_collection{'mixed'}{$rep}            = $$csps{'mixed'};
    $csp_collection{'sum of dissolved'}{$rep} = $$csps{'sum of dissolved'};
    $csp_collection{'cross'}{$rep}            = $$csps{'cross'};


    print $detail_fh "\n" if ($detail_fh);

  } # end of all repetitions
  #MPI_Finalize();

  # Average contig spectra and display final results
  my $csp_average = average_csps($metanames, \%csp_collection, $repnumber, $mixed,
    $cross, $detail_fh, $basic_fh, $minoverlap, $outstats);

  # Clean after ourselves
  $mco->destroy;
  print $detail_fh "Done!\n" if ($detail_fh);
  close $detail_fh if $detail_fh;

  return $csp_average, \%csp_collection;
}

#------------------------------------------------------------------------------#

sub average_csps {
  my ($metanames, $csp_collection, $repnumber, $mixed, $cross, $detail_fh,
    $basic_fh, $minoverlap, $outstats) = @_;

  if ($outstats) {
    $outstats = 1;
  } else {
    $outstats = 0;
  }

  # Result header
  print $basic_fh "# Circonspect ".$VERSION."\n";
  print $basic_fh "# Lines: sample name, contig spectrum, avg read length, min".
    " overlap specified, #repetitions\n";

  # List of the different contig spectrum names
  my @samples = @$metanames;
  push(@samples, 'mixed') if $mixed;
  if ($cross) {
    for my $sample (@$metanames) {
      push @samples, 'dissolved '.$sample;
    }
    push @samples, 'sum of dissolved', 'cross';
  }

  # Calculate an average for each contig spectrum
  my %avgs;
  print $detail_fh "Average\n" if ($detail_fh && $repnumber > 1);
  for my $sample (@samples) {

    # Get the final contig spectrum
    my $final_csp;
    if ($repnumber == 1) {
      $final_csp = $$csp_collection{$sample}{1};
    } elsif ($repnumber > 1) {
      # Average contig spectra if multiple repetitions
      my @csps;
      for (my $rep = 1 ; $rep <= $repnumber ; $rep++) {
        push @csps, $$csp_collection{$sample}{$rep};
      }
      $final_csp = Bio::Assembly::Tools::ContigSpectrum->new( -eff_asm_params => $outstats );
      $final_csp = $final_csp->average(\@csps);
      $avgs{$sample} = $final_csp;
    }

    # Display final results
    if ($sample eq 'mixed' || $sample eq 'cross') {
      print $basic_fh "$sample\n";
      print $basic_fh csp_result($final_csp, $minoverlap)."\n";
      print $detail_fh '  '.ucfirst($sample)." sample\n" if ($detail_fh && $repnumber > 1);
    } elsif ($sample eq 'sum of dissolved' || $sample =~ m/^dissolved /) {
      print $detail_fh '  '.ucfirst($sample)." sample\n" if ($detail_fh && $repnumber > 1);
    } else {
      print $basic_fh "$sample\n";
      print $basic_fh csp_result($final_csp, $minoverlap)."\n";
      print $detail_fh "  Sample $sample\n" if ($detail_fh && $repnumber > 1);
    }
    print $detail_fh csp_result_verbose($final_csp, $outstats) if ($detail_fh && $repnumber > 1);

  }
  print $detail_fh "\n" if ($detail_fh);

  return \%avgs;
}

#------------------------------------------------------------------------------#

sub process {

    my ( $fastafiles, $compolist, $metanames, $mco, $samplesize,$trim, $asmprog,
         $minidentity, $minoverlap, $mixed, $cross, $detail_fh, $outstats,
         $outace, $outdir, $outprefix, $rep ) = @_;

    my %csps; # save the different contig spectra here
    my @mixed_seqs;
    # Determine standard contig spectrum for each metagenome
    for (my $j = 0 ; $j < scalar @$fastafiles ; $j++) {

      my $mixed_nof_seq = $$compolist[$j];
      my $sample        = $$metanames[$j];
      print $detail_fh "  Sample $sample\n" if ($detail_fh);

      # Take a random sample (seq and score), trim it, put it in an object array
      my $seqs = $mco->get_rand_sample($sample, $samplesize, $trim);

      # Update mixed sample
      if ($mixed || $cross) {
        # Take random subsample of specified size and put it in mixed sample
        my $seq_subset = $mco->rand_subset($seqs, $mixed_nof_seq);
        push @mixed_seqs, @$seq_subset;
      }

      # Assemble it (singletons aren't counted as contigs)
      my $asm_io = assemble($seqs, $asmprog, $minidentity, $minoverlap);
      # Determine contig spectrum
      my $outace_base;
      $outace_base = File::Spec->catfile($outdir, "$outprefix-rep$rep-$sample") if $outace;
      my ($csp) = get_csp($asm_io, $seqs, $outstats, $outace_base, $metanames,
        0, $minidentity, $minoverlap, $detail_fh);

      # Clean sequences
      @$seqs = ();
      undef $seqs;

      # Store and output contig spectrum in a structure
      print $detail_fh csp_result_verbose($csp, $outstats) if ($detail_fh);
      $csps{$sample} = $csp;

    }

    # Take care of advanced contig spectrum types (mixed and cross)
    if ($mixed || $cross) {

      print $detail_fh "  Mixed sample\n" if ($detail_fh);

      # Assemble it (singletons aren't counted as contigs)
      my $mixed_asm_io = assemble(\@mixed_seqs, $asmprog, $minidentity,
        $minoverlap);
      # Calculate mixed contig spectrum
      my ($mixed_csp, $cross_csp, $dissolved_csps, $sum_dissolved_csp);
      if ( not $cross ) {
        # Calculate the mixed contig spectrum
        my $outace_base;
        $outace_base = File::Spec->catfile($outdir, "$outprefix-rep$rep-mixed") if $outace;
        ($mixed_csp) = get_csp($mixed_asm_io, \@mixed_seqs, $outstats,
          $outace_base, $metanames, 0, $minidentity, $minoverlap, $detail_fh);
      } else {
        # Calculate the cross-contig spectrum (mixed, dissolved)
        my $outace_base;
        $outace_base = File::Spec->catfile($outdir, "$outprefix-rep$rep") if $outace;
        ($mixed_csp, $cross_csp, $dissolved_csps, $sum_dissolved_csp) =
          get_csp($mixed_asm_io, \@mixed_seqs, $outstats, $outace_base,
          $metanames, 1, $minidentity, $minoverlap, $detail_fh);
      }

      # Clean mixed sample
      @mixed_seqs = ();
      undef @mixed_seqs;

      # Output and save mixed contig spectrum results
      print $detail_fh csp_result_verbose($mixed_csp, $outstats) if ($detail_fh);
      $csps{'mixed'} = $mixed_csp;

      # Store and display cross-contig spectrum results
      if ($cross) {
        # Dissolved contig spectra
        for my $sample (@$metanames) {
          my $dissolved_csp = $$dissolved_csps{$sample};
          print $detail_fh "  Dissolved $sample\n" if ($detail_fh);
          print $detail_fh csp_result_verbose($dissolved_csp, $outstats) if ($detail_fh);
          $csps{'dissolved '.$sample} = $dissolved_csp;
        }
        # Sum of dissolved contig spectra
        print $detail_fh "  Sum of dissolved samples\n" if ($detail_fh);
        print $detail_fh csp_result_verbose($sum_dissolved_csp, $outstats) if ($detail_fh);
        $csps{'sum of dissolved'} = $sum_dissolved_csp;
        # Cross-contig spectrum
        print $detail_fh "  Cross sample\n" if ($detail_fh);
        print $detail_fh csp_result_verbose($cross_csp, $outstats) if ($detail_fh);
        $csps{'cross'} = $cross_csp;
      }

    } # end of mixed and cross contig spectra

    return \%csps;
}

#------------------------------------------------------------------------------#

sub get_csp {
  # Get the contig spectrum for the current assembly
  # Input:
  #   asm_io:      Bio::Assembly::IO object initialized with the proper parser
  #   seqs:        arrayref of sequences used for the assembly
  #   outstats:    whether to calculated effective assembly parameters or not
  #   outace_base: if set, save assemblies in an ACE file with a name beginning
  #                like specified
  #   metanames:   arrayref of sample names
  #   cross:       do a cross contig spectrum? 1: yes, 0: no
  #   minidentity: minimum overlap identity
  #   minoverlap:  minimum overlap length
  #   details:      whether or not to create ACE files for dissolved and sum
  #                contig spectra
  my ($asm_io, $seqs, $outstats, $outace_base, $metanames, $cross, $minidentity,
    $minoverlap, $details) = @_;

  if ($outstats) {
    $outstats = 1;
  } else {
    $outstats = 0;
  }

  # Initialize contig spectrum objects
  my ($csp, $cross_csp, $dissolved_csps, $sum_dissolved_csp);
  $csp = Bio::Assembly::Tools::ContigSpectrum->new(
    -eff_asm_params => $outstats,
  );
  if ($cross) {
    for my $sample (@$metanames) {
      my $dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
        -eff_asm_params => $outstats
      );
      $dissolved_csps->{$sample} = $dissolved_csp;
    }
    $sum_dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -eff_asm_params => $outstats
    );
    $cross_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -eff_asm_params => $outstats
    );
  }

  # Initialize file handles for the ACE output
  my ($csp_io, $cross_csp_io, $dissolved_csps_ios, $sum_dissolved_csp_io) = 
    initialize_asm($outace_base, $metanames, $cross, $details);

  # Hash that contains ID of singlets: key=sequence ID, value=sequence object
  my %singlet_hash = map { $_->id => $_ } @$seqs;

  # One by one, add contigs first and then singlets to the contig spectra
  my $obj_num = 1;
  while ( my $obj = $asm_io->next_contig || (values %singlet_hash)[0] ) {

    # Convert sequences into singlets and rename contigs/singlets
    if ($obj->isa('Bio::PrimarySeqI') || $obj->isa('Bio::SeqI')) {
      $obj = Bio::Assembly::Singlet->new( -id => $obj_num, -seqref => $obj );
    } else {
      $obj->id($obj_num);
    }
    $obj_num++;
    my $contig = $obj; # singlets inherit from the contig object

    # Remove contig sequences from singlet hash
    for my $seq_id ( $contig->get_seq_ids ) {
      delete $singlet_hash{$seq_id};
    }

    # Add contig (or singlet) to the regular contig spectrum
    my $temp_csp = Bio::Assembly::Tools::ContigSpectrum->new(
      -assembly       => $contig,
      -eff_asm_params => $outstats,
    );
    $csp->add($temp_csp);
    write_asm($csp, $csp_io);
    $csp->drop_assembly();

    # Update advanced contig spectra
    if ($cross) {

      # Dissolved contig spectra
      for my $sample (@$metanames) {
        # Dissolved contig spectrum for 1 contig
        my $temp_dissolved_csp = Bio::Assembly::Tools::ContigSpectrum->new(
          -dissolve       => [$temp_csp, $sample],
          -eff_asm_params => $outstats,
          # specify minimum overlap and minimum identity values for dissolving
          -min_overlap    => $minoverlap,
          -min_identity   => $minidentity,
        );

        if ($details) {
          $dissolved_csps->{$sample}->add($temp_dissolved_csp);
          write_asm($dissolved_csps->{$sample}, $dissolved_csps_ios->{$sample});
          $dissolved_csps->{$sample}->drop_assembly();

          $sum_dissolved_csp->add($temp_dissolved_csp);
          write_asm($sum_dissolved_csp, $sum_dissolved_csp_io);
          $sum_dissolved_csp->drop_assembly();
        }
      }

      # Cross contig spectra
      my $temp_cross_csp = Bio::Assembly::Tools::ContigSpectrum->new(
        -cross          => $temp_csp,
        -eff_asm_params => $outstats,
        # specify minimum overlap and minimum identity values for dissolving
        -min_overlap    => $minoverlap,
        -min_identity   => $minidentity,
      );

      $cross_csp->add($temp_cross_csp);
      write_asm($cross_csp, $cross_csp_io);
      $cross_csp->drop_assembly();

    }

  }
  close_asm($csp_io, $cross_csp_io, $dissolved_csps_ios, $sum_dissolved_csp_io);

  # Number of repetitions is 1
  $csp->nof_rep(1);
  if ($cross) {
    $cross_csp->nof_rep(1);
    for my $sample (@$metanames) {
      $dissolved_csps->{$sample}->nof_rep(1);
    }
    $sum_dissolved_csp->nof_rep(1);
  }

  return $csp, $cross_csp, $dissolved_csps, $sum_dissolved_csp;
}

#------------------------------------------------------------------------------#

sub initialize_asm {
  # Initialize the assembly streams for the different contig spectra
  my ($outace_base, $metanames, $cross, $details) = @_;
  my ($csp_io, $cross_csp_io, $dissolved_csps_ios, $sum_dissolved_csp_io);
  my $fmt = 'ace';
  if ($outace_base) {
    my $csp_file;
    if (not $cross) {
      $csp_file = "$outace_base.$fmt";
    } else {
      $csp_file = "$outace_base-mixed.$fmt";
    }
    $csp_io = Bio::Assembly::IO->new( -format => $fmt, -file => ">$csp_file" );
    if ($cross) {
      if ($details) {
        for my $sample (@$metanames) {
          my $dissolved_csp_file = "$outace_base-dissolved-$sample.$fmt";
          my $dissolved_csp_io = Bio::Assembly::IO->new(-format => $fmt, -file => ">$dissolved_csp_file" );
          $dissolved_csps_ios->{$sample} = $dissolved_csp_io;
        }
        my $sum_dissolved_csp_file = "$outace_base-dissolved-sum.$fmt";
        $sum_dissolved_csp_io = Bio::Assembly::IO->new( -format => $fmt, -file => ">$sum_dissolved_csp_file" );
      }
      my $cross_csp_file = "$outace_base-cross.$fmt";
      $cross_csp_io = Bio::Assembly::IO->new( -format => $fmt, -file => ">$cross_csp_file" );
    }    
  }
  # else we do not initialize the io, so that the ACE files aren't produced
  return $csp_io, $cross_csp_io, $dissolved_csps_ios, $sum_dissolved_csp_io;
}

#------------------------------------------------------------------------------#

sub write_asm {
  # Given a contig spectrum object write its contigs on the provided initialized
  # assembly stream object unless the stream is undef
  my ($csp, $io) = @_;
  if (defined $io) {
    my @contigs = $csp->assembly;
    for my $contig( @contigs ) {
      $io->write_contig($contig);
    }
  }
  return 1;
}

#------------------------------------------------------------------------------#

sub close_asm {
  # For the different contig spectra, print assembly stream header, footer and
  # close io
  my ($csp_io, $cross_csp_io, $dissolved_csps_ios, $sum_dissolved_csp_io) = @_;
  for my $io ($csp_io, $cross_csp_io, values (%$dissolved_csps_ios), $sum_dissolved_csp_io) {
    if (defined $io) {
      $io->write_footer();
      $io->write_header();
      $io->close;
    }
  }
  return 1;
}

#------------------------------------------------------------------------------#

sub params_str {
  # Display command line passed parameters
  my ($fastafiles, $qualfiles, $discard, $trim, $outdir, $outprefix, $outdetails,
    $outstats, $outace, $samplesize, $metapercent, $samplecompo, $compolist,
    $repnumber, $mincoverage, $seed, $asmprog, $minidentity, $minoverlap,
    $mixed, $cross) = @_;

  my $size = scalar(@$fastafiles);
  my $s1 = '';
  my $s2 = '  ';
  my $str = '';

  # Input params
  $str .= "Input:\n";
  $str .= $s2."Fasta files: @$fastafiles\n";
  if ($qualfiles) {
    $str .= "Quality files: ";
    for (my $i = 0 ; $i < $size ; $i++) {
      if ($$qualfiles[$i]) { $str .= "$$qualfiles[$i] " };
    }
    $str .= "\n";
  }

  # Sequence preprocessing params
  if ($discard || $trim) {$str .= $s1."Sequence processing:\n"};
  if ($discard) {$str .= $s1.$s2."Discard size: $discard bp\n"};
  if ($trim) {$str .= $s1.$s2."Trim size: $trim bp\n"};
  
  # Random sampling params
  if ($samplesize || $metapercent || $repnumber || $mincoverage) {$str .= $s1.
    "Random sampling:\n"};
  if ($metapercent) {$str .= $s1.$s2."Sample size: $metapercent % of the ".
    "smallest metagenome\n"};
  if ($samplesize) {$str .= $s1.$s2."Sample size: $samplesize sequences\n"};
  if ($repnumber) {$str .= $s1.$s2."Number of repetitions: $repnumber\n"};
  if ($mincoverage) {$str .= $s1.$s2."Minimum coverage: $mincoverage x the ".
    "largest metagenome\n"};
  if ($seed) {$str .= $s1.$s2."Seed number: $seed\n"};
  
  # Assembly params
  if ($asmprog || $minidentity || $minoverlap) {$str .= $s1."Assembly:\n"};
  if ($asmprog) {$str .= $s1.$s2."Assembly program: $asmprog\n"};
  if ($minidentity) {$str .= $s1.$s2."Minimum sequence identity: $minidentity %".
    "\n"};
  if ($minoverlap) {$str .= $s1.$s2."Minimum sequence overlap: $minoverlap bp\n"};
  
  # Advanced contig spectra
  if ($mixed || $cross) {$str .= $s1."Advanced contig spectrum types:\n"};
  if ($mixed || $cross) {$str .= $s1.$s2."Mixed contig spectrum\n"};
  if ($cross) {$str .= $s1.$s2."Cross-contig spectrum\n"};
  if ($samplecompo) {$str .= $s1.$s2."Sample composition: @$samplecompo % from ".
    "each metagenome\n"};
  if ($samplecompo) {$str .= $s1.$s2."Sample composition: @$compolist sequences".
    "\n"};
  
  # Output params
  if ($outdir || $outprefix || $outdetails || $outstats || $outace) {
    $str .= $s1."Output:\n"
  };
  if ($outdir) {$str .= $s1.$s2."Directory: $outdir\n"};
  if ($outprefix) {$str .= $s1.$s2."Prefix: $outprefix\n"};
  if ($outdetails || $outstats || $outace) {$str .= $s1.$s2."Additional output:\n"};
  if ($outdetails) {$str .= $s1.$s2.$s2."Detailed results\n"};
  if ($outstats) {$str .= $s1.$s2.$s2."Contig statistics\n"};
  if ($outace) {$str .= $s1.$s2.$s2."Assembly ACE files\n"};

  # Function ends with success
  $str .= "\n";
  return $str;
}

#------------------------------------------------------------------------------#

sub basic_param_check {
  # Make sampling checks that don't require knowing the size of the metagenomes
  # Set the parameters that need to be set
  my ($discard, $trim, $asmprog, $samplecompo, $fastafiles, $mixed, $cross,
    $outdir, $outprefix, $outstats, $outace) = @_;

  # No mixed or cross-contigs if there is only one metagenome
  my $size = scalar(@{$fastafiles});
  if ($size == 1) {
    $mixed = 0;
    $cross = 0;
  }
 
  # Require a mixed contig spectrum when a cross-contig spectrum is required
  if ($cross) {
    $mixed = 1;
  }

  # Minimum sequence length supported by the assembler
  my $min_seq_length;
  if ($asmprog eq 'TigrAssembler') {
    $min_seq_length = 39;
  } elsif ($asmprog eq 'LigrAssembler') {
    $min_seq_length = 39;
  } elsif ($asmprog eq 'Cap3') {
  } elsif ($asmprog eq 'Newbler') {
  } elsif ($asmprog eq 'Minimo') {
  }
  if ( $min_seq_length && ( (not $discard) || ( $min_seq_length > $discard) ) ) {
    $discard = $min_seq_length;
    warn "Warning: the minimum sequence length supported by $asmprog is $min_seq_length";
  }

  # Check that trimming size is larger than discard size
  if ($trim && $discard && ($trim < $discard)) {
    die("Error: the specified trim size ($trim) is smaller than the discard ".
      "size ($discard)\n");
  }

  # Check the mixed/crossed sample composition (list of percent of sequences
  # from each metagenome)
  if ($mixed || $cross) {
    if ($samplecompo) {
      # Sample composition was specified
       # Verify that there are as many numbers as there are metagenomes to process
       unless (scalar(@{$samplecompo}) == $size) {
         die("Error: exactly one percent value must be specified for each input".
         " metagenome in the sample composition\n");
       }

       # Verify that the sum roughly adds up to 100%
       my $total = 0;
       for my $percent(@{$samplecompo}) {
         $total = $total + $percent;
       }

       unless ($total == 100) {
         # If the total is not 100
         if ($total >= 99 && $total <= 101) {
           # If the total is between 99 and 101%, we just adjust it to 100
           for (my $i = 0 ; $i < $size ; $i++) {
             $$samplecompo[$i] = $$samplecompo[$i]*100/$total; 
           }
         } else { 
           # The total is off, report an error
           die("Error: the total of the percent value for the sample composition".
             " should equal 100%\n");
         }
       }
     } else {
        # Sample composition not specified: set an equal number of sequences for
        # each metagenome
        my $value = 100 / $size;
        for (my $i = 0 ; $i < $size ; $i++) {
          $$samplecompo[$i] = $value; 
        }
      }
  }

  return ($mixed, $cross, $asmprog, $samplecompo, $outdir, $outprefix, $discard);
}

#------------------------------------------------------------------------------#

sub extended_param_check {
  # Some additional parameter validity check and compute missing parameters
  # necessary before computation begins
  my ($samplesize, $metapercent, $repnumber, $mincoverage, $samplecompo, $mixed,
    $cross, $fastafiles, $stats) = @_;

  my @compolist;
  # metagenome with smallest number of sequences
  my $minseqmeta = $$stats{'minseq'};
  # number of sequences in the metagenome with smallest number of sequences
  my $minseqsize = $$stats{$minseqmeta}{'numseq'};
  # metagenome with largest number of bp
  my $trimmedmaxbpmeta = $$stats{'trimmedmaxbp'};
  # number of bp in the metagenome with largest number of bp
  my $trimedmaxbpsize = $$stats{$trimmedmaxbpmeta}{'trimmednumbp'};

  # Case 1: repnumber and (metapercent or samplesize) given, compute mincoverage
  if ($repnumber && ($metapercent || $samplesize)) {

    # Compute samplesize based on metapercent
    if ($metapercent) {
      $samplesize = int($minseqsize * $metapercent / 100);
      # Set up samplesize to 1 if equal to 0
      if ($samplesize == 0) {
        $samplesize = 1;
        warn("Warning: since the sample size was but cannot be 0, it has been ".
          "set to 1 sequence. \n");
      }
    }
    # Compute metapercent (informative only)
    else {
      $metapercent = $samplesize / $minseqsize * 100;
    }

    # Compute compolist
    if ($mixed || $cross) {
      my $total = 0;
      for (my $i = 0 ; $i < scalar @$fastafiles ; $i++) {
        # Not rounding to prevent samplesize to exceed metagenome size when
        # doing a control cross- contig spectrum
        $compolist[$i] = int ( $samplesize * $$samplecompo[$i] / 100 );
        # Set up samplesize to 1 if equal to 0
        if ($compolist[$i] == 0) {
          $compolist[$i] = 1;
          warn("Warning: since the sample size in the sample composition was ".
            "but cannot be 0, it has been set to 1 sequence.\n");
        }
        # Set up samplesize to the metagenome size if it is over
        if ($compolist[$i] > $$stats{$$fastafiles[$i]}{'numseq'}) {
          $compolist[$i] = $$stats{$$fastafiles[$i]}{'numseq'};
        }
        $total += $compolist[$i];
      }
      $samplesize = $total;
    }

    # Compute mincoverage based on samplesize and repnumber (informative only)
    for (my $i = 0 ; $i < scalar @$fastafiles ; $i++) {
      my $fasta = $$fastafiles[$i];
      my $bpsize = $$stats{$fasta}{'numbp'};
      my $avgseq = $$stats{$fasta}{'trimmedseqbp'};
      my $newcoverage = $samplesize * $repnumber * $avgseq / $bpsize;
      if ($mixed || $cross) {
        my $percent = $$samplecompo[$i];
        $newcoverage = $newcoverage * $percent / 100;
      }
      if ( not($mincoverage) || $newcoverage < $mincoverage) {
        $mincoverage = $newcoverage;
      }
    }

  }

  # Case 2: repnumber and mincoverage given, compute samplesize
  elsif ($repnumber && $mincoverage) {
    # Compute samplesize
    $samplesize = 0;
    for (my $i = 0 ; $i < scalar @$fastafiles ; $i++) {
      my $fasta = $$fastafiles[$i];
      my $bpsize = $$stats{$fasta}{'numbp'};
      my $avgseq = $$stats{$fasta}{'trimmedseqbp'};
      my $newsamplesize = $mincoverage * $bpsize / $repnumber / $avgseq;
      if ($mixed || $cross) {
        my $percent = $$samplecompo[$i];
        $newsamplesize = $newsamplesize / $percent * 100;
      }
      if ($newsamplesize > $samplesize) {
        $samplesize = int( $newsamplesize + 0.5 );
      }
    }
    if ($mixed || $cross) {
      for (my $i = 0 ; $i < scalar @$samplecompo ; $i++) {
        my $percent = $$samplecompo[$i];
        $compolist[$i] = int( $samplesize * $percent /100 + 0.5);
      }
    }
    # Compute metapercent (informative only)
    $metapercent = $samplesize / $minseqsize * 100;
  }

  # Case 3: (metapercent or samplesize) and mincoverage given, compute repnumber
  elsif  (($samplesize || $metapercent) && $mincoverage) {
    # Compute samplesize based on metapercent
    if ($metapercent) {
      $samplesize = int($minseqsize * $metapercent / 100);
       # Set up samplesize to 1 if equal to 0
      if ($samplesize == 0) {
        $samplesize = 1;
        warn("Warning: since the sample size was $samplesize but cannot be 0, ".
          "it has been set to 1 sequence. \n");
      }
    }
    # Compute metapercent (informative only)
    else {
      $metapercent = $samplesize / $minseqsize * 100;
    }
    # Compute compolist
    if ($mixed || $cross) {
      my $total = 0;
      for (my $i = 0 ; $i < scalar @$fastafiles ; $i++) {
        $compolist[$i] = int ( $samplesize * $$samplecompo[$i] / 100 + 0.5 );
        if ($compolist[$i] == 0) {
          $compolist[$i] = 1;
          warn("Warning: since the sample size in the sample composition was ".
            "but cannot be 0, it has been set to 1 sequence. \n");
        }
        $total += $compolist[$i];
      }
      $samplesize = $total;
    }
    # Compute repnumber based on samplesize and mincoverage
    $repnumber = 0;
    for (my $i = 0 ; $i < scalar @$fastafiles ; $i++) {
      my $fasta = $$fastafiles[$i];
      my $bpsize = $$stats{$fasta}{'numbp'};
      my $avgseq = $$stats{$fasta}{'trimmedseqbp'};
      my $newrepnumber = $mincoverage * $bpsize / $samplesize / $avgseq;
      if ($mixed || $cross) {
        my $percent = $$samplecompo[$i];
        $newrepnumber = $newrepnumber * 100 / $percent;
      }
      if ($newrepnumber > $repnumber) {
        $repnumber = int( $newrepnumber + 0.5 );
      }
    }
  } else {
    die("Error: how did you get here?\n");
  }

  # Verify that samplesize does not exceed the size of the smallest metagenome
  if ($samplesize > $minseqsize) {
    die "Error: $samplesize sequences exceeds the size of the smallest ".
      "metagenome \"$minseqmeta\", $minseqsize sequences. Try decreasing the ".
      "sample size (-s or -p).\n";
  }

  return ($samplesize, $metapercent, $repnumber, $mincoverage, \@compolist);
}

#------------------------------------------------------------------------------#

sub assemble {
  # Assemble sequences. Input:
  #   seqs:        arrayref of sequences (can be Bio::Seq::Quality)
  #   asmprog:     the string of a valid bioperl-run wrapper to use for the assembler
  #   minidentity: the minimum overlap length to form contigs
  #   minoverlap:  the minimum overlap identity percentage to form contigs
  my ($seqs, $asmprog, $minidentity, $minoverlap) = @_;

  # Module name
  my $assembly_wrapper = 'Bio::Tools::Run::';
  if ($asmprog eq 'LigrAssembler') {
    # LIGR assembler runs through the TIGR assembler wrapper
    $assembly_wrapper .= 'TigrAssembler';
  } else {
    $assembly_wrapper .= $asmprog;
  }

  # Initialize module
  my $factory = Bio::Root::Root->new();
  eval { $factory->_load_module($assembly_wrapper); };
  if ($@) { die "Error: Could not load BioPerl wrapper for $asmprog.\n$@\n" };
  $factory = $assembly_wrapper->new();

  # Set the assembly parameters. Eval it because not all assemblers support
  # the minimum overlap and similarity parameters
  eval {
    $factory->minimum_overlap_length($minoverlap);
    $factory->minimum_overlap_similarity($minidentity);
  };
  if ($@) {
    die("Error: The $asmprog BioPerl wrapper does not support setting mininum".
      "overlap and similarity values.\n$@\n");
  }

  # Assembler-specific options
  if ($asmprog eq 'Minimo') {
    $factory->program_name('Minimo');
    $factory->good_qual(30); # quality score for good residues
    $factory->bad_qual(30);  # quality score for bad residues
    $factory->aln_wiggle(2); # consensus alignment wiggle value
  } elsif ($asmprog eq 'LigrAssembler') {
    $factory->program_name('LIGR_Assembler');
    $factory->use_tandem_32mers(1); # do not discard sequences with tandem 32mers
    $factory->maximum_end(1);       # maximum overhang
    $factory->incl_bad_seq(1);      # keep potential chimera and splice variant reads
    $factory->trimmed_seq(1);       # reads are high quality over their full length
    #$factory->max_err_32(32);      # max number + 1 of mismatches in a 32mer
  } elsif ($asmprog eq 'TigrAssembler') {
    $factory->program_name('TIGR_Assembler');
    $factory->use_tandem_32mers(1); # do not discard sequences with tandem 32mers
    $factory->maximum_end(1);       # maximum overhang
    #$factory->max_err_32(32);      # maximum number + 1 of mismatches in a 32 bp stretch
  } elsif ($asmprog eq 'Cap3') {
    $factory->program_name('cap3');
    $factory->clipping_quality_cutoff(6); # a quality of 6 is the minimum
    $factory->clipping_range(6);          # 6bp is the minimum
    $factory->max_gap_length(10000000);   # large gap length
    $factory->max_overhang_percent(3);    # 3% is the minimum
    $factory->overlap_score_cutoff(401);  # 401 is the minimum
  } elsif ($asmprog eq 'Newbler') {
    $factory->program_name('runAssembly');
    $factory->ace_raw(1);         # output the full "raw" read sequence
    $factory->expected_depth(0);  # 0 means to not use expected depth information
    $factory->no_duplicates(1);   # do not identify duplicate reads
    #$factory->large(1);          # use a faster but less precise heuristic
    #$factory->in_memory(1);      # keep sequence data in memory
    #$factory->no_auto_rescore();
    #$factory->seed_count();
    #$factory->seed_length();
    #$factory->seed_step();
  }

  # Quality scores
  my $quals = undef;
  if ($$seqs[0]->isa('Bio::Seq::Quality')) {
    # sequence and quality scores are contained in the same object
    $quals = $seqs;
  }

  # Run assembler
  ###############
  #$factory->verbose(2);
  #$factory->save_tempfiles(1);
  ###############
  $factory->out_type('Bio::Assembly::IO'); # output will be an assembly IO
  my $asm_io = $factory->run($seqs, $quals);

  return $asm_io;
}

#------------------------------------------------------------------------------#

sub max {
  # Takes the reference of an array containing positive numbers
  # Returns the largest number
  my ($arr) = @_;
  my $max = 0;
  for my $value (@$arr) {
    $max = $value if $value > $max;
  }
  return $max;
}

#------------------------------------------------------------------------------#

sub csp_result_verbose {
  # String representation of the results for a given contig_spectrum
  my ($csp_obj, $outstats) = @_;
  my $string = '    Contig spectrum         '.$csp_obj->to_string."\n";
  $string .=   '    Largest contig size     '.$csp_obj->max_size." sequences\n";
  $string .=   '    Number of sequences     '.$csp_obj->nof_seq." sequences\n";
  $string .=   '    Average sequence length '.$csp_obj->avg_seq_len." bp\n";
  if ($outstats) {
    $string .=   '    Minimum overlap length  ';
    if ($csp_obj->min_overlap) {
      $string .= $csp_obj->min_overlap;
    } else {
      $string .= '0';
    }
    $string .=   " bp\n";
    $string .=   '    Average overlap length  '.$csp_obj->avg_overlap." bp\n";
    $string .=   '    Minimum overlap match   ';
    if ($csp_obj->min_identity) {
      $string .= $csp_obj->min_identity;
    } else {
      $string .= '0';
    }
    $string .=   " %\n";
    $string .=   '    Average overlap match   '.$csp_obj->avg_identity." %\n";
  }
  return $string;
}

#------------------------------------------------------------------------------#

sub csp_result {
  # String representation of the results for a given contig_spectrum
  my ($csp_obj, $minoverlap) = @_;
  my $string = $csp_obj->to_string."\n";
  $string   .= $csp_obj->avg_seq_len."\n";
  $string   .= $minoverlap."\n";
  $string   .= $csp_obj->nof_rep."\n";
  return $string;
}

#------------------------------------------------------------------------------#

sub getqual {
  # Associate a QUAL file to each FASTA file
  # QUAL files are supposed to have the same basename as the corresponding FASTA
  # file. It is allowed for a FASTA file to have no QUAL file
  my ($fastafiles, $qualfiles) = @_;
  my $qualsize = scalar(@$qualfiles);
  my $fastasize = scalar(@$fastafiles);
  my @tmp = @$qualfiles;
  @$qualfiles = ();
  for (my $i = 0 ; $i < $qualsize ; $i++) {
    my @qualparse = fileparse($tmp[$i], '\..*');
    my $qualbasename = $qualparse[0];
    my $found = 0;
    for (my $j = 0 ; $j < $fastasize ; $j++) {
      my @fastaparse = fileparse(@$fastafiles[$j], '\..*');
      my $fastabasename = $fastaparse[0];
      if ($fastabasename eq $qualbasename) {
        $$qualfiles[$j] = $tmp[$i];
        $found = 1;
      }
    }
    if ($found == 0) {
      warn("Warning: no FASTA file seems to match the QUAL file \"$tmp[$i]\". ".
        "This QUAL file will be ignored\n");
    }
  }
  return $qualfiles;
}

#------------------------------------------------------------------------------#


1;
