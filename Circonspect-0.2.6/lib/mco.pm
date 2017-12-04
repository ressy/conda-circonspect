# Circonspect, copyright 2009 Florent Angly <florent.angly@gmail.com>, under
# the GNU GPLv3 license. Circonspect generates contig spectra by boostrapping
# assemblies of random subsets of metagenomes.

=head1  NAME

 The mco module implements a metagenome collection object.

=head1 NOTES

 Provides methods for the creation, manipulation, and destruction of a
 collection of metagenomes:
    * new
    * import_metainfo
    * import_sequences
    * get_size_stats
    * get_rand_sample
    * destroy

=head1 AUTHOR

 Florent Angly (florent.angly@gmail.com)

=cut

package mco;

use strict;
use warnings;
use Cwd;
use File::Spec;
use File::Basename;
# BioPerl modules
use Bio::SeqIO;
use Bio::DB::Fasta;
use Bio::DB::Qual;
use Bio::PrimarySeq;
use Bio::Seq::Quality;
use Bio::Root::Utilities;
use Math::Random::MT qw(srand rand);


=head2 new
  Title   : new
  Usage   : 
  Function: 
  Args    : 
  Returns : 
=cut

sub new {
  my ($class) = @_;
  my $self = bless{},$class;
  $self->{'metagenome'} = {}; # names and files of metagenomes
  #$self->{'metagenome'}{$metaname}{'fastafile'}
  #$self->{'metagenome'}{$metaname}{'qualfile'}
  $self->{'original_sequences'} = {}; #
  #$self->{'original_sequences'}{$fasta}{'seqdbfile'}
  #$self->{'original_sequences'}{$fasta}{'seq'}
  #$self->{'original_sequences'}{$fasta}{'qualdbfile'}
  #$self->{'original_sequences'}{$fasta}{'qual'}
  return $self;
}


=head2 import_metainfo
  Title   : import_metainfo
  Usage   : 
  Function: Give each metagenome a unique name, link FASTA and QUAL to a new
            location and register that for later use 
  Args    : 
  Returns : 
=cut

sub import_metainfo {
  my ($self, $fastafiles, $qualfiles, $outdir, $outprefix) = @_;
  my $metanames = $self->_create_unique_names($fastafiles);
  my %done;
  for (my $i = 0 ; $i < scalar(@$metanames) ; $i++) {
    my $metaname = $$metanames[$i];
    # Transform filenames to complete path
    my $fasta = Cwd::abs_path($$fastafiles[$i]);
    my $fastabasename = (fileparse($fasta, '\..*'))[0];
    my $new_fasta = File::Spec->catfile($outdir, $outprefix.'-'.$fastabasename.'.fa');
    my $qual;
    my $new_qual;
    if (defined $$qualfiles[$i]) {
      $qual = Cwd::abs_path($$qualfiles[$i]);
      my $qualbasename = (fileparse($qual, '\..*'))[0];
      $new_qual = File::Spec->catfile($outdir, $outprefix.'-'.$qualbasename.'.qual');
    }
    # Process non-already treated files
    if (not exists $done{$new_fasta}) {
      $self->link_uncompress_file($fasta, $new_fasta);
      if (defined $qual) {
        $self->link_uncompress_file($qual, $new_qual);
      }
      $done{$new_fasta} = undef;
    }
    # Save filenames for future use
    $$fastafiles[$i] = $new_fasta;
    $$qualfiles[$i] = $new_qual;
    $self->{'metagenome'}{$metaname}{'fastafile'} = $new_fasta;
    $self->{'metagenome'}{$metaname}{'qualfile'} = $new_qual;
  }
  return $metanames;
}


=head2 link_uncompress_file
  Title   : Uncompress a file to a given destination. If file is not in a
            compressed format, a link is done. Supports: .gz .Z .bz2 .zip
  Usage   : 
  Function: 
  Args    : 
  Returns : 
=cut

sub link_uncompress_file {
  my ($self, $in, $out) = @_;
  if ( $in =~ /^(.*)(\.gz|\.Z|\.bz2|\.zip)$/ ) {
    my $util = Bio::Root::Utilities->new();
    $util->uncompress(-file => $in, -outfile => $out, tmp => 1);
  } else {
    link $in, $out or die "Error: could not link $in to $out\n$!\n";
  }
  return 1;
}


=head2 import_sequences
  Title   : import_sequences
  Usage   : 
  Function: 
  Args    : 
  Returns : 
=cut

sub import_sequences {
  # Import original sequences
  my ($self, $fastafiles, $qualfiles, $discard, $cplxthres, $cplxwindow,
    $outdir, $outprefix) = @_;
  my $report = '';

  if ($discard) {
    $report .= "Keeping only sequences at least $discard bp long.\n";
  }

  my $nof_files = scalar @$fastafiles;
  my %stats;
  my %good_seqs;
  for (my $i = 0 ; $i < $nof_files ; $i++) {
    my $fasta = $$fastafiles[$i];

    # Skip already imported metagenomes
    next if exists $self->{'original_sequences'}{$fasta}{'seq'};

    my $qual = undef;
    $qual = $$qualfiles[$i] if ($qualfiles && $$qualfiles[$i]);

    $report .= "Importing $fasta";
    $report .= " and $qual" if $qual;
    $report .= "\n";

    # FASTA import
    my %sample_good_seqs;
    my $fastabasename = (fileparse($fasta, '\..*'))[0];
    my $fastadbfile = File::Spec->catfile($outdir, $outprefix.'-'.
      $fastabasename.'.fa');
    my $infasta  = Bio::SeqIO->new( -file => "$fasta", -format => 'fasta' );
    while ( my $seq_obj = $infasta->next_seq() ) {
      # Ignore sequences too small
      my $seqid  = $seq_obj->id;
      my $seqlen = $seq_obj->length;
      if ( $discard && ($seqlen < $discard) ) {
        next;
      }
      # Check for invalid sequences
      if (not $seq_obj->validate_seq()) {
        $seq_obj->throw("Cannot continue with invalid sequence $seqid\n");
      }
      # Filter out low complexity sequences/regions
      my $ranges;
      if ( not $cplxthres ) {
        # No complexity filtering
        $ranges = [ [ 1, $seqlen] ];
      } else {
        # Low complexity filtering
        if (not $cplxwindow) {
          # Filter out complete low complexity sequence
          $ranges = $self->_filter_low_cplx( $seq_obj, $cplxthres, $cplxwindow );
        } else {
          # Filter out low complexity regions, not full sequence
          $ranges = $self->_filter_low_cplx( $seq_obj, $cplxthres, $cplxwindow );
          # Remove subsequences that are not large enough
          next if not defined $ranges;
          for (my $i= 0; $i <= scalar @$ranges - 1; $i++) {
            my $range = $$ranges[$i];
            my $length = $$range[1] - $$range[0] + 1;
            if ( $discard && $range && ( $length < $discard) ) {
              splice @$ranges, $i, 1;
              $i--;
            }
          }
          if (scalar @$ranges == 0) {
            $ranges = undef;
          }
        }
      }
      if (defined $ranges) {
        $sample_good_seqs{$seqid} = $ranges;
      }
    }
    $infasta->close;
    $self->{'original_sequences'}{$fasta}{'seqdbfile'} = $fasta;
    $self->{'original_sequences'}{$fasta}{'seq'} = Bio::DB::Fasta->new(
      $fasta, -reindex => 1);

    # Exit if no sequences were imported
    my $nof_seqs = scalar keys %sample_good_seqs;
    if ($nof_seqs == 0) {
      die "Error: No sequences were imported. You can try to decrease the ".
        "discard size (-u) and low complexity filtering (-ht and -hw)\n";
    }

    # QUAL import
    if ($qual) {
      my $qualbasename = (fileparse($qual, '\..*'))[0];
      my $qualdbfile = File::Spec->catfile($outdir, $outprefix.'-'.
        $qualbasename.'.qual'); 
      my $inqual  = Bio::SeqIO->new( -file => "$qual", -format => 'qual' );
      my $nof_quals = 0;
      while ( my $qual_obj = $inqual->next_seq() ) {
        my $qualid = $qual_obj->id;
        # Ignore qual scores w/o sequence or with length too small
        next unless exists $sample_good_seqs{$qualid};
        # Check for invalid quality scores:
        # done automatically when creating quality score object. dies on failure
        $nof_quals++;
      }
      $inqual->close;
      my $nof_missing_quals = $nof_seqs - $nof_quals;
      if ($nof_missing_quals) {
        warn("Warning: $nof_missing_quals sequence(s) had missing or invalid ".
          "quality score information. This may confuse some assembly programs.".
          "\n");
      }
      $self->{'original_sequences'}{$fasta}{'qualdbfile'} = $qual;
      $self->{'original_sequences'}{$fasta}{'qual'} =
        Bio::DB::Qual->new( $qual, -reindex => 1);
    }

    # IDs of sequences that are ok to use
    $self->{'original_sequences'}{$fasta}{'goodseqids'} = \%sample_good_seqs;
  }

  return $report;
}

=head2 get_size_stats
  Title   : get_size_stats
  Usage   : 
  Function: 
  Args    : 
  Returns : 
=cut

sub get_size_stats {
  # Some statistics about all metagenomes
  #   - metagenome size in sequences
  #   - metagenome size in bp after discarding
  #   - average sequence length after discarding+trimming
  # Note about trimming/discarding:
  #   - the sequences smaller than the discarding size were not imported
  #   - the longer sequences than desired trimming size are not trimmed here
  #     but the reported size behaves as if the trimming had occured
  my ($self, $fastafiles, $trim) = @_;
  my $report = '';
  my %stats;

  my $trimmedmaxbp = $$fastafiles[0]; # name of largest metagenome (bp)
  my $minseq = $$fastafiles[0]; # name of smallest metagenome (#seq)

  # Get stats for each metagenome
  my %done;
  for ( my $i = 0 ; $i < scalar @$fastafiles ; $i++ ) {
    my $fasta = $$fastafiles[$i];
    # skip already done fasta files
    next if exists $done{$fasta};

    # get the stats
    my $numseq = 0;
    my $numbp  = 0;
    my $seqbp  = 0;
    my $trimmednumbp = 0;
    my $trimmedseqbp = 0;
    my $seqdb  = $self->{'original_sequences'}{$fasta}{'seq'};
    my %ranges = %{$self->{'original_sequences'}{$fasta}{'goodseqids'}};
    my @seqids = keys %ranges;
    $numseq = scalar(@seqids);
    for my $seqid (@seqids) {
      for my $range (@{$ranges{$seqid}}) {
        my $length = $$range[1] - $$range[0] + 1;
        $numbp += $length;
        if ($trim) {
          my $trimmedlength;
          if ($length <= $trim) {
            $trimmedlength = $length;
          } else {
            $trimmedlength = $trim;
          }
          $trimmednumbp += $trimmedlength;
        }
      }
    }
    if (not $trim) {
      $trimmednumbp = $numbp;
    }

    $seqbp = $numbp / $numseq;
    if ($trim) {
      $trimmedseqbp = $trimmednumbp / $numseq;
    } else {
      $trimmedseqbp = $seqbp;
    }

    $stats{$fasta}{'numseq'} = $numseq;
    $stats{$fasta}{'numbp'}  = $numbp;
    $stats{$fasta}{'seqbp'}  = $seqbp;
    $stats{$fasta}{'trimmednumbp'} = $trimmednumbp;
    $stats{$fasta}{'trimmedseqbp'} = $trimmedseqbp;

    # Update smallest/largest metagenome
    if ($stats{$fasta}{'trimmednumbp'} > $stats{$trimmedmaxbp}{'trimmednumbp'}) {
      $trimmedmaxbp = $fasta;
    }
    if ($stats{$fasta}{'numseq'} < $stats{$minseq}{'numseq'}) {
      $minseq = $fasta;
    }
    $done{$fasta} = undef; # create new hash entry (key without value)

    # report stats
    $report .= $fasta.", ".$stats{$fasta}{'numseq'}." sequences, ".
      $stats{$fasta}{'numbp'}." bp, ".$stats{$fasta}{'seqbp'}." bp/sequence\n";
  }

  # Update largest/smallest metagenome
  $stats{'minseq'} = $minseq;
  $stats{'trimmedmaxbp'} = $trimmedmaxbp;

  return \%stats, $report;
}


=head2 get_rand_sample
  Title   : get_rand_sample
  Usage   : 
  Function: 
  Args    : 
  Returns : 
=cut

sub get_rand_sample {
  # Get a random sample, trim it, and save it in an array of BioPerl sequences
  # and quality scores. Calling that function many times, so there are some time
  # optimizations. A list of the sequence ID to use can be passed. If none is
  # passed, then all the sequences in the database will be considered.
  my ($self, $samplename, $samplesize, $trim) = @_;
  my $fasta  = $self->{'metagenome'}{$samplename}{'fastafile'};
  my $seqdb  = $self->{'original_sequences'}{$fasta}{'seq'};
  my $qualdb = $self->{'original_sequences'}{$fasta}{'qual'}; # can be undef

  # Get random sequence IDs for this sample
  my $ranges = $self->{'original_sequences'}{$fasta}{'goodseqids'};
  my $idxs = [keys %{$ranges}];
  my $rand_idxs = $self->rand_subset($idxs, $samplesize);

  # Get sequence and score, trim them
  my @seqs;
  for my $seqname (@$rand_idxs) {

    # Retrieve sequence, qual, and ID from database
    my $newid = $samplename.'|'.$seqname;

    # A sequence can have several regions of high complexity. Pick one randomly
    my $seq_regions = $$ranges{$seqname};
    my $nof_regions = scalar @$seq_regions;
    my $rand_region_idx = int(rand($nof_regions));
    my $rand_region_range = $$seq_regions[$rand_region_idx];

    # Now pick a subset of that region if trimming was required
    my $start_pos = $$rand_region_range[0];
    my $stop_pos  = $$rand_region_range[1];

    if ( $trim ) {
      my $length = $stop_pos - $start_pos + 1;
      $start_pos += int( rand( $length - $trim + 1 ) );
      $stop_pos  = $start_pos + $trim - 1;
    }

    my $seqstr  = $seqdb->seq($seqname, $start_pos => $stop_pos);
    my $qualarr = $qualdb->qual($seqname, $start_pos => $stop_pos) if $qualdb;

    # Create BioPerl object
    my $seq_obj;
    if (not $qualdb) {
      $seq_obj = Bio::PrimarySeq->new( -seq => $seqstr,
                                       -id  => $newid  );
    } else {
      $seq_obj = Bio::Seq::Quality->new( -seq  => $seqstr,
                                         -id   => $newid,
                                         -qual => $qualarr );
    }

    # Save sequence
    push @seqs, $seq_obj;

  }
  # Return arrayref of sequences and quals
  return \@seqs; #success
}


=head2 destroy
  Title   : destroy
  Usage   : 
  Function: 
  Args    : 
  Returns : 
=cut

sub destroy {
  # Destroy circonspect object
  my ($self) = @_;
  # Remove metagenome sequences
  for my $fasta (keys %{$self->{'original_sequences'}}) {
    my $seqdbfile =  $self->{'original_sequences'}{$fasta}{'seqdbfile'};
    unlink $seqdbfile || die("Error: could not delete sequence database file ".
      "'$seqdbfile'\n");
    my $seqindex  = ($self->{'original_sequences'}{$fasta}{'seq'})->index_name(
      $seqdbfile, 0);
    unlink $seqindex  || die("Error: could not delete sequence database index ".
      "file \"$seqindex\"\n");
    if (defined $self->{'original_sequences'}{$fasta}{'qualdbfile'}) {
      my $qualdbfile =  $self->{'original_sequences'}{$fasta}{'qualdbfile'};
      unlink $qualdbfile || die("Error: could not delete quality database file ".
        "'$qualdbfile'\n");
      my $qualindex   = ($self->{'original_sequences'}{$fasta}{'qual'})->index_name($qualdbfile, 0);
      unlink $qualindex  || die("Error: could not delete quality database ".
        "index file \"$qualindex\"\n");
    }
  }
  delete $self->{'original_sequences'};
  # Remove metagenome information
  delete $self->{'metagenome'};
  return;
}


=head2 rand_subset
  Title   : rand_subset
  Usage   : $subset_arrayref = $self->rand_subset(\@array, 10);
  Function: Take a given number of elements from a given array at random. It
            uses the Fisher-Yates algorithm to shuffle the tail of the array and
            returns only that shuffled tail.
  Args    : arrayref
            number of elements (integer)
  Returns : arrayref
=cut

sub rand_subset {
  my ($self, $array, $n) = @_;
  my $max_offset = $#$array;
  my $stop_offset = $max_offset - $n + 1;
  for (my $i = $max_offset; $i >= $stop_offset ; $i--) {
    my $j = int rand($i+1);
    next if $i == $j;
    @$array[$i,$j] = @$array[$j,$i];
  }
  my @rand_subset = @$array[$stop_offset .. $max_offset];
  return \@rand_subset;
}


#
# Internal functions
#

=head2 _trim_obj
  Title   : _trim_obj
  Usage   : $trimmed_seq = $self->_trim($trim_l, $seq);
  Function: Trim a DNA sequence and quality scores to a specified length
            starting at random position in the sequence
  Args    : - trimming length
            - Bio::PrimarySeqI or Bio::SeqI sequence
  Returns : - truncated Bio::PrimarySeqI or Bio::SeqI sequence
=cut

sub _trim_obj { 
  my ($self, $trim, $seq_obj) = @_;
  my $length = $seq_obj->length;
  if ($trim && $length > $trim) {
    my $start = int( rand( $length - $trim + 1 ) ) + 1;
    my $stop  = $start + $trim - 1;
    $seq_obj  = $seq_obj->trunc($start, $stop);
  }
  return $seq_obj;
}


=head2 _trim_str
  Title   : _trim_str
  Usage   : ($trimmed_seq, $trimmed_qual) = $self->_trim_str($trim_l, $seq_str, $qual_str);
  Function: Trim a DNA sequence string and quality scores arrayred to a specified length
            starting at random position in the sequence
  Args    : - trimming length
            - sequence string
            - quality score arrayref
  Returns : - truncated sequence string
            - truncated quality score arrayref
=cut

sub _trim_str { 
  my ($self, $trim, $seq_str, $qual_arr) = @_;
  my $length = length($seq_str);
  if ($trim && $length > $trim) {
    my $start_offset = int( rand( $length - $trim + 1 ) );
    my $stop_offset  = $start_offset + $trim - 1;
    $seq_str   = substr $seq_str, $start_offset, $trim;
    @$qual_arr = @$qual_arr[$start_offset .. $stop_offset] if $qual_arr;
  }
  return $seq_str, $qual_arr;
}


=head2 _create_unique_names
  Title   : _create_unique_names
  Usage   : @names = @{$self->_create_unique_names(\@fastafiles)};
  Function: produces unique metagenome names
  Args    : fasta files arrayref
  Returns : unique names arrayref
=cut

sub _create_unique_names {
  my ($self, $fastafile) = @_;
  # Get names from file basename and count how many times each basename occurs
  my @names;
  my %counts;
  for my $fasta (@$fastafile) {
    my @parse = fileparse($fasta, '\..*');
    my $name = $parse[0]; # filename basename
    push @names, $name;
    if (exists $counts{$name}) {
      $counts{$name}++;
    } else {
      $counts{$name} = 1
    }
  }
  # Generate new names for duplicates
  my $separator = '_';
  for (my $i = $#names; $i >= 0; $i--) {
    my $name = $names[$i];
    if ($counts{$name} > 1) {
      my $rep = $counts{$name};
      my $newname = $name.$separator.$rep;
      $names[$i] = $newname;
      $counts{$name}--;
    }
  }
  return \@names;
}


=head2 _filter_low_cplx
  Title   : _filter_low_cplx
  Usage   : my @seqs = _filter_low_cplx( $seq_obj, $cplxthres, $cplxwindow );
  Function: Filter out low complexity sequences by specifying a minimum
            dinucleotide entropy (Shannon-Wiener index h). If a sliding window
            length is specified, regions of low complexity are removed (instead
            of discarding the entire sequence), thereby splitting the sequence
            in several pieces.
  Args    : - Bio::Seq object
            - minimum entropy threshold
            - sliding window length
  Returns : - ranges of high complexity regions in the sequence
=cut

sub _filter_low_cplx {
  my ($self, $seq, $thres, $window) = @_;
  my $seqstr = $seq->seq();
  my $seqlen = length($seqstr);
  my $half   = int($window / 2);
  my $ranges; # regions of the sequence with high complexity
  if ( defined $window && $window == 0 ) {
    $window = undef;
  }
  if ( $window && $window <= $seqlen ) {
    # Excise regions of low complexity from a sequence
    my ($entropy, $dn_freqs, $nof_dns);
    my $range = undef;
    my $kept_prev = 0;
    my $last_good_bp;
    my $last_pos = $seqlen - $window + 1;
    # Entropy of a bp = average entropy of all sliding windows covering this bp
    my @window_entropies;
    my $sum = 0;
    my $nof_windows = 0;
    for my $bp (1 .. $seqlen) {
      my $pos_to_add = $bp;
      # Calculate entropy of the new window, add it to the sum
      if ($pos_to_add <= $last_pos) {
        ($entropy, $dn_freqs, $nof_dns) = $self->_dinucleotide_entropy( $seqstr,
          $pos_to_add, $window, $entropy, $dn_freqs, $nof_dns, $last_pos);
        $window_entropies[$pos_to_add-1] = $entropy;
        $sum += $window_entropies[$pos_to_add-1];
        $nof_windows++;
      }
      # Remove old non-covering window from sum
      my $pos_to_rem = $bp - $window;
      if ($pos_to_rem > 0) {
        $sum -= $window_entropies[$pos_to_rem-1];
        $nof_windows--;
      }
      # Average entropy of the bp
      my $bp_entropy = $sum / $nof_windows;
      # Update
      ($ranges, $last_good_bp) = update_ranges($ranges, $bp, $bp_entropy, $thres, $last_good_bp, $half, $seqlen);
    }
    @window_entropies = undef;

    sub update_ranges {
      my ($ranges, $bp, $bp_entropy, $thres, $last_good_bp, $half, $seqlen) = @_;
      my $win_left  = $bp - $half;
      my $win_right = $bp + $half;
      if ($bp == 1) {
        if ($bp_entropy < $thres) {
          $last_good_bp = undef;
        } else {
          $last_good_bp = $bp;
        }
      } elsif ($bp == $seqlen) {
        if ($bp_entropy < $thres) {
        } else {
          if ($last_good_bp <= $seqlen) {
            my $range_left  = $last_good_bp;
            my $range_right = $win_right;
            $range_right = $seqlen if $range_right > $seqlen;
            push @$ranges, [$range_left, $range_right];
          }
        }
      } else {
        if ($bp_entropy < $thres) {
          if (defined $last_good_bp) {
            my $range_left  = $last_good_bp;
            my $range_right = $win_left - 1;
            if ( $range_right > $range_left) {
              $range_right = $seqlen if $range_right > $seqlen;
              push @$ranges, [$range_left, $range_right];
            }
          }
          $last_good_bp = $win_right + 1;
        } else {
          if (not defined $last_good_bp) {
            $last_good_bp = $bp;
          }
        }
      }
      return $ranges, $last_good_bp;
    }

  } else {
    # Filter out entire sequence if it has low complexity
    my $position = 1;
    my $window = $seqlen;
    my ($entropy, $dn_counts, $nof_dinucs) =
      $self->_dinucleotide_entropy($seqstr, $position, $window);
    if ( $entropy >= $thres ) {
      my $range = [1, $seqlen];
      push @$ranges, $range;
    }
  }

  return $ranges;
}


=head2 _dinucleotide_entropy
  Title   : _dinucleotide_entropy
  Usage   : ($H, $dn_counts, $nof_dinucs) = $self->_dinucleotide_entropy($seqstr,
            $start_pos, $window, $prev_H, $prev_dn_counts, $prev_nof_dinucs)
  Function: Calculate the entropy (Shannon index) of the dinucleotides in a
            sequence fragment 
  Args    : 
  Returns : 
=cut

sub _dinucleotide_entropy {
  my ($self, $seqstr, $start_pos, $window, $prev_H, $prev_dn_counts, $prev_nof_dinucs) = @_;
  my $H;
  my $dn_counts;
  my $nof_dinucs;
  my $stop_pos = $start_pos + $window - 1;
  $seqstr = uc($seqstr);
  return if ( $stop_pos > length($seqstr) );
  # Calculate entropy only if the window does not stick out of the sequence
  $H = $prev_H;
  my $update_possible = 1;
  if ( defined $prev_dn_counts) {
    # Update entropy value
    # Remove previous first dinucleotide
    $dn_counts = $prev_dn_counts;
    $nof_dinucs = $prev_nof_dinucs;
    my $prev_start_off = $start_pos - 2;
    my $prev_dn = substr $seqstr, $prev_start_off, 2;
    if ( exists $$dn_counts{$prev_dn} ) {
      # Dinucleotide is valid
      if ( $$prev_dn_counts{$prev_dn} != 0 ) {
        my $prev_p = $$prev_dn_counts{$prev_dn} / $prev_nof_dinucs;
        $H += $prev_p * log($prev_p);
      }
      $$dn_counts{$prev_dn}--;
      $nof_dinucs--;
      if ( $$dn_counts{$prev_dn} != 0 ) {
        my $p = $$dn_counts{$prev_dn} / $prev_nof_dinucs;
        $H -= $p * log($p);
      }
    }
    # Add next last dinucleotide
    my $next_start_off = $stop_pos - 2;
    my $next_dn = substr $seqstr, $next_start_off, 2;
    if ( exists $$dn_counts{$next_dn} ) {
      if ( $$dn_counts{$next_dn} != 0 ) {
        my $prev_p = $$prev_dn_counts{$next_dn} / $prev_nof_dinucs;
        $H += $prev_p * log($prev_p);
      }
      $nof_dinucs++;
      $$dn_counts{$next_dn}++;
      if ( $$dn_counts{$next_dn} != 0 ) {
        my $p = $$dn_counts{$next_dn} / $nof_dinucs;
        $H -= $p * log($p);
      }
    }
    # If the number of dinucleotides has changed, updating the entropy like I
    # just did is wrong. Fall back to calculating entropy of full window
    if ( $nof_dinucs != $prev_nof_dinucs ) {
      $update_possible = 0;
    }
  }
  if ( not (defined $prev_dn_counts) || $update_possible == 0 ) {
    # Counts for full window
    $nof_dinucs = 0;
    $dn_counts = { 'AA' => 0, 'AC' => 0, 'AG' => 0, 'AT' => 0,
                   'CA' => 0, 'CC' => 0, 'CG' => 0, 'CT' => 0,
                   'GA' => 0, 'GC' => 0, 'GG' => 0, 'GT' => 0,
                   'TA' => 0, 'TC' => 0, 'TG' => 0, 'TT' => 0  };
    for my $pos ($start_pos+1 .. $stop_pos) {
      my $start_off = $pos - 2; # dinucleotide start offset
      my $dinuc = substr $seqstr, $start_off, 2;
      next if not exists $$dn_counts{$dinuc};
      $$dn_counts{$dinuc}++;
      $nof_dinucs++;
    }
    # Entropy of full window
    $H = 0;
    for my $dinuc ( keys %$dn_counts ) {
      next if $$dn_counts{$dinuc} == 0;
      my $p = $$dn_counts{$dinuc} / $nof_dinucs; # dinucleotide probability
      $H -= $p * log($p);
    }
  }
  return $H, $dn_counts, $nof_dinucs;
}


=head2 _valid_dinucleotide
  Title   : _valid_dinucleotide
  Usage   : $is_valid = $self->_dinucleotide_entropy($dinuc, $dinuc_hashref)
  Function: Check if a dinucleotide is a valid sequence (containing only 'A's,
            'C's, 'G's and 'T's
  Args    : - dinucleotide to test
            - hashref of all valid dinucleotides
  Returns : 1 (valid) or 0 (not valid)
=cut

sub _valid_dinucleotide {
  my ($self, $dinuc, $possible_dinuc_hashref) = @_;
  if (exists $$possible_dinuc_hashref{$dinuc}) {
    return 1; # valid dinucleotide
  } else {
    return 0; # invalid
  }
}

1;
