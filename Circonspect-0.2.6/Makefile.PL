use inc::Module::Install;

name     'Circonspect';
all_from 'lib/Circonspect.pm';

build_requires 'Test::More';

# Perl module dependencies for Circonspect
requires 'Math::Random::MT'                     => '1.16';
requires 'Getopt::Declare'                      => 0;

requires 'Graph'                                => 0;
# optional dependency of BioPerl, but necessary for Bio::Assembly::ContigSpectrum

# Require a specific bioperl version: 1.6.2 (1.006002) or development (1.0069)
requires 'Bio::Root::Version'                   => '1.006002';
# Require specific bioperl-live modules
requires 'Bio::PrimarySeq'                      => 0;
requires 'Bio::SeqIO'                           => 0;
requires 'Bio::Seq::Quality'                    => 0;
requires 'Bio::Root::Utilities'                 => 0;
requires 'Bio::DB::Fasta'                       => 0;
requires 'Bio::DB::Qual'                        => 0;
requires 'Bio::Assembly::Singlet'               => 0;
requires 'Bio::Assembly::Tools::ContigSpectrum' => 0;
# Require specific bioperl-run modules
requires 'Bio::Tools::Run::TigrAssembler'       => 0;
requires 'Bio::Tools::Run::Cap3'                => 0;
requires 'Bio::Tools::Run::Minimo'              => 0;
requires 'Bio::Tools::Run::Newbler'             => 0;

# External Circonspect dependency: Minimo
requires_external_bin 'Minimo';
# Other assemblers are supported, by order of preference: LIGR Assembler, CAP3,
# TIGR Assembler and Newbler. I made Minimo a requirement since it is Circonspect's
# default and because there is no way to require here one assembler amongst the
# 5 supported.

install_script 'script/Circonspect';
install_script 'utils/csp2csv';
install_script 'utils/batch_cross_contigs.sh';
install_script 'utils/reformat_quality_file';
install_script 'utils/average_contig_spectra';

auto_install;

WriteAll;

