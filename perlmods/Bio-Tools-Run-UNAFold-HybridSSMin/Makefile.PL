use 5.008008;
use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    NAME              => 'Bio::Tools::Run::UNAFold::HybridSSMin',
    VERSION_FROM      => 'lib/Bio/Tools/Run/UNAFold/HybridSSMin.pm', # finds $VERSION
    PREREQ_PM         => {}, # e.g., Module::Name => 1.1
    ($] >= 5.005 ?     ## Add these new keywords supported since 5.005
      (ABSTRACT       =>  'Bioperl wrapper for the UNAFold hybrid-ss-min tool',
       AUTHOR         => 'Caroline Johnston <cass@>') : ()),
);
