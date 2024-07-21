use ExtUtils::MakeMaker;
my %options = %{       
{
  'TYPEMAPS' => [
    '/usr/lib/perl5/5.8.8/ExtUtils/typemap',
    '/usr/lib/perl5/vendor_perl/5.8.8/i386-linux/PDL/Core/typemap.pdl'
  ],
  'NAME' => 'foo_pl_596c',
  'INC' => '-I/usr/lib/perl5/vendor_perl/5.8.8/i386-linux/PDL/Core',
  'VERSION' => '0.00'
}
};
WriteMakefile(%options);

# Remove the Makefile dependency. Causes problems on a few systems.
sub MY::makefile { '' }
