/usr/bin/perl5.8.8 /usr/lib/perl5/5.8.8/ExtUtils/xsubpp  -typemap /usr/lib/perl5/5.8.8/ExtUtils/typemap -typemap /usr/lib/perl5/vendor_perl/5.8.8/i386-linux/PDL/Core/typemap.pdl   foo_pl_d362.xs > foo_pl_d362.xsc && mv foo_pl_d362.xsc foo_pl_d362.c
gcc -c  -I/usr/lib/perl5/vendor_perl/5.8.8/i386-linux/PDL/Core -fno-strict-aliasing -pipe -Wdeclaration-after-statement -I/usr/local/include -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -I/usr/include/gdbm -O2 -pipe -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fomit-frame-pointer -march=i586 -mtune=pentiumpro -fasynchronous-unwind-tables   -DVERSION=\"0.00\" -DXS_VERSION=\"0.00\" -fPIC "-I/usr/lib/perl5/5.8.8/i386-linux/CORE"   foo_pl_d362.c
In file included from /usr/lib/perl5/vendor_perl/5.8.8/i386-linux/PDL/Core/pdlcore.h:10,
                 from foo_pl_d362.xs:9:
/usr/lib/perl5/vendor_perl/5.8.8/i386-linux/PDL/Core/ppport.h:227:1: warning: "PERL_UNUSED_DECL" redefined
In file included from foo_pl_d362.xs:6:
/usr/lib/perl5/5.8.8/i386-linux/CORE/perl.h:163:1: warning: this is the location of the previous definition
foo_pl_d362.xs: In function ‘pdl_shift_negs_readdata’:
foo_pl_d362.xs:173: warning: comparison is always false due to limited range of data type
foo_pl_d362.xs:186: error: expected ‘while’ before ‘break’
foo_pl_d362.xs:219: error: expected ‘while’ before ‘break’
foo_pl_d362.xs:239: warning: comparison is always false due to limited range of data type
foo_pl_d362.xs:252: error: expected ‘while’ before ‘break’
foo_pl_d362.xs:285: error: expected ‘while’ before ‘break’
foo_pl_d362.xs:318: error: expected ‘while’ before ‘break’
foo_pl_d362.xs:351: error: expected ‘while’ before ‘break’
foo_pl_d362.xs:385: error: expected ‘while’ before ‘default’
foo_pl_d362.xs:393: warning: ISO C90 forbids mixed declarations and code
foo_pl_d362.c:1404: error: static declaration of ‘XS_foo_pl_d362_set_debugging’ follows non-static declaration
foo_pl_d362.c:1402: error: previous declaration of ‘XS_foo_pl_d362_set_debugging’ was here
foo_pl_d362.c:1424: error: static declaration of ‘XS_foo_pl_d362_set_boundscheck’ follows non-static declaration
foo_pl_d362.c:1422: error: previous declaration of ‘XS_foo_pl_d362_set_boundscheck’ was here
foo_pl_d362.c:1446: error: static declaration of ‘XS_PDL_shift_negs’ follows non-static declaration
foo_pl_d362.c:1444: error: previous declaration of ‘XS_PDL_shift_negs’ was here
foo_pl_d362.c:1521: error: static declaration of ‘XS_PDL_inc’ follows non-static declaration
foo_pl_d362.c:1519: error: previous declaration of ‘XS_PDL_inc’ was here
foo_pl_d362.c:1630: error: static declaration of ‘XS_PDL_tcumul’ follows non-static declaration
foo_pl_d362.c:1628: error: previous declaration of ‘XS_PDL_tcumul’ was here
foo_pl_d362.c:1735: error: static declaration of ‘boot_foo_pl_d362’ follows non-static declaration
foo_pl_d362.c:1733: error: previous declaration of ‘boot_foo_pl_d362’ was here
foo_pl_d362.c:1769: error: expected declaration or statement at end of input
make: *** [foo_pl_d362.o] Error 1
