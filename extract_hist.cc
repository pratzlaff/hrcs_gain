#include <cstdlib>
#include <cstdio>
#include <getopt.h>
#include <iostream>
#include <vector>
#include <iterator>
#include <cpputil/ss_cast.hh>
#include "lab.hh"

using std::vector;
using std::cout;
using std::cerr;

namespace {
  char* version_string = "0.1";
  int hflag = 0;
  int vflag = 0;
  int iflag = 0;
  option lopts[] = {
    { "help",    no_argument, &hflag, 1 },
    { "version", no_argument, &vflag, 1 },
    { "info",    no_argument, &iflag, 1 },
    { 0, 0, 0, 0 }
  };
  int help();
  int version();
}

int main(int argc, char** argv) {

  int c;
  while ((c=getopt_long_only(argc, argv, "", lopts, 0))!=-1) {
    switch (c) {
    // a flag was set/unset on our behalf, nothing more to do
    case 0:
      break;
    // problem occurred
    case '?':
    case ':':
      cerr << "Try `--help' for more information.\n";
      return EXIT_FAILURE;
    // didn't handle all of our specified options
    default:
      cerr << "programmer error, unhandled option = "; cerr.put(c); cerr << '\n';
      return EXIT_FAILURE;
    }
  }

  if (hflag) return help();
  if (vflag) return version();

  if (
      ( !iflag && (argc-optind != 5) ) ||
      (  iflag && (argc-optind != 1) )     ) {
    cerr << "Usage: " << argv[0] << " [options] (--info | ytap ysubtap xtap xsubtap) binfile\n";
    return EXIT_FAILURE;
  }

  int ytap, ysubtap, xtap, xsubtap, x1, x2, y1, y2;

  // just print information about the subtaps in the file
  if (iflag) {

    lab::binfile_input f(argv[optind++]);

    cout << "ytap\tysubtap\txtap\txsubtap\ty1\ty2\tx1\tx2\n";
    cout << "N\tN\tN\tN\tN\tN\tN\tN\n";
    while (f.next_subtap(ytap, ysubtap, xtap, xsubtap, y1, y2, x1, x2)) {
      cout <<
	ytap << '\t' << ysubtap << '\t' <<
	xtap << '\t' << xsubtap << '\t' <<
	y1 << '\t' << y2 << '\t' << x1 << '\t' << x2 << '\n';
    }
    return 0;
  }

  // print histogram of a single subtap
  else {

    int yytap    = util::ss_cast<int>(argv[optind++]);
    int yysubtap = util::ss_cast<int>(argv[optind++]);
    int xxtap    = util::ss_cast<int>(argv[optind++]);
    int xxsubtap = util::ss_cast<int>(argv[optind++]);

    lab::binfile_input f(argv[optind++]);

    vector<int> x, y;
    while (f.next_subtap(ytap, ysubtap, xtap, xsubtap, y1, y2, x1, x2, x, y)) {

      if (yytap == ytap &&
	  yysubtap == ysubtap &&
	  xxtap == xtap &&
	  xxsubtap == xsubtap) {
	cout << "pha\tn\n";
	cout << "N\tN\n";
	for (vector<int>::size_type i=0; i!=x.size(); ++i)
	  cout << x[i] << "\t" << y[i] << "\n";

	return 0;
      }
    }
  }

} // main

namespace {

  int version() {
    cout << version_string << '\n';
    return 0;
  }

  int help() {
    const char* help_text = "\
=head1 NAME\n\
\n\
extract_hist - extract PHA histogram from a binfile\n\
\n\
=head1 SYNOPSIS\n\
\n\
extract_hist [options] (--info | ytap ysubtap xtap xsubtap) binfile\n\
\n\
=head1 DESCRIPTION\n\
\n\
Extracts PHA histogram values for a given subtap from a binfile (c.f.\n\
F<genstats>) containing the PHA information.\n\
\n\
=head1 OPTIONS\n\
\n\
=over 4\n\
\n\
=item --help\n\
\n\
Print this help text and exit.\n\
\n\
=item --version\n\
\n\
Print the program version and exit.\n\
\n\
=back\n\
\n\
=head1 AUTHOR\n\
\n\
Peter Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> August 2007\n\
\n\
=head1 SEE ALSO\n\
\n\
nothing to see here\n\
\n\
=cut\n\
";

    const char* pager = std::getenv("PAGER");
    if (!pager) pager = "more";

    FILE* pd = popen((std::string("pod2text -c | ")+pager).c_str(), "w");
    if (!pd) {
      std::perror("error starting pod2text");
      return EXIT_FAILURE;
    }

    int n = 0;
    int len = std::strlen(help_text);
    while (n < len) {
      int written = std::fwrite(help_text, 1, len-n, pd);
      if (!written) {
	std::perror("error writing help");
	return EXIT_FAILURE;
      }
      n+=written;
    }

    if (pclose(pd) == -1) {
      std::perror("error writing help");
      return EXIT_FAILURE;
    }

    return 0;
  }

}
