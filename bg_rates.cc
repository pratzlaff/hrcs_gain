#include <vector>
#include <string>
#include <cstdlib>
#include <iterator>
#include <cstdio>
#include <getopt.h>
#include <iostream>
#include "lab.hh"

using std::cout;
using std::cerr;
using std::vector;
using std::string;
using std::ostream_iterator;

namespace {

  namespace opts {
    const char* config = lab::testfile;
    const char* evtdir = lab::evtdir;
    const char* bgdir = lab::bgdir;
    int filter = 1;
  
    char* version_string = "0.1";
    int help = 0;
    int version = 0;
    option lopts[] = {
      { "help",    no_argument, &help, 1 },
      { "version", no_argument, &version, 1 },
      { "nofilter", no_argument, &filter, 0 },
      { 0, 0, 0, 0 }
    };
  }

  int help();
  int version();
}

int main(int argc, char** argv) {

  int c;
  while ((c=getopt_long_only(argc, argv, "", opts::lopts, 0))!=-1) {
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

  if (opts::help) return help();
  if (opts::version) return version();

  if ( argc-optind != 0) {
    cerr << "Usage: " << argv[0] << " [options] file\n";
    return EXIT_FAILURE;
  }

  vector<string> line;
  vector<int> energy;
  vector<int> mcp;
  vector<int> time;
  vector<int> bg_time;
  vector<string> hrc_file;
  vector<string> bg_hrc_file;

  lab::test_data(vector<string>(), line, energy, mcp, time, hrc_file, bg_time, bg_hrc_file);

  copy(line.begin(), line.end(), ostream_iterator<string>(cout, " "));
  cout << '\n';

  copy(hrc_file.begin(), hrc_file.end(), ostream_iterator<string>(cout, " "));
  cout << '\n';

  copy(bg_hrc_file.begin(), bg_hrc_file.end(), ostream_iterator<string>(cout, " "));
  cout << '\n';

  copy(mcp.begin(), mcp.end(), ostream_iterator<int>(cout, " "));
  cout << '\n';

  copy(energy.begin(), energy.end(), ostream_iterator<int>(cout, " "));
  cout << '\n';

  copy(time.begin(), time.end(), ostream_iterator<int>(cout, " "));
  cout << '\n';

  copy(bg_time.begin(), bg_time.end(), ostream_iterator<int>(cout, " "));
  cout << '\n';

} // main

namespace {

  int version() {
    cout << opts::version_string << '\n';
    return 0;
  }

  int help() {
    const char* help_text = "\
=head1 NAME\n\
\n\
template - a template C++ program\n\
\n\
=head1 SYNOPSIS\n\
\n\
template [options] arg1 arg2\n\
\n\
=head1 DESCRIPTION\n\
\n\
A template for a C++ program with option processing and --help\n\
text placeholder.\n\
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
