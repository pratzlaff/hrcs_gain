#include <iterator>
#include <stdexcept>
#include <map>
#include <algorithm>
#include <numeric>
#include <vector>
#include <cpputil/byte.hh>
#include <cpputil/sequence.hh>
#include <cpputil/rdb.hh>
#include <cpputil/ss_cast.hh>
#include "lab.hh"

namespace lab {

  using util::isbigendian;
  using util::bswap;
  using util::sequence;
  using util::rdb_read;

  using std::vector;
  using std::string;
  using std::map;

  bool binfile_input::read_header(int &ytap, int &ysubtap,
				  int &xtap, int &xsubtap,
				  int &y1, int &y2,
				  int &x1, int &x2)
  {
    util::uint32 yytap, yysubtap, xxtap, xxsubtap, yy1, yy2, xx1, xx2;

    in.read(reinterpret_cast<char*>(&yytap), 4);
    if (in.eof())
      return false;
    in.read(reinterpret_cast<char*>(&yysubtap), 4);
    in.read(reinterpret_cast<char*>(&xxtap), 4);
    in.read(reinterpret_cast<char*>(&xxsubtap), 4);
    in.read(reinterpret_cast<char*>(&yy1), 4);
    in.read(reinterpret_cast<char*>(&yy2), 4);
    in.read(reinterpret_cast<char*>(&xx1), 4);
    in.read(reinterpret_cast<char*>(&xx2), 4);

    if (!in)
      throw binfile_error("file truncated");

    if (!isbigendian()) {
      bswap(yytap);
      bswap(yysubtap);
      bswap(xxtap);
      bswap(xxsubtap);
      bswap(yy1);
      bswap(yy2);
      bswap(xx1);
      bswap(xx2);
    }

    ytap = yytap;
    ysubtap = yysubtap;
    xtap = xxtap;
    xsubtap = xxsubtap;
    y1 = yy1;
    y2 = yy2;
    x1 = xx1;
    x2 = xx2;

    return true;
  }

  void binfile_input::read_data(vector<int> &x, vector<int> &y)
  {
    util::uint32 yy[256];

    in.read(reinterpret_cast<char*>(&yy[0]), 256 * 4);

    if (!in)
      throw binfile_error("file truncated");

    if (!isbigendian())
      bswap(yy, yy+sizeof(yy)/sizeof(yy[0]));

    x.resize(256);
    generate(x.begin(), x.end(), sequence<int>(0,1));
    y = vector<int>(yy, yy+sizeof(yy)/sizeof(yy[0]));
  }

  bool binfile_input::next_subtap(int &ytap, int &ysubtap,
				  int &xtap, int &xsubtap,
				  int &y1, int &y2, int &x1, int &x2)
  {
    if (!read_header(ytap, ysubtap, xtap, xsubtap, y1, y2, x1, x2))
      return false;

    in.seekg(256 * 4, std::ios_base::cur);

    if (!in)
      throw binfile_error("file truncated");

    return true;
  }


  bool binfile_input::next_subtap(int &ytap, int &ysubtap,
				  int &xtap, int &xsubtap,
				  int &y1, int &y2, int &x1, int &x2,
				  vector<int> &x, vector<int> &y)
  {

    if (!read_header(ytap, ysubtap, xtap, xsubtap, y1, y2, x1, x2))
      return false;

    read_data(x, y);

    if (!in)
      throw binfile_error("file truncated");

    return true;
  }

  namespace {
    int str2int(const string& s) {
      return util::ss_cast<int, string>(s);
    }
  }

  void test_data(const vector<string>& anodes,
		 vector<string>& line,
		 vector<int>&         energy,
		 vector<int>&         mcp,
		 vector<int>&         time,
		 vector<string>& hrc_file,
		 vector<int>&         bg_time,
		 vector<string>& bg_hrc_file
		 )
  {
    map<string, vector<string> > cols;

    const char* names[] = { "line",
			    "energy",
			    "MCP",
			    "time",
			    "HRC_file",
			    "b_time",
			    "b_HRC_file",
    };

    rdb_read(lab::testfile, cols, vector<string>(names, names+sizeof(names)/sizeof(names[0])));

    line        = cols["line"];
    hrc_file    = cols["HRC_file"];
    bg_hrc_file = cols["b_HRC_file"];

    //    int util::ss_cast(const int&, string);
    //    int (*func)(string&) = &util::ss_cast<int, string>;

    // convert strings to output types
    //    transform(cols["energy"].begin(), cols["energy"].end(), back_inserter(energy), &util::ss_cast<int, string>);
    transform(cols["energy"].begin(), cols["energy"].end(),
	      back_inserter(energy), str2int);
    transform(cols["MCP"].begin(), cols["MCP"].end(),
	      back_inserter(mcp), str2int);
    transform(cols["time"].begin(), cols["time"].end(),
	      back_inserter(time), str2int);
    transform(cols["b_time"].begin(), cols["b_time"].end(),
	      back_inserter(bg_time), str2int);

  }

} // namespace lab
