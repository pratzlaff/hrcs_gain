#ifndef LAB_HH
#define LAB_HH

#include <vector>
#include <string>
#include <stdexcept>
#include <fstream>

namespace lab {

  const char* const bgdir = "/data/legs/rpete/data/hrcs_lab/bg";
  const char* const mergedbg = "/data/legs/rpete/data/hrcs_lab/bg/merged_bg.bin";
  const char* const evtdir = "/data/legs/rpete/data/hrcs_lab/evt1";
  const char* const analdir = "/data/legs/rpete/data/hrcs_lab/analysis";
  const char* const testfile = "/data/legs/rpete/cal/hrcs_gain/hrcs_lab.rdb";

  const std::size_t chipx_min = 1;
  const std::size_t chipx_max = 4096;
  const std::size_t chipy_min = 1;
  const std::size_t chipy_max = 16384;

  const std::size_t rawx_min = chipx_min;
  const std::size_t rawx_max = chipx_max;
  const std::size_t rawy_min = chipy_min;
  const std::size_t rawy_max = 3 * chipy_max;

  const std::size_t tapsize = 256;
  const std::size_t subtaps = 3;

  void test_data(const std::vector<std::string>& anodes,
		 std::vector<std::string>& line,
		 std::vector<int>&         energy,
		 std::vector<int>&         mcp,
		 std::vector<int>&         time,
		 std::vector<std::string>& hrc_file,
		 std::vector<int>&         bg_time,
		 std::vector<std::string>& bg_hrc_file);

  class binfile_error : public std::runtime_error
  {
  public:
    binfile_error(const std::string& s = "unidentified error")
      : std::runtime_error(s)
    { }
  };

  class binfile_input {
  private:

    bool read_header(int&, int&, int&, int&, int&, int&, int&, int&);
    void read_data(std::vector<int>&, std::vector<int>&);
    std::fstream in;

  public:

    binfile_input ( const std::string& s )
      : in(s.c_str(), std::ios_base::binary | std::ios_base::in)
    {
      if (!in)
	throw binfile_error("unable to open file "+s);
    
    }

    bool next_subtap( int &ytap, int &ysubtap,
		      int &xtap, int &xsubtap,
		      int &y1, int &y2,
		      int &x1, int &x2);

    bool next_subtap( int &ytap, int &ysubtap,
		      int &xtap, int &xsubtap,
		      int &y1, int &y2,
		      int &x1, int &x2,
		      std::vector<int> &x,
		      std::vector<int> &y );

  };


} // namespace lab

#endif
