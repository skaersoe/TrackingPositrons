#ifndef NA63_FILEHANDLER  
#define NA63_FILEHANDLER
#include <string>
#include <vector>

namespace na63{

  std::vector<hit> LoadHits(std::string filename); // Loads hits from a file. Returns 1 if everything goes okay, 0 otherwise.

}

#endif