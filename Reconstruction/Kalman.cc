#include "FileHandler.hh"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>

using namespace na63;

  int main() {

    std::cout << "Path to hits file: ";

    std::string filename;

    std::getline(std::cin, filename);

    std::cout << "\n" << "Loading hits..." << "\n";

    std::vector<hit> hits = LoadHits(filename);

    for (int i = 0; i < hits.size(); i++) {
      std::cout << hits[i].x_pixel << std::endl;
    }

    return 1;
  }