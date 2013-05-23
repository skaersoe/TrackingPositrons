#include "FileHandler.hh"
#include <stdlib.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <cstdio>

namespace na63{

  std::vector<hit> LoadHits(const std::string filename) {

    std::vector<hit> hitarray;

    std::ifstream hitfile;

    hitfile.open(filename.c_str(), std::ifstream::in);
    if (hitfile.fail()) {
      std::cout << "The file could not be opened" << std::endl;
    }

    std::string filebuffer;
    int newpos;

    while (hitfile.good()) {
      getline(hitfile, filebuffer);
      if (filebuffer[0] == '[') {
        hit newhit;
        sscanf(filebuffer.c_str(), "[%lf, %lf, %lf, %i][%i, %i, %i", &newhit.x, &newhit.y, &newhit.z, &newhit.det,
                                    &newhit.x_pixel, &newhit.y_pixel, &newhit.z_pixel);
        hitarray.push_back(newhit);
      }
    }

//    for (int i = 0; i < hitarray.size(); i++) {
//
//      std::cout << hitarray[i].det << std::endl;
//
//    }

    hitfile.close();

    return hitarray;
  }

}