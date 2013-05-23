#include "TrackSeeder.hh"

namespace na63 {
  
  std::vector<seed> Seed(std::vector<hit> hits) {

    std::vector<seed> tracks;

    for (int i = 0; i < hits.size(); i++) {
      if (hits[i].det != 1)
        break;
      seed newSeed;
      newSeed.hitIndex = i;
      newSeed.x = 0.0;
      newSeed.y = 0.0;
      newSeed.z = 0.0;
      newSeed.dx = (hits[i].x - newSeed.x) / (hits[i].z - newSeed.z);
      newSeed.dy = (hits[i].y - newSeed.y) / (hits[i].z - newSeed.z);
      seeds.push_back(newSeed);
    }

    return seeds;

  }

}