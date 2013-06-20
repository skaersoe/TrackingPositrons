#ifndef NA63_SIMULATION_ENERGYLOSS_H
#define NA63_SIMULATION_ENERGYLOSS_H

#ifndef Float
typedef float Float;
#endif
#include <TH2.h>
#include <TH3.h>

namespace na63 {

class EnergyLoss {

public:
  EnergyLoss(
    const unsigned xbins,
    const Float xmin,
    const Float xmax,
    const unsigned ybins,
    const Float ymin,
    const Float ymax,
    const unsigned zbins,
    const Float zmin,
    const Float zmax
  );
  ~EnergyLoss();
  void Fill(const Float x, const Float y, const Float z, const Float E);
  Float Get(const Float x, const Float y, const Float z) const;
  Float Get(const unsigned xbin, const unsigned ybin, const unsigned zbin) const;
  Float GetSum() const;
  TH1F* BuildXHistogram() const;
  TH2F* BuildXYHistogram() const;
  TH2F* BuildXZHistogram() const;
  TH2F* BuildYZHistogram() const;
  TH3F* Build3DHistogram() const;

private:
  const Float xmin, xmax, ymin, ymax, zmin, zmax;
  const unsigned xbins, ybins, zbins;
  const unsigned size;
  Float* grid, sum;

  static unsigned TranslateToIndex(const Float f, const Float min,
      const Float max, const unsigned bins);
  unsigned GetRawIndex(const unsigned x, const unsigned y,
      const unsigned z) const;
  unsigned GetXBin(const Float f) const;
  unsigned GetYBin(const Float f) const;
  unsigned GetZBin(const Float f) const;

};

} // End namespace na63

#endif /* NA63_SIMULATION_ENERGYLOSS_HH */