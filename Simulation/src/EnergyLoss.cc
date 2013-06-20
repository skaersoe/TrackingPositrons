#include "Simulation/EnergyLoss.hh"
#include <cassert>

namespace na63 {

unsigned EnergyLoss::TranslateToIndex(const Float f, const Float min,
    const Float max, const unsigned bins) {

  Float percentage = (f - min) / (max - min);
  unsigned index = unsigned(percentage * (Float)(bins - 1) + 0.5f);

  return index;
}

unsigned EnergyLoss::GetRawIndex(const unsigned x, const unsigned y,
    const unsigned z) const {
  unsigned index = x + y*xbins + z*xbins*ybins;
  // assert(index < size);
  return index;
}

unsigned EnergyLoss::GetXBin(const Float f) const {
  return TranslateToIndex(f,xmin,xmax,xbins);
}

unsigned EnergyLoss::GetYBin(const Float f) const {
  return TranslateToIndex(f,ymin,ymax,ybins);
}

unsigned EnergyLoss::GetZBin(const Float f) const {
  return TranslateToIndex(f,zmin,zmax,zbins);
}
  
EnergyLoss::EnergyLoss(
    const unsigned xbins_,
    const Float xmin_,
    const Float xmax_,
    const unsigned ybins_,
    const Float ymin_,
    const Float ymax_,
    const unsigned zbins_,
    const Float zmin_,
    const Float zmax_)
    : xbins(xbins_), xmin(xmin_), xmax(xmax_),
      ybins(ybins_), ymin(ymin_), ymax(ymax_),
      zbins(zbins_), zmin(zmin_), zmax(zmax_),
      size(xbins_*ybins_*zbins_), sum(0) {

  grid = new Float[size];
  for (int i=0;i<size;i++) grid[i] = 0;

}

EnergyLoss::~EnergyLoss() {
  delete grid;
}

void EnergyLoss::Fill(const Float x, const Float y, const Float z,
    const Float E) {
  unsigned xbin = GetXBin(x);
  unsigned ybin = GetYBin(y);
  unsigned zbin = GetZBin(z);
  assert(E >= 0);
  assert(E == E);
  unsigned index = GetRawIndex(xbin,ybin,zbin);
  if (index >= size) return;
  sum += E;
  grid[index] += E;
}

Float EnergyLoss::Get(const unsigned xbin, const unsigned ybin,
    const unsigned zbin) const {
  assert(xbin < xbins && ybin < ybins && zbin < zbins);
  return grid[GetRawIndex(xbin,ybin,zbin)];
}

Float EnergyLoss::Get(const Float x, const Float y, const Float z) const {
  unsigned xbin = GetXBin(x);
  unsigned ybin = GetYBin(y);
  unsigned zbin = GetZBin(z);
  return Get(xbin,ybin,zbin);
}

Float EnergyLoss::GetSum() const {
  return sum;
}

TH1F* EnergyLoss::BuildXHistogram() const {
  TH1F *hist = new TH1F(
    "Energy loss",
    "Energy loss;x",
    xbins,xmin,xmax
  );
  for (unsigned i=0;i<xbins;i++) {
    Float sum = 0;
    for (unsigned j=0;j<ybins;j++) {
      for (unsigned k=0;k<zbins;k++) {
        sum += Get(i,j,k);
      }
    }
    hist->SetBinContent(i,sum);
  }
  return hist;
}

TH2F* EnergyLoss::BuildXYHistogram() const {
  TH2F *hist = new TH2F(
    "Energy loss",
    "Energy loss;x;y",
    xbins,xmin,xmax,
    ybins,ymin,ymax
  );
  for (unsigned i=0;i<xbins;i++) {
    for (unsigned j=0;j<ybins;j++) {
      Float sum = 0;
      for (unsigned k=0;k<zbins;k++) {
        sum += Get(i,j,k);
      }
      hist->SetBinContent(i,j,sum);
    }
  }
  return hist;
}

TH2F* EnergyLoss::BuildXZHistogram() const {
  TH2F *hist = new TH2F(
    "Energy loss",
    "Energy loss;x;z",
    xbins,xmin,xmax,
    zbins,zmin,zmax
  );
  for (unsigned i=0;i<xbins;i++) {
    for (unsigned j=0;j<zbins;j++) {
      Float sum = 0;
      for (unsigned k=0;k<ybins;k++) {
        sum += Get(i,k,j);
      }
      hist->SetBinContent(i,j,sum);
    }
  }
  return hist;
}

TH2F* EnergyLoss::BuildYZHistogram() const {
  TH2F *hist = new TH2F(
    "Energy loss",
    "Energy loss;y;z",
    ybins,ymin,ymax,
    zbins,zmin,zmax
  );
  for (unsigned i=0;i<ybins;i++) {
    for (unsigned j=0;j<zbins;j++) {
      Float sum = 0;
      for (unsigned k=0;k<xbins;k++) {
        sum += Get(k,i,j);
      }
      hist->SetBinContent(i,j,sum);
    }
  }
  return hist;
}

TH3F* EnergyLoss::Build3DHistogram() const {
  TH3F *hist = new TH3F(
    "Energy loss",
    "Energy loss;x;y;z",
    xbins,xmin,xmax,
    ybins,ymin,ymax,
    zbins,zmin,zmax
  );
  for (unsigned i=0;i<xbins;i++) {
    for (unsigned j=0;j<ybins;j++) {
      for (unsigned k=0;k<zbins;k++) {
        hist->SetBinContent(i,j,k,Get(i,j,k));
      }
    }
  }
  return hist;
}

} // End namespace na63