#pragma once
#include <array>
#include <openbabel/mol.h>
#include <vector>
using XYZ = std::array<float, 3>;
using USRTRIPPLE = std::array<float, 3>;
using USRCATDESCRIPTORS = std::array<float, 60>;

class Usrcat {
public:
  Usrcat();
  ~Usrcat();
  float score(const USRCATDESCRIPTORS &a, const USRCATDESCRIPTORS &b) const;
  USRCATDESCRIPTORS GetUsrcatValues(OpenBabel::OBMol &m) const; // Mol cannot be
                                                                // const because
                                                                // OBMolAtomIter
                                                                // returns
                                                                // non-const
                                                                // iterators

private:
	const std::array<std::string,5> smartspatterns{ {
		"[!#1]",                                                   // All atoms
		"[#6+0!$(*~[#7,#8,F]),SH0+0v2,s+0,S^3,Cl+0,Br+0,I+0]", // Hydrophobic
		"[a]",                                                 // Aromatic
		"[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N&v3;H1,"
		"H2]-[!$(*=[O,N,P,S])]),$([N;v3;H0]),$([n,o,s;+0]),F]", // Acceptor
		"[N!H0v3,N!H0+v4,OH+0,SH+0,nH+0]"} };                     // Donor

  inline USRTRIPPLE GetTripple(const std::vector<XYZ> &coords,
                               const XYZ &point) const;
  inline float Euclid3DDistSquared(const XYZ &point1, const XYZ &point2) const;
  inline XYZ GetCentroid(const std::vector<XYZ> &coords) const;
  inline std::vector<XYZ>::const_iterator
  WhichClosestToPoint(const std::vector<XYZ> &coords, const XYZ &point) const;
  inline std::vector<XYZ>::const_iterator
  WhichFarthestFrom(const std::vector<XYZ> &coords, const XYZ &point) const;
};

Usrcat::Usrcat() {}

Usrcat::~Usrcat() {}

USRCATDESCRIPTORS Usrcat::GetUsrcatValues(OpenBabel::OBMol &m) const {
  USRCATDESCRIPTORS descriptors;
  std::vector<XYZ> coords;
  coords.reserve(m.NumAtoms());
  for (OpenBabel::OBMolAtomIter a(m); a; ++a) {
    coords.push_back({{static_cast<float>(a->x()), static_cast<float>(a->y()),
                       static_cast<float>(a->z())}});
  }

  XYZ centroid = GetCentroid(coords);
  auto closesttocentroid = WhichClosestToPoint(coords, centroid);
  auto farthestfromcentroid = WhichFarthestFrom(coords, centroid);
  auto ffromfarthestfromcentroid =
      WhichFarthestFrom(coords, *farthestfromcentroid);

  unsigned int valuescounter = 0;
  std::vector<std::vector<int>> p4maps(smartspatterns.size(),
                                       std::vector<int>());
  for (unsigned int p = 0; p < smartspatterns.size(); p++) {
    OpenBabel::OBSmartsPattern smarts;
    smarts.Init(smartspatterns[p]);
    std::vector<std::vector<int>> maplist;

    if (smarts.Match(m, maplist)) {
      if (maplist.size() > 0) {
        for (unsigned int i = 0; i < maplist.size(); i++) {
          p4maps[p].push_back(maplist[i][0]);
        }
        std::vector<XYZ> subset(p4maps[p].size());
        for (unsigned int i = 0; i < p4maps[p].size(); i++) {
          subset[i] = coords[(p4maps[p][i] - 1)];
        }
        USRTRIPPLE temptripple;
        temptripple = GetTripple(subset, centroid);
        descriptors[valuescounter++] = temptripple[0];
        descriptors[valuescounter++] = temptripple[1];
        descriptors[valuescounter++] = temptripple[2];
        temptripple = GetTripple(subset, *closesttocentroid);
        descriptors[valuescounter++] = temptripple[0];
        descriptors[valuescounter++] = temptripple[1];
        descriptors[valuescounter++] = temptripple[2];
        temptripple = GetTripple(subset, *farthestfromcentroid);
        descriptors[valuescounter++] = temptripple[0];
        descriptors[valuescounter++] = temptripple[1];
        descriptors[valuescounter++] = temptripple[2];
        temptripple = GetTripple(subset, *ffromfarthestfromcentroid);
        descriptors[valuescounter++] = temptripple[0];
        descriptors[valuescounter++] = temptripple[1];
        descriptors[valuescounter++] = temptripple[2];
      }
    }
  }
  return descriptors;
}

inline USRTRIPPLE Usrcat::GetTripple(const std::vector<XYZ> &coords,
                                     const XYZ &point) const {
  USRTRIPPLE mvs = {{0.0, 0.0, 0.0}};

  if (coords.size() == 0)
    return mvs;
  if (coords.size() == 1) {
    mvs[0] = sqrtf(Euclid3DDistSquared(coords[0], point));
    return mvs;
  }

  float fnumpoints = static_cast<float>(coords.size());
  std::vector<float> dfrompoint(coords.size(), 0.0f);

  // Populate dfrompoint + determine mean;
  for (unsigned i = 0; i < coords.size(); ++i) {
    dfrompoint[i] = sqrtf(Euclid3DDistSquared(coords[i], point));
    mvs[0] += dfrompoint[i];
  }
  mvs[0] /= fnumpoints;

  // Variance AND Skew
  float tmp;
  float topsum = 0;
  for (auto val : dfrompoint) {
    tmp = val - mvs[0];
    mvs[1] += tmp * tmp;
    topsum += tmp * tmp * tmp;
  }

  // If mvs[1] is zero at this stage, skew (mvs[2]) will be nan.  Check mvs[1]
  // after its correction.
  mvs[2] = (topsum / (powf(mvs[1], 1.5))) * sqrtf(fnumpoints);
  mvs[1] *= (1.0 / (fnumpoints - 1.0f));
  if (mvs[1] == 0)
    mvs[2] = 0;
  return mvs;
}

inline float Usrcat::Euclid3DDistSquared(const XYZ &point1,
                                         const XYZ &point2) const {
  return (((point1[0] - point2[0]) * (point1[0] - point2[0]))) +
         (((point1[1] - point2[1]) * (point1[1] - point2[1]))) +
         (((point1[2] - point2[2]) * (point1[2] - point2[2])));
}

inline XYZ Usrcat::GetCentroid(const std::vector<XYZ> &coords) const {
  // Very modern C++ way to get molecular centroid using standard library and
  // lambdas
  std::vector<XYZ>::const_iterator minx, miny, minz, maxx, maxy, maxz;
  std::tie(minx, maxx) = std::minmax_element(
      coords.cbegin(), coords.cend(),
      [](const XYZ &a, const XYZ &b) { return a[0] < b[0]; });
  std::tie(miny, maxy) = std::minmax_element(
      coords.cbegin(), coords.cend(),
      [](const XYZ &a, const XYZ &b) { return a[1] < b[1]; });
  std::tie(minz, maxz) = std::minmax_element(
      coords.cbegin(), coords.cend(),
      [](const XYZ &a, const XYZ &b) { return a[2] < b[2]; });
  return {{(((*maxx)[0] - (*minx)[0]) / 2.0f) + (*minx)[0],
           (((*maxy)[1] - (*miny)[1]) / 2.0f) + (*miny)[1],
           (((*maxz)[2] - (*minz)[2]) / 2.0f) + (*minz)[2]}};
}

inline std::vector<XYZ>::const_iterator
Usrcat::WhichClosestToPoint(const std::vector<XYZ> &coords,
                            const XYZ &point) const {
  float distsquared, mindistsquared = std::numeric_limits<float>::max();
  auto mindistindex = coords.cbegin();
  for (auto i = coords.cbegin(); i != coords.cend(); ++i) {
    distsquared = Euclid3DDistSquared(*i, point);
    if (distsquared < mindistsquared) {
      mindistsquared = distsquared;
      mindistindex = i;
    }
  }
  return mindistindex;
}
inline std::vector<XYZ>::const_iterator
Usrcat::WhichFarthestFrom(const std::vector<XYZ> &coords,
                          const XYZ &point) const {
  auto maxdistindex = coords.cbegin();
  float distsquared, maxdistsquared = 0;
  for (auto i = coords.cbegin(); i != coords.cend(); ++i) {
    distsquared = Euclid3DDistSquared(*i, point);
    if (distsquared > maxdistsquared) {
      maxdistsquared = distsquared;
      maxdistindex = i;
    }
  }
  return maxdistindex;
}
inline float Usrcat::score(const USRCATDESCRIPTORS &a,
                           const USRCATDESCRIPTORS &b) const {
  float tot{0.0f};
  for (unsigned i = 0; i < a.size(); ++i) {
    tot += ((a[i] - b[i]) * (a[i] - b[i]));
  }
  return tot;
}