#include <fstream>
#include <sstream>

#include "Point.hh"
#include "read_Mesh2D.hh"

namespace Amanzi {
namespace AmanziGeometry {

Mesh2D
readFile(const std::string& filename,
                 std::vector<int>& soil_types,
                 std::vector<int>& bedrock_types,
                 std::vector<double>& depths_to_bedrock) {
  std::ifstream fin(filename);
  bedrock_types.clear();
  soil_types.clear();
  depths_to_bedrock.clear();
  
  std::string name;
  Point p1(3),p2(3),p3(3), c(2);
  std::vector<double> depths(3);
  int tri_index;
  std::vector<int> veg_types;
  int veg_type, bedrock_type, soil_type;
  double depth_to_bedrock;

  std::string line;
  std::getline(fin, line);

  // data structures for mesh
  std::vector<std::vector<int> > conn;
  PointFactory fac;
  
  while (std::getline(fin, line)) {
    std::istringstream ss(line);
    ss >> c[0] >> c[1]
       >> p1[0] >> p1[1] >> p1[2] >> depths[0]
       >> p2[0] >> p2[1] >> p2[2] >> depths[1]
       >> p3[0] >> p3[1] >> p3[2] >> depths[2]
       >> veg_type >> soil_type >> bedrock_type;

    veg_types.push_back(veg_type);
    soil_types.push_back(soil_type);
    bedrock_types.push_back(bedrock_type);

    std::vector<int> tri_conn(3);
    bool isnew = fac.addPoint(p1, tri_conn[0]);
    if (isnew) depths_to_bedrock.push_back(depths[0]);

    isnew = fac.addPoint(p2, tri_conn[1]);
    if (isnew) depths_to_bedrock.push_back(depths[1]);

    isnew = fac.addPoint(p3, tri_conn[2]);
    if (isnew) depths_to_bedrock.push_back(depths[2]);

    conn.emplace_back(tri_conn);
  }

  std::vector<std::vector<int> > sets(1);
  sets[0] = std::move(veg_types);
  Mesh2D m(std::move(fac.points), std::move(conn), std::move(sets));
  return m;
}

}
}
