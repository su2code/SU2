#include "catch.hpp"
#include <sstream>
#include <iomanip>
#include "../../../Common/include/geometry/dual_grid/CEdge.hpp"
#include "../../../Common/include/geometry/dual_grid/CVertex.hpp"

TEST_CASE("Volume Computation", "[Dual Grid]") {

  const int nDim = 3;
  
  su2double Coord_FaceiPoint[3];
  su2double Coord_FaceElem_CG[3];
  su2double Coord_Elem_CG[3];
  su2double Coord_Edge_CG[3];
  su2double scaling = 100;

  Coord_FaceiPoint[0]  = scaling*0.664995; Coord_FaceiPoint[1]  = scaling*1.1462;  Coord_FaceiPoint[2]  = scaling*0.00926223;
  Coord_FaceElem_CG[0] = scaling*0.655997; Coord_FaceElem_CG[1] = scaling*1.13054; Coord_FaceElem_CG[2] = scaling*0.00945181;
  Coord_Elem_CG[0]     = scaling*0.653846; Coord_Elem_CG[1]     = scaling*1.12927; Coord_Elem_CG[2]     = scaling*0.00835789;
  Coord_Edge_CG[0]     = scaling*0.664943; Coord_Edge_CG[1]     = scaling*1.14623; Coord_Edge_CG[2]     = scaling*0.00935524;

  CEdge edge(0, 1, nDim);
  SECTION("2D Edge"){
    su2double volume = edge.GetVolume(Coord_FaceiPoint, Coord_Edge_CG, Coord_Elem_CG);    
    REQUIRE(volume == Approx(0.00607415));
  }
  
  SECTION("3D Edge"){
    su2double volume = edge.GetVolume(Coord_FaceiPoint, Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);    
    REQUIRE(volume == Approx(0.000546832));
  }
  
  CVertex vertex(0, nDim);
  SECTION("2D Vertex"){
    vertex.SetNodes_Coord(Coord_Edge_CG, Coord_Elem_CG);
    REQUIRE(vertex.GetNormal()[0] == Approx(-1.696));
    REQUIRE(vertex.GetNormal()[1] == Approx(1.1097));
  }
  SECTION("3D Vertex"){
    vertex.SetNodes_Coord(Coord_Edge_CG, Coord_FaceElem_CG, Coord_Elem_CG);
    REQUIRE(vertex.GetNormal()[0] == Approx(-0.0864312));
    REQUIRE(vertex.GetNormal()[1] == Approx(0.0499696));
    REQUIRE(vertex.GetNormal()[2] == Approx(0.111938));
  }
  
}
