#ifndef STANDARD_FEM_DATA
#define STANDARD_FEM_DATA

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <unordered_map>

using namespace std;

enum EnumElementsForm
{
    TRGLE3  = 0,
    TRGLE6  = 1,
    QUAD4   = 2,
    TETRA4  = 3,
    TETRA10 = 4,
    PRISM6  = 5,
    PRISM8  = 6,
    QUAD8   = 7,
    PRISM15 = 8,
    PRISM20 = 9
};


class ShapeIndexer
{
 public:
  static int get_vtk_index(const int shape_index) {return vtk_indices[shape_index];}

  typedef std::unordered_map<int,int> MapVTK;
  static MapVTK vtk_indices;
};


struct StdElement
{
    int vtkIndex;
    int nodesInElement;
    int facesInElement;

    vector<vector<int> > vvFacesNodes;
};

class StandardElements
{
public:
    StandardElements();
   ~StandardElements();
    vector<StdElement> elementProps;

};
#endif
