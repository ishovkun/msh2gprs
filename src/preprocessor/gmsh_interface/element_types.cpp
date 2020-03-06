#include "element_types.hpp"

namespace gprs_data {

// map gmsh element id to vtk id
std::vector<int> msh_id_to_vtk_id = {
    VTK_ID::InvalidElementID,     // [0]  does not exist
    VTK_ID::LineID,               // [1]  2-node edge
    VTK_ID::TriangleID,           // [2]  3-node triangle
    VTK_ID::QuadrangleID,         // [3]  4-node quadrangle
    VTK_ID::TetrahedronID,        // [4]  4-node tetrahedron
    VTK_ID::HexahedronID,         // [5]  8-node hexahedron
    VTK_ID::WedgeID,              // [6]  6-node prism
    VTK_ID::PyramidID,            // [7]  5-node pyramid
    VTK_ID::QuadraticEdgeID,      // [8]  3-node second order line
    VTK_ID::QuadraticTriangleID,  // [9]  6-node second order triangle
    VTK_ID::QuadraticQuadID,      // [10] 9-node second orderer quadrangle
    VTK_ID::QuadraticTetraID,     // [11] 10-node second order tetrahedron
    VTK_ID::InvalidElementID,     // [12] 27-node second order hexahedron
    VTK_ID::InvalidElementID,     // [13] 18-node second order prism
    VTK_ID::InvalidElementID,     // [14] 14-node second order pyramid
    VTK_ID::VertexID,             // [15] 1-node point
    VTK_ID::InvalidElementID,     // [16] 8-node second order quadrangle
    VTK_ID::QudraticHexahedronID, // [17] 20-node second order hexahedron
    VTK_ID::QuadraticWedgeID,     // [18] 15-node second order prism
    VTK_ID::QuadraticPyramidID,   // [19] 13-node second order pyramid
    VTK_ID::InvalidElementID,     // [20] 9-node third order incomplete triangle
    VTK_ID::InvalidElementID,     // [21] 10-node third order triangle
    VTK_ID::InvalidElementID,     // [22] 12-node fourth order incomplete triangle
    VTK_ID::InvalidElementID,     // [23] 15-node fourth order triangle
    VTK_ID::InvalidElementID,     // [24] 15-node fifth order incomplete triangle
    VTK_ID::InvalidElementID,     // [25] 21-node fifth order complete triangle
    VTK_ID::InvalidElementID,     // [26] 4-node third order edge
    VTK_ID::InvalidElementID,     // [27] 5-node fourth order edge
    VTK_ID::InvalidElementID,     // [28] 6-node fifth order edge
    VTK_ID::InvalidElementID,     // [29] 20-node third order tetrahedron
    VTK_ID::InvalidElementID,     // [30] 35-node fourth order tetrahedron
    VTK_ID::InvalidElementID      // [31] 56-node fifth order tetrahedron
};

std::vector<GmshElementType> gmsh_element_types =
{
 GmshElementType::invalid_element, // [0]  does not exist
 GmshElementType::edge,            // [1]  2-node edge
 GmshElementType::face,            // [2]  3-node triangle
 GmshElementType::face,            // [3]  4-node quadrangle
 GmshElementType::cell,            // [4]  4-node tetrahedron
 GmshElementType::cell,            // [5]  8-node hexahedron
 GmshElementType::cell,            // [6]  6-node prism
 GmshElementType::cell,            // [7]  5-node pyramid
 GmshElementType::invalid_element, // [8]  3-node second order line
 GmshElementType::invalid_element, // [9]  6-node second order triangle
 GmshElementType::invalid_element, // [10] 9-node second orderer quadrangle
 GmshElementType::invalid_element, // [11] 10-node second order tetrahedron
 GmshElementType::invalid_element, // [12] 27-node second order hexahedron
 GmshElementType::invalid_element, // [13] 18-node second order prism
 GmshElementType::invalid_element, // [14] 14-node second order pyramid
 GmshElementType::node,            // [15] 1-node point
 GmshElementType::invalid_element, // [16] 8-node second order quadrangle
 GmshElementType::invalid_element, // [17] 20-node second order hexahedron
 GmshElementType::invalid_element, // [18] 15-node second order prism
 GmshElementType::invalid_element, // [19] 13-node second order pyramid
 GmshElementType::invalid_element, // [20] 9-node third order incomplete triangle
 GmshElementType::invalid_element, // [21] 10-node third order triangle
 GmshElementType::invalid_element, // [22] 12-node fourth order incomplete triangle
 GmshElementType::invalid_element, // [23] 15-node fourth order triangle
 GmshElementType::invalid_element, // [24] 15-node fifth order incomplete triangle
 GmshElementType::invalid_element, // [25] 21-node fifth order complete triangle
 GmshElementType::invalid_element, // [26] 4-node third order edge
 GmshElementType::invalid_element, // [27] 5-node fourth order edge
 GmshElementType::invalid_element, // [28] 6-node fifth order edge
 GmshElementType::invalid_element, // [29] 20-node third order tetrahedron
 GmshElementType::invalid_element, // [30] 35-node fourth order tetrahedron
 GmshElementType::invalid_element  // [31] 56-node fifth order tetrahedron
};

std::vector<size_t> gmsh_element_nvertices =
{
 0,  // [0]  does not exist
 2,  // [1]  2-node edge
 3,  // [2]  3-node triangle
 4,  // [3]  4-node quadrangle
 4,  // [4]  4-node tetrahedron
 8,  // [5]  8-node hexahedron
 6,  // [6]  6-node prism
 5,  // [7]  5-node pyramid
 3,  // [8]  3-node second order line
 6,  // [9]  6-node second order triangle
 9,  // [10] 9-node second orderer quadrangle
 10, // [11] 10-node second order tetrahedron
 27, // [12] 27-node second order hexahedron
 18, // [13] 18-node second order prism
 14, // [14] 14-node second order pyramid
 1,  // [15] 1-node point
 8,  // [16] 8-node second order quadrangle
 20, // [17] 20-node second order hexahedron
 15, // [18] 15-node second order prism
 13, // [19] 13-node second order pyramid
 20, // [20] 9-node third order incomplete triangle
 10, // [21] 10-node third order triangle
 22, // [22] 12-node fourth order incomplete triangle
 15, // [23] 15-node fourth order triangle
 15, // [24] 15-node fifth order incomplete triangle
 21, // [25] 21-node fifth order complete triangle
 4,  // [26] 4-node third order edge
 5,  // [27] 5-node fourth order edge
 6,  // [28] 6-node fifth order edge
 20, // [29] 20-node third order tetrahedron
 35, // [30] 35-node fourth order tetrahedron
 56  // [31] 56-node fifth order tetrahedron
};

}  // end namespace gprs_data
