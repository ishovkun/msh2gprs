#include "element.hpp"

ShapeIndexer::MapVTK ShapeIndexer::vtk_indices = {
    {TRGLE3,  5},
    {QUAD4,   9},
    {TETRA4,  10},
    {PRISM8,  12},
    {PRISM6,  13},
    {TRGLE6,  22},
    {QUAD8,   23},
    {PRISM20, 25},
    {PRISM15, 26}
};


StandardElements::StandardElements()
{

  elementProps.resize(10);

  /********************************************/
  /*          Triangle elements props         */
  /********************************************/
  elementProps[TRGLE3].vtkIndex = 5;
  elementProps[TRGLE3].nodesInElement = 3;
  elementProps[TRGLE3].facesInElement = 1;

  elementProps[TRGLE3].vvFacesNodes.resize(1);
  elementProps[TRGLE3].vvFacesNodes[0].push_back(0);
  elementProps[TRGLE3].vvFacesNodes[0].push_back(1);
  elementProps[TRGLE3].vvFacesNodes[0].push_back(2);

  /********************************************/
  /*          Triangle elements props         */
  /********************************************/
  elementProps[TRGLE6].vtkIndex = 22;
  elementProps[TRGLE6].nodesInElement = 6;
  elementProps[TRGLE6].facesInElement = 1;

  elementProps[TRGLE6].vvFacesNodes.resize(1);
  elementProps[TRGLE6].vvFacesNodes[0].push_back(0);
  elementProps[TRGLE6].vvFacesNodes[0].push_back(3);
  elementProps[TRGLE6].vvFacesNodes[0].push_back(1);
  elementProps[TRGLE6].vvFacesNodes[0].push_back(4);
  elementProps[TRGLE6].vvFacesNodes[0].push_back(2);
  elementProps[TRGLE6].vvFacesNodes[0].push_back(5);

  /********************************************/
  /*          Quad elements props             */
  /********************************************/
  elementProps[QUAD4].vtkIndex = 9;
  elementProps[QUAD4].nodesInElement = 4;
  elementProps[QUAD4].facesInElement = 1;

  elementProps[QUAD4].vvFacesNodes.resize(1);
  elementProps[QUAD4].vvFacesNodes[0].push_back(0);
  elementProps[QUAD4].vvFacesNodes[0].push_back(1);
  elementProps[QUAD4].vvFacesNodes[0].push_back(2);
  elementProps[QUAD4].vvFacesNodes[0].push_back(3);

  elementProps[QUAD8].vtkIndex = 23;
  elementProps[QUAD8].nodesInElement = 8;
  elementProps[QUAD8].facesInElement = 1;

  elementProps[QUAD8].vvFacesNodes.resize(1);
  elementProps[QUAD8].vvFacesNodes[0].push_back(0);
  elementProps[QUAD8].vvFacesNodes[0].push_back(4);
  elementProps[QUAD8].vvFacesNodes[0].push_back(1);
  elementProps[QUAD8].vvFacesNodes[0].push_back(5);
  elementProps[QUAD8].vvFacesNodes[0].push_back(2);
  elementProps[QUAD8].vvFacesNodes[0].push_back(6);
  elementProps[QUAD8].vvFacesNodes[0].push_back(3);
  elementProps[QUAD8].vvFacesNodes[0].push_back(7);

  /********************************************/
  /*   Tetrahedral elements props             */
  /********************************************/
  elementProps[TETRA4].vtkIndex = 10;
  elementProps[TETRA4].nodesInElement = 4;
  elementProps[TETRA4].facesInElement = 4;

  elementProps[TETRA4].vvFacesNodes.resize(4);
  elementProps[TETRA4].vvFacesNodes[0].push_back(0);
  elementProps[TETRA4].vvFacesNodes[0].push_back(1);
  elementProps[TETRA4].vvFacesNodes[0].push_back(2);

  elementProps[TETRA4].vvFacesNodes[1].push_back(0);
  elementProps[TETRA4].vvFacesNodes[1].push_back(1);
  elementProps[TETRA4].vvFacesNodes[1].push_back(3);

  elementProps[TETRA4].vvFacesNodes[2].push_back(1);
  elementProps[TETRA4].vvFacesNodes[2].push_back(2);
  elementProps[TETRA4].vvFacesNodes[2].push_back(3);

  elementProps[TETRA4].vvFacesNodes[3].push_back(2);
  elementProps[TETRA4].vvFacesNodes[3].push_back(0);
  elementProps[TETRA4].vvFacesNodes[3].push_back(3);

  /********************************************/
  /*   Tetrahedral elements props             */
  /********************************************/
  elementProps[TETRA10].vtkIndex = 24;
  elementProps[TETRA10].nodesInElement = 10;
  elementProps[TETRA10].facesInElement = 4;

  elementProps[TETRA10].vvFacesNodes.resize(4);
  elementProps[TETRA10].vvFacesNodes[0].push_back(0);
  elementProps[TETRA10].vvFacesNodes[0].push_back(4);
  elementProps[TETRA10].vvFacesNodes[0].push_back(1);
  elementProps[TETRA10].vvFacesNodes[0].push_back(5);
  elementProps[TETRA10].vvFacesNodes[0].push_back(2);
  elementProps[TETRA10].vvFacesNodes[0].push_back(6);

  elementProps[TETRA10].vvFacesNodes[1].push_back(0);
  elementProps[TETRA10].vvFacesNodes[1].push_back(4);
  elementProps[TETRA10].vvFacesNodes[1].push_back(1);
  elementProps[TETRA10].vvFacesNodes[1].push_back(8);
  elementProps[TETRA10].vvFacesNodes[1].push_back(3);
  elementProps[TETRA10].vvFacesNodes[1].push_back(7);

  elementProps[TETRA10].vvFacesNodes[2].push_back(1);
  elementProps[TETRA10].vvFacesNodes[2].push_back(5);
  elementProps[TETRA10].vvFacesNodes[2].push_back(2);
  elementProps[TETRA10].vvFacesNodes[2].push_back(9);
  elementProps[TETRA10].vvFacesNodes[2].push_back(3);
  elementProps[TETRA10].vvFacesNodes[2].push_back(8);

  elementProps[TETRA10].vvFacesNodes[3].push_back(2);
  elementProps[TETRA10].vvFacesNodes[3].push_back(6);
  elementProps[TETRA10].vvFacesNodes[3].push_back(0);
  elementProps[TETRA10].vvFacesNodes[3].push_back(7);
  elementProps[TETRA10].vvFacesNodes[3].push_back(3);
  elementProps[TETRA10].vvFacesNodes[3].push_back(9);

  /********************************************/
  /*   Hexahedral elements props              */
  /********************************************/
  elementProps[PRISM8].vtkIndex = 12;
  elementProps[PRISM8].nodesInElement = 8;
  elementProps[PRISM8].facesInElement = 6;

  elementProps[PRISM8].vvFacesNodes.resize(6);
  elementProps[PRISM8].vvFacesNodes[0].push_back(0);
  elementProps[PRISM8].vvFacesNodes[0].push_back(1);
  elementProps[PRISM8].vvFacesNodes[0].push_back(2);
  elementProps[PRISM8].vvFacesNodes[0].push_back(3);

  elementProps[PRISM8].vvFacesNodes[1].push_back(4);
  elementProps[PRISM8].vvFacesNodes[1].push_back(5);
  elementProps[PRISM8].vvFacesNodes[1].push_back(6);
  elementProps[PRISM8].vvFacesNodes[1].push_back(7);

  elementProps[PRISM8].vvFacesNodes[2].push_back(0);
  elementProps[PRISM8].vvFacesNodes[2].push_back(1);
  elementProps[PRISM8].vvFacesNodes[2].push_back(5);
  elementProps[PRISM8].vvFacesNodes[2].push_back(4);

  elementProps[PRISM8].vvFacesNodes[3].push_back(1);
  elementProps[PRISM8].vvFacesNodes[3].push_back(2);
  elementProps[PRISM8].vvFacesNodes[3].push_back(6);
  elementProps[PRISM8].vvFacesNodes[3].push_back(5);

  elementProps[PRISM8].vvFacesNodes[4].push_back(2);
  elementProps[PRISM8].vvFacesNodes[4].push_back(3);
  elementProps[PRISM8].vvFacesNodes[4].push_back(7);
  elementProps[PRISM8].vvFacesNodes[4].push_back(6);

  elementProps[PRISM8].vvFacesNodes[5].push_back(3);
  elementProps[PRISM8].vvFacesNodes[5].push_back(0);
  elementProps[PRISM8].vvFacesNodes[5].push_back(4);
  elementProps[PRISM8].vvFacesNodes[5].push_back(7);

  /**********************************************/
  /*        wedge elements props                */
  /*                                            */
  /* wedge = hexahedron with agglutinated nodes */
  /**********************************************/
  /*
                w
                ^
                |
                3                       3                      3
              ,/|`\                   ,/|`\                  ,/|`\
            ,/  |  `\               12  |  13              12  |  13
          ,/    |    `\           ,/    |    `\          ,/    |    `\
         4------+------5         4------14-----5        4------14-----5
         |      |      |         |      8      |        |      8      |
         |    ,/|`\    |         |      |      |        |    ,/|`\    |
         |  ,/  |  `\  |         |      |      |        |  15  |  16  |
         |,/    |    `\|         |      |      |        |,/    |    `\|
        ,|      |      |\        10     |      11       10-----17-----11
      ,/ |      0      | `\      |      0      |        |      0      |
     u   |    ,/ `\    |    v    |    ,/ `\    |        |    ,/ `\    |
         |  ,/     `\  |         |  ,6     `7  |        |  ,6     `7  |
         |,/         `\|         |,/         `\|        |,/         `\|
         1-------------2         1------9------2        1------9------2
  */

  elementProps[PRISM6].vtkIndex = 13;
  elementProps[PRISM6].nodesInElement = 6;
  elementProps[PRISM6].facesInElement = 5;

  elementProps[PRISM6].vvFacesNodes.resize(5);
  elementProps[PRISM6].vvFacesNodes[0].push_back(2);
  elementProps[PRISM6].vvFacesNodes[0].push_back(1);
  elementProps[PRISM6].vvFacesNodes[0].push_back(0);

  elementProps[PRISM6].vvFacesNodes[1].push_back(3);
  elementProps[PRISM6].vvFacesNodes[1].push_back(4);
  elementProps[PRISM6].vvFacesNodes[1].push_back(5);

  elementProps[PRISM6].vvFacesNodes[2].push_back(0);
  elementProps[PRISM6].vvFacesNodes[2].push_back(1);
  elementProps[PRISM6].vvFacesNodes[2].push_back(4);
  elementProps[PRISM6].vvFacesNodes[2].push_back(3);

  elementProps[PRISM6].vvFacesNodes[3].push_back(1);
  elementProps[PRISM6].vvFacesNodes[3].push_back(2);
  elementProps[PRISM6].vvFacesNodes[3].push_back(5);
  elementProps[PRISM6].vvFacesNodes[3].push_back(4);

  elementProps[PRISM6].vvFacesNodes[4].push_back(2);
  elementProps[PRISM6].vvFacesNodes[4].push_back(0);
  elementProps[PRISM6].vvFacesNodes[4].push_back(3);
  elementProps[PRISM6].vvFacesNodes[4].push_back(5);


  elementProps[PRISM15].vtkIndex = 26;
  elementProps[PRISM15].nodesInElement = 15;
  elementProps[PRISM15].facesInElement = 5;

  elementProps[PRISM15].vvFacesNodes.resize(5);
  elementProps[PRISM15].vvFacesNodes[0].push_back(2); elementProps[PRISM15].vvFacesNodes[0].push_back(9); elementProps[PRISM15].vvFacesNodes[0].push_back(1);
  elementProps[PRISM15].vvFacesNodes[0].push_back(6); elementProps[PRISM15].vvFacesNodes[0].push_back(0); elementProps[PRISM15].vvFacesNodes[0].push_back(7);

  elementProps[PRISM15].vvFacesNodes[1].push_back(5); elementProps[PRISM15].vvFacesNodes[1].push_back(14); elementProps[PRISM15].vvFacesNodes[1].push_back(4);
  elementProps[PRISM15].vvFacesNodes[1].push_back(12); elementProps[PRISM15].vvFacesNodes[1].push_back(3); elementProps[PRISM15].vvFacesNodes[1].push_back(13);

  elementProps[PRISM15].vvFacesNodes[2].push_back(0); elementProps[PRISM15].vvFacesNodes[2].push_back(6);
  elementProps[PRISM15].vvFacesNodes[2].push_back(1); elementProps[PRISM15].vvFacesNodes[2].push_back(10);
  elementProps[PRISM15].vvFacesNodes[2].push_back(4); elementProps[PRISM15].vvFacesNodes[2].push_back(12);
  elementProps[PRISM15].vvFacesNodes[2].push_back(3); elementProps[PRISM15].vvFacesNodes[2].push_back(8);

  elementProps[PRISM15].vvFacesNodes[3].push_back(1); elementProps[PRISM15].vvFacesNodes[3].push_back(9);
  elementProps[PRISM15].vvFacesNodes[3].push_back(2); elementProps[PRISM15].vvFacesNodes[3].push_back(11);
  elementProps[PRISM15].vvFacesNodes[3].push_back(5); elementProps[PRISM15].vvFacesNodes[3].push_back(14);
  elementProps[PRISM15].vvFacesNodes[3].push_back(4); elementProps[PRISM15].vvFacesNodes[3].push_back(10);

  elementProps[PRISM15].vvFacesNodes[4].push_back(2); elementProps[PRISM15].vvFacesNodes[4].push_back(7);
  elementProps[PRISM15].vvFacesNodes[4].push_back(0); elementProps[PRISM15].vvFacesNodes[4].push_back(8);
  elementProps[PRISM15].vvFacesNodes[4].push_back(3); elementProps[PRISM15].vvFacesNodes[4].push_back(13);
  elementProps[PRISM15].vvFacesNodes[4].push_back(5); elementProps[PRISM15].vvFacesNodes[4].push_back(11);

  /**********************************************/
  /*        prism 20		                */
  /**********************************************/
  /*
Hexahedron:             Hexahedron20:          Hexahedron27:

       v
3----------2            3----13----2           3----13----2
|\     ^   |\           |\         |\          |\         |\
| \    |   | \          | 15       | 14        |15    24  | 14
|  \   |   |  \         9  \       11 \        9  \ 20    11 \
|   7------+---6        |   7----19+---6       |   7----19+---6
|   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23|
0---+---\--1   |        0---+-8----1   |       0---+-8----1   |
 \  |    \  \  |         \  17      \  18       \ 17    25 \  18
  \ |     \  \ |         10 |        12|        10 |  21    12|
   \|      w  \|           \|         \|          \|         \|
    4----------5            4----16----5           4----16----5
*/
  elementProps[PRISM20].vtkIndex = 25;
  elementProps[PRISM20].nodesInElement = 20;
  elementProps[PRISM20].facesInElement = 6;

  elementProps[PRISM20].vvFacesNodes.resize(6);
  elementProps[PRISM20].vvFacesNodes[0].push_back(0);
  elementProps[PRISM20].vvFacesNodes[0].push_back(8);
  elementProps[PRISM20].vvFacesNodes[0].push_back(1);
  elementProps[PRISM20].vvFacesNodes[0].push_back(11);
  elementProps[PRISM20].vvFacesNodes[0].push_back(2);
  elementProps[PRISM20].vvFacesNodes[0].push_back(13);
  elementProps[PRISM20].vvFacesNodes[0].push_back(3);
  elementProps[PRISM20].vvFacesNodes[0].push_back(9);

  elementProps[PRISM20].vvFacesNodes[1].push_back(4);
  elementProps[PRISM20].vvFacesNodes[1].push_back(16);
  elementProps[PRISM20].vvFacesNodes[1].push_back(5);
  elementProps[PRISM20].vvFacesNodes[1].push_back(18);
  elementProps[PRISM20].vvFacesNodes[1].push_back(6);
  elementProps[PRISM20].vvFacesNodes[1].push_back(19);
  elementProps[PRISM20].vvFacesNodes[1].push_back(7);
  elementProps[PRISM20].vvFacesNodes[1].push_back(17);

  elementProps[PRISM20].vvFacesNodes[2].push_back(0);
  elementProps[PRISM20].vvFacesNodes[2].push_back(8);
  elementProps[PRISM20].vvFacesNodes[2].push_back(1);
  elementProps[PRISM20].vvFacesNodes[2].push_back(12);
  elementProps[PRISM20].vvFacesNodes[2].push_back(5);
  elementProps[PRISM20].vvFacesNodes[2].push_back(16);
  elementProps[PRISM20].vvFacesNodes[2].push_back(4);
  elementProps[PRISM20].vvFacesNodes[2].push_back(10);

  elementProps[PRISM20].vvFacesNodes[3].push_back(1);
  elementProps[PRISM20].vvFacesNodes[3].push_back(11);
  elementProps[PRISM20].vvFacesNodes[3].push_back(2);
  elementProps[PRISM20].vvFacesNodes[3].push_back(14);
  elementProps[PRISM20].vvFacesNodes[3].push_back(6);
  elementProps[PRISM20].vvFacesNodes[3].push_back(18);
  elementProps[PRISM20].vvFacesNodes[3].push_back(5);
  elementProps[PRISM20].vvFacesNodes[3].push_back(12);

  elementProps[PRISM20].vvFacesNodes[4].push_back(2);
  elementProps[PRISM20].vvFacesNodes[4].push_back(13);
  elementProps[PRISM20].vvFacesNodes[4].push_back(3);
  elementProps[PRISM20].vvFacesNodes[4].push_back(15);
  elementProps[PRISM20].vvFacesNodes[4].push_back(7);
  elementProps[PRISM20].vvFacesNodes[4].push_back(19);
  elementProps[PRISM20].vvFacesNodes[4].push_back(6);
  elementProps[PRISM20].vvFacesNodes[4].push_back(14);

  elementProps[PRISM20].vvFacesNodes[5].push_back(3);
  elementProps[PRISM20].vvFacesNodes[5].push_back(9);
  elementProps[PRISM20].vvFacesNodes[5].push_back(0);
  elementProps[PRISM20].vvFacesNodes[5].push_back(10);
  elementProps[PRISM20].vvFacesNodes[5].push_back(4);
  elementProps[PRISM20].vvFacesNodes[5].push_back(17);
  elementProps[PRISM20].vvFacesNodes[5].push_back(7);
  elementProps[PRISM20].vvFacesNodes[5].push_back(15);

}

StandardElements::~StandardElements()
{
}
