#include "simdata.hpp"

// library for analytical geometry
#include <angem/Point.hpp>
#include <angem/Rectangle.hpp>
#include <angem/PointSet.hpp>
#include "angem/CollisionGJK.hpp"
#include <angem/Collisions.hpp>
#include <angem/utils.hpp>
#include <mesh/utils.hpp> // to remesh embedded fractures
#include <mesh/Mesh.hpp> // 3D mesh format
// parser for user-defined expressions for reservoir data
#include <muparser/muParser.h>

// profiling
// #include <ctime>

#define SPECIAL_CELL = 999
#include <algorithm>
#include <exception>
#include <unordered_set>

using Point = angem::Point<3, double>;
const int MARKER_BELOW_FRAC = 0;
const int MARKER_ABOVE_FRAC = 1;
const int MARKER_FRAC = 2;


// SimData::SimData(const string & inputstream, const SimdataConfig & config)
SimData::SimData(mesh::Mesh & grid, const SimdataConfig & config)
    :
    grid(grid),
    config(config)
{
//   nNodes = 0;
//   nCells = 0;
// //
//   dNotNumber = -999.999;

//   vPointPass.push_back(0);
//   vPointPass.push_back(2);
//   vPointPass.push_back(1);

//   vPointCoord.resize(3, 0.0);

//   vvPlate.resize(3);
//   for (int i = 0; i < 3; i++) vvPlate[i].resize(3, 0.0);

//   // Boundary conditions and initial state
//   Sxx = 0.0;// 399.0= 0.748 * Szz bar [default: 0]
//   Syy = 0.0; // 424.0 = 0.795 * Szz bar [default: 300.0]
//   Szz = 0.0; // 533 = 2172*9.81*2500*1e-5 (rho*g*h) bar [default: 600]
//   Syz = 0.0;
//   Sxz = 0.0;
//   Sxy = 0.0;

//   //wells
//   nWells = 2;
//   vsWell.resize(nWells);

//   // well 1 at the footwall block
//   vsWell[0].vWellCoordinate.clear();
//   vsWell[0].vWellCoordinate.push_back(-250.0); // x
//   vsWell[0].vWellCoordinate.push_back(-0.0); // 45-30-90 model
//   vsWell[0].vWellCoordinate.push_back(-1450.0); // z0
//   vsWell[0].vWellCoordinate.push_back(-1550.0); // z1
//   vsWell[0].Type = "WCONPROD";
//   vsWell[0].radius_poisk = 6.0; // m

//   // well 2 at the hangingwall block
//   vsWell[1].vWellCoordinate.clear();
//   vsWell[1].vWellCoordinate.push_back(250.0); // x
//   vsWell[1].vWellCoordinate.push_back(-0.0); // 45-30-90 model
//   vsWell[1].vWellCoordinate.push_back(-1450.0); // z0
//   vsWell[1].vWellCoordinate.push_back(-1550.0); // z1
//   vsWell[1].Type = "WCONPROD";
//   vsWell[1].radius_poisk = 6.0; // m

  // Kirill's renumbering
  // pRenum = new renum();
}

SimData::~SimData()
{
}


void SimData::defineEmbeddedFractureProperties()
{
  std::size_t ef_ind = 0;
  // class that checks if shapes collide
  angem::CollisionGJK<double> collision;
  // non-const since fracture is adjusted to avoid collision with vertices
  for (auto & frac_conf : config.fractures)
  {
    vEfrac.emplace_back();
    auto & frac = vEfrac.back();
    Point total_shift = {0, 0, 0};

 redo_collision:
    std::cout << "fracture polygon" << std::endl;
    std::cout << *(frac_conf.body) << std::endl;

    // find cells intersected by the fracture
    frac.cells.clear();
    for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
    {
      const auto poly_cell = cell.polyhedron();

      if (collision.check(*frac_conf.body, poly_cell))
      {
        frac.cells.push_back(cell.index());

        // check if some vertices are too close to the fracture
        // and if so move a fracture a little bit
        // for (const auto & ivertex : cell.vVertices)
        for (const auto & vertex : poly_cell.get_points())
        {
          // const auto & vert = vvVrtxCoords[ivertex];
          const auto vc = poly_cell.center() - vertex;
          if ( fabs(frac_conf.body->plane.distance(vertex)/vc.norm()) < 1e-4 )
          {
            // shift in the direction perpendicular to fracture
            const double h = (poly_cell.get_points()[1] -
                              poly_cell.get_points()[0]).norm();
            const Point shift = h/5 * frac_conf.body->plane.normal();
            total_shift += shift;
            std::cout << "shifting fracture: " << shift ;
            std::cout << " due to collision with vertex: " << vertex;
            std::cout << std::endl;
            frac_conf.body->move(shift);
            goto redo_collision;
          }
        }  // end adjusting fracture
      }  // end collision processing
    }    // end cell loop

    std::cout << "Total shift = " << total_shift << std::endl;
    const std::size_t n_efrac_cells = frac.cells.size();
    std::cout << "fracture " << ef_ind
              << " occupies " << n_efrac_cells << " cells"
              << std::endl;
    if (n_efrac_cells == 0)
    {
      vEfrac.pop_back();
      continue;
    }

    // fill out properties
    std::cout << "n_efrac_cells = " << n_efrac_cells << std::endl;

    vEfrac[ef_ind].points.assign(n_efrac_cells,
                                 frac_conf.body->center());
    vEfrac[ef_ind].dip.assign(n_efrac_cells,
                              frac_conf.body->plane.dip_angle());
    vEfrac[ef_ind].strike.assign(n_efrac_cells,
                                 frac_conf.body->plane.strike_angle());

    vEfrac[ef_ind].cohesion       = frac_conf.cohesion;
    vEfrac[ef_ind].friction_angle = frac_conf.friction_angle;
    vEfrac[ef_ind].dilation_angle = frac_conf.dilation_angle;
    vEfrac[ef_ind].aperture       = frac_conf.aperture;
    vEfrac[ef_ind].conductivity   = frac_conf.conductivity;

    ef_ind++;
  }  // end efracs loop

}


void SimData::computeCellClipping()
{
  // determine points of intersection of embedded fractures with
  // the mesh

  // criterion for point residing on the plane
  const double tol = 1e-8;

  // too lazy to account for fractures not collided with any cells
  assert(config.fractures.size() == vEfrac.size());

  for (std::size_t ifrac=0; ifrac<config.fractures.size(); ++ifrac)
  {
    const auto & frac_cells = vEfrac[ifrac].cells;
    const auto & frac_plane = config.fractures[ifrac].body->plane;

    std::vector<std::vector<angem::Point<3,double>>> vvSection;
    vvSection.resize(frac_cells.size());
    std::vector<angem::PolyGroup<double>> splits(frac_cells.size());

    /* loop faces:
     * if any neighbor cell is in collision list,
     * determine the intersection points of the face with the fracture plane
     * add these points to the point set for cells
     */
    for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
    {
      // vector of cells containing efrac and neighboring face
      std::vector<std::size_t> v_neighbors;
      for (const std::size_t & ineighbor : face.neighbors())
      {
        const std::size_t frac_cell_local_ind = find(ineighbor, frac_cells);
        if (frac_cell_local_ind != frac_cells.size())
          v_neighbors.push_back(frac_cell_local_ind);
      }

      if (v_neighbors.size() > 0)
      {
        // construct polygon and determine intersection points
        angem::Polygon<double> poly_face(face.vertices());
        std::vector<Point> section;
        angem::collision(poly_face, frac_plane, section);

        // no intersection
        // we still need to add polygon into splits for transmissibility
        if (section.size() < 2)
        {
          angem::PolyGroup<double> split;
          angem::split(poly_face, frac_plane, split,
                       MARKER_BELOW_FRAC, MARKER_ABOVE_FRAC);
          angem::Polygon<double>::reorder_indices(split.vertices.points,
                                                  split.polygons[0]);
          for (const auto & ineighbor : v_neighbors)
            splits[ineighbor].add(split);
          continue;
        }

        // save intersection data into neighbor fracture cells
        for (const auto & ineighbor : v_neighbors)
          for (const auto & p : section)
            vvSection[ineighbor].push_back(p);

        // build polygons from intersection and save to scratch
        angem::PolyGroup<double> split;
        angem::split(poly_face, frac_plane, split,
                     MARKER_BELOW_FRAC, MARKER_ABOVE_FRAC);
        angem::Polygon<double>::reorder_indices(split.vertices.points,
                                                split.polygons[0]);
        angem::Polygon<double>::reorder_indices(split.vertices.points,
                                                split.polygons[1]);
        // add split to neighbor cess splits
        for (const auto & ineighbor : v_neighbors)
          splits[ineighbor].add(split);

      }  // end if has ef cells neighbors

    }  // end face loop

    angem::PointSet<3,double> setVert(tol);
    mesh::SurfaceMesh<double> frac_msh(1e-6);
    for (std::size_t i=0; i<vEfrac[ifrac].cells.size(); ++i)
    {
      // loop through sda cells
      auto & section_points = vvSection[i];

      // some point among those we obtain in the previous part of code
      // are duplicated since two adjacent faces intersecting a plane
      // have one point in common
      std::vector<Point> set_points;
      angem::remove_duplicates(section_points, set_points, tol);

      // correct ordering for quads
      // if (set_points.size() > 3)
      // {
      //   angem::Polygon<double>::reorder(set_points);
      // }
      angem::Polygon<double> poly_section(set_points);
      vvSection[i] = poly_section.get_points();
      frac_msh.insert(poly_section);

      // remove cell if number of points < 3 <=> area = 0
      if (set_points.size() < 3)
      {
        std::cout << "erasing fracture cell" << vEfrac[ifrac].cells[i] << std::endl;
        vEfrac[ifrac].cells.erase(vEfrac[ifrac].cells.begin() + i);
        vvSection.erase(vvSection.begin() + i);
        vEfrac[ifrac].points.erase(vEfrac[ifrac].points.begin() + i);
        vEfrac[ifrac].strike.erase(vEfrac[ifrac].strike.begin() + i);
        vEfrac[ifrac].dip.erase(vEfrac[ifrac].dip.begin() + i);
        i--;
        continue;
      }

      // add fracture polygon to splits to compute transes
      splits[i].add(angem::Polygon<double>(set_points), MARKER_FRAC);

      // write points into a global set so we have an ordered set
      // of vertices and we can retreive indices
      for (const Point & p : set_points)
        setVert.insert(p);

    }  // end sda cells loop

    vEfrac[ifrac].mesh = std::move(frac_msh);

    // get indices of frac vertices
    // vEfrac[ifrac].vIndices.resize(vvSection.size());
    // std::size_t icell = 0;
    // for (const auto & cell_section : vvSection)
    // {
    //   for (const Point & p : cell_section)
    //   {
    //     const std::size_t ind = setVert.find(p);
    //     vEfrac[ifrac].vIndices[icell].push_back(ind);
    //   }
    //   icell++;
    // }

    // convert set of vertices to std::vector
    // vEfrac[ifrac].vVertices.resize(setVert.size());
    // std::size_t ivert = 0;
    // for (const Point & p : setVert)
    // {
    //   // vEfrac[ifrac].vVertices[ivert] = p;
    //   ivert++;
    // }

    // std::cout << "computing edfm transes" << std::endl;
    computeEDFMTransmissibilities(splits, ifrac);


//     for(int iface = 0; iface < vsFaceCustom.size(); iface++)
//     {
//       // vector of cells containing efrac and neighboring face
//       std::vector<std::size_t> v_neighbors;
//       for (const std::size_t & ineighbor : vsFaceCustom[iface].vNeighbors)
//       {
//         const std::size_t frac_cell_local_ind = find(ineighbor, frac_cells);
//         if (frac_cell_local_ind != frac_cells.size())
//           v_neighbors.push_back(frac_cell_local_ind);
//       }

//       if (v_neighbors.size() > 0)
//       {
//         // construct polygon and determine intersection points
//         angem::Polygon<double> poly_face(vvVrtxCoords,
//                                          vsFaceCustom[iface].vVertices);
//         std::vector<Point> section;
//         angem::collision(poly_face, frac_plane, section);

//         // no intersection
//         // we still need to add polygon into splits for transmissibility
//         if (section.size() < 2)
//         {
//           angem::PolyGroup<double> split;
//           angem::split(poly_face, frac_plane, split,
//                        MARKER_BELOW_FRAC, MARKER_ABOVE_FRAC);
//           angem::Polygon<double>::reorder_indices(split.vertices.points,
//                                                   split.polygons[0]);
//           for (const auto & ineighbor : v_neighbors)
//             splits[ineighbor].add(split);
//           continue;
//         }

//         // save intersection data into neighbor fracture cells
//         for (const auto & ineighbor : v_neighbors)
//           for (const auto & p : section)
//             vvSection[ineighbor].push_back(p);

//         // build polygons from intersection and save to scratch
//         angem::PolyGroup<double> split;
//         angem::split(poly_face, frac_plane, split,
//                      MARKER_BELOW_FRAC, MARKER_ABOVE_FRAC);
//         angem::Polygon<double>::reorder_indices(split.vertices.points,
//                                                 split.polygons[0]);
//         angem::Polygon<double>::reorder_indices(split.vertices.points,
//                                                 split.polygons[1]);
//         // add split to neighbor cess splits
//         for (const auto & ineighbor : v_neighbors)
//           splits[ineighbor].add(split);

//       }  // end if has ef cells neighbors
//     }    // end face loop

//     angem::PointSet<3,double> setVert(tol);
//     mesh::SurfaceMesh<double> frac_msh(1e-6, /* max_edges = */ nCells);
//     for (std::size_t i=0; i<vEfrac[ifrac].cells.size(); ++i)
//     {
//       // loop through sda cells
//       auto & section_points = vvSection[i];

//       // some point among those we obtain in the previous part of code
//       // are duplicated since two adjacent faces intersecting a plane
//       // have one point in common
//       std::vector<Point> set_points;
//       angem::remove_duplicates(section_points, set_points, tol);

//       // correct ordering for quads
//       // if (set_points.size() > 3)
//       // {
//       //   angem::Polygon<double>::reorder(set_points);
//       // }
//       angem::Polygon<double> poly_section(set_points);
//       vvSection[i] = poly_section.get_points();
//       frac_msh.insert(poly_section);

//       // remove cell if number of points < 3 <=> area = 0
//       if (set_points.size() < 3)
//       {
//         std::cout << "erasing fracture cell" << vEfrac[ifrac].cells[i] << std::endl;
//         vEfrac[ifrac].cells.erase(vEfrac[ifrac].cells.begin() + i);
//         vvSection.erase(vvSection.begin() + i);
//         vEfrac[ifrac].points.erase(vEfrac[ifrac].points.begin() + i);
//         vEfrac[ifrac].strike.erase(vEfrac[ifrac].strike.begin() + i);
//         vEfrac[ifrac].dip.erase(vEfrac[ifrac].dip.begin() + i);
//         i--;
//         continue;
//       }

//       // add fracture polygon to splits to compute transes
//       // @TODO: add marker to consider this an active poly
//       splits[i].add(angem::Polygon<double>(set_points), MARKER_FRAC);

//       // write points into a global set so we have an ordered set
//       // of vertices and we can retreive indices
//       for (const Point & p : set_points)
//         setVert.insert(p);

//     }

//     vEfrac[ifrac].mesh = std::move(frac_msh);

//     // get indices of frac vertices
//     vEfrac[ifrac].vIndices.resize(vvSection.size());
//     std::size_t icell = 0;
//     for (const auto & cell_section : vvSection)
//     {
//       for (const Point & p : cell_section)
//       {
//         const std::size_t ind = setVert.find(p);
//         vEfrac[ifrac].vIndices[icell].push_back(ind);
//       }
//       icell++;
//     }

//     // convert set of vertices to std::vector
//     vEfrac[ifrac].vVertices.resize(setVert.size());
//     std::size_t ivert = 0;
//     for (const Point & p : setVert)
//     {
//       vEfrac[ifrac].vVertices[ivert] = p;
//       ivert++;
//     }

//     // std::cout << "computing edfm transes" << std::endl;
//     computeEDFMTransmissibilities(splits, ifrac);
  }  // end efrac loop


//   // clear cell connections
//   // for (std::size_t f=0; f<vEfrac.size(); ++f)
//   // {
//   //   const auto & efrac = vEfrac[f];
//   //   for (std::size_t i=0; i<efrac.cells.size(); ++i)
//   //     for (std::size_t j=i; j<efrac.cells.size(); ++j)
//   //       if (flow_data.connection_exists(efrac.cells[i], efrac.cells[j]))
//   //       {
//   //         const std::size_t dead_conn = flow_data.connection_index(efrac.cells[i], efrac.cells[j]);
//   //         // const double T_ij = flow_data.trans_ij[dead_conn];
//   //         flow_data.clear_connection(efrac.cells[i], efrac.cells[j]);
//   //         const std::size_t con_fi =
//   //             flow_data.connection_index(efrac.cells[i],
//   //                                        get_flow_element_index(f, i) );
//   //         const double T_fi = flow_data.trans_ij[con_fi];

//   //         flow_data.insert_connection(get_flow_element_index(f, i),
//   //                                     efrac.cells[j]);
//   //         flow_data.trans_ij.push_back(T_fi);

//   //         const std::size_t con_fj =
//   //             flow_data.connection_index(efrac.cells[j],
//   //                                        get_flow_element_index(f, j) );
//   //         const double T_fj = flow_data.trans_ij[con_fj];

//   //         flow_data.insert_connection(get_flow_element_index(f, j),
//   //                                     efrac.cells[i]);
//   //         flow_data.trans_ij.push_back(T_fj);
//   //       }
//   // }

//   // for (std::size_t f=0; f<vEfrac.size(); ++f)
//   // {
//   //   const auto & efrac = vEfrac[f];
//   //   for (std::size_t i=0; i<efrac.cells.size(); ++i)
//   //   {
//   //     // const auto neighbors = flow_data.v_neighbors[get_flow_element_index( f, efrac.cells[i] )];
//   //     const auto neighbors = flow_data.v_neighbors[efrac.cells[i]];
//   //     for (const auto & neighbor : neighbors)
//   //     {
//   //       if (neighbor < nCells)
//   //       {
//   //         const std::size_t dead_conn = flow_data.connection_index(neighbor, efrac.cells[i]);
//   //         const double T_ij = flow_data.trans_ij[dead_conn];
//   //         flow_data.clear_connection(neighbor, efrac.cells[i]);
//   //         flow_data.insert_connection(neighbor, get_flow_element_index(f, i));
//   //         flow_data.trans_ij.push_back(T_ij);

//   //       }
//   //     }
//   //   }
//   // }

}


// void SimData::mergeSmallFracCells()
// {
//   for (std::size_t ifrac=0; ifrac<config.fractures.size(); ++ifrac)
//   {
//     auto & msh = vEfrac[ifrac].mesh;

//     double max_area = 0;
//     for (const auto & element : msh.polygons)
//     {
//       angem::Polygon<double> poly(msh.vertices.points, element);
//       const double area = poly.area();
//       max_area = std::max(area, max_area);
//     }

//     // merge tiny cells
//     std::size_t ielement = 0;
//     std::size_t n_frac_elements = msh.polygons.size();

//     // loop is with variable upper limit since elements can be
//     // merged and deleted
//     while (ielement < n_frac_elements)
//     {
//       angem::Polygon<double> poly(msh.vertices.points,
//                                   msh.polygons[ielement]);
//       const double area_factor = poly.area() / max_area;

//       if (area_factor < config.frac_cell_elinination_factor)
//       {
//         const std::size_t global_ielement = get_flow_element_index(ifrac, ielement);
//         const std::size_t new_element = msh.merge_element(ielement);
//         // update flow data
//         flow_data.merge_elements(get_flow_element_index(ifrac, new_element),
//                                  global_ielement);

//         n_frac_elements = msh.polygons.size();
//         if (ielement >= n_frac_elements)
//           break;
//         continue;
//       }
//       ielement++;
//     }

//   }
// }


void SimData::computeReservoirTransmissibilities()
{
  // init tran
  CalcTranses calc;
  calc.NbNodes     = grid.n_vertices();
  calc.NbPolyhedra = grid.n_cells();
  calc.NbPolygons  = grid.n_faces();
  calc.NbFracs     = n_flow_dfm_faces;
  calc.NbZones     = n_flow_dfm_faces + grid.n_cells();
  calc.NbOptions   = 1;
  calc.fracporo    = 1.0;
  calc.init();

  // fill data
  // nodes
  for ( std::size_t i = 0; i < grid.n_vertices(); i++ )
  {
    Point vertex = grid.vertices[i];
    calc.vCoordinatesX[i] = vertex.x();
    calc.vCoordinatesY[i] = vertex.y();
    calc.vCoordinatesZ[i] = vertex.z();
  }

  // polygons (2d elements)
  int code_polygon = 0;
  calc.vvFNodes.resize(grid.n_faces());

  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    const std::size_t ipoly = face.index();
    calc.vvFNodes[ipoly] = face.vertex_indices();

    if (is_fracture(face.marker()))
    {
      calc.vCodePolygon[ipoly] = code_polygon;
      code_polygon++;
    }
    else  // non-frac faces
      calc.vCodePolygon[ipoly] = -1;
  }

  // polyhedra (3d elements)
  calc.vvVFaces.resize(grid.n_faces());
  calc.vCodePolyhedron.resize(grid.n_cells());
  for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
  {
    const auto faces = cell.faces();
    const std::size_t icell = cell.index();
    calc.vvVFaces[icell].reserve(faces.size());
    for (const mesh::face_iterator & face : faces)
    {
      // TODO: get face indices some other way
      // if (dfm_faces.find(hash) == dfm_faces.end())
      //   throw std::out_of_range("bad cell hash");
      const std::size_t ipolygon = face.index();
      calc.vvVFaces[icell].push_back(ipolygon);
    }
    calc.vCodePolyhedron[icell] = dfm_faces.size() + icell;
  }

  // Properties
  calc.vZPermeability.assign(calc.NbZones * 3, 0.0);
  calc.vZConduction.assign( (calc.NbPolyhedra + n_flow_dfm_faces) * 3, 0.0);

  // DFM fractures
  for (auto it : dfm_faces)
  {
    const int i = it.second.nfluid;
    if (i > 0)  // if active
    {
      calc.vZoneCode[i] = i;
      calc.vZVolumeFactor[i] = it.second.aperture;
      calc.vZPorosity[i] = 1.0;
      calc.vZPermCode[i] = 1;

      //@HACK default permeability for all fractures
      const double perm = it.second.conductivity / it.second.aperture;
      calc.vZPermeability[i*3+0] = perm;
      calc.vZPermeability[i*3+1] = perm;
      calc.vZPermeability[i*3+2] = perm;

      calc.vZConduction[i*3+0] = 1;
      calc.vZConduction[i*3+1] = 1;
      calc.vZConduction[i*3+2] = 1;
      calc.vTimurConnectionFactor[i] = 1.0;
    }
  }

  // properties regular cells
  for ( std::size_t i = 0; i < calc.NbPolyhedra; i++ )
  {
    const std::size_t n = i + n_flow_dfm_faces;
    calc.vZoneCode[n] = calc.vCodePolyhedron[i];
    calc.vZPorosity[n] = get_property(i, "PORO");
    calc.vZPermCode[n] = 1;

    const angem::Point<3,double> perm = get_permeability(i);
    calc.vZPermeability[n*3+0] = perm[0];
    calc.vZPermeability[n*3+1] = perm[1];
    calc.vZPermeability[n*3+2] = perm[2];

    double thc = 0;
    try
    {
      thc = get_property(i, "THCROCK");
    }
    catch (const std::out_of_range& e)
    {
      calc.vZConduction[n*3+0] = thc;
      calc.vZConduction[n*3+1] = thc;
      calc.vZConduction[n*3+2] = thc;
    }

    calc.vZVolumeFactor[n] = get_volume_factor(i);
    calc.vTimurConnectionFactor[n] = 1.0;
  }


  FlowData matrix_flow_data;
  std::cout << "compute karimi" << std::endl;
  calc.compute_flow_data();
  std::cout << "end compute karimi" << std::endl;
  calc.extractData(matrix_flow_data);

  // copy to global
  flow_data.volumes.reserve(matrix_flow_data.volumes.size());
  flow_data.poro.reserve(matrix_flow_data.volumes.size());
  flow_data.depth.reserve(matrix_flow_data.volumes.size());
  for (std::size_t i=0; i<matrix_flow_data.volumes.size(); ++i)
  {
    flow_data.volumes.push_back(matrix_flow_data.volumes[i]);
    flow_data.poro.push_back(matrix_flow_data.poro[i]);
    flow_data.depth.push_back(matrix_flow_data.depth[i]);
  }

  for (const auto & conn : matrix_flow_data.map_connection)
  {
    const std::size_t iconn = conn.second;
    const auto element_pair = matrix_flow_data.invert_hash(conn.first);
    flow_data.insert_connection(element_pair.first, element_pair.second);
    flow_data.trans_ij.push_back(matrix_flow_data.trans_ij[iconn]);
  }


  // save custom user-defined cell data for flow output
  const std::size_t n_vars = rockPropNames.size();

  // save flow variable names
  flow_data.custom_names.clear();
  for (std::size_t j=0; j<n_vars; ++j)
    if (config.expression_type[j] == 0)
      flow_data.custom_names.push_back(rockPropNames[j]);

  // save values
  flow_data.custom_data.resize(n_flow_dfm_faces + grid.n_cells());
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
    if (is_fracture(face.marker()))
    {
      const std::size_t icell = face.neighbors()[0];
      for (std::size_t j=0; j<n_vars; ++j)
        if (config.expression_type[j] == 0)
        {
          const std::size_t ielement = dfm_faces[face.index()].nfluid;
          flow_data.custom_data[ielement].push_back(vsCellRockProps[icell].v_props[j]);
        }
    }

  for (std::size_t i=0; i<grid.n_cells(); ++i)
  {
    // flow_data.custom_data.emplace_back();
    const std::size_t ielement = n_flow_dfm_faces + i;
    for (std::size_t j=0; j<n_vars; ++j)
      if (config.expression_type[j] == 0)
        flow_data.custom_data[ielement].push_back( vsCellRockProps[i].v_props[j] );
  }

  // new_flow_data = flow_data;
}


std::size_t SimData::get_flow_element_index(const std::size_t ifrac,
                                            const std::size_t ielement) const
{
  std::size_t result = n_flow_dfm_faces + grid.n_cells();
  for (std::size_t i=0; i<ifrac; ++i)
  {
    result += vEfrac[i].mesh.polygons.size();
  }

  result += ielement;
  return result;
}


void SimData::computeFracFracTran(const std::size_t                 frac,
                                  const EmbeddedFracture          & efrac,
                                  const mesh::SurfaceMesh<double> & mesh,
                                  FlowData                        & frac_flow_data)
{
  CalcTranses calc;

  calc.NbNodes     = mesh.vertices.size();
  calc.NbPolyhedra = 0;
  calc.NbPolygons  = mesh.polygons.size();
  calc.NbZones     = calc.NbPolygons;  // 2 block + 1 frac
  calc.NbOptions   = 1;  // when 2 runs volume correction procedure
  calc.fracporo    = 1.0;
  calc.init();
  // -------------------- geometry -----------------------
  for (std::size_t i=0; i<mesh.vertices.size(); ++i)
  {
    calc.vCoordinatesX[i] = mesh.vertices[i][0];
    calc.vCoordinatesY[i] = mesh.vertices[i][1];
    calc.vCoordinatesZ[i] = mesh.vertices[i][2];
  }

  // polygons (2d elements)
  const std::size_t n_poly = mesh.polygons.size();
  int code_polygon = 0;
  for(std::size_t ipoly = 0; ipoly < n_poly; ++ipoly)
  {
    calc.vvFNodes[ipoly] = mesh.polygons[ipoly];
    calc.vCodePolygon[ipoly] = code_polygon;
    code_polygon++;
  }

  // --------------- Properties ------------------------
  for ( std::size_t ipoly = 0; ipoly < n_poly; ++ipoly )
  {
    calc.vZoneCode[ipoly] = ipoly;
    calc.vZVolumeFactor[ipoly] = efrac.aperture;
    calc.vZPorosity[ipoly] = 1.0;
    calc.vZPermCode[ipoly] = 1;

    const double f_perm = efrac.conductivity / efrac.aperture;
    calc.vZPermeability[ipoly*3 + 0] = f_perm;
    calc.vZPermeability[ipoly*3 + 1] = f_perm;
    calc.vZPermeability[ipoly*3 + 2] = f_perm;

    calc.vZConduction[ipoly*3 + 0] = 1;
    calc.vZConduction[ipoly*3 + 1] = 1;
    calc.vZConduction[ipoly*3 + 2] = 1;
    calc.vTimurConnectionFactor[ipoly] = 1.0;
  }

  calc.compute_flow_data();
  calc.extractData(frac_flow_data);
}


void SimData::computeEDFMTransmissibilities(const std::vector<angem::PolyGroup<double>> & splits,
                                            const int frac_ind)
{
  const auto & efrac = vEfrac[frac_ind];
  // compute transmissibilities between one embedded fracture and cells
  // run karimi class once per cell
  std::cout << "computing frac-matrix transes" << std::endl;
  std::size_t ecell = 0; // index of cell in efrac cell array
  for (const auto & split : splits)  // loop split edfm cells
  {
    // global cell index
    const std::size_t icell = efrac.cells[ecell];

    // init tran
    CalcTranses tran;
    tran.NbNodes     = split.vertices.size();
    tran.NbPolyhedra = 2;  // frac splits a cell into two polehedra
    tran.NbPolygons  = split.polygons.size();
    tran.NbFracs     = 1;
    tran.NbZones     = 3;  // 2 block + 1 frac
    tran.NbOptions   = 1;  // when 2 runs volume correction procedure
    tran.fracporo    = 1.0;
    tran.init();

    // --------------- fill geometry -------------------
    // vertices
    for ( std::size_t i = 0; i < split.vertices.size(); i++ )
    {
      tran.vCoordinatesX[i] = split.vertices[i][0];
      tran.vCoordinatesY[i] = split.vertices[i][1];
      tran.vCoordinatesZ[i] = split.vertices[i][2];
    }

    // polygons (2d elements)
    int code_polygon = 0;
    tran.vvFNodes.clear();
    tran.vCodePolygon.clear();
    std::vector<double> vConductivity, vAperture;

    for(std::size_t ipoly = 0; ipoly < split.polygons.size(); ipoly++)
    {
      tran.vvFNodes.push_back(split.polygons[ipoly]);

      if(split.markers[ipoly] == MARKER_FRAC)  // fracture polygons
      {
        tran.vCodePolygon.push_back( code_polygon );
        code_polygon++;
        vConductivity.push_back(efrac.conductivity);
        vAperture.push_back(efrac.aperture);
      }
      else  // cell polygons (faces)
      {
        tran.vCodePolygon.push_back( -1 );
      }
    }

    // polyhedra (3d elements)
    tran.vvVFaces.resize(2);
    tran.vCodePolyhedron.clear();

    // polyhedra resulting from split, on both sides of frac
    // only two polyhedra result from split of a cell by a fracture
    // int n_poly_above = 0, n_poly_below = 0;
    for(int ipoly = 0; ipoly < split.polygons.size(); ipoly++)
    {
      if (split.markers[ipoly] == MARKER_BELOW_FRAC)  // below frac
        tran.vvVFaces[0].push_back( ipoly );

      else if (split.markers[ipoly] == MARKER_ABOVE_FRAC)  // above frac
        tran.vvVFaces[1].push_back( ipoly );

      else if (split.markers[ipoly] == MARKER_FRAC)  // fracture itself included in each polyhedera
      {
        tran.vvVFaces[0].push_back( ipoly );
        tran.vvVFaces[1].push_back( ipoly );
      }
      else
      {
        std::cout << "unknown split market " << split.markers[ipoly] << std::endl;
        exit(0);
      }
    }

    tran.vCodePolyhedron.push_back( 1 + 0 );
    tran.vCodePolyhedron.push_back( 1 + 1 );

    // --------------- Properties ------------------------
    tran.vZPermeability.assign(tran.NbZones * 3, 0.0);
    tran.vZConduction.assign( (tran.NbPolyhedra + tran.NbFracs) * 3, 0.0);

    // 1 (one) edfm frac
    int efrac_zone = 0;
    tran.vZoneCode[efrac_zone] = efrac_zone;
    tran.vZVolumeFactor[efrac_zone] = efrac.aperture;
    tran.vZPorosity[efrac_zone] = 1.0;
    tran.vZPermCode[efrac_zone] = 1;

    const double f_perm = efrac.conductivity / efrac.aperture;
    tran.vZPermeability[0] = f_perm;
    tran.vZPermeability[1] = f_perm;
    tran.vZPermeability[2] = f_perm;

    tran.vZConduction[0] = 1;
    tran.vZConduction[1] = 1;
    tran.vZConduction[2] = 1;
    tran.vTimurConnectionFactor[0] = 1.0;

    // properties of polyhedra in cell
    const angem::Point<3,double> perm = get_permeability(icell);
    for ( std::size_t i = 0; i < tran.NbPolyhedra; i++ )
    {
      const std::size_t n = 1 + i;
      tran.vZoneCode[n] = tran.vCodePolyhedron[i];
      tran.vZPorosity[n] = get_property(icell, "PORO");
      assert(tran.vZPorosity[n] > 1e-16);
      tran.vZPermCode[n] = n;

      tran.vZPermeability[n*3+0] = perm[0];
      tran.vZPermeability[n*3+1] = perm[1];
      tran.vZPermeability[n*3+2] = perm[2];
      double thc = 0;
      try
      {
        thc = get_property(icell, "THCROCK");
      }
      catch (const std::out_of_range& e)
      {
        tran.vZConduction[n*3+0] = thc;
        tran.vZConduction[n*3+1] = thc;
        tran.vZConduction[n*3+2] = thc;
      }

      tran.vZVolumeFactor[n] = 1;
      tran.vTimurConnectionFactor[n] = 1.0;
    }

    tran.compute_flow_data();

    FlowData matrix_fracture_flow_data;
    tran.extractData(matrix_fracture_flow_data);

    // fill global flow data
    double f_m_tran = 0;
    { // compute frac-matrix trans as a sum of two frac-block trances weighted by volume
      const double t1 = matrix_fracture_flow_data.trans_ij[0];
      const double t2 = matrix_fracture_flow_data.trans_ij[1];
      const double v1 = matrix_fracture_flow_data.volumes[1];
      const double v2 = matrix_fracture_flow_data.volumes[2];
      f_m_tran = (t1*v1 + t2*v2) / (v1 + v2);
      // f_m_tran = (v1 + v2) / (v1/t1 + v2/t2);
      // f_m_tran = std::min(t1, t2);
    }

    flow_data.trans_ij.push_back(f_m_tran);
    flow_data.insert_connection(get_flow_element_index(frac_ind, ecell),
                                n_flow_dfm_faces + icell);

    ecell++;
  }  // end splits loop

  // compute transmissibilities between EDFM segments
  std::cout << "frac-frac whithin one frac approximations" << std::endl;
  FlowData frac_flow_data;
  computeFracFracTran(frac_ind, efrac, efrac.mesh, frac_flow_data);

  // fill global data
  for (std::size_t i=0; i<efrac.mesh.polygons.size(); ++i)
  {
    flow_data.volumes.push_back(frac_flow_data.volumes[i]);
    flow_data.poro.push_back(frac_flow_data.poro[i]);
    flow_data.depth.push_back(frac_flow_data.depth[i]);
  }

  for (const auto & conn : frac_flow_data.map_connection)
  {
    const std::size_t iconn = conn.second;
    const auto element_pair = frac_flow_data.invert_hash(conn.first);
    flow_data.insert_connection(get_flow_element_index(frac_ind, element_pair.first),
                                get_flow_element_index(frac_ind, element_pair.second));
    flow_data.trans_ij.push_back(frac_flow_data.trans_ij[iconn]);
  }

  // save custom cell data
  const std::size_t n_vars = rockPropNames.size();
  for (std::size_t i=0; i<efrac.cells.size(); ++i)
  {
    flow_data.custom_data.emplace_back();
    for (std::size_t j=0; j<n_vars; ++j)
      if (config.expression_type[j] == 0)
        flow_data.custom_data.back().push_back( vsCellRockProps[efrac.cells[i]].v_props[j] );
  }
}


std::size_t SimData::n_default_vars() const
{
  SimdataConfig dummy;
  return dummy.all_vars.size();
}


void SimData::
compute_frac_frac_intersection_transes(const std::vector<Point>                    & verts,
                                       const std::vector<std::vector<std::size_t>> & polys,
                                       const std::vector<int>                      & markers,
                                       FlowData                                    & flow_data) const
{
  CalcTranses tran;

  tran.NbNodes     = verts.size();
  tran.NbPolyhedra = 0;
  tran.NbPolygons  = polys.size();
  tran.NbZones     = tran.NbPolygons;  // 2 block + 1 frac
  tran.NbOptions   = 1;  // when 2 runs volume correction procedure
  tran.fracporo    = 1.0;
  tran.init();

  // -------------------- geometry -----------------------
  for (std::size_t i=0; i<verts.size(); ++i)
  {
    tran.vCoordinatesX[i] = verts[i][0];
    tran.vCoordinatesY[i] = verts[i][1];
    tran.vCoordinatesZ[i] = verts[i][2];
  }

  // polygons (2d elements)
  int code_polygon = 0;
  // tran.vNbFNodes.clear();
  tran.vvFNodes.clear();
  tran.vCodePolygon.clear();
  vector<double> vConductivity, vAperture;
  for(int ipoly = 0; ipoly < polys.size(); ipoly++)
  {
    tran.vvFNodes.push_back( polys[ipoly] );
    // tran.vNbFNodes.push_back( polys[ipoly].size() );

    tran.vCodePolygon.push_back( code_polygon );
    code_polygon++;
    vConductivity.push_back(vEfrac[markers[ipoly]].conductivity);
    vAperture.push_back(vEfrac[markers[ipoly]].aperture);
  }
  // --------------- Properties ------------------------
  tran.vZPermeability.assign(tran.NbZones * 3, 0.0);
  tran.vZConduction.assign(tran.NbPolygons * 3, 0.0);

  for ( int i = 0; i < polys.size(); i++ )
  {
    tran.vZoneCode[i] = i;
    tran.vZVolumeFactor[i] = vAperture[i];
    tran.vZPorosity[i] = 1.0;
    tran.vZPermCode[i] = 1;

    //@HACK default permeability for all fractures
    const double cond = vConductivity[i];
    const double w = vAperture[i];
    tran.vZPermeability[i*3 + 0] = cond / w;
    tran.vZPermeability[i*3 + 1] = cond / w;
    tran.vZPermeability[i*3 + 2] = cond / w;

    tran.vZConduction[i*3 + 0] = 1;
    tran.vZConduction[i*3 + 1] = 1;
    tran.vZConduction[i*3 + 2] = 1;
    tran.vTimurConnectionFactor[i] = 1.0;
  }

  tran.compute_flow_data();
  tran.extractData(flow_data);
}


void SimData::defineRockProperties()
{
  // print header
  std::cout << "function parsers setup" << std::endl;
  std::cout << "Variables:" << std::endl;
  for (std::size_t i=0; i<config.all_vars.size(); ++i)
  {
    std::cout << config.all_vars[i] << "\t";
    if ((i+1)%10 == 0)
      std::cout << std::endl;
  }
  std::cout << std::endl;

  // get number of variables in an empty simdataconfig - should be 3=x+y+z
  const std::size_t shift = n_default_vars();

  // resize rock properties
  vsCellRockProps.resize(grid.n_cells());

  const std::size_t n_variables = config.all_vars.size();
  std::vector<double> vars(n_variables);

  // save variables name for output
  rockPropNames.resize(n_variables - shift);
  for (std::size_t i=shift; i<config.all_vars.size(); ++i)
    rockPropNames[i - shift] = config.all_vars[i];

  // loop various domain configs:
  // they may have different number of variables and expressions
  for (const auto & conf: config.domains)
  {
    const std::size_t n_expressions = conf.expressions.size();
    std::vector<mu::Parser> parsers(n_expressions);

    // assign expressions and variables
    for (std::size_t i=0; i<n_expressions; ++i)
    {
      for (std::size_t j=0; j<n_variables; ++j)
      {
        try {
          parsers[i].DefineVar(config.all_vars[j], &vars[j]);
        }
        catch(mu::Parser::exception_type & e)
        {
          std::cout << _T("Initialization error:  ") << e.GetMsg() << endl;
          std::cout << "when setting variable '"
                    << config.all_vars[j]
                    << "'" << std::endl;
          exit(-1);
        }

      }
      try {
        parsers[i].SetExpr(conf.expressions[i]);
      }
      catch(mu::Parser::exception_type & e)
      {
        std::cout << _T("Initialization error:  ") << e.GetMsg() << endl;
        std::cout << "when setting expression '"
                  << conf.expressions[i]
                  << "'" << std::endl;
        exit(-1);
      }
    }

    // loop cells and evaluate expressions
    for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
    {
      if ( cell.marker() == conf.label ) // cells
      {
        std::fill(vars.begin(), vars.end(), 0);
        Point center = cell.center();
        vars[0] = center[0];  // x
        vars[1] = center[1];  // y
        vars[2] = center[2];  // z

        // Evaluate expression -> write into variable
        for (std::size_t i=0; i<n_expressions; ++i)
        {
          try {
          vars[conf.local_to_global_vars.at(i)] = parsers[i].Eval();
          }
          catch(mu::Parser::exception_type & e)
          {
            std::cout << _T("Evaluation error:  ") << e.GetMsg() << endl;
            std::cout << "when evaluating expression '"
                      << conf.expressions[i]
                      << "'" << std::endl;
            exit(-1);
          }
        }

        // copy vars to cell properties
        vsCellRockProps[cell.index()].v_props.resize(n_variables - shift);
        // start from 3 to skip x,y,z
        for (std::size_t j=shift; j<n_variables; ++j)
        {
          try
          {
            vsCellRockProps[cell.index()].v_props[j - shift] = vars[j];
          }
          catch (std::out_of_range & e)
          {
            vsCellRockProps[cell.index()].v_props[j - shift] = 0;
          }
        }
      }  // end match label
    }    // end cell loop
  }  // end domain loop

}


void SimData::splitInternalFaces()
{
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++ face)
  {
    if (is_fracture(face.marker()))
      grid.mark_for_split(face);
  }

  grid.split_faces();
//   int counter_;
//   cout << endl << "Create set of stick vertices (slow)" << endl;
//   vector<set<int> > vsetGlueVerticies;
//   counter_ = 0;
//   for ( int icell = 0; icell < nCells; icell++ )
//   {
//     counter_ += vsCellCustom[icell].nVertices;
//   }
//   vsetGlueVerticies.resize(counter_);
//   for(int i = 0; i < counter_; i++) vsetGlueVerticies[i].insert(i);

//   // TODO
//   // vector<bool> exlude_some_verticies(nNodes,false);
//   // for ( int iface = 0; iface < nFaces; iface++ )
//   // {
//   //     // vector<int>::iterator it_face_vrtx;
//   //     // for ( it_face_vrtx = vsFaceCustom[ iface ].vVertices.begin();  it_face_vrtx != vsFaceCustom[ iface ].vVertices.end(); ++it_face_vrtx)
//   //   for (const auto & vert : vsFaceCustom[ iface ].vVertices)
//   //     {
//   //       if(vsFaceCustom[ iface ].nMarker == -3333332 || vsFaceCustom[ iface ].nMarker == -3333331)
//   //         // exlude_some_verticies[*it_face_vrtx] = true;
//   //         exlude_some_verticies[vert] = true;
//   //     }
//   // }

//   vector<double> vVerticesPair; vVerticesPair.resize(2, 0);
//   for ( int iface = 0; iface < nFaces; iface++ )
//   {
//     int job_percent = int ( ( 100. * iface ) / ( nFaces ) );
//     cout << "\r    " << job_percent << "%";
//     /// non physical face
//     if ( vsFaceCustom[ iface ].nMarker == 0 && vsFaceCustom[ iface ].nNeighbors == 2 )
//     {
//       /// loop by all face vertices
//       // for ( it_face_vrtx = vsFaceCustom[ iface ].vVertices.begin();  it_face_vrtx != vsFaceCustom[ iface ].vVertices.end(); ++it_face_vrtx)
//       for (const auto & face_vertex : vsFaceCustom[ iface ].vVertices)
//       {
//         /// loop by neighbor cells
//         std::size_t n_polyhedron = 0;
//         // set<int>::iterator it_polyhedron;
//         // for ( it_polyhedron = vsetPolygonPolyhedron[iface].begin(); it_polyhedron != vsetPolygonPolyhedron[iface].end(); ++it_polyhedron, n_polyhedron++)
//         for (const auto & polyhedron : vsetPolygonPolyhedron[iface] )
//         {
//           /// loop by all cell vertices
//           // int ncell = *it_polyhedron;
//           std::size_t ncell = polyhedron;
//           std::size_t ivrtx = 0;

//           // vector<int>::iterator it_cell_vrtx;
//           // for ( it_cell_vrtx = vsCellCustom[ncell].vVertices.begin() ; it_cell_vrtx < vsCellCustom[ncell].vVertices.end(); ++it_cell_vrtx, ivrtx++ )
//           for (const auto & cell_vertex : vsCellCustom[ncell].vVertices)
//           {
//             // if ( *it_cell_vrtx == *it_face_vrtx )
//             if ( cell_vertex == face_vertex )
//               break;
//             ivrtx++;
//           }
//           vVerticesPair[n_polyhedron] = vsCellCustom[ncell].vVerticesNewnum[ivrtx];
//           n_polyhedron++;
//         }
//        vsetGlueVerticies[ vVerticesPair[0] ].insert ( vVerticesPair[1] );
//        vsetGlueVerticies[ vVerticesPair[1] ].insert ( vVerticesPair[0] );
//      }
//     }
//     /*
//     if ( vsFaceCustom[ iface ].nMarker > 0 && vsFaceCustom[ iface ].nNeighbors == 2 )
//     {
//       vector<int>::iterator it_face_vrtx;

//       /// loop by all face vertices
//       for ( it_face_vrtx = vsFaceCustom[ iface ].vVertices.begin();  it_face_vrtx != vsFaceCustom[ iface ].vVertices.end(); ++it_face_vrtx )
//       {
//         if ( exlude_some_verticies[*it_face_vrtx] == true)
//         {
//           /// loop by neighbor cells
//           int n_polyhedron = 0;
//           set<int>::iterator it_polyhedron;
//           for ( it_polyhedron = vsetPolygonPolyhedron[iface].begin(); it_polyhedron != vsetPolygonPolyhedron[iface].end(); ++it_polyhedron, n_polyhedron++ )
//           {
//             /// loop by all cell vertices
//             int ncell = *it_polyhedron;
//             int ivrtx = 0;

//             vector<int>::iterator it_cell_vrtx;
//             for ( it_cell_vrtx = vsCellCustom[ncell].vVertices.begin() ; it_cell_vrtx < vsCellCustom[ncell].vVertices.end(); ++it_cell_vrtx, ivrtx++ )
//             {
//               if ( *it_cell_vrtx == *it_face_vrtx )
//                 break;
//             }
//             vVerticesPair[n_polyhedron] = vsCellCustom[ncell].vVerticesNewnum[ivrtx];
//           }
//           vsetGlueVerticies[ vVerticesPair[0] ].insert ( vVerticesPair[1] );
//           vsetGlueVerticies[ vVerticesPair[1] ].insert ( vVerticesPair[0] );
//         }
//      }
//     } */
//   }

//   cout << endl << "Distinguish authenic vertices (might be very slow)" << endl;
//   int n_possible_verticies = vsetGlueVerticies.size();
//   vector<int> v_buf_storage;
//   int total_size_old = 0;
//   int total_size_new = 1;
//   int cycle = 0;
//   while ( total_size_old != total_size_new )
//   {
//     total_size_old = total_size_new;
//     total_size_new = 0;
//     cout << endl << "cycle   :" << cycle << "\t \t Hopefully < 10"; cycle++;

//     for ( int ivrtx = vsCellCustom[0].nVertices; ivrtx < n_possible_verticies; ivrtx++ )
//     {
//       int job_percent = int ( ( 100. * ivrtx ) / ( n_possible_verticies ) );
//       cout << "\r    " << job_percent << "%";
//       v_buf_storage.clear();
//       set<int>::iterator it_set;

//       for ( it_set = vsetGlueVerticies[ivrtx].begin(); it_set != vsetGlueVerticies[ivrtx].end(); ++it_set )
//       {
//         set<int>::iterator it_set_down;

//         for ( it_set_down = vsetGlueVerticies[ *it_set ].begin(); it_set_down != vsetGlueVerticies[ *it_set ].end(); ++it_set_down )
//         {
//           v_buf_storage.push_back ( *it_set_down );
//         }
//       }

//       vector<int>::iterator it_vec;

//       for ( it_vec = v_buf_storage.begin(); it_vec != v_buf_storage.end(); it_vec++ )
//       {
//         vsetGlueVerticies[ivrtx].insert ( *it_vec );
//       }
//       total_size_new += vsetGlueVerticies[ivrtx].size();
//     }
//   }

//   cout << endl << "Renumber vector of stick vertices" << endl;
//   vector<int> vRenumVerticies;
//   vRenumVerticies.resize(n_possible_verticies);

//   /// take first cell
//   for(int ivrtx = 0; ivrtx < vsCellCustom[0].nVertices; ivrtx++)
//   {
//     vRenumVerticies[ivrtx] = ivrtx;
//   }

//   counter_ = vsCellCustom[0].nVertices;
//   for(int ivrtx = vsCellCustom[0].nVertices; ivrtx < n_possible_verticies; ivrtx++)
//   {
//     if( *vsetGlueVerticies[ivrtx].begin() == ivrtx )
//     {
//       // this vertex is not stick
//       vRenumVerticies[ivrtx] = counter_;
//       counter_++;
//     }
//     else
//     {
//       // this vertex is stick
//       // set is sorted by c++ defaults, and we take minimum value
//       vRenumVerticies[ivrtx] = vRenumVerticies[ *vsetGlueVerticies[ivrtx].begin() ];
//     }
//   }

//   for(int icell = 0; icell < nCells; icell++)
//   {
//     for(int ivrtx = 0; ivrtx < vsCellCustom[icell].nVertices; ivrtx++)
//     {
//       vsCellCustom[icell].vVerticesNewnum[ ivrtx ] = vRenumVerticies[ vsCellCustom[icell].vVerticesNewnum[ ivrtx ] ];
//     }
//   }
//   int nTotalNodes = counter_;

//   cout << "\t check renumbering consistency " << endl;
//   vector<double> vStatus;
//   vStatus.resize(counter_, false);
//   for ( int icell = 0; icell < nCells; icell++ )
//   {
//     for ( int ivrtx = 0; ivrtx < vsCellCustom[icell].nVertices; ivrtx++ )
//       vStatus[ vsCellCustom[icell].vVerticesNewnum[ivrtx] ] = true;
//   }

//   cout << "\t change face numbering" << endl;
//   for(int iface = 0; iface < nFaces; iface++)
//   {
//     vsFaceCustom[iface].vVerticesNewnum.resize( vsFaceCustom[iface].nVertices, -1);
//     // we take always [0] - support cell
//     int icell = vsFaceCustom[iface].vNeighbors[0];
//     for(int ivrtx = 0; ivrtx < vsFaceCustom[iface].nVertices; ivrtx++)
//     {
//       for(int inode = 0; inode < vsCellCustom[icell].nVertices; inode++)
//       {
//         if( vsFaceCustom[iface].vVertices[ivrtx] == vsCellCustom[icell].vVertices[inode] )
//         {
//           vsFaceCustom[iface].vVerticesNewnum[ivrtx] = vsCellCustom[icell].vVerticesNewnum[inode];
//           break;
//         }
//       }
//     }
//   }

//   cout << "\t create atoms list" << endl;
//   set<int>::iterator itintset;
//   vector<int> vTempAtoms;
//   vTempAtoms.resize( nTotalNodes, -999 );
//   for (itintset = setIdenticalInternalMarker.begin(); itintset != setIdenticalInternalMarker.end(); ++itintset)
//   {
//     for(int iface = 0; iface < nFaces; iface++)
//     {
//       int ncell = vsFaceCustom[iface].vNeighbors[1];
//       if( vsFaceCustom[iface].nMarker == *itintset)
//       {
//         for(int ivrtx = 0; ivrtx < vsFaceCustom[iface].nVertices; ivrtx++)
//         {
//           for(int inode = 0; inode < vsCellCustom[ncell].nVertices; inode++)
//           {
//             if( vsFaceCustom[iface].vVertices[ivrtx] == vsCellCustom[ncell].vVertices[inode] )
//             {
//               vTempAtoms[ vsCellCustom[ncell].vVerticesNewnum[inode] ] = vsFaceCustom[iface].vVerticesNewnum[ivrtx];
//             }
//           }
//         }
//       }
//     }
//   }

//   cout << endl << "Change coordanates vector" << endl;
//   // vector<vector<double> > vvNewCoordinates;
//   vector<Point> vvNewCoordinates;
//   vvNewCoordinates.resize( nTotalNodes, vector<double>(3,0) );
//   for(std::size_t icell = 0; icell < nCells; icell++)
//   {
//     vector<std::size_t>::iterator it_old, it_new;
//     it_new = vsCellCustom[icell].vVerticesNewnum.begin();
//     for(it_old = vsCellCustom[icell].vVertices.begin(); it_old != vsCellCustom[icell].vVertices.end(); ++it_old, ++it_new)
//     {
//       vvNewCoordinates[ *it_new ] = vvVrtxCoords[ *it_old ];
//     }
//   }

//   vvVrtxCoords.resize(nTotalNodes, vector<double>(3,0));
//   for(int ivrtx = 0; ivrtx < nTotalNodes; ivrtx++)
//   {
//     vvVrtxCoords[ivrtx] = vvNewCoordinates[ivrtx];
//   }

//   cout << endl << "Unify previous and splitted data" << endl;
//   cout << "Verticies : " << nNodes << "\t \t After splitting : " << nTotalNodes << endl;
//   nNodes = nTotalNodes;
//   for(int icell = 0; icell < nCells; icell++)
//   {
//     vsCellCustom[icell].vVertices = vsCellCustom[icell].vVerticesNewnum;
//   }

//   for(int iface = 0; iface < nFaces; iface++)
//   {
//     vsFaceCustom[iface].vVertices = vsFaceCustom[iface].vVerticesNewnum;
//   }

//   vvAtoms.resize(nNodes, vector<int>(2,0) );
//   nAtoms = 0;
//   for(int iatom = 0; iatom < nNodes; iatom++)
//   {
//     if(vTempAtoms[iatom] >= 0)
//     {
//       vvAtoms[nAtoms][0] = iatom;
//       vvAtoms[nAtoms][1] = vTempAtoms[iatom];
//       nAtoms++;
//     }
//   }

//   return;
//   cout << endl << "Smart verticies renumbering" << endl;
//   vector<int> vNodesID;
//   vector<int> vIA;
//   vector<int> vJA;
//   vector<int> vRCM;

//   vIA.push_back(0);
//   for(int ic = 0; ic < nCells; ++ic)
//   {
//     int n = vsCellCustom[ic].vVertices.size();
//     for(int iv = 0; iv < n; ++iv)
//     {
//       vJA.push_back( vsCellCustom[ic].vVertices[iv] );
//     }
//     n += vIA[ic];
//     vIA.push_back(n);
//   }

//   vRCM.assign(nNodes,-1);
//   pRenum->convert(nCells, nNodes, vIA, vJA, vRCM);

//   for(int icell = 0; icell < nCells; icell++)
//   {
//     for(int iv = 0; iv < vsCellCustom[icell].vVertices.size(); iv++)
//       vsCellCustom[icell].vVertices[iv] = vRCM[vsCellCustom[icell].vVertices[iv]];
//   }

//   for(int iface = 0; iface < nFaces; iface++)
//   {
//     for(int iv = 0; iv < vsFaceCustom[iface].vVertices.size(); iv++)
//       vsFaceCustom[iface].vVertices[iv] = vRCM[vsFaceCustom[iface].vVertices[iv]];
//   }

//   for(int iatom = 0; iatom < nNodes; iatom++)
//   {
//     if(vTempAtoms[iatom] >= 0)
//     {
//       vvAtoms[nAtoms][0] = vRCM[vvAtoms[nAtoms][0]];
//       vvAtoms[nAtoms][1] = vRCM[vvAtoms[nAtoms][1]];
//     }
//   }

//   for(int ivrtx = 0; ivrtx < nTotalNodes; ivrtx++)
//     vvNewCoordinates[vRCM[ivrtx]] = vvVrtxCoords[ivrtx];

//   for(int ivrtx = 0; ivrtx < nTotalNodes; ivrtx++)
//     vvVrtxCoords[ivrtx] = vvNewCoordinates[ivrtx];
}


void SimData::handleConnections()
{
  // cells
  for (auto cell = grid.begin_cells(); cell!=grid.end_cells(); ++cell)
    for (const auto & conf : config.domains)
      if ( cell.marker() == conf.label and conf.coupled) // cells
        gm_cell_to_flow_cell[cell.index()].push_back(n_flow_dfm_faces + cell.index());

  // finally embedded fractures
  for (std::size_t ifrac=0; ifrac<vEfrac.size(); ++ifrac)
  {
    const auto & efrac = vEfrac[ifrac];
    for (std::size_t i=0; i<efrac.cells.size(); ++i)
    {
      const std::size_t icell = efrac.cells[i];
      for (const auto & conf: config.domains)
        if (grid.cell_markers[icell] == conf.label and conf.coupled)
          gm_cell_to_flow_cell[icell].push_back(get_flow_element_index(ifrac, i));
    }
  }
}


void SimData::definePhysicalFacets()
{
  gm_cell_to_flow_cell.resize(grid.n_cells());
  std::size_t n_facets = 0;
  nNeumannFaces = 0;
  nDirichletFaces = 0;
  nDirichletNodes = 0;
  int nfluid = 0;

  std::size_t iface = 0;
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face, ++iface)
  {
    const bool is_boundary = (face.neighbors().size() < 2);
    const int marker = face.marker();
    for (const auto & conf : config.bc_faces)
      if (marker == conf.label)
      {
        PhysicalFace facet;
        facet.nface = iface;
        facet.ntype = conf.type;
        facet.nmark = conf.label;
        facet.condition = conf.value;
        facet.nfluid = -1;
        boundary_faces.insert({face.index(), facet});
        boundary_face_markers.insert(marker);
        break;
      }

    if (!is_boundary and marker != 0)
    {
      PhysicalFace facet;
      facet.nface = iface;
      facet.ntype = 0;
      facet.nmark = marker;
      fracture_face_markers.insert(marker);
      bool coupled = false;
      const auto neighbors = face.neighbors();
      for (const auto & neighbor : neighbors)
        for (const auto & conf: config.domains)
          if (grid.cell_markers[neighbor] == conf.label and conf.coupled)
            coupled = true;

      if (coupled)
        facet.nfluid = nfluid;
      else
        facet.nfluid = -2;  // just a negative number (-2 + 1 < 0)

      bool found_label = false;
      for (std::size_t ifrac=0; ifrac<config.discrete_fractures.size(); ++ifrac)
        if( marker == config.discrete_fractures[ifrac].label)
        {
          facet.aperture = config.discrete_fractures[ifrac].aperture; //m
          facet.conductivity = config.discrete_fractures[ifrac].conductivity; //mD.m
          found_label = true;
        }
      if (!found_label)
      {
        std::cout << "No properties for DFM label "
                  << marker
                  << " found! Setting sealed fault."
                  << std::endl;
        facet.aperture = 1; //m
        facet.conductivity = 0; //mD.m
      }

      dfm_faces.insert({face.index(), facet});

      // save to gm_cell_to_flow_cell array
      for (const auto ineighbor : face.neighbors())
        for (const auto & conf : config.domains)
          if ( grid.cell_markers[ineighbor] == conf.label ) // cells
            if (conf.coupled)
              gm_cell_to_flow_cell[ineighbor].push_back(facet.nfluid);

      if (coupled)
        nfluid++;
    }
  }
  n_flow_dfm_faces = static_cast<std::size_t>(nfluid);
}


// void SimData::createSimpleWells()
// {
//   vector<double> center(3,0);
//   /// choise support cell for internal facets
//   // FAULT/FRACTURE PART
//   for ( int iwell = 0; iwell < nWells; iwell++ )
//   {
//     int icell = 0;
//     for ( int iface = 0; iface < nFaces; iface++ )
//     {
//       center[0] = vsFaceCustom[iface].center[0];
//       center[1] = vsFaceCustom[iface].center[1];
//       center[2] = vsFaceCustom[iface].center[2];

//       if ( vsFaceCustom[iface].nMarker > 0 )
//       {
//         if ( abs ( center[0] - vsWell[iwell].vWellCoordinate[0] ) < vsWell[iwell].radius_poisk &&
//              abs ( center[1] - vsWell[iwell].vWellCoordinate[1] ) < vsWell[iwell].radius_poisk )
//         {
//           if(center[0] < -250 || center[0] > 250)
// 	  {
// 	    vsWell[iwell].vRadiusPoisk.push_back ( abs ( center[0] - vsWell[iwell].vWellCoordinate[0] ) );
// 	    vsWell[iwell].vID.push_back ( icell );
// 	    vsWell[iwell].vWi.push_back ( 100.0 );
// 	  }
//         }
//         ++icell;
//       }
//     }
//   }
//   // MATRIX PART
//   /*
//   int n = 0;
//   for(int iface = 0; iface < nFaces; iface++)
//   {
//     if(vsFaceCustom[iface].nMarker > 0) n++;
//   }
//   for ( int iwell = 0; iwell < nWells; iwell++ )
//   {
//     int icell = 0;
//     for ( int ic = 0; ic < nCells; ic++ )
//     {
//       center[0] = vsCellCustom[ic].center[0];
//       center[1] = vsCellCustom[ic].center[1];
//       center[2] = vsCellCustom[ic].center[2];
//       if ( abs ( center[0] - vsWell[iwell].vWellCoordinate[0] ) < vsWell[iwell].radius_poisk
// 	 && (center[2] > vsWell[iwell].vWellCoordinate[2]) && (center[2] < vsWell[iwell].vWellCoordinate[3]) )
//       {
//         vsWell[iwell].vRadiusPoisk.push_back ( abs ( center[0] - vsWell[iwell].vWellCoordinate[0] ) )  ; //&& (center[2] > vsWell[iwell].vWellCoordinate[2]) && (center[2] < vsWell[iwell].vWellCoordinate[3])
//         vsWell[iwell].vID.push_back ( ic + n);
//         vsWell[iwell].vWi.push_back ( 10.0 ); //well index [default: 10.0]
//       }
//     }
//     if(vsWell[iwell].vWi.size() < 1)
//     {
//       cout << "Well definition is wrong" << endl;
//     }
//   }*/

//   for ( int iwell = 0; iwell < nWells; iwell++ )
//   {
//     vsWell[iwell].datum = 1e16;
//     for(int i = 0; i < vsWell[iwell].vID.size(); ++i)
//     {
//       int icell = vsWell[iwell].vID.size();
//       if( fabs(vsCellCustom[icell].center[2]) < vsWell[iwell].datum )
// 	      vsWell[iwell].datum = abs(vsCellCustom[icell].center[2]);
//     }
//   }

// #if 0
//   // limit number of perforations
//   for ( int iwell = 0; iwell < nWells; iwell++ )
//   {
//     vector<double> vrad_;
//     vrad_.assign(vsWell[iwell].vRadiusPoisk.size(),0);
//     vrad_ = vsWell[iwell].vRadiusPoisk;
//     sort(vrad_.begin(), vrad_.end());

//     int n_ = std::min(2, int(vrad_.size()));
//     for(int j = 0; j < n_; ++j)
//     {
//       int m_ = 0;
//       for(int k = j; k < vsWell[iwell].vRadiusPoisk.size(); ++k)
//       {
//        cout << vrad_[j] << "\t" << vsWell[iwell].vRadiusPoisk[k] << endl;
//        if(abs(vrad_[j] - vsWell[iwell].vRadiusPoisk[k]) < 1e-8 )
//        {
// 	  m_ = k;
//           break;
//        }
//       }
//       double rad_ = vsWell[iwell].vRadiusPoisk[j];
//       vsWell[iwell].vRadiusPoisk[j] = vsWell[iwell].vRadiusPoisk[m_];
//       vsWell[iwell].vRadiusPoisk[m_] = rad_;
//       int l_ = vsWell[iwell].vID[j];
//       vsWell[iwell].vID[j] = vsWell[iwell].vID[m_];
//       vsWell[iwell].vID[m_] = l_;
//       vsWell[iwell].vWi[j] = vsWell[iwell].vWi[m_];
//     }
//     vsWell[iwell].vRadiusPoisk.resize(n_);
//     vsWell[iwell].vID.resize(n_);
//     vsWell[iwell].vWi.resize(n_);

//     for(int k_ = 0; k_ < vsWell[iwell].vRadiusPoisk.size(); ++k_)
//     {
//       cout << "WELL(" << iwell << "), Perf(" << k_ << ") = "  <<  vsWell[iwell].vID[k_] << endl;
//     }
//   }
// #endif
// }


double SimData::get_property(const std::size_t cell,
                             const std::string & key) const
{
  // // query property by key
  const std::size_t ikey = find(key, rockPropNames);

  if (ikey == rockPropNames.size())
      throw std::out_of_range(key);

  if (ikey >= vsCellRockProps[cell].v_props.size())
    throw std::out_of_range("You most probably haven't specified props for this part of domain");

  if (cell >= vsCellRockProps.size())
    throw std::out_of_range(std::to_string(cell));

  return vsCellRockProps[cell].v_props[ikey];
}


angem::Point<3,double> SimData::get_permeability(const std::size_t cell) const
{
  try
  {
    const double permx = get_property(cell, "PERMX");
    const double permy = get_property(cell, "PERMY");
    const double permz = get_property(cell, "PERMZ");
    return angem::Point<3,double>(permx, permy, permz);
  }
  catch (const std::out_of_range& e)
  {
    try
    {
      const double perm = get_property(cell, "PERM");
      return angem::Point<3,double>(perm, perm, perm);
    }
    catch (const std::out_of_range& e)
    {
      return angem::Point<3,double>(config.default_permeability,
                                    config.default_permeability,
                                    config.default_permeability);
    }
  }
}


double SimData::get_volume_factor(const std::size_t cell) const
{
  try
  {
    const double vf = get_property(cell, "VFACTOR");
    if (vf < 1e-16)
      return config.default_volume_factor;
    return vf;
  }
  catch (const std::out_of_range& e)
  {
    return config.default_volume_factor;
  }
}


void SimData::computeTransBetweenDifferentEfracs()
{
  for (std::size_t i=0; i<vEfrac.size(); ++i)
    for (std::size_t j=i+1; j<vEfrac.size(); ++j)
    {
      const auto & ifrac = vEfrac[i];
      const auto & jfrac = vEfrac[j];

      angem::CollisionGJK<double> collision;
      const auto & iShape = config.fractures[i].body;
      const auto & jShape = config.fractures[j].body;

      const double tol = std::min(ifrac.mesh.minimum_edge_size() / 3,
                                  jfrac.mesh.minimum_edge_size() / 3);
      // fast check
      if (collision.check(*iShape, *jShape))
      {
        for (std::size_t ielement=0; ielement<ifrac.mesh.polygons.size(); ++ielement)
          for (std::size_t jelement=0; jelement<jfrac.mesh.polygons.size(); ++jelement)
          {
            const angem::Polygon<double> poly_i(ifrac.mesh.vertices,
                                                ifrac.mesh.polygons[ielement]);
            const angem::Polygon<double> poly_j(jfrac.mesh.vertices,
                                                jfrac.mesh.polygons[jelement]);
            std::vector<Point> section;
            if (angem::collision(poly_i, poly_j, section, tol))
            {
              std::cout << "found intersection between edfm fracs" << std::endl;
              std::cout << section << std::endl;
              if (!section.empty())
              {
                angem::PolyGroup<double> splits(1e-8);
                angem::split(poly_i, poly_j.plane, splits, i, i);
                angem::split(poly_j, poly_i.plane, splits, j, j);

                FlowData frac_frac_flow_data;
                compute_frac_frac_intersection_transes(splits.vertices.points,
                                                       splits.polygons,
                                                       splits.markers,
                                                       frac_frac_flow_data);
                std::cout << "frac-frac intersection tran data" << std::endl;
                double trans = 0;
                for (const auto & conn : frac_frac_flow_data.map_connection)
                {
                  const std::size_t iconn = conn.second;
                  const auto element_pair = frac_frac_flow_data.invert_hash(conn.first);
                  if (splits.markers[element_pair.first] != splits.markers[element_pair.second])
                    trans += frac_frac_flow_data.trans_ij[iconn];
                }

                flow_data.trans_ij.push_back(trans);
                flow_data.insert_connection(get_flow_element_index(i, ielement),
                                            get_flow_element_index(j, jelement));
              }

            }


          }
      }
    }
}


// void SimData::meshFractures()
// {
//   bool should_do_remeshing = false;
//   for (std::size_t f=0; f<vEfrac.size(); ++f)
//     if (config.fractures[f].n1 > 0)
//       should_do_remeshing = true;
//   if (!should_do_remeshing)
//     return;

//   std::vector<mesh::SurfaceMesh<double>> new_frac_meshes(vEfrac.size());
//   std::size_t old_shift = nDFMFracs + nCells;
//   std::size_t new_shift = nDFMFracs + nCells;

//   for (std::size_t f=0; f<vEfrac.size(); ++f)
//   {
//     auto & efrac = vEfrac[f];
//     const auto & frac_rect = *(config.fractures[f].body);
//     const angem::Basis<3,double> frac_basis = frac_rect.plane.get_basis();
//     Point t1 = frac_basis(0);
//     Point t2 = frac_basis(1);

//     const auto & points = frac_rect.get_points();
//     const double length = (points[1] - points[0]).norm();
//     const double width = (points[2] - points[1]).norm();

//     t1 *= length;
//     t2 *= width;

//     const std::size_t n1 = config.fractures[f].n1;
//     const std::size_t n2 = config.fractures[f].n2;

//     if ( n1 > 0 and n2 > 0)  // do remeshing
//     {
//       mesh::SurfaceMesh<double> new_frac_mesh =
//           mesh::make_surface_mesh(t1, t2, points[0], n1, n2);

//       const double tol = std::min(new_frac_mesh.minimum_edge_size() / 3,
//                                   efrac.mesh.minimum_edge_size() / 3);

//       for (std::size_t i=0; i<new_frac_mesh.polygons.size(); ++i)
//         for (std::size_t j=0; j<efrac.mesh.polygons.size(); ++j)
//         {
//           const angem::Polygon<double> poly_i(new_frac_mesh.vertices,
//                                               new_frac_mesh.polygons[i]);
//           const angem::Polygon<double> poly_j(efrac.mesh.vertices,
//                                               efrac.mesh.polygons[j]);
//           std::vector<Point> section;
//           if (angem::collision(poly_i, poly_j, section, tol))
//           {
//             if (section.size() == 2) // only touching sides
//               continue;

//             std::cout << "collision " << i << " " << j << std::endl;
//             const angem::Polygon<double> poly_section(section);

//             // const std::size_t old_element = get_flow_element_index(f, j);
//             const std::size_t old_element = old_shift + j;
//             const auto neighbors = flow_data.v_neighbors[old_element];
//             for (const auto & neighbor : neighbors)
//             {
//               if (neighbor < nDFMFracs + nCells and neighbor > nDFMFracs)
//               {
//                 std::size_t new_conn;
//                 if (new_flow_data.connection_exists(new_shift + i, neighbor))
//                   new_conn = new_flow_data.connection_index(new_shift + i, neighbor);
//                 else
//                   new_conn = new_flow_data.insert_connection(new_shift + i, neighbor);

//                 const std::size_t old_conn = flow_data.connection_index(old_element, neighbor);
//                 const double factor = poly_section.area() / poly_j.area();
//                 const double T_ij = flow_data.trans_ij[old_conn];
//                 if (new_conn >= new_flow_data.trans_ij.size())
//                   new_flow_data.trans_ij.push_back(T_ij * factor);
//                 else
//                   new_flow_data.trans_ij[new_conn] += T_ij * factor;
//               }
//             }
//           }
//           else
//           {
//             // std::cout << "no collision " << i << " " << j << std::endl;
//           }
//         }  // end frac element loop

//       FlowData frac_flow_data;
//       computeFracFracTran(f, efrac, new_frac_mesh, frac_flow_data);

//       for (std::size_t i=0; i<new_frac_mesh.polygons.size(); ++i)
//       {
//         new_flow_data.volumes.push_back(frac_flow_data.volumes[i]);
//         new_flow_data.poro.push_back(frac_flow_data.poro[i]);
//         new_flow_data.depth.push_back(frac_flow_data.depth[i]);
//       }

//       for (const auto & conn : frac_flow_data.map_connection)
//       {
//         const std::size_t iconn = conn.second;
//         const auto element_pair = frac_flow_data.invert_hash(conn.first);
//         const std::size_t i = new_shift + element_pair.first;
//         const std::size_t j = new_shift + element_pair.second;
//         new_flow_data.insert_connection(i, j);
//         new_flow_data.trans_ij.push_back(frac_flow_data.trans_ij[iconn]);
//       }

//       // save custom cell data
//       const std::size_t n_vars = rockPropNames.size();
//       std::size_t n_flow_vars = 0;
//       for (std::size_t j=0; j<n_vars; ++j)
//         if (config.expression_type[j] == 0)
//           n_flow_vars++;

//       for (std::size_t i=0; i<new_frac_mesh.polygons.size(); ++i)
//       {
//         std::vector<double> new_custom_data(n_flow_vars);
//         const std::size_t ielement = new_shift + i;
//         const auto & neighbors = new_flow_data.v_neighbors[ielement];
//         std::vector<std::size_t> rock_cell_neighbors;
//         for (const auto & neighbor : neighbors)
//           if (neighbor < nCells)
//             rock_cell_neighbors.push_back(neighbor);

//         for (const auto & neighbor : rock_cell_neighbors)
//         {
//           std::size_t counter = 0;
//           for (std::size_t j=0; j<n_vars; ++j)
//             if (config.expression_type[j] == 0)
//             {
//               new_custom_data[counter] += vsCellRockProps[neighbor].v_props[j];
//               counter++;
//             }
//         }
//         // divide by number of neighbors
//         for (double & value : new_custom_data)
//           value /= static_cast<double>(rock_cell_neighbors.size());

//         new_flow_data.custom_data.push_back(new_custom_data);
//       }

//       new_frac_meshes[f] = std::move(new_frac_mesh);
//     }
//     else  // copy old efrac data
//     {
//       std::pair<std::size_t,std::size_t> range =
//           {old_shift, old_shift + vEfrac[f].mesh.polygons.size()};

//       for (const auto & conn : flow_data.map_connection)
//       {
//         const auto elements = flow_data.invert_hash(conn.first);
//         if (elements.second >= range.first and elements.second < range.second)
//         {
//           const double Tij = flow_data.trans_ij[conn.second];
//           std::size_t ielement = elements.first;
//           std::size_t jelement = elements.second;
//           if (ielement >= nCells)
//             ielement = ielement - old_shift + new_shift;
//           if (jelement >= nCells)
//             jelement = jelement - old_shift + new_shift;
//           new_flow_data.insert_connection(ielement, jelement);
//           new_flow_data.trans_ij.push_back(Tij);
//         }
//       }

//       for (std::size_t i=0; i<efrac.mesh.polygons.size(); ++i)
//       {
//         new_flow_data.volumes.push_back(flow_data.volumes[old_shift + i]);
//         new_flow_data.poro.push_back(flow_data.poro[old_shift + i]);
//         new_flow_data.depth.push_back(flow_data.depth[old_shift + i]);
//         new_flow_data.custom_data.push_back(flow_data.custom_data[old_shift + i]);
//       }
//     }

//     old_shift += vEfrac[f].mesh.polygons.size();
//     if (new_frac_meshes[f].empty())
//       new_shift = old_shift;
//     else
//       new_shift += new_frac_meshes[f].polygons.size();
//   }  // end efrac loop

//   // const std::string vtk_file2 = "./ababa.vtk";
//   // IO::VTKWriter::write_vtk(new_frac_mesh.vertices.points,
//   //                          new_frac_mesh.polygons,
//   //                          vtk_file2);

//   std::cout << "saving frac meshes" << std::endl;

//   bool any_new = false;
//   for (std::size_t f=0; f<vEfrac.size(); ++f)
//   {
//     if (!new_frac_meshes[f].empty())
//     {
//       vEfrac[f].mesh = std::move(new_frac_meshes[f]);
//       any_new = true;
//     }
//   }

//   if (any_new)
//     flow_data = std::move(new_flow_data);
// }
