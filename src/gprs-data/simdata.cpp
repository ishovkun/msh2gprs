#include <simdata.hpp>

// library for analytical geometry
#include <angem/Point.hpp>
#include <angem/Rectangle.hpp>
#include <angem/PointSet.hpp>
#include <angem/CollisionGJK.hpp>
#include <angem/Collisions.hpp>
#include <angem/PolyGroup.hpp>
#include <angem/utils.hpp>
#include <mesh/utils.hpp> // to remesh embedded fractures
#include <mesh/Mesh.hpp> // 3D mesh format
// parser for user-defined expressions for reservoir data
#include <muparser/muParser.h>

#include <algorithm>
#include <exception>
#include <unordered_set>

using Point = angem::Point<3, double>;
const int MARKER_BELOW_FRAC = 0;
const int MARKER_ABOVE_FRAC = 1;
const int MARKER_FRAC = 2;

namespace gprs_data
{

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

    std::cout << "computing EDFM collisions for mechanics" << std::endl;
 redo_collision:
    // find cells intersected by the fracture
    frac.cells.clear();
    for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
    {
      const std::unique_ptr<angem::Polyhedron<double>> p_poly_cell = cell.polyhedron();
      const auto & poly_cell = *p_poly_cell;

      if (collision.check(*frac_conf.body, poly_cell))
      {
        frac.cells.push_back(cell.index());

        // check if some vertices are too close to the fracture
        // and if so move a fracture a little bit
        // for (const auto & ivertex : cell.vVertices)
        for (const auto & vertex : poly_cell.get_points())
        {
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


void SimData::handleEmbeddedFractures()
{
  cout << "Compute cell clipping and EDFM transmissibilities" << endl;
  computeCellClipping();

  // // cout << "Merge small edfm cells" << endl;
  // // pSimData->mergeSmallFracCells();

  // // std::cout << "mesh fractures" << std::endl;
  // // pSimData->meshFractures();

  std::cout << "Compute transmissibilities between edfm fracs" << std::endl;
  computeTransBetweenDifferentEfracs();
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
      angem::Polygon<double> poly_section(set_points);
      vvSection[i] = poly_section.get_points();

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

      frac_msh.insert(poly_section);

      // add fracture polygon to splits to compute transes
      splits[i].add(angem::Polygon<double>(set_points), MARKER_FRAC);

      // write points into a global set so we have an ordered set
      // of vertices and we can retreive indices
      for (const Point & p : set_points)
        setVert.insert(p);

    }  // end sda cells loop

    vEfrac[ifrac].mesh = std::move(frac_msh);

    computeEDFMTransmissibilities(splits, ifrac);

  }  // end efrac loop

}


void SimData::mergeSmallFracCells()
{
  for (std::size_t ifrac=0; ifrac<config.fractures.size(); ++ifrac)
  {
    auto & msh = vEfrac[ifrac].mesh;

    double max_area = 0;
    for (const auto & element : msh.polygons)
    {
      angem::Polygon<double> poly(msh.vertices.points, element);
      const double area = poly.area();
      max_area = std::max(area, max_area);
    }

    // merge tiny cells
    std::size_t ielement = 0;
    std::size_t n_frac_elements = msh.polygons.size();

    // loop is with variable upper limit since elements can be
    // merged and deleted
    while (ielement < n_frac_elements)
    {
      angem::Polygon<double> poly(msh.vertices.points,
                                  msh.polygons[ielement]);
      const double area_factor = poly.area() / max_area;

      if (area_factor < config.frac_cell_elinination_factor)
      {
        const std::size_t global_ielement = efrac_flow_index(ifrac, ielement);
        const std::size_t new_element = msh.merge_element(ielement);
        // update flow data
        flow_data.merge_elements(efrac_flow_index(ifrac, new_element),
                                 global_ielement);

        n_frac_elements = msh.polygons.size();
        if (ielement >= n_frac_elements)
          break;
        continue;
      }
      ielement++;
    }

  }
}


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

  // FILL DATA
  // nodes
  for ( std::size_t i = 0; i < grid.n_vertices(); i++ )
  {
    Point vertex = grid.vertices[i];
    calc.X[i] = vertex.x();
    calc.Y[i] = vertex.y();
    calc.Z[i] = vertex.z();
  }

  // polygons (2d elements)
  calc.vvFNodes.resize(grid.n_faces());

  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    const std::size_t ipoly = face.index();
    calc.vvFNodes[ipoly] = face.vertex_indices();

    if (is_fracture(face.marker()))
      calc.vCodePolygon[ipoly] = dfm_faces.find(face.index())->second.nfluid;
    else  // non-frac faces
      calc.vCodePolygon[ipoly] = -1;
  }

  // polyhedra (3d elements)
  calc.vvVFaces.resize(grid.n_cells());
  calc.vCodePolyhedron.resize(grid.n_cells());
  for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
  {
    // std::cout << "cell.volume() = " << cell.volume() << std::endl;
    const auto faces = cell.faces();
    const std::size_t icell = cell.index();
    calc.vvVFaces[icell].reserve(faces.size());
    for (const mesh::face_iterator & face : faces)
    {
      const std::size_t ipolygon = face.index();
      calc.vvVFaces[icell].push_back(ipolygon);
    }
    calc.vCodePolyhedron[icell] = n_flow_dfm_faces + icell;
  }

  // Properties
  // DFM fractures
  for (const auto & it : dfm_faces)
  {
    if (it.second.coupled)  // if active
    {
      const int i = static_cast<int>(it.second.nfluid);
      calc.vZoneCode[i] = i;
      calc.vZVolumeFactor[i] = it.second.aperture;
      calc.vZPorosity[i] = 1.0;
      calc.vZPermCode[i] = 1;

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

  // avoid catching exception each time
  bool thermal_conductivity_available = true;
  try
  {
    get_property(0, "THCROCK");
  }
  catch (const std::out_of_range& e)
  {
    thermal_conductivity_available = false;
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
    if (thermal_conductivity_available)
      thc = get_property(i, "THCROCK");

    calc.vZConduction[n*3+0] = thc;
    calc.vZConduction[n*3+1] = thc;
    calc.vZConduction[n*3+2] = thc;

    calc.vZVolumeFactor[n] = get_volume_factor(i);
    calc.vTimurConnectionFactor[n] = 1.0;
  }

  std::cout << "Compute karimi" << std::endl;
  FlowData matrix_flow_data;
  calc.compute_flow_data();
  std::cout << "end compute karimi" << std::endl;
  calc.extractData(matrix_flow_data);

  // copy to global
  const std::size_t n_volumes = matrix_flow_data.volumes.size();
  std::cout << "n_volumes = " << n_volumes << std::endl;
  std::cout << "grid.n_cells() = " << grid.n_cells() << std::endl;
  flow_data.volumes.reserve(n_volumes);
  flow_data.poro.reserve(n_volumes);
  flow_data.depth.reserve(n_volumes);
  for (std::size_t i=0; i<n_volumes; ++i)
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
      if(dfm_faces.find(face.index())->second.coupled)
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


std::size_t SimData::efrac_flow_index(const std::size_t ifrac,
                                      const std::size_t ielement) const
{
  std::size_t result = n_flow_dfm_faces + grid.n_cells();
  for (std::size_t i=0; i<ifrac; ++i)
    result += vEfrac[i].mesh.polygons.size();

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
    calc.X[i] = mesh.vertices[i][0];
    calc.Y[i] = mesh.vertices[i][1];
    calc.Z[i] = mesh.vertices[i][2];
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
      tran.X[i] = split.vertices[i][0];
      tran.Y[i] = split.vertices[i][1];
      tran.Z[i] = split.vertices[i][2];
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
        exit(-1);
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
      f_m_tran = (t1*v1 + t2*v2) / (v1 + v2);  // arithmetic
      // f_m_tran = (v1 + v2) / (v1/t1 + v2/t2);  // harmonic
      // f_m_tran = std::min(t1, t2);
    }

    flow_data.trans_ij.push_back(f_m_tran);
    flow_data.insert_connection(efrac_flow_index(frac_ind, ecell),
                                res_cell_flow_index(icell));

    if (config.edfm_method == EDFMMethod::projection)
      apply_projection_edfm(frac_ind, ecell, icell);

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
    flow_data.insert_connection(efrac_flow_index(frac_ind, element_pair.first),
                                efrac_flow_index(frac_ind, element_pair.second));
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


void SimData::apply_projection_edfm(const std::size_t ifrac,     // embedded frac index ndex of embedded fracture
                                    const std::size_t ielement,  // frac element index  / fracture element index
                                    const std::size_t icell)     // reservoir cell index
{
  const mesh::cell_iterator cell = grid.create_cell_iterator(icell);
  const std::vector<Point> frac_element_vertices =
      vEfrac[ifrac].mesh.create_poly_iterator(ielement).vertices();
  const Point frac_normal = angem::Polygon(frac_element_vertices).plane.normal();

  for (const auto & face : cell.faces())
  {
    // don't connect to cells that are perpendicular to the fracture
    if (fabs(frac_normal.dot(face.normal())) < 1e-10)
      continue;

    // compute projection
    const auto face_poly = face.polygon();
    const std::vector<Point> projected_frac_vertices =
        face_poly.plane.project_points(frac_element_vertices);
    const double projection_area = angem::Polygon(projected_frac_vertices).area();

    // modify neighbor map
    const auto neighbor = cell.neighbor_by_face(face);

    const double k_cell_n = get_permeability(cell.index()) * face.normal();
    const double k_neighbor_n = get_permeability(neighbor.index()) * face.normal();
    const double volume_cell = cell.volume();
    const double volume_neighbor = neighbor.volume();
    const double k_face = (volume_cell + volume_neighbor) /
                          (volume_cell/k_cell_n + volume_neighbor/k_neighbor_n);
    // new face trans
    const double T_face_mm = (face_poly.area() - projection_area) /
                             (neighbor.center() - cell.center()).norm() *
                             CalcTranses::transmissibility_conversion_factor;
    // std::cout << "T_face_mm = " << T_face_mm << std::endl;

    // old matrix-matrix transmissibility
    std::size_t con = flow_data.connection_index(res_cell_flow_index(icell),
                                                 res_cell_flow_index(neighbor.index()));
    const double T_face_mm_old = flow_data.trans_ij[con];

    if (T_face_mm / T_face_mm_old < 1e-4)
      flow_data.remove_connection(con);
    else
      flow_data.trans_ij[con] = T_face_mm;

    // compute projection connection

    // double T_mm = flow_data.trans_ij[con];
    // double T_if =
    //     face_poly.area() / (neighbor.center() - cell.center()).norm() * 10;
    // std::cout << "trans = " << trans*0.0085267146719160104986876640419948 << std::endl;
    // std::size_t con = flow_data.connection_index(res_cell_flow_index(icell),
    //                                              res_cell_flow_index(neighbor.index()));
    // const double t1 = flow_data.trans_ij[con];
    // std::cout << "t1 = " << t1 << std::endl;
  }
  abort();
}


void SimData::
compute_frac_frac_intersection_transes(const std::vector<angem::Point<3,double>>   & verts,
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
    tran.X[i] = verts[i][0];
    tran.Y[i] = verts[i][1];
    tran.Z[i] = verts[i][2];
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

  dfm_master_grid = grid.split_faces();
}


void SimData::handleConnections()
{
  std::cout << "handle connections" << std::endl;
  // gm_cell_to_flow_cell.resize(flow_data.volumes.size());
  gm_cell_to_flow_cell.resize(grid.n_cells(), std::vector<std::size_t>());

  // cells
  for (auto cell = grid.begin_cells(); cell!=grid.end_cells(); ++cell)
    for (const auto & conf : config.domains)
    {
      if (cell.marker() == conf.label) // cells
        if (conf.coupled)
        {
          gm_cell_to_flow_cell[cell.index()].push_back(n_flow_dfm_faces + cell.index());
          break;
        }
    }

  // finally embedded fractures
  for (std::size_t ifrac=0; ifrac<vEfrac.size(); ++ifrac)
  {
    const auto & efrac = vEfrac[ifrac];
    for (std::size_t i=0; i<efrac.cells.size(); ++i)
    {
      const std::size_t icell = efrac.cells[i];
      for (const auto & conf: config.domains)
        if (grid.cell_markers[icell] == conf.label and conf.coupled)
          gm_cell_to_flow_cell[icell].push_back(efrac_flow_index(ifrac, i));
    }
  }

  // update dfm face indices in map
  std::unordered_set<std::size_t> face_touched;
  std::size_t counter = 0;
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++ face)
    if (is_fracture(face.marker()))
    {
      counter++;
      auto flow_face_it = dfm_faces.find(face.master_index());
      auto & facet = flow_face_it->second;
      if (face_touched.insert(face.master_index()).second)
      {
        facet.ifracture = find(face.marker(), fracture_face_markers);
        facet.nface = face.index();
      }
    }

  // std::cout << "n new dfm face = " << counter << std::endl;
}


void SimData::definePhysicalFacets()
{
  std::size_t n_facets = 0;
  int nfluid = 0;
  n_neumann_faces = 0;
  n_dirichlet_faces = 0;

  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    const bool is_boundary = (face.neighbors().size() < 2);
    const int marker = face.marker();
    for (const auto & conf : config.bc_faces)
      if (marker == conf.label)
      {
        PhysicalFace facet;
        facet.nface = face.index();
        facet.ntype = conf.type;
        facet.nmark = conf.label;
        facet.condition = conf.value;
        facet.coupled = false;
        boundary_faces.insert({face.index(), facet});
        boundary_face_markers.insert(marker);
        if (conf.type == 1)
          n_dirichlet_faces++;
        else if (conf.type == 2)
          n_neumann_faces++;

        break;
      }

    if (!is_boundary and marker != 0)
    {
      PhysicalFace facet;
      facet.nface = face.index();
      facet.ntype = 0;
      facet.nmark = marker;
      facet.neighbor_cells = face.neighbors();
      fracture_face_markers.insert(marker);
      bool coupled = false;
      const auto neighbors = face.neighbors();
      for (const auto & neighbor : neighbors)
        for (const auto & conf: config.domains)
          if (grid.cell_markers[neighbor] == conf.label and conf.coupled)
            coupled = true;

      if (coupled)
      {
        facet.nfluid = nfluid;
        facet.coupled = true;
      }
      else
      {
        facet.nfluid = -2;  // just a negative number (-2 + 1 < 0)
        facet.coupled = false;
      }

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

      if (coupled)
        nfluid++;
    }  //  end DFM fracture face case
  }
  n_flow_dfm_faces = static_cast<std::size_t>(nfluid);

  std::cout << "Number of Neumann faces = " << n_neumann_faces << std::endl;
  std::cout << "Number of Dirichlet faces = " << n_dirichlet_faces << std::endl;
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
                flow_data.insert_connection(efrac_flow_index(i, ielement),
                                            efrac_flow_index(j, jelement));
              }

            }


          }
      }
    }
}


void SimData::meshFractures()
{
  bool should_do_remeshing = false;
  for (std::size_t f=0; f<vEfrac.size(); ++f)
    if (config.fractures[f].n1 > 0)
      should_do_remeshing = true;
  if (!should_do_remeshing)
    return;

  std::vector<mesh::SurfaceMesh<double>> new_frac_meshes(vEfrac.size());
  std::size_t old_shift = n_flow_dfm_faces + grid.n_cells();
  std::size_t new_shift = n_flow_dfm_faces + grid.n_cells();

  for (std::size_t f=0; f<vEfrac.size(); ++f)
  {
    auto & efrac = vEfrac[f];
    const auto & frac_rect = *(config.fractures[f].body);
    const angem::Basis<3,double> frac_basis = frac_rect.plane.get_basis();
    Point t1 = frac_basis(0);
    Point t2 = frac_basis(1);

    const auto & points = frac_rect.get_points();
    const double length = (points[1] - points[0]).norm();
    const double width = (points[2] - points[1]).norm();

    t1 *= length;
    t2 *= width;

    const std::size_t n1 = config.fractures[f].n1;
    const std::size_t n2 = config.fractures[f].n2;

    if ( n1 > 0 and n2 > 0)  // do remeshing
    {
      mesh::SurfaceMesh<double> new_frac_mesh =
          mesh::make_surface_mesh(t1, t2, points[0], n1, n2);

      const double tol = std::min(new_frac_mesh.minimum_edge_size() / 3,
                                  efrac.mesh.minimum_edge_size() / 3);

      for (std::size_t i=0; i<new_frac_mesh.polygons.size(); ++i)
        for (std::size_t j=0; j<efrac.mesh.polygons.size(); ++j)
        {
          const angem::Polygon<double> poly_i(new_frac_mesh.vertices,
                                              new_frac_mesh.polygons[i]);
          const angem::Polygon<double> poly_j(efrac.mesh.vertices,
                                              efrac.mesh.polygons[j]);
          std::vector<Point> section;
          if (angem::collision(poly_i, poly_j, section, tol))
          {
            if (section.size() == 2) // only touching sides
              continue;

            std::cout << "collision " << i << " " << j << std::endl;
            const angem::Polygon<double> poly_section(section);

            const std::size_t old_element = old_shift + j;
            const auto neighbors = flow_data.v_neighbors[old_element];
            for (const auto & neighbor : neighbors)
            {
              if (neighbor < n_flow_dfm_faces + grid.n_cells() and neighbor > n_flow_dfm_faces)
              {
                std::size_t new_conn;
                if (new_flow_data.connection_exists(new_shift + i, neighbor))
                  new_conn = new_flow_data.connection_index(new_shift + i, neighbor);
                else
                  new_conn = new_flow_data.insert_connection(new_shift + i, neighbor);

                const std::size_t old_conn = flow_data.connection_index(old_element, neighbor);
                const double factor = poly_section.area() / poly_j.area();
                const double T_ij = flow_data.trans_ij[old_conn];
                if (new_conn >= new_flow_data.trans_ij.size())
                  new_flow_data.trans_ij.push_back(T_ij * factor);
                else
                  new_flow_data.trans_ij[new_conn] += T_ij * factor;
              }
            }
          }
          else
          {
            // std::cout << "no collision " << i << " " << j << std::endl;
          }
        }  // end frac element loop

      FlowData frac_flow_data;
      computeFracFracTran(f, efrac, new_frac_mesh, frac_flow_data);

      for (std::size_t i=0; i<new_frac_mesh.polygons.size(); ++i)
      {
        new_flow_data.volumes.push_back(frac_flow_data.volumes[i]);
        new_flow_data.poro.push_back(frac_flow_data.poro[i]);
        new_flow_data.depth.push_back(frac_flow_data.depth[i]);
      }

      for (const auto & conn : frac_flow_data.map_connection)
      {
        const std::size_t iconn = conn.second;
        const auto element_pair = frac_flow_data.invert_hash(conn.first);
        const std::size_t i = new_shift + element_pair.first;
        const std::size_t j = new_shift + element_pair.second;
        new_flow_data.insert_connection(i, j);
        new_flow_data.trans_ij.push_back(frac_flow_data.trans_ij[iconn]);
      }

      // save custom cell data
      const std::size_t n_vars = rockPropNames.size();
      std::size_t n_flow_vars = 0;
      for (std::size_t j=0; j<n_vars; ++j)
        if (config.expression_type[j] == 0)
          n_flow_vars++;

      for (std::size_t i=0; i<new_frac_mesh.polygons.size(); ++i)
      {
        std::vector<double> new_custom_data(n_flow_vars);
        const std::size_t ielement = new_shift + i;
        const auto & neighbors = new_flow_data.v_neighbors[ielement];
        std::vector<std::size_t> rock_cell_neighbors;
        for (const auto & neighbor : neighbors)
          if (neighbor < grid.n_cells())
            rock_cell_neighbors.push_back(neighbor);

        for (const auto & neighbor : rock_cell_neighbors)
        {
          std::size_t counter = 0;
          for (std::size_t j=0; j<n_vars; ++j)
            if (config.expression_type[j] == 0)
            {
              new_custom_data[counter] += vsCellRockProps[neighbor].v_props[j];
              counter++;
            }
        }
        // divide by number of neighbors
        for (double & value : new_custom_data)
          value /= static_cast<double>(rock_cell_neighbors.size());

        new_flow_data.custom_data.push_back(new_custom_data);
      }

      new_frac_meshes[f] = std::move(new_frac_mesh);
    }
    else  // copy old efrac data
    {
      std::pair<std::size_t,std::size_t> range =
          {old_shift, old_shift + vEfrac[f].mesh.polygons.size()};

      for (const auto & conn : flow_data.map_connection)
      {
        const auto elements = flow_data.invert_hash(conn.first);
        if (elements.second >= range.first and elements.second < range.second)
        {
          const double Tij = flow_data.trans_ij[conn.second];
          std::size_t ielement = elements.first;
          std::size_t jelement = elements.second;
          if (ielement >= grid.n_cells())
            ielement = ielement - old_shift + new_shift;
          if (jelement >= grid.n_cells())
            jelement = jelement - old_shift + new_shift;
          new_flow_data.insert_connection(ielement, jelement);
          new_flow_data.trans_ij.push_back(Tij);
        }
      }

      for (std::size_t i=0; i<efrac.mesh.polygons.size(); ++i)
      {
        new_flow_data.volumes.push_back(flow_data.volumes[old_shift + i]);
        new_flow_data.poro.push_back(flow_data.poro[old_shift + i]);
        new_flow_data.depth.push_back(flow_data.depth[old_shift + i]);
        new_flow_data.custom_data.push_back(flow_data.custom_data[old_shift + i]);
      }
    }

    old_shift += vEfrac[f].mesh.polygons.size();
    if (new_frac_meshes[f].empty())
      new_shift = old_shift;
    else
      new_shift += new_frac_meshes[f].polygons.size();
  }  // end efrac loop

//   // const std::string vtk_file2 = "./ababa.vtk";
//   // IO::VTKWriter::write_vtk(new_frac_mesh.vertices.points,
//   //                          new_frac_mesh.polygons,
//   //                          vtk_file2);

  std::cout << "saving frac meshes" << std::endl;

  bool any_new = false;
  for (std::size_t f=0; f<vEfrac.size(); ++f)
  {
    if (!new_frac_meshes[f].empty())
    {
      vEfrac[f].mesh = std::move(new_frac_meshes[f]);
      any_new = true;
    }
  }

  if (any_new)
    flow_data = std::move(new_flow_data);

  std::cout << "need to restore. aborting" << std::endl;
  abort();
}


bool SimData::is_embedded_fracture(const std::size_t flow_element_index) const
{
  if (flow_element_index > grid.n_cells() + n_flow_dfm_faces)
    return true;
  return false;
}


bool SimData::is_discrete_fracture(const std::size_t flow_element_index) const
{
  if (flow_element_index < n_flow_dfm_faces)
    return true;
  return false;
}


bool SimData::is_reservoir_element(const std::size_t flow_element_index) const
{
  if (!is_embedded_fracture(flow_element_index) and !is_discrete_fracture(flow_element_index))
    return true;
  return false;
}



void SimData::setupWells()
{
  for (const auto & conf : config.wells)
  {
    Well well(conf);
    std::cout << "setting up " << well.name << std::endl;
    if (well.simple())
      setupSimpleWell(well);
    else
      setupComplexWell(well);

    std::cout << "computing WI" << std::endl;
    computeWellIndex(well);
    wells.push_back(std::move(well));
  }
}


void SimData::setupSimpleWell(Well & well)
{
  std::cout << "simple well " << well.name << std::endl;
  const Point direction = {0, 0, -1};
  // well assigned with a single coordinate
  for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
  {
    const std::unique_ptr<angem::Polyhedron<double>> p_poly_cell = cell.polyhedron();
    if (p_poly_cell->point_inside(well.coordinate))
    {
      // define sufficiantly-large artificial segment to compute intersection
      // with cell
      const double h = (p_poly_cell->center() -
                        p_poly_cell->get_points()[0]).norm();
      const Point p1 = well.coordinate - direction*h*10;
      const Point p2 = well.coordinate + direction*h*10;
      std::vector<Point> section_data;
      if (angem::collision(p1, p2, *p_poly_cell, section_data, 1e-6))
      {
        well.connected_volumes.push_back(n_flow_dfm_faces + cell.index());
        well.segment_length.push_back(section_data[0].distance(section_data[1]));
        well.directions.push_back(direction);
        break;
      }
    }
  }
}


void SimData::setupComplexWell(Well & well)
{
  // setup well with segments
  std::cout << "complex well " << well.name << std::endl;
  for (std::size_t isegment=0; isegment<well.segments.size(); ++isegment)
  {
    std::cout << "isegment = " << isegment << ":\t";
    auto segment = well.segments[isegment];
    std::cout << segment.first << " | " << segment.second << std::endl;
    std::vector<Point> section_data;
    for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
    {
      const std::unique_ptr<angem::Polyhedron<double>> p_poly_cell = cell.polyhedron();
      if (angem::collision(segment.first, segment.second,
                           *p_poly_cell, section_data, 1e-6))
      {
        if (section_data.size() != 2)
        {
          std::cout << "something is wrong with segment" << std::endl;
          section_data.clear();
          continue;
        }
        well.connected_volumes.push_back(n_flow_dfm_faces + cell.index());

        well.segment_length.push_back(section_data[0].distance(section_data[1]));
        well.directions.push_back(segment.second - segment.first);
        section_data.clear();
      }
    }
  }
}


angem::Point<3,double> SimData::get_dx_dy_dz(const std::size_t icell) const
{
  const auto cell_poly = grid.create_cell_iterator(icell).polyhedron();
  angem::Point<3,double> dir, neg_dir, result;
  dir = {1, 0, 0}; neg_dir = - dir;
  result[0] = (cell_poly->support(dir) - cell_poly->support(neg_dir)).norm();
  dir = {0, 1, 0}; neg_dir = - dir;
  result[1] = (cell_poly->support(dir) - cell_poly->support(neg_dir)).norm();
  dir = {0, 0, 1}; neg_dir = - dir;
  result[2] = (cell_poly->support(dir) - cell_poly->support(neg_dir)).norm();
  return result;
}


double compute_productivity(const double k1, const double k2,
                            const double dx1, const double dx2,
                            const double length, const double radius,
                            const double skin = 0)
{
  // pieceman radius
  const double r =
      0.28*std::sqrt(std::sqrt(k2/k1)*dx1*dx1 + std::sqrt(k1/k2)*dx2*dx2) /
      (std::pow(k2/k1, 0.25) + std::pow(k1/k2, 0.25));
  const double j_ind = 2*M_PI*std::sqrt(k1*k2)*length/(std::log(r/radius) + skin);
  // std::cout << "pieceman, rwell " << r << "\t" << radius << std::endl << std::flush;
  // std::cout << "length " << length << std::endl << std::flush;
  // std::cout << "other "<< 2*M_PI*std::sqrt(k1*k2)*length << std::endl;
  assert(j_ind >= 0);
  return j_ind;


}


void SimData::computeWellIndex(Well & well)
{
  well.indices.resize(well.connected_volumes.size());
  for (std::size_t i=0; i<well.connected_volumes.size(); ++i)
  {
    const std::size_t icell = well.connected_volumes[i] - n_flow_dfm_faces;
    const angem::Point<3,double> perm = get_permeability(icell);
    angem::Point<3,double> dx_dy_dz = get_dx_dy_dz(icell);
    // std::cout << "dx_dy_dz = " << dx_dy_dz << std::endl;
    // std::cout << "perm = " << perm << std::endl;
    std::cout << "well.segment_length[i] = " << well.segment_length[i] << std::endl;
    angem::Point<3,double> productivity;
    productivity[0] =
        compute_productivity(perm[1], perm[2], dx_dy_dz[1], dx_dy_dz[2],
                             well.segment_length[i]*fabs(well.directions[i][0]),
                             well.radius);
    productivity[1] =
        compute_productivity(perm[0], perm[2], dx_dy_dz[0], dx_dy_dz[2],
                             well.segment_length[i]*fabs(well.directions[i][1]),
                             well.radius);
    productivity[2] =
        compute_productivity(perm[0], perm[1], dx_dy_dz[0], dx_dy_dz[1],
                             well.segment_length[i]*fabs(well.directions[i][2]),
                             well.radius);
    well.indices[i] = productivity.norm();
    std::cout << "productivity.norm() " << productivity.norm() << std::endl;
  }
}

}  // end namespace
