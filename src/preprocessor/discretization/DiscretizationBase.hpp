#pragma once

#include "mesh/Mesh.hpp"
#include "PhysicalFace.hpp"
#include "angem/Tensor2.hpp"

namespace discretization
{

using PhysicalFace = gprs_data::FractureFace;


enum ControlVolumeType
{
  cell,
  face
};


enum ConnectionType
{
  matrix_matrix = 1,
  matrix_fracture = 2,
  fracture_fracture = 3
};


struct ConnectionData
{
  ConnectionType type;
  std::vector<double> coefficients;
  std::vector<size_t> elements;
  double area;
  angem::Point<3,double> normal;
  angem::Point<3,double> center;
};


struct ControlVolumeData
{
  ControlVolumeType type;
  // map control volume to cell/face index
  std::size_t       master;
  angem::Point<3,double> center;
  double volume;
  double porosity;
  angem::Tensor2<3,double> permeability;
  std::vector<double> custom;
};


/* This is an abstract base class for
 * all discretization classes out there. */
class DiscretizationBase
{
 public:
  DiscretizationBase(const mesh::Mesh                                    & grid,
                     const std::set<int>                                 & dfm_markers,
                     const std::vector<std::vector<double>>              & props,
                     const std::vector<std::string>                      & keys);

  // get a reference to the face_data vector
  std::vector<ConnectionData> & get_face_data();
  // get a reference to the cell_data vector
  std::vector<ControlVolumeData> & get_cell_data();
  // get vector of strings of custom data keys
  std::vector<std::string> get_custom_keys() const;
  // main method. build the discretization
  virtual void build() = 0;

 protected:
  // build control volumes data
  virtual void build_cell_data();
  // is a face a dfm fracture
  bool is_fracture (const int marker);
  angem::Tensor2<3,double> get_permeability(const std::size_t cell) const;
  double get_porosity(const std::size_t cell) const;
  // calculate the positions of the perm keys (in keys array)
  void infer_perm_assignment();
  // calculate the positions of the porosity keys (in keys array)
  void infer_poro_assignment();
  // calculate keys that are not reserved (poro, perm)
  void infer_custom_keys();

  // ATTRIBUTES
  // input
  const mesh::Mesh & grid;
  const std::set<int> & dfm_markers;
  const std::vector<std::vector<double>> & props;
  const std::vector<std::string> & keys;

  //  computed
  std::array<int, 9> perm_keys = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
  int poro_key = -1;
  std::vector<size_t> custom_keys;
  std::vector<ConnectionData> con_data;
  std::vector<ControlVolumeData> cv_data;
};

}
