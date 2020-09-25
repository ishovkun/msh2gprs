#pragma once
#include "DiscretizationEDFM.hpp"
#include "../../EmbeddedFractureManager.hpp"

namespace discretization {

class DiscretizationPEDFM : public DiscretizationEDFM {
 public:
  DiscretizationPEDFM(const DoFNumbering & split_dof_numbering,
                      const DoFNumbering & combined_dof_numbering,
                      gprs_data::SimData & data,
                      std::vector<ControlVolumeData> & cv_data,
                      std::vector<ConnectionData> & connection_data,
                      const gprs_data::EmbeddedFractureManager & edfm_mgr);
  virtual ~DiscretizationPEDFM() = default;

  void build() override;

 private:
  void build_non_neighboring_connections_();
  std::list<const mesh::Face*> select_faces_(const mesh::Face & frac_face) const;
  bool non_branching_face_(const mesh::Face & frac_face, const mesh::Face & face) const;
  const mesh::Cell* smallest_neighbor_(const mesh::Face & face) const;
  bool smaller_cut_part_above_(const mesh::Cell* neighbor, const angem::Plane<double> &frac_plane) const;
  bool have_common_vertices_(const mesh::Face & face1, const mesh::Face & face2) const;
  size_t other_cell_(const mesh::Face & frac_face, const mesh::Face & isolating_face) const;
  void build_matrix_matrix_(size_t dof1, size_t dof2, double projection_area);
  void build_fracture_matrix_(size_t dof1, size_t dof2, double projection_area,
                              const angem::Point<3,double> & normal);
  void finalize_connections_();
  size_t host_cell_index_(const mesh::Face & frac_face) const;

  std::pair<bool, std::unique_ptr<angem::Polygon<double>>>
  project_(const angem::Polygon<double> & poly1, const angem::Polygon<double> & poly2) const;

  const gprs_data::EmbeddedFractureManager & _edfm_mgr;
  PureConnectionMap _cleared_connections;
};

}  // end namespace discretization
