#pragma once
#include <string>

struct OutputConfigGPRS
{
  std::string geometry_file          = "gm_geometry.txt";
  std::string mechanics_kwd_file     = "gm_keywords.txt";
  std::string efrac_file             = "gm_SDA.txt";
  std::string discrete_frac_file     = "gm_DFM.txt";
  std::string bcond_file             = "gm_bcond.txt";
  std::string wells_file             = "wells.txt";
  std::string mech_ms_file           = "ms_mech.txt";
  std::string flow_ms_file           = "ms_flow.txt";
  std::string flow_cv_file           = "fl_cell_data.txt";
  std::string flow_connection_file   = "fl_face_data.txt";
  std::string mech_trans_update_file = "gm_update_trans.txt";
  std::string fem_file               = "gm_fem.txt";
  double transmissibility_mult = 0.0085267146719160104986876640419948;
};
