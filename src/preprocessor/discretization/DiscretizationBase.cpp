#include "DiscretizationBase.hpp"

namespace discretization
{

DiscretizationBase::
DiscretizationBase(const mesh::Mesh                                    & grid,
                   const std::set<int>                                 & dfm_markers,
                   const std::vector<std::vector<double>>              & props,
                   const std::vector<std::string>                      & keys)
    : grid(grid),
      dfm_markers(dfm_markers),
      props(props),
      keys(keys)
{
  infer_perm_assignment();
  infer_poro_assignment();
  infer_custom_keys();
}


bool DiscretizationBase::is_fracture(const int marker)
{
  const auto it = dfm_markers.find(marker);
  if (it != dfm_markers.end())
    return true;
  else return false;
}


angem::Tensor2<3,double>
DiscretizationBase::get_permeability(const std::size_t cell) const
{
  assert(cell < props.size());
  angem::Tensor2<3,double> K;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      K(i, j) = (perm_keys[3*i+j] >= 0) ? props[cell][perm_keys[3*i+j]] : 0;
  return K;
}


double DiscretizationBase::get_porosity(const std::size_t cell) const
{
  assert(cell < props.size());
  assert(poro_key > 0 && poro_key < keys.size());
  return props[cell][poro_key];
}


void DiscretizationBase::infer_perm_assignment()
{
  bool found_perm_x = false;
  bool found_perm_y = false;
  bool found_perm_z = false;

  for (std::size_t i = 0; i < keys.size(); i++)
  {
    const auto & key = keys[i];
    if (key == "PERMX")
    {
      found_perm_x = true;
      perm_keys[0] = i;
    }
    else if (key == "PERMY")
    {
      found_perm_y = true;
      perm_keys[1*3 + 1] = i;
    }
    else if (key == "PERMZ")
    {
      found_perm_z = true;
      perm_keys[2*3 + 2] = i;
    }
    else if (key == "PERM")
    {
      perm_keys[0] = i;
      perm_keys[1*3 + 1] = i;
      perm_keys[2*3 + 2] = i;
      found_perm_x = true;
      found_perm_y = true;
      found_perm_z = true;
    }
  }

  if (found_perm_x)
    if (found_perm_y)
      if (found_perm_z)
        return;

  throw std::invalid_argument("permebility is undefined");
  return;
}


void DiscretizationBase::infer_custom_keys()
{
  for (size_t j = 0; j < keys.size(); j++)
  {
    if (std::find(perm_keys.begin(), perm_keys.end(), j) == perm_keys.end())
        if (j != poro_key)
      custom_keys.push_back(j);
  }
}



void DiscretizationBase::infer_poro_assignment()
{
  for (std::size_t i = 0; i < keys.size(); i++)
  {
    const auto & key = keys[i];
    if (key == "PORO")
    {
      poro_key = i;
      return;
    }
  }

  throw std::invalid_argument("porosity is undefined");
}


void DiscretizationBase::build_cell_data()
{
  cv_data.resize(grid.n_cells());
  for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
  {
    const std::size_t i = cell.index();
    auto & data = cv_data[i];
    data.type = ControlVolumeType::cell;
    data.master = i;
    data.porosity = get_porosity(i);
    data.permeability = get_permeability(i);
    data.center = cell.center();
    data.volume = cell.volume() * data.porosity;

    data.custom.resize(custom_keys.size());
    for (size_t j = 0; j < custom_keys.size(); ++j)
      data.custom[j] = props[i][custom_keys[j]];
  }
}


std::vector<ControlVolumeData> & DiscretizationBase::get_cell_data()
{
  return cv_data;
}


std::vector<ConnectionData> & DiscretizationBase::get_face_data()
{
  return con_data;
}


std::vector<std::string> DiscretizationBase::get_custom_keys() const
{
  std::vector<std::string> custom_key_names;
  for (const size_t i : custom_keys)
    custom_key_names.push_back(keys[i]);
  return custom_key_names;
}

}
