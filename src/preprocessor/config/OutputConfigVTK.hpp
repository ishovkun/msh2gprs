#pragma once
#include <string>

enum VTKOutputFlags {
  save_fractures  = 0x0001,
  save_wells      = 0x0002,
  save_flow_graph = 0x0003,
};

struct OutputConfigVTK
{
  std::string flow_reservoir_grid_file      = "flow_reservoir_mesh.vtk";
  std::string mechanics_reservoir_grid_file = "mechanics_reservoir_mesh.vtk";
  std::string fracture_grid_file            = "fractures.vtk";
  std::string wells_file                    = "wells.vtk";
  std::string flow_graph_file               = "flow_graph.vtk";
  VTKOutputFlags flags                      = VTKOutputFlags::save_wells;
};

inline VTKOutputFlags operator|(VTKOutputFlags f1, VTKOutputFlags f2)
{
  return static_cast<VTKOutputFlags>(static_cast<unsigned int>(f1) |
                                     static_cast<unsigned int>(f2));
}

inline VTKOutputFlags& operator|=(VTKOutputFlags &f1, VTKOutputFlags f2)
{
  f1 = f1 | f2;
  return f1;
}


inline VTKOutputFlags operator&(VTKOutputFlags f1, VTKOutputFlags f2)
{
  return static_cast<VTKOutputFlags>(static_cast<unsigned int>(f1) &
                                     static_cast<unsigned int>(f2));
}

inline VTKOutputFlags& operator&=(VTKOutputFlags& f1, VTKOutputFlags f2)
{
  f1 = f1 & f2;
  return f1;
}

template <class StreamType>
inline StreamType &
operator<<(StreamType &s, VTKOutputFlags flags)
{
  if ( flags & save_fractures )
    s << "fractures|";
  if ( flags & save_wells)
    s << "wells|";
  if ( flags & save_flow_graph )
    s << "flow_graph|";
  return s;
}
