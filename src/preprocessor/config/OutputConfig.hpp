#pragma once
#include "OutputConfigINSIM.hpp"
#include "OutputConfigGPRS.hpp"
#include "OutputConfigVTK.hpp"
#include "OutputConfigPostprocessor.hpp"
#include <vector>

enum class OutputFormat
{
  gprs, vtk, postprocessor, insim
};

struct OutputConfig
{
  std::vector<OutputFormat> output_formats = {OutputFormat::gprs,
                                              OutputFormat::vtk,
                                              OutputFormat::postprocessor};

  OutputConfigGPRS gprs;
  OutputConfigINSIM insim;
  OutputConfigVTK vtk;
  OutputConfigPostprocessor postprocessor;

  std::string dir = "output";
};
