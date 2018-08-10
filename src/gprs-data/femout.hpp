#ifndef __OUTDATA
#define __OUTDATA
#include "simdata.hpp"
#include "transes.hpp"
class SimData;
class CalcTranses;

class OutputData
{
public:
   OutputData(SimData * pSimData);
  ~OutputData();

  void writeGeomechDataNewKeywords(const std::string & output_path);

protected:
  SimData * pSim;
};
#endif
