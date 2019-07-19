/************************************************************************/
/* DISCRETE FEATURE MODEL                                               */
/* by Mohammad Karimi-Fard (karimi@stanford.edu)                        */
/* March 2007, Stanford, CA.                                            */
/*                                                                      */
/* December 2012. Modifications by Timur Garipov                        */
/* Geomechanical interface. Convert Tetgen data to Karimi data          */
/* September 2018 - porting into C++ by Igor Shovkun                    */
/************************************************************************/
#pragma once

#include <FlowData.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include <algorithm>
#include <math.h>
#include <iterator>
#include <vector>
#include <set>
#include <time.h>


namespace flow
{


/* This is a class for computing flow properties.
 * It compute transmissibilities and volumes.
 * It features fracture transmissibilities and
 * fracture-fracture transmissibilities through the
 * start-delta transformation.*/
class CalcTranses
{
public:
   CalcTranses();
  ~CalcTranses();
  void compute_flow_data();
  static void save_output(const FlowData    & data,
                          const std::string & output_dir);
  // this guy writes text output in a series of files
  void writeOutputFiles(const std::string & output_path) const;
  void extractData(FlowData & data) const;
  void init();


public:
  int NbNodes;
  int NbPolygons;
  int NbPolyhedra;
  int NbZones;
  int NbFracs;
  int NbOptions;
  // coordinates
  std::vector<double>	X,Y,Z;
  // faces
  std::vector<std::vector<std::size_t> > vvVFaces;
  std::vector<int> vCodePolygon;
  // elements
  std::vector<std::vector<std::size_t> > vvFNodes;
  // polyhedron code is the index of a cell
  std::vector<int> vCodePolyhedron;
  // -------
  //  PROPERTIES
  // region from which the properties are taken
  std::vector<int> vZoneCode;
  // volume factor for each zone
  std::vector<double> vZVolumeFactor;
  // modify conenction area
  std::vector<double> vTimurConnectionFactor;
  // porosity
  std::vector<double> vZPorosity;
  // not actually used as a matter of fact i think
  std::vector<int> vZPermCode;
  // permeability
  std::vector<double> vZPermeability;
  // thermal conductivity
  std::vector<double> vZConduction;

protected:
  void ProjectionA(double mx,double my,double mz,
                   double px,double py,double pz,
                   double ix,double iy,double iz,
                   double nx,double ny,double nz,
                   double *hx,double *hy,double *hz);
  void ProjectionB(double mx,double my,double mz,
                   double ix,double iy,double iz,
                   double ux,double uy,double uz,
                   double *hx,double *hy,double *hz);

  void ComputeBasicGeometry();
  void ComputeControlVolumeList();
  void PrepareConnectionList();
  void ConstructConnectionList();
  void VolumeCorrection();
  void ComputeContinuityNode();
  void ComputeDirectionalPermeability();
  void ComputeTransmissibilityPart();
  void ComputeTransmissibilityList();

protected:
  std::size_t	NbCVs,NbVolumes,NbInterfaces,NbEquations,NbFeatures,NbCF,NbConnections,NbIntersections;
  std::size_t	NbTransmissibility,NbMetric;
  int		OptionVC,OptionGO,OptionMC;

  ///// Grid information /////
  std::vector<std::vector<std::size_t>>		FNodes;

  std::vector<int>		CodePolygon;
  std::vector<int>		CodePolyhedron;

  ///// Control volume numbers /////
  std::vector<int>		EQF;
  std::vector<int>		EQV;

  ///// Additional polygon information /////
  std::vector<double>		FArea;
  std::vector<double>		FXG,FYG,FZG;
  std::vector<double>		Fnx,Fny,Fnz;		// unit vector

  ///// Additional polyhedron information /////
  std::vector<double>		VVolume;
  std::vector<double>		VXG,VYG,VZG;

  ///// Definition of the control volumes /////
  double		Tolerance;

  std::vector<int>		CVType;
  std::vector<int>		CVZone;
  std::vector<double>	CVx, CVy, CVz;
  std::vector<double>		CVVolume;

  std::vector<int>		ZoneCode;
  std::vector<double>	ZVolumeFactor;
  std::vector<double>	ZPorosity;
  std::vector<int>		ZPermCode;
  std::vector<std::vector<double>>	ZPermeability;
  std::vector<std::vector<double>>	ZConduction;
  double		K1,K2,K3,K4,K5,K6;
  double		Kx,Ky,Kz;

///// Definition of the connections /////

  std::vector<int>		ConType;
  std::vector<int>		ConN;
  std::vector<std::vector<int>>		ConCV;
  std::vector<std::vector<double>>		ConTr;

  std::vector<std::vector<double>>		ConGeom;
  std::vector<std::vector<double>>		ConMult;

  std::vector<std::vector<double>>	ConArea;
  std::vector<std::vector<double>>	ConPerm;
  std::vector<double>		ConP1x, ConP1y, ConP1z;
  std::vector<double>		ConP2x, ConP2y, ConP2z;
  std::vector<double>		ConIx, ConIy, ConIz;
  std::vector<double>		ConVx, ConVy, ConVz;
  std::vector<double>		Conhx, Conhy, Conhz;

///// Transmissibility List /////

  std::vector<int>		iTr, jTr;
  std::vector<double>	Tij,SumTr, TConductionIJ;

  int		NbEdges;	// Actives
  std::vector<int>		ListV1, ListV2;
  std::vector<int>		ListE1, ListE2, ListF;

  double		TotalVolume,FaceArea,delta;
  double		LocalDistance,t,xt,yt,zt;

  clock_t		t1,t2,t3,t4,t5,t6,t7,t8,t9;
  clock_t		Deb_Computing,Fin_Computing;

 public:
  double fracporo;
  static constexpr double transmissibility_conversion_factor =
      0.0085267146719160104986876640419948;
};

}
