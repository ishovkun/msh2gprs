#define _CRT_SECURE_NO_DEPRECATE
#include "transes.hpp"
#include "simdata.hpp"
#include <random>
class SimData;

CalcTranses::CalcTranses()
{}


void CalcTranses::init()
{
  //coordinates
  vCoordinatesX.resize(NbNodes);
  vCoordinatesY.resize(NbNodes);
  vCoordinatesZ.resize(NbNodes);
  //faces
  vNbVFaces.resize(NbPolyhedra);
  vvVFaces.resize(NbPolyhedra);
  //elements
  vNbFNodes.resize(NbPolygons, 3);
  vvFNodes.resize(NbPolygons);
  vCodePolyhedron.resize(NbPolyhedra);
  //properties
  vZoneCode.resize(NbZones);
  vZVolumeFactor.resize(NbZones);
  vTimurConnectionFactor.resize(NbZones);
  vZPorosity.resize(NbZones);
  vZPermCode.resize(NbZones);
  vZPermeability.resize(NbZones);

  OptionVC = 0;
  OptionGO = 0;
  OptionMC = 0;
  Tolerance = 0.05;

  vCoordinatesX.resize( NbNodes );
  vCoordinatesY.resize( NbNodes );
  vCoordinatesZ.resize( NbNodes );
}

CalcTranses::~CalcTranses()
{
  // for (std::size_t)
  // free(vNbFNodes);
  // // **FNodes -- need to remove subelements but i'm lazy
  // for ()
  // free(FNodes);
  // free(NbVFaces);
  // free (VFaces);
}

void CalcTranses::createKarimiData()
{}

/********************************************************************/
double CalcTranses::ABS(double v)
{
    if (v>=0.) return(v);
    return(-v);
}
/********************************************************************/
void CalcTranses::CheckIt(int TEST)
{
  if (TEST == 1)
  {
    printf("Memory Allocation Problem.\n");
    exit(0);
  }
}
/********************************************************************/
void CalcTranses::ProjectionA( double mx,double my,double mz,
                               double px,double py,double pz,
                               double ix,double iy,double iz,
                               double nx,double ny,double nz,
                               double *hx,double *hy,double *hz)
{
  double  t;

  t=((ix-mx)*nx+(iy-my)*ny+(iz-mz)*nz)/(px*nx+py*ny+pz*nz);
  *hx = mx + t*px;
  *hy = my + t*py;
  *hz = mz + t*pz;
}
/********************************************************************/
void CalcTranses::ProjectionB( double mx,  double my,  double mz,
                               double ix,  double iy,  double iz,
                               double ux,  double uy,  double uz,
                               double *hx, double *hy, double *hz)
{
    double  t;

    t = (ux*(mx-ix) + uy*(my-iy) + uz*(mz-iz)) / (ux*ux + uy*uy + uz*uz);

    *hx = ix + t*ux;
    *hy = iy + t*uy;
    *hz = iz + t*uz;
}
/********************************************************************/
void CalcTranses::ComputeBasicGeometry()
{
    // int i, j, k;
    int k;
    double  xi, yi, zi, areatmp, volumetmp;
    double  ux, uy, uz,
        vx, vy, vz,
        nx,ny,nz,
        nl,h;

//////////////////////////////////////////////////////////////////
///// Polygon Area, Center of Mass, and Normal (Unit vector) /////
//////////////////////////////////////////////////////////////////

    // FArea = (double*)malloc(NbPolygons*(sizeof(double)));
    // FXG   = (double*)malloc(NbPolygons*(sizeof(double)));
    // FYG   = (double*)malloc(NbPolygons*(sizeof(double)));
    // FZG   = (double*)malloc(NbPolygons*(sizeof(double)));
    // Fnx   = (double*)malloc(NbPolygons*(sizeof(double)));
    // Fny   = (double*)malloc(NbPolygons*(sizeof(double)));
    // Fnz   = (double*)malloc(NbPolygons*(sizeof(double)));
    FArea.resize(NbPolygons);
    FXG.resize(NbPolygons);
    FYG.resize(NbPolygons);
    FZG.resize(NbPolygons);
    Fnx.resize(NbPolygons);
    Fny.resize(NbPolygons);
    Fnz.resize(NbPolygons);

    for (std::size_t i=0; i<NbPolygons; i++)
    {
      FArea[i] = FXG[i] = FYG[i] = FZG[i] = Fnx[i] = Fny[i] = Fnz[i] = 0;
      for (std::size_t j=1; j<NbFNodes[i]-1; j++)
      {
        ux = X[FNodes[i][j]] - X[FNodes[i][0]];
        uy = Y[FNodes[i][j]] - Y[FNodes[i][0]];
        uz = Z[FNodes[i][j]] - Z[FNodes[i][0]];
        vx = X[FNodes[i][j+1]] - X[FNodes[i][0]];
        vy = Y[FNodes[i][j+1]] - Y[FNodes[i][0]];
        vz = Z[FNodes[i][j+1]] - Z[FNodes[i][0]];
        nx = (uy*vz - vy*uz);
        ny = (vx*uz - ux*vz);
        nz = (ux*vy - uy*vx);
        areatmp = .5*sqrt(nx*nx + ny*ny + nz*nz);

        FArea[i] += areatmp;
        FXG[i] += areatmp*(X[FNodes[i][0]] + X[FNodes[i][j]] + X[FNodes[i][j+1]])/3.;
        FYG[i] += areatmp*(Y[FNodes[i][0]] + Y[FNodes[i][j]] + Y[FNodes[i][j+1]])/3.;
        FZG[i] += areatmp*(Z[FNodes[i][0]] + Z[FNodes[i][j]] + Z[FNodes[i][j+1]])/3.;

        Fnx[i] +=.5*nx;
        Fny[i] +=.5*ny;
        Fnz[i] +=.5*nz;
      }

        FXG[i] = FXG[i] / FArea[i];
        FYG[i] = FYG[i] / FArea[i];
        FZG[i] = FZG[i] / FArea[i];

        Fnx[i] = Fnx[i] / FArea[i];
        Fny[i] = Fny[i] / FArea[i];
        Fnz[i] = Fnz[i] / FArea[i];

        nl = sqrt(Fnx[i]*Fnx[i] + Fny[i]*Fny[i] + Fnz[i]*Fnz[i]);

        Fnx[i] = Fnx[i] / nl;
        Fny[i] = Fny[i] / nl;
        Fnz[i] = Fnz[i] / nl;
    }

////////////////////////////////////////////////
///// Polyhedron Volume and Center of Mass /////
////////////////////////////////////////////////

    // VVolume = (double*)malloc(NbPolyhedra*(sizeof(double)));
    // VXG     = (double*)malloc(NbPolyhedra*(sizeof(double)));
    // VYG     = (double*)malloc(NbPolyhedra*(sizeof(double)));
    // VZG     = (double*)malloc(NbPolyhedra*(sizeof(double)));

    VVolume.resize(NbPolyhedra);
    VXG.resize(NbPolyhedra);
    VYG.resize(NbPolyhedra);
    VZG.resize(NbPolyhedra);

    for (std::size_t i=0; i<NbPolyhedra; i++)
    {
        VVolume[i] = VXG[i] = VYG[i] = VZG[i] = 0;
        xi = yi = zi = 0;
        for (std::size_t j=0; j<NbVFaces[i]; j++) // Defining a node inside the polyhedron
        {
            xi += FXG[VFaces[i][j]];
            yi += FYG[VFaces[i][j]];
            zi += FZG[VFaces[i][j]];
        }

        xi = xi / NbVFaces[i];
        yi = yi / NbVFaces[i];
        zi = zi / NbVFaces[i];

        for (std::size_t j=0; j<NbVFaces[i]; j++)
        {
          k = VFaces[i][j];
          h = Fnx[k]*(FXG[k]-xi) +
              Fny[k]*(FYG[k]-yi) +
              Fnz[k]*(FZG[k]-zi);

          volumetmp = ABS(h*FArea[k]) / 3.;

          VXG[i] += (FXG[k] + .25*(xi-FXG[k])) * volumetmp;
          VYG[i] += (FYG[k] + .25*(yi-FYG[k])) * volumetmp;
          VZG[i] += (FZG[k] + .25*(zi-FZG[k])) * volumetmp;
          VVolume[i] += volumetmp;
        }

        VXG[i] = VXG[i] / VVolume[i];
        VYG[i] = VYG[i] / VVolume[i];
        VZG[i] = VZG[i] / VVolume[i];

    }
}
/********************************************************************/
void CalcTranses::ComputeControlVolumeList()
{
  // int j;

  // CVType   =  (int*)malloc(NbCVs*(sizeof(int)));
  // CVZone   =  (int*)malloc(NbCVs*(sizeof(int)));
  // CVVolume =  (double*)malloc(NbCVs*(sizeof(double)));
  // CVx      =  (double*)malloc(NbCVs*(sizeof(double)));
  // CVy      =  (double*)malloc(NbCVs*(sizeof(double)));
  // CVz      =  (double*)malloc(NbCVs*(sizeof(double)));

  CVType.resize(NbCVs);
  CVZone.resize(NbCVs);
  CVVolume.resize(NbCVs);
  CVx.resize(NbCVs);
  CVy.resize(NbCVs);
  CVz.resize(NbCVs);

  for (std::size_t i=0; i<NbPolygons; i++)
    if (EQF[i] != -1) // Active polygon
    {
      std::size_t j = EQF[i];
      CVZone[j] = CodePolygon[i];
      CVType[j] = 2;  // Feature
      CVx[j] = FXG[i];
      CVy[j] = FYG[i];
      CVz[j] = FZG[i];
      CVVolume[j] = FArea[i] * ZVolumeFactor[CVZone[j]];
    }

  for (int i=0;i<NbPolyhedra;i++)
  {
    if (EQV[i] != -1) // Active polyhedron
    {
      std::size_t j = EQV[i];
      CVZone[j] = CodePolyhedron[i];
      CVType[j] = 1;  // Volume
      CVx[j] = VXG[i];
      CVy[j] = VYG[i];
      CVz[j] = VZG[i];
      CVVolume[j] = VVolume[i] * ZVolumeFactor[CVZone[j]];
    }
  }
}
/********************************************************************/
void CalcTranses::PrepareConnectionList()
{
    // int i,j,k,iswap;
    int iswap;

///// Direct construction of M-M and M-F connections /////

    // ListV1 = (int*)malloc(NbPolygons*(sizeof(int)));
    // ListV2 = (int*)malloc(NbPolygons*(sizeof(int)));

    ListV1.resize(NbPolygons);
    ListV2.resize(NbPolygons);

    for (std::size_t i=0; i<NbPolygons; i++)
      ListV1[i] = ListV2[i] = -1;

    for (std::size_t i=0; i<NbPolyhedra; i++)
      for (std::size_t j=0; j<NbVFaces[i]; j++)
      {
        if (ListV1[VFaces[i][j]] == -1)
          ListV1[VFaces[i][j]] = i;
        else
          ListV2[VFaces[i][j]] = i;
      }

///// Construction of F-F needs sorting... less efficient /////

    // ListE1 = (int*)malloc(NbEdges*(sizeof(int)));
    // ListE2 = (int*)malloc(NbEdges*(sizeof(int)));
    // ListF  = (int*)malloc(NbEdges*(sizeof(int)));

    ListE1.resize(NbEdges);
    ListE2.resize(NbEdges);
    ListF.resize(NbEdges);

    std::size_t k=0;
    for (std::size_t i=0; i<NbPolygons; i++)
    {
      if (CodePolygon[i] >= 0)
      {
        for (std::size_t j=0; j<NbFNodes[i]-1; j++)
        {
          if (FNodes[i][j] < FNodes[i][j+1])
          {
            ListE1[k] = FNodes[i][j];
            ListE2[k] = FNodes[i][j+1];
            ListF [k] = i;
          }
          else
          {
            ListE1[k] = FNodes[i][j+1];
            ListE2[k] = FNodes[i][j];
            ListF [k] = i;
          }
          k++;
        }
        if (FNodes[i][0] < FNodes[i][NbFNodes[i]-1])
        {
          ListE1[k] = FNodes[i][0];
          ListE2[k] = FNodes[i][NbFNodes[i]-1];
          ListF [k] = i;
        }
        else
        {
          ListE1[k] = FNodes[i][NbFNodes[i]-1];
          ListE2[k] = FNodes[i][0];
          ListF [k] = i;
        }
        k++;
      }
    }

///// Two step sorting /////

// First

    for (std::size_t i=0; i<NbEdges-1; i++)
      for (std::size_t j=i+1; j<NbEdges; j++)
      {
        if (ListE1[i] > ListE1[j])
        {
          iswap     = ListE1[i];
          ListE1[i] = ListE1[j];
          ListE1[j] = iswap;

          iswap     = ListE2[i];
          ListE2[i] = ListE2[j];
          ListE2[j] = iswap;

          iswap     = ListF[i];
          ListF[i]  = ListF[j];
          ListF[j]  = iswap;
        }
      }

// Second

    for (std::size_t i=0; i<NbEdges-1; i++)
      for (std::size_t j=i+1; j<NbEdges; j++)
      {
        if (ListE1[i] == ListE1[j] && ListE2[i] > ListE2[j])
        {
          iswap     = ListE2[i];
          ListE2[i] = ListE2[j];
          ListE2[j] = iswap;

          iswap    = ListF[i];
          ListF[i] = ListF[j];
          ListF[j] = iswap;
        }
      }

///// Evaluation of the number of connections /////

// M-M and M-F

    NbConnections = 0;

    for (std::size_t i=0; i<NbPolygons; i++)
    {
      if (CodePolygon[i] < 0)
      {
        if (ListV2[i] != -1)
          NbConnections++;
      }
      else
      {
        if (ListV2[i] == -1)
          NbConnections++;
        else
          NbConnections += 2;
      }
    }

// F-F

    std::size_t i = 0;
    while (i < NbEdges-1)
    {
      std::size_t j = i + 1;
      while ((j < NbEdges) && (ListE1[i] == ListE1[j]) && (ListE2[i] == ListE2[j]))
        j++;

      if ((j-i) >= 2)
        NbConnections++;

      i = j;
    }

}
/********************************************************************/
void CalcTranses::ConstructConnectionList()
{
    int i,j,k,n;

    // ConType = (int*)malloc(NbConnections*(sizeof(int)));
    // CheckIt(ConType==NULL);
    // ConN  = (int*)malloc(NbConnections*(sizeof(int)));
    // CheckIt(ConN==NULL);
    // ConCV = (int**)malloc(NbConnections*(sizeof(int*)));
    // CheckIt(ConCV==NULL);
    // ConTr = (double**)malloc(NbConnections*(sizeof(double*)));
    // CheckIt(ConTr==NULL);
    // ConArea = (double**)malloc(NbConnections*(sizeof(double*)));
    // CheckIt(ConArea==NULL);
    ConType.resize(NbConnections);
    ConN.resize(NbConnections);
    ConCV.resize(NbConnections);
    ConTr.resize(NbConnections);
    ConArea.resize(NbConnections);

    ConGeom.resize(NbConnections);
    ConMult.resize(NbConnections);

    ConPerm.resize(NbConnections);
    ConP1x.resize(NbConnections);
    ConP1y.resize(NbConnections);
    ConP1z.resize(NbConnections);

    ConP2x.resize(NbConnections);
    ConP2y.resize(NbConnections);
    ConP2z.resize(NbConnections);

    ConIx.resize(NbConnections);
    ConIy.resize(NbConnections);
    ConIz.resize(NbConnections);

    ConVx.resize(NbConnections);
    ConVy.resize(NbConnections);
    ConVz.resize(NbConnections);


    NbTransmissibility = 0;

    k=0;
    std::cout << "NbPolygons = " << NbPolygons << std::endl;
    for (i=0; i<NbPolygons; i++)
    {
      std::cout << "polygon " << i << std::endl;
      if (CodePolygon[i] < 0 && ListV2[i] >= 0)  // M-M
      {
        std::cout << "M-M" << std::endl;
        ConType[k] = 1;
        ConN[k] = 2;

        ConCV[k].resize(2);
        ConCV[k][0] = EQV[ListV1[i]];
        ConCV[k][1] = EQV[ListV2[i]];

        // ConArea[k] = (double*)malloc(2*(sizeof(double)));
        // CheckIt(ConArea[k] == NULL);
        ConArea[k].resize(2);
        ConArea[k][0] = ConArea[k][1] = FArea[i];

        ConP1x[k] = Fnx[i];
        ConP1y[k] = Fny[i];
        ConP1z[k] = Fnz[i];

        ConP2x[k] = Fnx[i];
        ConP2y[k] = Fny[i];
        ConP2z[k] = Fnz[i];

        ConIx[k] = FXG[i];
        ConIy[k] = FYG[i];
        ConIz[k] = FZG[i];

        ConVx[k] = Fnx[i];
        ConVy[k] = Fny[i];
        ConVz[k] = Fnz[i];

        ConTr[k].resize(2);
        ConPerm[k].resize(2);

        ConGeom[k].resize(2);
        ConMult[k].resize(2);

        k++;
        NbTransmissibility++;
      }

      if (CodePolygon[i] >= 0 && ListV2[i] < 0) // M-F
      {
        std::cout << "M-F" << std::endl;
        std::cout << "CodePolygon[i] = " << CodePolygon[i] << std::endl;
        std::cout << "ListV2[i] = " << ListV2[i] << std::endl;
        ConType[k] = 2;
        ConN[k] = 2;

        std::cout << "pops" << std::endl;
        ConCV[k].resize(2);
        // std::cout << ListV1[i] << std::endl;
        if (ListV1[i] < 0)
        {
          std::cout << "polygon " << i << std::endl;
          std::cout << "wrong ListV2 index " << ListV1[i] << std::endl;
          exit(0);
        }
        ConCV[k][0] = EQV[ListV1[i]];
        ConCV[k][1] = EQF[i];

        ConArea[k].resize(2);
        ConArea[k][0] = ConArea[k][1] = FArea[i];

        ConP1x[k] = Fnx[i];
        ConP1y[k] = Fny[i];
        ConP1z[k] = Fnz[i];
        std::cout << "popun" << std::endl;

        ConIx[k] = FXG[i];
        ConIy[k] = FYG[i];
        ConIz[k] = FZG[i];
        ConVx[k] = Fnx[i];
        ConVy[k] = Fny[i];
        ConVz[k] = Fnz[i];

        ConTr[k].resize(2);
        ConPerm[k].resize(2);

        ConGeom[k].resize(2);
        ConMult[k].resize(2);

        k++;
        NbTransmissibility++;
      }
      std::cout << "wheeps" << std::endl;

        if (CodePolygon[i] >= 0 && ListV2[i] >= 0)  // M-F and F-M
        {
            ConType[k] = 2;
            ConN[k] = 2;
            // ConCV[k] = (int*)malloc(2*(sizeof(int)));
            ConCV[k].resize(2);
            ConCV[k][0] = EQV[ListV1[i]];
            ConCV[k][1] = EQF[i];

            ConArea[k].resize(2);
            ConArea[k][0] = ConArea[k][1] = FArea[i];

            ConP1x[k] = Fnx[i];
            ConP1y[k] = Fny[i];
            ConP1z[k] = Fnz[i];

            ConIx[k] = FXG[i];
            ConIy[k] = FYG[i];
            ConIz[k] = FZG[i];
            ConVx[k] = Fnx[i];
            ConVy[k] = Fny[i];
            ConVz[k] = Fnz[i];

            ConTr[k] .resize(2);
            ConPerm[k].resize(2);

            ConGeom[k].resize(2);
            ConMult[k].resize(2);

            k++;

            ConType[k] = 2;
            ConN[k] = 2;
            ConCV[k].resize(2);
            ConCV[k][0] = EQV[ListV2[i]];
            ConCV[k][1] = EQF[i];

            ConArea[k].resize(2);
            ConArea[k][0] = ConArea[k][1] = FArea[i];

            ConP1x[k] = Fnx[i];
            ConP1y[k] = Fny[i];
            ConP1z[k] = Fnz[i];

            ConIx[k] = FXG[i];
            ConIy[k] = FYG[i];
            ConIz[k] = FZG[i];
            ConVx[k] = Fnx[i];
            ConVy[k] = Fny[i];
            ConVz[k] = Fnz[i];

            ConTr[k] .resize(2);
            ConPerm[k].resize(2);

            ConGeom[k].resize(2);
            ConMult[k].resize(2);
            k++;
            NbTransmissibility += 2;
        }
    }


  // F-F Connections
  i=0;
  while ( i < NbEdges-1 )
  {
    j = i+1;
    while ( ( j<NbEdges ) && ( ListE1[i]==ListE1[j] ) && ( ListE2[i]==ListE2[j] ) )
    {
      j++;
    }
    if ( ( j-i ) >= 2 )
    {
      ConType[k] = 3;
      ConN[k] = ( j-i );
      // ConCV[k] = ( int* ) malloc ( ( j-i ) * ( sizeof ( int ) ) );
      ConCV[k].resize(j-i);
      for ( n=i; n<j; n++ )
        ConCV[k][n-i] = EQF[ListF[n]];

      ConIx[k] = X[ListE1[i]];
      ConIy[k] = Y[ListE1[i]];
      ConIz[k] = Z[ListE1[i]];
      ConVx[k] = X[ListE1[i]] - X[ListE2[i]];
      ConVy[k] = Y[ListE1[i]] - Y[ListE2[i]];
      ConVz[k] = Z[ListE1[i]] - Z[ListE2[i]];

      // ConArea[k] = ( double* ) malloc ( ( j-i ) * ( sizeof ( double ) ) );
      ConArea[k].resize(j-1);

      for ( n = i; n < j; n++ ) // Double check the formula
      {
        ConArea[k][n - i] = ZVolumeFactor[CVZone[ConCV[k][n - i]]] *
            sqrt ( ConVx[k] * ConVx[k] + ConVy[k] * ConVy[k] + ConVz[k] * ConVz[k] );
        // TODO TIMUR (F-F connection)
        ConArea[k][n - i] *= vTimurConnectionFactor[CVZone[ConCV[k][n - i]]];
      }
      // ConTr[k] = ( double* ) malloc ( ( j-i ) * ( sizeof ( double ) ) );
      ConTr[k].resize(j-i);
      ConPerm[k].resize(j-i);

      ConGeom[k].resize(j-i);
      ConMult[k].resize(j-i);

      NbTransmissibility += ( ConN[k]* ( ConN[k]-1 ) ) /2;
      k++;
    }
    i = j;
  }

}
/********************************************************************/
void CalcTranses::VolumeCorrection()  // Volume should be at least twice bigger...
{
    int i;

    for (i=0; i<NbPolygons; i++)
    {
      if (CodePolygon[i] >= 0 && ListV2[i] < 0) // M-F ///
      {
        if (CVVolume[EQV[ListV1[i]]] > 2*CVVolume[EQF[i]])
          CVVolume[EQV[ListV1[i]]] -= CVVolume[EQF[i]];
      }

      if (CodePolygon[i] >= 0 && ListV2[i] >= 0)  // M-F and F-M ///
      {
        if (CVVolume[EQV[ListV1[i]]] > 2*.5*CVVolume[EQF[i]])
          CVVolume[EQV[ListV1[i]]] -= .5*CVVolume[EQF[i]];

        if (CVVolume[EQV[ListV2[i]]] > 2*.5*CVVolume[EQF[i]])
          CVVolume[EQV[ListV2[i]]] -= .5*CVVolume[EQF[i]];
      }
    }

}
/********************************************************************/
void CalcTranses::ComputeContinuityNode()
{
    int i,j,k;
    double  hx,hy,hz,px,py,pz;

    Conhx.resize(NbConnections);
    Conhy.resize(NbConnections);
    Conhz.resize(NbConnections);

    for (i=0;i<NbConnections;i++)
    {
      Conhx[i] = Conhy[i] = Conhz[i] = 0;

      if (ConType[i] == 1)  // M-M ///
      {
        px = CVx[ConCV[i][1]] - CVx[ConCV[i][0]];
        py = CVy[ConCV[i][1]] - CVy[ConCV[i][0]];
        pz = CVz[ConCV[i][1]] - CVz[ConCV[i][0]];

        ProjectionA(  CVx[ConCV[i][0]], CVy[ConCV[i][0]], CVz[ConCV[i][0]],
                      px, py, pz,
                      ConIx[i], ConIy[i], ConIz[i],
                      ConVx[i], ConVy[i], ConVz[i],
                      &Conhx[i], &Conhy[i], &Conhz[i]);
      }
      else if (ConType[i] == 2)  // M-F ///
      {
        Conhx[i] = CVx[ConCV[i][1]];
        Conhy[i] = CVy[ConCV[i][1]];
        Conhz[i] = CVz[ConCV[i][1]];
      }
      else if (ConType[i] == 3)  // F-F ///
      {
        for (j=0;j<ConN[i];j++)
        {
          ProjectionB( CVx[ConCV[i][j]], CVy[ConCV[i][j]], CVz[ConCV[i][j]],
                       ConIx[i], ConIy[i], ConIz[i],
                       ConVx[i], ConVy[i], ConVz[i],
                       &hx, &hy, &hz);
          Conhx[i] += hx;
          Conhy[i] += hy;
          Conhz[i] += hz;
        }

        Conhx[i] = Conhx[i] / ConN[i];
        Conhy[i] = Conhy[i] / ConN[i];
        Conhz[i] = Conhz[i] / ConN[i];
      }
    }
}
/********************************************************************/
void CalcTranses::ComputeDirectionalPermeability()
{
    int j,k;
    double  fx,fy,fz,fl;

    for (std::size_t i=0; i<NbConnections; i++)
    {
      if (ConType[i]==1)  // M-M /////////////////////////////////////////////////
      {
        k = ConCV[i][0];
        fx = CVx[k] - Conhx[i];
        fy = CVy[k] - Conhy[i];
        fz = CVz[k] - Conhz[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);
        fx = fx/fl;
        fy = fy/fl;
        fz = fz/fl;

        Kx = ZPermeability[CVZone[k]][0]*fx +
             ZPermeability[CVZone[k]][3]*fy+
             ZPermeability[CVZone[k]][4]*fz;

        Ky = ZPermeability[CVZone[k]][3]*fx+
             ZPermeability[CVZone[k]][1]*fy +
             ZPermeability[CVZone[k]][5]*fz;

        Kz = ZPermeability[CVZone[k]][4]*fx +
             ZPermeability[CVZone[k]][5]*fy +
             ZPermeability[CVZone[k]][2]*fz;

        ConPerm[i][0] = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);

        // @HACK
        // We dont need a REAL value of conductivity
        // We just calculate geometric part (@HACK)
        Kx=1.0; Ky=1.0; Kz=1.0;
        if(ZConduction[CVZone[k]][0] != 0.0)
        {
          Kx = (ZConduction[CVZone[k]][0]*fx +
                ZConduction[CVZone[k]][3]*fy +
                ZConduction[CVZone[k]][4]*fz) / ZConduction[CVZone[k]][0];

          Ky = (ZConduction[CVZone[k]][3]*fx +
                ZConduction[CVZone[k]][1]*fy +
                ZConduction[CVZone[k]][5]*fz) / ZConduction[CVZone[k]][0];

          Kz = (ZConduction[CVZone[k]][4]*fx +
                ZConduction[CVZone[k]][5]*fy +
                ZConduction[CVZone[k]][2]*fz) / ZConduction[CVZone[k]][0];
        }

        ConMult[i][0] = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);

        k = ConCV[i][1];
        fx = CVx[k] - Conhx[i];
        fy = CVy[k] - Conhy[i];
        fz = CVz[k] - Conhz[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);
        fx = fx/fl;
        fy = fy/fl;
        fz = fz/fl;

        Kx = ZPermeability[CVZone[k]][0]*fx +
             ZPermeability[CVZone[k]][3]*fy +
             ZPermeability[CVZone[k]][4]*fz;

        Ky = ZPermeability[CVZone[k]][3]*fx +
             ZPermeability[CVZone[k]][1]*fy +
             ZPermeability[CVZone[k]][5]*fz;

        Kz = ZPermeability[CVZone[k]][4]*fx +
             ZPermeability[CVZone[k]][5]*fy +
             ZPermeability[CVZone[k]][2]*fz;

        ConPerm[i][1] = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);

        // @HACK
        // We dont need a REAL value of conductivity
        // We just calculate geometric part
        Kx=1.0; Ky=1.0; Kz=1.0;
        if(ZConduction[CVZone[k]][0] != 0.0)
        {
          Kx = (ZConduction[CVZone[k]][0]*fx +
                ZConduction[CVZone[k]][3]*fy +
                ZConduction[CVZone[k]][4]*fz) / ZConduction[CVZone[k]][0];

          Ky = (ZConduction[CVZone[k]][3]*fx +
                ZConduction[CVZone[k]][1]*fy +
                ZConduction[CVZone[k]][5]*fz) / ZConduction[CVZone[k]][0];

          Kz = (ZConduction[CVZone[k]][4]*fx +
                ZConduction[CVZone[k]][5]*fy +
                ZConduction[CVZone[k]][2]*fz) / ZConduction[CVZone[k]][0];
        }
        ConMult[i][1]=sqrt(Kx*Kx+Ky*Ky+Kz*Kz);
      }
      if (ConType[i]==2)  // M-F /////////////////////////////////////////////////
      {
        k = ConCV[i][0];
        fx = CVx[k] - Conhx[i];
        fy = CVy[k] - Conhy[i];
        fz = CVz[k] - Conhz[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);
        fx = fx/fl;
        fy = fy/fl;
        fz = fz/fl;

        Kx = ZPermeability[CVZone[k]][0]*fx +
             ZPermeability[CVZone[k]][3]*fy +
             ZPermeability[CVZone[k]][4]*fz;

        Ky = ZPermeability[CVZone[k]][3]*fx +
             ZPermeability[CVZone[k]][1]*fy +
             ZPermeability[CVZone[k]][5]*fz;

        Kz = ZPermeability[CVZone[k]][4]*fx +
             ZPermeability[CVZone[k]][5]*fy +
             ZPermeability[CVZone[k]][2]*fz;

        ConPerm[i][0] = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);

        k = ConCV[i][1];
        ConPerm[i][1] = ZPermeability[CVZone[k]][0];    // Kn (0)

        // @HACK
        // We dont need a REAL value of conductivity
        // We just calculate geometric part
        Kx=1.0; Ky=1.0; Kz=1.0;
        if(ZConduction[CVZone[k]][0] != 0.0)
        {
          Kx = (ZConduction[CVZone[k]][0]*fx +
                ZConduction[CVZone[k]][3]*fy +
                ZConduction[CVZone[k]][4]*fz) / ZConduction[CVZone[k]][0];

          Ky = (ZConduction[CVZone[k]][3]*fx +
                ZConduction[CVZone[k]][1]*fy +
                ZConduction[CVZone[k]][5]*fz) / ZConduction[CVZone[k]][0];

          Kz = (ZConduction[CVZone[k]][4]*fx +
                ZConduction[CVZone[k]][5]*fy +
                ZConduction[CVZone[k]][2]*fz) / ZConduction[CVZone[k]][0];
        }

        ConMult[i][0] = sqrt(Kx*Kx + Ky*Ky + Kz*Kz);
        ConMult[i][1] = 1.0;
      }
      if (ConType[i] == 3)  // F-F /////////////////////////////////////////////////
      {
        for (std::size_t j=0; j<ConN[i]; j++)
        {
          k = ConCV[i][j];
          ConPerm[i][j] = ZPermeability[CVZone[k]][0];  // Kp (1)

          // @HACK
          // for fractures all directions are identiacal
          // and multiplicator is UNIT
          ConMult[i][j] = 1.0;
        }
      }
    }
}
/********************************************************************/
void CalcTranses::ComputeTransmissibilityPart()
{
    // int i,j,k;
    double  nx,ny,nz,nl,fx,fy,fz,fl;

    for (std::size_t i=0; i<NbConnections; i++)
    {
      if (ConType[i] == 1)  // M-M ///
      {
        std::size_t k  = ConCV[i][0];
        nx = ConP1x[i];
        ny = ConP1y[i];
        nz = ConP1z[i];

        fx = CVx[k] - Conhx[i];
        fy = CVy[k] - Conhy[i];
        fz = CVz[k] - Conhz[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);
        fx = fx/fl;
        fy = fy/fl;
        fz = fz/fl;

        //    ConTr[i][0]=ConArea[i][0]*ConPerm[i][0]*ABS(nx*fx+ny*fy+nz*fz)/fl;
        ConTr[i][0]   = ConArea[i][0]*ConPerm[i][0] * 1./fl;
        ConGeom[i][0] = ConArea[i][0]*ConMult[i][0] * 1./fl;

        k  = ConCV[i][1];
        nx = ConP2x[i];
        ny = ConP2y[i];
        nz = ConP2z[i];

        fx = CVx[k] - Conhx[i];
        fy = CVy[k] - Conhy[i];
        fz = CVz[k] - Conhz[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);
        fx = fx/fl;
        fy = fy/fl;
        fz = fz/fl;

        //    ConTr[i][1]=ConArea[i][1]*ConPerm[i][1]*ABS(nx*fx+ny*fy+nz*fz)/fl;
        ConTr[i][1]   = ConArea[i][1] * ConPerm[i][1] * 1./fl;
        ConGeom[i][1] = ConArea[i][1] * ConMult[i][1] * 1./fl;
      }
      else if (ConType[i] == 2)  // M-F //
      {
        std::size_t k = ConCV[i][0];
        nx = ConP1x[i];
        ny = ConP1y[i];
        nz = ConP1z[i];

        fx = CVx[k] - Conhx[i];
        fy = CVy[k] - Conhy[i];
        fz = CVz[k] - Conhz[i];
        fl = sqrt(fx*fx + fy*fy + fz*fz);
        fx = fx/fl;
        fy = fy/fl;
        fz = fz/fl;

        ConTr[i][0]   = ConArea[i][0] * ConPerm[i][0] * 1./fl;
        ConGeom[i][0] = ConArea[i][0] * ConMult[i][0] * 1./fl;

        k = ConCV[i][1];

        ConTr[i][1] = ConArea[i][1] * ConPerm[i][1] *
            1./(.5*ZVolumeFactor[CVZone[k]]);

        ConGeom[i][1] = ConArea[i][1] * ConMult[i][1] *
            1./(.5*ZVolumeFactor[CVZone[k]]);
      }
      else if (ConType[i] == 3)  // F-F //
      {
        for (std::size_t j=0; j<ConN[i]; j++)
        {
          std::size_t k  = ConCV[i][j];
          fx = CVx[k] - Conhx[i];
          fy = CVy[k] - Conhy[i];
          fz = CVz[k] - Conhz[i];
          fl = sqrt(fx*fx + fy*fy + fz*fz);
          fx = fx/fl;
          fy = fy/fl;
          fz = fz/fl;

          ConTr[i][j] = ConArea[i][j] * ConPerm[i][j] * 1./fl;
          ConGeom[i][j] = ConArea[i][j] * ConMult[i][j] * 1./fl;
          }
      }
    }
}
/********************************************************************/
void CalcTranses::ComputeTransmissibilityList()
{
    // int i,j,k,n,m;
    double  SumTr;
    FILE * poutfile;
    FILE * pFile;

    // iTr           = (int*)malloc(NbTransmissibility*(sizeof(int)));
    // jTr           = (int*)malloc(NbTransmissibility*(sizeof(int)));
    // Tij           = (double*)malloc(NbTransmissibility*(sizeof(double)));
    // TConductionIJ = (double*)malloc(NbTransmissibility*(sizeof(double)));
    iTr.resize(NbTransmissibility);
    jTr.resize(NbTransmissibility);
    Tij.resize(NbTransmissibility);
    TConductionIJ.resize(NbTransmissibility);

    std::size_t k=0;
    for ( std::size_t i=0; i<NbConnections; i++ )
    {
      if ( ConType[i]==1 || ConType[i]==2 ) // M-M, M-F ////////////////////////////
      {
        iTr[k] = ConCV[i][0];
        jTr[k] = ConCV[i][1];

        if ( ConTr[i][0] + ConTr[i][1] != 0.0 )
        {
          Tij[k] = ( ConTr[i][0]*ConTr[i][1] ) / ( ConTr[i][0] + ConTr[i][1] );
        }

        if ( ConGeom[i][0] + ConGeom[i][1] != 0.0 )
          TConductionIJ[k] = ( ConGeom[i][0] * ConGeom[i][1] ) /
                             ( ConGeom[i][0] + ConGeom[i][1] );

        if ( TConductionIJ[k] < 0.0 )
        {
          cout << "M-M or M-F : Conduction is negative: check values" << endl;
          exit ( 0 );
        }
        k++;
      }
      else if ( ConType[i] == 3 ) // F-F ///
      {
        SumTr = 0;
        double SumTr2 = 0;
        for ( std::size_t j=0; j<ConN[i]; j++ )
          SumTr+=ConTr[i][j];

        for ( std::size_t j=0; j<ConN[i]; j++ )
          SumTr2+=ConGeom[i][j];

        for ( std::size_t j=0; j<ConN[i]-1; j++ )
          for ( std::size_t n=j+1; n<ConN[i]; n++ )
          {
            iTr[k] = ConCV[i][j];
            jTr[k] = ConCV[i][n];
            Tij[k] = ( ConTr[i][j]*ConTr[i][n] ) / SumTr;

            TConductionIJ[k] = ( ConGeom[i][j]*ConGeom[i][n] ) / SumTr2;
            if ( TConductionIJ[k] < 0.0 )
            {
              cout << "M-M or M-F : Conduction is negative: check values" << endl;
              exit ( 0 );
            }
            k++;
          }
      }
      else
      {
        cout << "Wrong connection type" << endl;
        exit ( 0 );
      }

    }

}

/********************************************************************/
void CalcTranses::createKarimiApproximation()
{
    // pSim->outstream = "adgprs";
    // FILE  *poutfile;
    // int   i,j,k,n,m;

    double m1x,m1y,m1z,p1x,p1y,p1z,m2x,m2y,m2z,p2x,p2y,p2z,ix,iy,iz,inx,iny,inz;
    double mx,my,mz,ux,uy,uz,vx,vy,vz,nx,ny,nz,xi,yi,zi,nl,h,fx,fy,fz,fl;
    int   iswap;
    double    areatmp,volumetmp;
    double    hx,hy,hz;
    int   NbActivePolygon;
    int   NbFeatureCode;

    // printf("/*************************************************/\n");
    // printf("/* DISCRETE FEATURE MODEL                        */\n");
    // printf("/* by Mohammad Karimi-Fard (karimi@stanford.edu) */\n");
    // printf("/* March 2007, Stanford, CA.                     */\n");
    // printf("/*************************************************/\n");
    // printf("\n");

    // X = (double*)malloc(NbNodes*(sizeof(double)));
    // Y = (double*)malloc(NbNodes*(sizeof(double)));
    // Z = (double*)malloc(NbNodes*(sizeof(double)));
    X.resize(NbNodes);
    Y.resize(NbNodes);
    Z.resize(NbNodes);

    for (std::size_t i = 0; i < NbNodes; i++)
    {
      X[i] = vCoordinatesX[i];
      Y[i] = vCoordinatesY[i];
      Z[i] = vCoordinatesZ[i];
    }

    std::cout << 1 << std::endl;

    NbFNodes.resize(NbPolygons);
    EQF.resize(NbPolygons);
    CodePolygon.resize(NbPolygons);
    NbFNodes.resize(NbPolygons);
    FNodes.resize(NbPolygons);

    NbEdges = 0;
    NbCVs = 0;
    for (std::size_t i=0; i<NbPolygons; i++)
    {
      NbFNodes[i] = vNbFNodes[i];
      FNodes[i].resize(NbFNodes[i]);

      for (std::size_t j=0; j<NbFNodes[i]; j++)
        FNodes[i][j] = vvFNodes[i][j];

      CodePolygon[i] = vCodePolygon[i];
      if (CodePolygon[i] >= 0)
      {
        EQF[i] = NbCVs++;
        NbEdges += NbFNodes[i];
      }
      else EQF[i] = -1;
    }

    std::cout << 2 << std::endl;

    NbActivePolygon = NbCVs;

    NbVFaces.resize(NbPolyhedra);
    VFaces.resize(NbPolyhedra);
    EQV.resize(NbPolyhedra);
    CodePolyhedron.resize(NbPolyhedra);

    for (std::size_t i=0; i<NbPolyhedra; i++)
    {
      NbVFaces[i] = vNbVFaces[i];
      VFaces[i].resize(vNbVFaces[i]);
      for (std::size_t j=0; j<NbVFaces[i]; j++)
        VFaces[i][j] = vvVFaces[i][j];

      CodePolyhedron[i] = vCodePolyhedron[i];
      if (CodePolyhedron[i] >= 0)
      {
        EQV[i] = NbCVs++;
      }
      else EQV[i] = -1;
    }

    std::cout << 3 << std::endl;

    ZoneCode.resize(NbZones);
    ZVolumeFactor.resize(NbZones);
    ZPorosity.resize(NbZones);
    ZPermCode.resize(NbZones);

    ZPermeability.resize(NbZones);
    ZConduction.resize(NbZones);

    for (std::size_t i=0; i<NbZones; i++)
    {
      ZoneCode[i] = vZoneCode[i];
      std::size_t j = ZoneCode[i];
      ZVolumeFactor[j] = vZVolumeFactor[j];

      ZPorosity[j] = vZPorosity[j];
      ZPermCode[j] = vZPermCode[j];

      ZPermeability[j].resize(6);
      ZConduction[j].resize(6);

      ZPermeability[j][0] = vZPermeability[j * 3 + 0];
      ZPermeability[j][1] = vZPermeability[j * 3 + 1];
      ZPermeability[j][2] = vZPermeability[j * 3 + 2];

      ZPermeability[j][3] = 0.0;
      ZPermeability[j][4] = 0.0;
      ZPermeability[j][5] = 0.0;

      ZConduction[j][0] = vZConduction[j * 3 + 0];
      ZConduction[j][1] = vZConduction[j * 3 + 1];
      ZConduction[j][2] = vZConduction[j * 3 + 2];

      ZConduction[j][3] = 0.0;
      ZConduction[j][4] = 0.0;
      ZConduction[j][5] = 0.0;
   }
    std::cout << 4 << std::endl;

    ComputeBasicGeometry();
    std::cout << 41 << std::endl;
    ComputeControlVolumeList();
    std::cout << 42 << std::endl;
    PrepareConnectionList();
    std::cout << 43 << std::endl;
    ConstructConnectionList();

    std::cout << 5 << std::endl;

    if (NbOptions == 1) VolumeCorrection();

    // printf("NbCVs = %d\n",NbCVs);
    // printf("NbTransmissibility = %d\n",NbTransmissibility);

    ComputeContinuityNode();
    ComputeDirectionalPermeability();
    ComputeTransmissibilityPart();
    ComputeTransmissibilityList();

    std::cout << 6 << std::endl;

//////////////////////////////////
///// Computing Total Volume /////
//////////////////////////////////

    TotalVolume = 0;
    for (std::size_t i=0; i<NbCVs; i++)
      TotalVolume += CVVolume[i];

    // printf("TotalVolume=%e\n",TotalVolume);
    std::cout << 7 << std::endl;


/* OUTPUT MAPPING FOR GEOMECHANICS  */

// ASSUMING FIRST THE FEATURE CODES ARE NUMBERED
  // printf("Prepare reservoir-geomechanical mapping...\n");
  NbFeatureCode = 0;
  for ( std::size_t i = 0; i < NbPolygons; i++ )
  {
    if ( CodePolygon[i] > NbFeatureCode )
    {
      NbFeatureCode++;
    }
  }

  // Features
  vector<int> vPriority; vPriority.resize(NbPolygons);
  int kk = 1;
  for (std::size_t i = 0; i < NbPolygons; i++)
    if(CodePolygon[i] > -1)
    {
      vPriority[i] = kk;
      kk++;
    }
}


void CalcTranses::writeOutputFiles(const std::string & output_path) const
{
  /* Creates files:
   * fl_vol.txt
   * fl_poro.txt
   * fl_depth.txt
   * fl_tran.txt
   * fl_tran_n.txt
   * fl_tran_u.txt
   */

  stringstream out;
  out << "model/";
  out << "/";
  string outputPath_ = output_path;// out.str();
 //////////////////////////
///// OUTPUT Volumes /////
//////////////////////////
  int i;
  printf ( "Output Volumes...\n" );
  string outfile = outputPath_ + "fl_vol.txt";
  FILE * poutfile;
  poutfile=fopen ( outfile.c_str(),"w");
  //fprintf(out,"%d\n",NbCVs);
  fprintf ( poutfile,"%s\n","VOLUME" );
  for ( i=0; i<NbCVs; i++ )
    fprintf ( poutfile,"%e\n", CVVolume[i] );
  fprintf(poutfile,"%s\n","/\n");
  fclose(poutfile);
  ///////////////////////////
  ///// OUTPUT Porosity /////
  ///////////////////////////

    printf("Output Porosity...\n");
    outfile = outputPath_ + "fl_poro.txt";
    poutfile=fopen( outfile.c_str(),"w");
    fprintf(poutfile,"%s\n","PORO");
    for (i=0;i<NbCVs;i++) fprintf(poutfile,"%e\n",ZPorosity[CVZone[i]]);

    fprintf(poutfile,"%s\n","/\n");
    fclose(poutfile);

  ///////////////////////////
  ///// OUTPUT Depth  /////
  ///////////////////////////
    printf("Output Depth (z)...\n");
    outfile = outputPath_ + "fl_depth.txt";
    poutfile=fopen( outfile.c_str(),"w");
    //fprintf(out,"%d\n",NbCVs);
    fprintf(poutfile,"%s\n","DEPTH");
    for (i=0;i<NbCVs;i++)
      fprintf(poutfile,"%e\n",-CVz[i]);

    fprintf(poutfile,"%s\n","/\n");
    fclose(poutfile);

    /* OUTPUT Transmissibility */
    printf("Output Transmissibility...\n");
    outfile = outputPath_ + "fl_tran.txt";
    poutfile=fopen( outfile.c_str(),"w");
    fprintf(poutfile,"%s\n","TPFACONNS");
    fprintf(poutfile,"%d\n",NbTransmissibility);

    for (i=0;i<NbTransmissibility;i++)
      fprintf(poutfile,
              "%d\t%d\t%e\n",
              iTr[i],jTr[i],Tij[i] * transmissibility_conversion_factor);

    fprintf(poutfile,"%s\n","/\n");
    fclose(poutfile);

    printf("Output Transmissibility and Geometrical part...\n");
    outfile = outputPath_ + "fl_tran_n.txt";
    poutfile=fopen( outfile.c_str(),"w");
    fprintf(poutfile,"%s\n","TPFACONNSN");
    fprintf(poutfile,"%d\n",NbTransmissibility);

    for (i=0;i<NbTransmissibility;i++)
      fprintf(poutfile,
              "%d\t%d\t%e\t%e\n",
              iTr[i],jTr[i],Tij[i] * transmissibility_conversion_factor,
              TConductionIJ[i]);

    fprintf(poutfile,"%s\n","/\n");
    fclose(poutfile);

    string outstring =  outputPath_ + "fl_tran_u.txt";
    poutfile = fopen(outstring.c_str(),"w");
    fprintf(poutfile,"GMUPDATETRANS\n");

    int k=0;
    for(i=0;i<NbConnections;i++)
    {
      if(ConType[i]==1 || ConType[i]==2)  // M-M, M-F /
      {
      if(ConType[i]==1)   // M-M ////////////////////////////
      {
        // Con# ConType i ai j aj -> Tij=ai*aj/(ai+aj)
        fprintf(poutfile,"%d\t%d\t%d\t%e\t%d\t%e\n", k,
                ConType[i], ConCV[i][0], ConTr[i][0], ConCV[i][1], ConTr[i][1]);
      }
      if(ConType[i]==2)   // M-F ////////////////////////////
      {
        //Con# ConType m am i ci ei ki ai=ci*ki/ei-> Tmi=am*ai/(am+ai)
        fprintf(poutfile,"%d\t%d\t%d\t%e\t%d\t%e\t%e\t%e\n", k,
                ConType[i], ConCV[i][0], ConTr[i][0], ConCV[i][1],
                2.*ConArea[i][1], ZVolumeFactor[CVZone[ConCV[i][1]]], ConPerm[i][1]);
      }
      k++;
      }
      if(ConType[i]==3) // F-F /////////////////////////////////////////////////
      {
        //Con# ConType i ci ei ki j cj ej kj N n cn en kn
        //ai=ci*ki*ei aj=cj*kj*ej an=cn*kn*en
        //-> Tij=ai*aj/(SUM an)
        //
        for(int j=0;j<ConN[i]-1;j++)
          for(int n=j+1;n<ConN[i];n++)
          {
            fprintf(poutfile,"%d\t%d\t%d\t%e\t%e\t%e\t%d\t%e\t%e\t%e\t",
                    k,
                    ConType[i],

                    ConCV[i][j],
                    ConTr[i][j]/(ConPerm[i][j]*ZVolumeFactor[CVZone[ConCV[i][j]]]),
                    ZVolumeFactor[CVZone[ConCV[i][j]]],
                    ConPerm[i][j],

                    ConCV[i][n],
                    ConTr[i][n]/(ConPerm[i][n]*ZVolumeFactor[CVZone[ConCV[i][n]]]),
                    ZVolumeFactor[CVZone[ConCV[i][n]]],
                    ConPerm[i][n]);

            fprintf(poutfile,"%d\t",ConN[i]);
            for(int m=0;m<ConN[i];m++)
              fprintf(poutfile,"%d\t%e\t%e\t%e\t",
                      ConCV[i][m],
                      ConTr[i][m]/(ConPerm[i][m]*ZVolumeFactor[CVZone[ConCV[i][m]]]),
                      ZVolumeFactor[CVZone[ConCV[i][m]]],
                      ConPerm[i][m]);

            fprintf(poutfile,"\n");
            k++;
          }
    }
    }
    fprintf(poutfile,"/\n");
    fclose(poutfile);
}


void CalcTranses::extractData(FlowData & data) const
{
  // std::cout << "extracting volume data" << std::endl;
  // Extract Volumes, porosity, depth
  data.volumes.resize(NbCVs);
  data.poro.resize(NbCVs);
  data.depth.resize(NbCVs);
  for (std::size_t i=0; i<NbCVs; i++ )
  {
    data.volumes[i] = CVVolume[i];
    data.poro[i]    = ZPorosity[CVZone[i]];
    data.depth[i]   = -CVz[i];
  }

  // Transmissibility
  // std::cout << "extracting transes" << std::endl;
  data.ielement.resize(NbTransmissibility);
  data.jelement.resize(NbTransmissibility);
  data.trans_ij.resize(NbTransmissibility);
  data.conduct_ij.resize(NbTransmissibility);
  for (std::size_t i=0;i<NbTransmissibility; i++)
  {
    data.ielement[i]   = iTr[i];
    data.jelement[i]   = jTr[i];
    data.trans_ij[i]   = Tij[i];
    data.conduct_ij[i] = TConductionIJ[i];
  }

  // std::cout << "extracting geomechanics-related stuff" << std::endl;
  // Geomechanics
  data.connection_type.resize(NbConnections);
  for(std::size_t i=0;i<NbConnections;i++)
  {
    data.connection_type[i] = ConType[i];
    // data.connection_n[i] =
  }


}
