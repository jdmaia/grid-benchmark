#include "Checks.hpp"

#include <Grid/Grid.h>
#include "Benchmark-Grid/Roctx.hpp"

using namespace Grid;


void Benchmark::checkWilson(void)
{
  GRID_ROCTX_RANGE("checkWilson");
  int L=8;
  std::cout << GridLogMessage << "Validating Wilson propagator" << std::endl;

  Grid::Coordinate mpi = GridDefaultMpi();
  assert(mpi.size() == 4);
  Coordinate local({L, L, L, L});
  Coordinate latt4(
      {local[0] * mpi[0], local[1] * mpi[1], local[2] * mpi[2], local[3] * mpi[3]});

  ///////// Lattice Init ////////////
  GridCartesian *FGrid = SpaceTimeGrid::makeFourDimGrid(
      latt4, GridDefaultSimd(Nd, vComplexD::Nsimd()), GridDefaultMpi());
  GridRedBlackCartesian *FrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(FGrid);

  ///////// RNG Init ////////////
  std::vector<int> seeds({1, 2, 3, 4});
  GridParallelRNG pRNG(FGrid);
  pRNG.SeedFixedIntegers(seeds);
  GridSerialRNG sRNG;
  sRNG.SeedFixedIntegers(seeds);
  std::cout << GridLogMessage << "Initialised RNGs" << std::endl;

  ///////// Init Fields /////////
  LatticeGaugeFieldD Umu(FGrid);
  SU<Nc>::ColdConfiguration(pRNG,Umu); // Unit gauge

  LatticeFermionD    src(FGrid);
  gaussian(pRNG,src);
  LatticeFermionD    tmp(FGrid);
  LatticeFermionD    ref(FGrid);

  Coordinate point(4,0);
  src=Zero();
  SpinColourVectorD ferm;
  gaussian(sRNG,ferm);
  pokeSite(ferm,src,point);

  RealD mass = 0.01;
  WilsonFermionD Dw(Umu,*FGrid,*FrbGrid,mass);

  /////////////////////////////////////
  // Momentum space propagator + FFT //
  /////////////////////////////////////
  std::cout << GridLogMessage <<  " Solving by FFT and Feynman rules" <<std::endl;
  Dw.FreePropagator(src,ref,mass);

  LatticeFermionD    result(FGrid);

  ////////////////////////////////////////////////////////////////////////
  // Conjugate gradient on normal equations system
  ////////////////////////////////////////////////////////////////////////
  std::cout << GridLogMessage << " Solving by Conjugate Gradient (CGNE)" <<std::endl;
  Dw.Mdag(src,tmp);
  src=tmp;
  MdagMLinearOperator<WilsonFermionD,LatticeFermionD> HermOp(Dw);
  ConjugateGradient<LatticeFermionD> CG(1.0e-8,10000);
  {
    GRID_ROCTX_RANGE("CG solve");
    CG(HermOp,src,result);
  }

  std::cout << GridLogMessage << " Taking difference" <<std::endl;
  std::cout << GridLogMessage << "Wilson result "<<norm2(result)<<std::endl;
  std::cout << GridLogMessage << "Wilson ref    "<<norm2(ref)<<std::endl;

  RealD diff = norm2(ref - result);
  RealD sum  = norm2(ref + result);
  RealD normratio = diff/sum;
  std::cout << GridLogMessage << "||result - ref||                      "<< diff <<std::endl;
  std::cout << GridLogMessage << "||result - ref||/||result + ref||     "<< normratio <<std::endl;
  if (normratio >= 1e-8)
  {
    std::cout << GridLogError << "Failed to validate free Wilson propagator: ||(result - ref)||/||(result + ref)|| >= 1e-8" << std::endl;
    exit(EXIT_FAILURE);
  }
}
