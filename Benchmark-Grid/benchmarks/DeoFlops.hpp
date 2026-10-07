#pragma once

#include <Grid/Grid.h>
#include "instantiation/instantiations.hpp"
#include "Common.hpp"
#include "Benchmark-Grid/ActionInfo.hpp"
#include "Benchmark-Grid/Utils.hpp"
#include "Benchmark-Grid/Roctx.hpp"


using namespace Grid;

  
struct controls
{
  int Opt;
  int CommsOverlap;
  Grid::CartesianCommunicator::CommunicatorPolicy_t CommsAsynch;
};

/* Reimplement projection to higher representation because Grid's built-in one is hardcoded to configure-time Nc */
template<typename Action>
void update_representation(
  typename ActionFactory<Action>::RepresentedGaugeField &Uas,
  const typename ActionFactory<Action>::FundamentalGaugeField &Uin)
{
  static_assert(
    !std::is_same_v<Action,Action>,
    "update_representation not implemented for provided Action"
  );
}

typedef ActionFactory<Grid::DomainWallFermionSp4TwoIndexAntiSymmetricD> DWFSp4TwoIndexASInfo;

template<>
void update_representation<DWFSp4TwoIndexASInfo::Action>(
  typename DWFSp4TwoIndexASInfo::RepresentedGaugeField &Uas,
  const typename DWFSp4TwoIndexASInfo::FundamentalGaugeField &Uin);


namespace Benchmark
{
  template<typename Action>
  static double DeoFlops(int Ls, int L, Grid::Coordinate pattern)
  {
    std::string rangeName = std::string("DeoFlops L=") + std::to_string(L);
    GRID_ROCTX_RANGE(rangeName.c_str());
    double gflops;
    double gflops_best = 0;
    double gflops_worst = 0;
    std::vector<double> gflops_all;

    ///////////////////////////////////////////////////////
    // Set/Get the layout & grid size
    ///////////////////////////////////////////////////////
    int threads = GridThread::GetThreads();
    Coordinate mpi = GridDefaultMpi();
    assert(mpi.size() == 4);
    Coordinate local({L*pattern[0], L*pattern[1], L*pattern[2], L*pattern[3]});
    Coordinate latt4(
        {local[0] * mpi[0], local[1] * mpi[1], local[2] * mpi[2], local[3] * mpi[3]});

    GridCartesian *TmpGrid = SpaceTimeGrid::makeFourDimGrid(
        latt4, GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
    uint64_t NP = TmpGrid->RankCount();
    uint64_t NN = TmpGrid->NodeCount();
    uint64_t SHM = NP / NN;
    delete TmpGrid;

    ///////// Welcome message ////////////
    grid_big_sep();
    typedef ActionFactory<Action> ActionInfo;
    typedef typename ActionInfo::Representation Representation;
    typedef typename ActionInfo::GroupInfo GroupInfo;
    typedef typename GroupInfo::group Group;
    constexpr const char* groupFamily = GroupInfo::groupName;
    constexpr int Nc = ActionInfo::Nc;

    std::cout << GridLogMessage << "Benchmark " << groupFamily << "(" << Nc << ") "
              << ActionInfo::name() << " on " << L << "^4 local volume "
              << std::endl;
    if (pattern[0]>1 || pattern[1]>1 || pattern[2]>1 || pattern[3]>1)
    {
      std::cout << GridLogMessage << "Using pattern "
                << pattern[0] << "."
                << pattern[1] << "."
                << pattern[2] << "."
                << pattern[3] << ": local volume scaled to "
                << local[0] << "."
                << local[1] << "."
                << local[2] << "."
                << local[3] << std::endl;
    }
    std::cout << GridLogMessage << "* Group family   : " << groupFamily << std::endl;
    std::cout << GridLogMessage << "* Nc             : " << Nc << std::endl;
    std::cout << GridLogMessage << "* Representation : " << GroupInfo::getRepresentationName()
              << std::endl;
    std::cout << GridLogMessage
              << "* Global volume  : " << GridCmdVectorIntToString(latt4) << std::endl;
    if (Ls > 0) std::cout << GridLogMessage << "* Ls             : " << Ls << std::endl;
    std::cout << GridLogMessage << "* ranks          : " << NP << std::endl;
    std::cout << GridLogMessage << "* nodes          : " << NN << std::endl;
    std::cout << GridLogMessage << "* ranks/node     : " << SHM << std::endl;
    std::cout << GridLogMessage << "* ranks geom     : " << GridCmdVectorIntToString(mpi)
              << std::endl;
    std::cout << GridLogMessage << "* Using " << threads << " threads" << std::endl;
    grid_big_sep();

    ///////// Lattice Init ////////////
    GRID_ROCTX_PUSH("Setup");
    GridCartesian *UGrid = SpaceTimeGrid::makeFourDimGrid(
        latt4, GridDefaultSimd(Nd, Action::Simd::Nsimd()), GridDefaultMpi());
    GridRedBlackCartesian *UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);
    GridCartesian *FGrid;
    GridRedBlackCartesian *FrbGrid;

    if (Ls > 0)
    {
      FGrid   = SpaceTimeGrid::makeFiveDimGrid(Ls, UGrid);
      FrbGrid = SpaceTimeGrid::makeFiveDimRedBlackGrid(Ls, UGrid);
    }
    else
    {
      FGrid   = UGrid;
      FrbGrid = UrbGrid;
    }

    ///////// RNG Init ////////////
    std::vector<int> seeds4({1, 2, 3, 4});
    std::vector<int> seeds5({5, 6, 7, 8});
    GridParallelRNG RNG4(UGrid);
    RNG4.SeedFixedIntegers(seeds4);
    GridParallelRNG RNG5(FGrid);
    GridParallelRNG* FRNG;
    if (Ls > 0)
    {
      FRNG = &RNG4;
    }
    else
    {
      RNG5.SeedFixedIntegers(seeds5);
      FRNG = &RNG5;
    }
    std::cout << GridLogMessage << "Initialised RNGs" << std::endl;

    typedef typename Action::FermionField Fermion;
    typedef typename ActionInfo::RepresentedGaugeField RepresentedGaugeField;
    typedef typename ActionInfo::FundamentalGaugeField FundamentalGaugeField;

    ///////// Source preparation ////////////

    // The reference Grid commit does not support Group::HotConfiguration
    // for FP32 Sp(4).
    // We'll get around this by generating the gauge in double-precision
    // and casting to single-precision.
    RepresentedGaugeField Umu(UGrid);
    if constexpr (std::is_same_v<typename GroupInfo::label, GroupName::Sp> &&
                  std::is_same_v<typename Action::Simd::Real, float>)
    {
      // Set up Grid + RNG for double-precision.
      // Assume DWF for simplicity: this will just cause a compilation
      // failure otherwise.
      typedef DomainWallFermion<WilsonImpl<Grid::vComplexD,
                   Representation,
                   Grid::CoeffReal>> ActionD;
      typedef ActionFactory<ActionD> ActionInfoD;

      GridCartesian *UGridD = SpaceTimeGrid::makeFourDimGrid(
          latt4, GridDefaultSimd(Nd, ActionD::Simd::Nsimd()), GridDefaultMpi());
      GridParallelRNG RNG4D(UGridD);
      RNG4D.SeedFixedIntegers(seeds4);

      // Gauge field generation block -- same as top-level else
      // statement.
      typename ActionInfoD::RepresentedGaugeField Umu_fp64(UGridD);
      if constexpr (Representation::isFundamental)
      {
        Group::HotConfiguration(RNG4D, Umu_fp64);
      }
      else
      {
        typename ActionInfoD::FundamentalGaugeField Umu_fund_fp64(UGridD);
        Group::HotConfiguration(RNG4D, Umu_fund_fp64);
        update_representation<ActionD>(Umu_fp64, Umu_fund_fp64);
      }

      // Convert FP64 -> FP32
      precisionChange(Umu, Umu_fp64);
      delete UGridD;
    }
    else
    {
      if constexpr (Representation::isFundamental)
      {
        Group::HotConfiguration(RNG4, Umu);
      }
      else
      {
        FundamentalGaugeField Umu_fund(UGrid);
        Group::HotConfiguration(RNG4, Umu_fund);
        update_representation<Action>(Umu, Umu_fund);
      }
    }

    Fermion src(FGrid);
    random(*FRNG, src);
    Fermion src_e(FrbGrid);
    Fermion src_o(FrbGrid);
    Fermion r_e(FrbGrid);
    Fermion r_o(FrbGrid);
    Fermion r_eo(FGrid);
    Action action = ActionInfo::create(Umu, UGrid, UrbGrid, FGrid, FrbGrid);
    GRID_ROCTX_POP(); // Setup

    {

      pickCheckerboard(Even, src_e, src);
      pickCheckerboard(Odd, src_o, src);

      const int num_cases = 4;
      std::string fmt("G/S/C ; G/O/C ; G/S/S ; G/O/S ");

      controls Cases[] = {
          {Action::Kernels::OptGeneric, Action::Kernels::CommsThenCompute,
           CartesianCommunicator::CommunicatorPolicyConcurrent},
          {Action::Kernels::OptGeneric, Action::Kernels::CommsAndCompute,
           CartesianCommunicator::CommunicatorPolicyConcurrent},
          {Action::Kernels::OptGeneric, Action::Kernels::CommsThenCompute,
           CartesianCommunicator::CommunicatorPolicySequential},
          {Action::Kernels::OptGeneric, Action::Kernels::CommsAndCompute,
           CartesianCommunicator::CommunicatorPolicySequential}};

      for (int c = 0; c < num_cases; c++)
      {
        const char* commsDesc =
            (Cases[c].CommsOverlap == Action::Kernels::CommsAndCompute)
            ? "overlapped-comms-compute"
            : "sequential-comms-then-compute";
        const char* policyDesc =
            (Cases[c].CommsAsynch == CartesianCommunicator::CommunicatorPolicyConcurrent)
            ? "concurrent-comms-policy"
            : "sequential-comms-policy";
        std::string caseLabel = std::string("Case ") + std::to_string(c) + ": " +
                                 commsDesc + ", " + policyDesc;
        GRID_ROCTX_RANGE(caseLabel.c_str());

        Action::Kernels::Comms = Cases[c].CommsOverlap;
        Action::Kernels::Opt = Cases[c].Opt;
        CartesianCommunicator::SetCommunicatorPolicy(Cases[c].CommsAsynch);

        grid_small_sep();
        if (Action::Kernels::Opt == Action::Kernels::OptGeneric)
          std::cout << GridLogMessage << "* Using GENERIC Nc " << getClassName<typename Action::Kernels>() << std::endl;
        if (Action::Kernels::Comms == Action::Kernels::CommsAndCompute)
          std::cout << GridLogMessage << "* Using Overlapped Comms/Compute" << std::endl;
        if (Action::Kernels::Comms == Action::Kernels::CommsThenCompute)
          std::cout << GridLogMessage << "* Using sequential Comms/Compute" << std::endl;
        std::cout << GridLogMessage << "* " << actionPrec<Action>() << " precision " << std::endl;
        grid_small_sep();

        int nwarm = 10;
        double t0 = usecond();
        action.FermionGrid()->Barrier();
        {
          GRID_ROCTX_RANGE("Warmup");
          for (int i = 0; i < nwarm; i++)
          {
            action.DhopEO(src_o, r_e, DaggerNo);
          }
        }
        action.FermionGrid()->Barrier();
        double t1 = usecond();
        uint64_t ncall = 500;

        action.FermionGrid()->Broadcast(0, &ncall, sizeof(ncall));

        time_statistics timestat;
        std::vector<double> t_time(ncall);
        {
          GRID_ROCTX_RANGE("DhopEO timed loop");
          for (uint64_t i = 0; i < ncall; i++)
          {
            t0 = usecond();
            action.DhopEO(src_o, r_e, DaggerNo);
            t1 = usecond();
            t_time[i] = t1 - t0;
          }
        }
        action.FermionGrid()->Barrier();

        double volume = Ls > 0? Ls : 1;
        for (int mu = 0; mu < Nd; mu++)
          volume = volume * latt4[mu];

        double fps   = ActionInfo::fps();
        double flops = (fps * volume) / 2.;
        double gf_hi, gf_lo, gf_err;

        timestat.statistics(t_time);
        gf_hi = flops / timestat.min / 1000.;
        gf_lo = flops / timestat.max / 1000.;
        gf_err = flops / timestat.min * timestat.err / timestat.mean / 1000.;

        gflops = flops / timestat.mean / 1000.;
        gflops_all.push_back(gflops);
        if (gflops_best == 0)
          gflops_best = gflops;
        if (gflops_worst == 0)
          gflops_worst = gflops;
        if (gflops > gflops_best)
          gflops_best = gflops;
        if (gflops < gflops_worst)
          gflops_worst = gflops;

        std::cout << GridLogMessage << "Deo FlopsPerSite is " << fps << std::endl;
        std::cout << GridLogMessage << std::fixed << std::setprecision(1)
                  << "Deo Gflop/s =   " << gflops << " (" << gf_err << ") " << gf_lo
                  << "-" << gf_hi << std::endl;
        std::cout << GridLogMessage << std::fixed << std::setprecision(1)
                  << "Deo Gflop/s per rank   " << gflops / NP << std::endl;
        std::cout << GridLogMessage << std::fixed << std::setprecision(1)
                  << "Deo Gflop/s per node   " << gflops / NN << std::endl;
      }

      grid_small_sep();
      std::cout << GridLogMessage << L << "^4"
                << " Deo Best  Gflop/s        =   " << gflops_best << " ; "
                << gflops_best / NN << " per node " << std::endl;
      std::cout << GridLogMessage << L << "^4"
                << " Deo Worst Gflop/s        =   " << gflops_worst << " ; "
                << gflops_worst / NN << " per node " << std::endl;
      std::cout << GridLogMessage << fmt << std::endl;
      std::cout << GridLogMessage;

      for (int i = 0; i < gflops_all.size(); i++)
      {
        std::cout << gflops_all[i] / NN << " ; ";
      }
      std::cout << std::endl;
    }

    // Clean up Grids -- prevents too many MPI Communicators existing + crashing
    delete UGrid;
    delete UrbGrid;
    if (Ls > 0)
    {
      delete FGrid;
      delete FrbGrid;
    }

    return gflops_best;
  }
}
