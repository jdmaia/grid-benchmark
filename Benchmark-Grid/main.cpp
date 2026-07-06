/*
 * Copyright © 2015 Peter Boyle <paboyle@ph.ed.ac.uk>
 * Copyright © 2024 Simon Buerger <simon.buerger@rwth-aachen.de>
 * Copyright © 2022-2025 Antonin Portelli <antonin.portelli@me.com>
 *
 * This is a fork of Benchmark_ITT.cpp from Grid
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "Common.hpp"
#include "json.hpp"
#include <Grid/Grid.h>
#include "instantiation/instantiations.hpp"
#include "Benchmark-Grid/benchmarks/Benchmarks.hpp"
#include "Benchmark-Grid/Utils.hpp"
#include "Benchmark-Grid/checks/Checks.hpp"


using namespace Grid;


std::vector<std::string> get_mpi_hostnames()
{
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  char hostname[MPI_MAX_PROCESSOR_NAME];
  int name_len = 0;
  MPI_Get_processor_name(hostname, &name_len);

  // Allocate buffer to gather all hostnames
  std::vector<char> all_hostnames(world_size * MPI_MAX_PROCESSOR_NAME);

  // Use MPI_Allgather to gather all hostnames on all ranks
  MPI_Allgather(hostname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, all_hostnames.data(),
                MPI_MAX_PROCESSOR_NAME, MPI_CHAR, MPI_COMM_WORLD);

  // Convert the gathered hostnames back into a vector of std::string
  std::vector<std::string> hostname_list(world_size);
  for (int i = 0; i < world_size; ++i)
  {
    hostname_list[i] = std::string(&all_hostnames[i * MPI_MAX_PROCESSOR_NAME]);
  }

  return hostname_list;
}


namespace Benchmark
{
  static void Decomposition(int& NNout, Grid::Coordinate& pattern, nlohmann::json& json_results)
  {
    nlohmann::json tmp;
    int threads = GridThread::GetThreads();
    Grid::Coordinate mpi = GridDefaultMpi();
    assert(mpi.size() == 4);
    Coordinate local({8, 8, 8, 8});
    Coordinate latt4(
        {local[0] * mpi[0], local[1] * mpi[1], local[2] * mpi[2], local[3] * mpi[3]});
    GridCartesian *TmpGrid = SpaceTimeGrid::makeFourDimGrid(
        latt4, GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
    Grid::Coordinate shm(4, 1);
    GlobalSharedMemory::GetShmDims(mpi, shm);

    uint64_t NP = TmpGrid->RankCount();
    uint64_t NN = TmpGrid->NodeCount();
    NNout = NN;
    uint64_t SHM = NP / NN;
    delete TmpGrid;

    grid_big_sep();
    std::cout << GridLogMessage << "Grid Default Decomposition patterns\n";
    grid_small_sep();
    std::cout << GridLogMessage << "* OpenMP threads : " << GridThread::GetThreads()
              << std::endl;

    std::cout << GridLogMessage << "* MPI layout     : " << GridCmdVectorIntToString(mpi)
              << std::endl;
    std::cout << GridLogMessage << "* Shm layout     : " << GridCmdVectorIntToString(shm)
              << std::endl;

    std::cout << GridLogMessage << "* vReal          : " << sizeof(vReal) * 8 << "bits ; "
              << GridCmdVectorIntToString(GridDefaultSimd(4, vReal::Nsimd()))
              << std::endl;
    std::cout << GridLogMessage << "* vRealF         : " << sizeof(vRealF) * 8
              << "bits ; "
              << GridCmdVectorIntToString(GridDefaultSimd(4, vRealF::Nsimd()))
              << std::endl;
    std::cout << GridLogMessage << "* vRealD         : " << sizeof(vRealD) * 8
              << "bits ; "
              << GridCmdVectorIntToString(GridDefaultSimd(4, vRealD::Nsimd()))
              << std::endl;
    std::cout << GridLogMessage << "* vComplex       : " << sizeof(vComplex) * 8
              << "bits ; "
              << GridCmdVectorIntToString(GridDefaultSimd(4, vComplex::Nsimd()))
              << std::endl;
    std::cout << GridLogMessage << "* vComplexF      : " << sizeof(vComplexF) * 8
              << "bits ; "
              << GridCmdVectorIntToString(GridDefaultSimd(4, vComplexF::Nsimd()))
              << std::endl;
    std::cout << GridLogMessage << "* vComplexD      : " << sizeof(vComplexD) * 8
              << "bits ; "
              << GridCmdVectorIntToString(GridDefaultSimd(4, vComplexD::Nsimd()))
              << std::endl;
    std::cout << GridLogMessage << "* ranks          : " << NP << std::endl;
    std::cout << GridLogMessage << "* nodes          : " << NN << std::endl;
    std::cout << GridLogMessage << "* ranks/node     : " << SHM << std::endl;

    for (unsigned int i = 0; i < mpi.size(); ++i)
    {
      tmp["mpi"].push_back(mpi[i]);
      tmp["shm"].push_back(shm[i]);
    }
    tmp["ranks"] = NP;
    tmp["nodes"] = NN;
    tmp["pattern"] = pattern;
    json_results["geometry"] = tmp;
  }
};

void printUsage()
{
  std::cout<<GridLogMessage<<"Usage: Benchmark_Grid <options>"<<std::endl;
  std::cout<<GridLogMessage<<"If conflicting or repeated options are given, the final option takes precedence."<<std::endl;
  std::cout<<GridLogMessage<<"Options:"<<std::endl;
  std::cout<<GridLogMessage<<"  --help                       : This message"<<std::endl;
  std::cout<<GridLogMessage<<"  --json-out <path>            : Export results to a JSON at <path>."<<std::endl;
  std::cout<<GridLogMessage<<"  --benchmark-memory           : Enable axpy memory benchmark (default=on)."<<std::endl;
  std::cout<<GridLogMessage<<"  --no-benchmark-memory        : Disable axpy memory benchmark."<<std::endl;
  std::cout<<GridLogMessage<<"  --benchmark-su4              : Enable SU(4) memory benchmark (default=on)."<<std::endl;
  std::cout<<GridLogMessage<<"  --no-benchmark-su4           : Disable SU(4) memory benchmark."<<std::endl;
  std::cout<<GridLogMessage<<"  --benchmark-comms            : Enable communications benchmark (default=on)."<<std::endl;
  std::cout<<GridLogMessage<<"  --no-benchmark-comms         : Disable communications benchmark."<<std::endl;
  std::cout<<GridLogMessage<<"  --benchmark-flops            : Enable all Dirac Matrix Flops benchmarks (default=on)."<<std::endl;
  std::cout<<GridLogMessage<<"  --no-benchmark-flops         : Disable all Dirac Matrix Flops benchmarks."<<std::endl;
  std::cout<<GridLogMessage<<"  --benchmark-flops-su4        : Enable SU(4) Dirac Matrix Flops benchmark (default=on)."<<std::endl;
  std::cout<<GridLogMessage<<"  --no-benchmark-flops-su4     : Disable SU(4) Dirac Matrix Flops benchmark."<<std::endl;
  std::cout<<GridLogMessage<<"  --benchmark-flops-sp4-f      : Enable Sp(4) Dirac Matrix Flops benchmark (default=on)."<<std::endl;
  std::cout<<GridLogMessage<<"  --no-benchmark-flops-sp4-f   : Disable Sp(4) Dirac Matrix Flops benchmark."<<std::endl;
  std::cout<<GridLogMessage<<"  --benchmark-flops-sp4-2as    : Enable Sp(4) Two-Index AntiSymmetric Dirac Matrix Flops benchmark (default=on)."<<std::endl;
  std::cout<<GridLogMessage<<"  --no-benchmark-flops-sp4-2as : Disable Sp(4) Two-Index AntiSymmetric Dirac Matrix Flops benchmark."<<std::endl;
  std::cout<<GridLogMessage<<"  --benchmark-flops-fp64       : Enable FP64 Dirac Matrix Flops benchmarks (default=on)."<<std::endl;
  std::cout<<GridLogMessage<<"  --no-benchmark-flops-fp64    : Disable FP64 Dirac Matrix Flops benchmarks."<<std::endl;
  std::cout<<GridLogMessage<<"  --benchmark-latency          : Enable point-to-point communications latency benchmark (default=off)."<<std::endl;
  std::cout<<GridLogMessage<<"  --no-benchmark-latency       : Disable point-to-point communications latency benchmark."<<std::endl;
  std::cout<<GridLogMessage<<"  --benchmark-p2p              : Enable point-to-point communications bandwidth benchmark (default=off)."<<std::endl;
  std::cout<<GridLogMessage<<"  --no-benchmark-p2p           : Disable point-to-point communications bandwidth benchmark."<<std::endl;
  std::cout<<GridLogMessage<<"  --check-wilson               : Enable Wilson Fermion correctness check (default=on)."<<std::endl;
  std::cout<<GridLogMessage<<"  --no-check-wilson            : Disable Wilson Fermion correctness check."<<std::endl;
  std::cout<<GridLogMessage<<"  --pattern <x.y.z.t>          : Scales the local lattice dimensions by the factors in the string x.y.z.t."<<std::endl;
  std::cout<<GridLogMessage<<"  --max-L                      : Sets the maximum lattice size for the flops benchmarks. This must decompose to 2^n 3^m, for n>0 and m={0,1}."<<std::endl;
  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<"See below for Grid usage."<<std::endl;
  std::cout<<GridLogMessage<<std::endl;
  std::cout<<GridLogMessage<<std::endl;
}

int main(int argc, char **argv)
{
  for (int i = 0; i < argc; i++)
  {
    auto arg = std::string(argv[i]);
    if (arg == "--help")
    {
      printUsage();
      break;
    }
  }
  Grid_init(&argc, &argv);

  bool do_su4 = true;
  bool do_memory = true;
  bool do_comms = true;
  bool do_flops = true;
  bool do_flops_su4  = true;
  bool do_flops_sp4_f = true;
  bool do_flops_sp4_2as = true;
  bool do_flops_fp64 = true;
  bool do_check_wilson = true;

  // NOTE: these two take O((number of ranks)^2) time, which might be a lot, so they are
  // off by default
  bool do_latency = false;
  bool do_p2p = false;

  Grid::Coordinate pattern({1,1,1,1});
  int maxL = 32;
  std::string json_filename = ""; // empty indicates no json output
  
  for (int i = 0; i < argc; i++)
  {
    auto arg = std::string(argv[i]);
    if (arg == "--json-out")
    {
      if ((i+1) < argc)
      {
        json_filename = argv[i + 1];
        ++i;
      }
      else
      {
        std::cerr << GridLogError << "--json-out provided without an output filepath." << std::endl;
        exit(EXIT_FAILURE);
      }
    }
    if (arg == "--benchmark-su4")
      do_su4 = true;
    if (arg == "--benchmark-memory")
      do_memory = true;
    if (arg == "--benchmark-comms")
      do_comms = true;
    if (arg == "--benchmark-flops")
      do_flops = true;
    if (arg == "--benchmark-latency")
      do_latency = true;
    if (arg == "--benchmark-p2p")
      do_p2p = true;
    if (arg == "--benchmark-flops-su4")
      do_flops_su4 = true;
    if (arg == "--benchmark-flops-sp4-f")
      do_flops_sp4_f = true;
    if (arg == "--benchmark-flops-sp4-2as")
      do_flops_sp4_2as = true;
    if (arg == "--benchmark-flops-fp64")
      do_flops_fp64 = true;
    if (arg == "--check-wilson")
      do_check_wilson = true;
    if (arg == "--no-benchmark-su4")
      do_su4 = false;
    if (arg == "--no-benchmark-memory")
      do_memory = false;
    if (arg == "--no-benchmark-comms")
      do_comms = false;
    if (arg == "--no-benchmark-flops")
      do_flops = false;
    if (arg == "--no-benchmark-latency")
      do_latency = false;
    if (arg == "--no-benchmark-p2p")
      do_p2p = false;
    if (arg == "--no-benchmark-flops-su4")
      do_flops_su4 = false;
    if (arg == "--no-benchmark-flops-sp4-f")
      do_flops_sp4_f = false;
    if (arg == "--no-benchmark-flops-sp4-2as")
      do_flops_sp4_2as = false;
    if (arg == "--no-benchmark-flops-fp64")
      do_flops_fp64 = false;
    if (arg == "--no-check-wilson")
      do_check_wilson = false;
    if (arg == "--pattern")
    {
      // Make sure there's another argument to parse
      if (i == (argc - 1))
      {
        std::cout << GridLogError << "No argument provided for --max-L" << std::endl;
        exit(EXIT_FAILURE);
      }
      ++i;
      GridCmdOptionIntVector(argv[i], pattern);
    }
    if (arg == "--max-L")
    {
      // Make sure there's another argument to parse
      if (i == (argc - 1))
      {
        std::cout << GridLogError << "No argument provided for --max-L" << std::endl;
        exit(EXIT_FAILURE);
      }

      // Check that the argument is convertible to an int
      ++i;
      std::string arg2 = std::string(argv[i]);
      try
      {
        maxL = stoi(arg2);
      }
      catch (const std::invalid_argument &e)
      {
        std::cerr << GridLogError << "Invalid argument for --max-L: " << arg2 << std::endl;
        exit(EXIT_FAILURE);
      }
      catch (const std::out_of_range &e)
      {
        std::cerr << GridLogError << "Argument for --max-L is out of range: " << arg2
                  << std::endl;
        exit(EXIT_FAILURE);
      }

      // Check limits on L
      if (maxL < 32)
      {
        std::cerr << GridLogError << "--max-L must be at least 32" << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  CartesianCommunicator::SetCommunicatorPolicy(
      CartesianCommunicator::CommunicatorPolicySequential);

  int NN;
  nlohmann::json json_results;
  Benchmark::Decomposition(NN, pattern, json_results);

  // Generate DeoFlops local volumes
  int sel = 4;
  int selm1 = sel - 1;
  std::vector<int> L_list({8, 12, 16, 24});
  for (int curL = 16; curL < maxL; curL *= 2)
  {
    if (curL * 2 <= maxL)
      L_list.push_back(curL * 2);
    if (curL * 3 <= maxL)
      L_list.push_back(curL * 3);
  }

  // Error check: make sure input maxL is a member of the local
  // volume list.
  if (L_list.back() != maxL)
  {
    int lowerL, higherL;
    if (L_list.size() % 2) // 2^n*2 < L < 2^n*3
    {
      int refL = L_list[L_list.size() - 3];
      lowerL = refL * 2;
      higherL = refL * 3;
    }
    else // 2^n*3 < L < 2^(n+1)*2
    {
      int refL = L_list[L_list.size() - 4];
      lowerL = refL * 3;
      higherL = refL * 4;
    }

    std::cout
        << GridLogError
        << "maxL must decompose to 2^n 3^m, for n>0 and m={0,1}. Nearest valid maxLs are "
        << lowerL << " and " << higherL << "." << std::endl;
    exit(EXIT_FAILURE);
  }

  std::vector<double> wilson_fp32;
  std::vector<double> wilson_fp64;
  std::vector<double> wilson_sp4_fp32;
  std::vector<double> wilson_sp4_fp64;
  std::vector<double> wilson_sp4_2as_fp32;
  std::vector<double> wilson_sp4_2as_fp64;
  std::vector<double> dwf4_fp32;
  std::vector<double> dwf4_fp64;
  std::vector<double> dwf4_su4_fp32;
  std::vector<double> dwf4_su4_fp64;
  std::vector<double> dwf4_sp4_fp32;
  std::vector<double> dwf4_sp4_fp64;
  std::vector<double> dwf4_sp4_2as_fp32;
  std::vector<double> dwf4_sp4_2as_fp64;
  std::vector<double> staggered_fp32;
  std::vector<double> staggered_fp64;

  auto runBenchmark = [&json_results](const std::string& name, const std::function<void(nlohmann::json&)>& fn)
  {
    grid_big_sep();
    std::cout << GridLogMessage << " " << name << " benchmark " << std::endl;
    grid_big_sep();
    fn(json_results);
  };

  if (do_memory)  runBenchmark("Memory",         &Benchmark::Memory);
  if (do_su4)     runBenchmark("SU(4)",          &Benchmark::SU4);
  if (do_comms)   runBenchmark("Communications", &Benchmark::Comms);
  if (do_latency) runBenchmark("Latency",        &Benchmark::Latency);
  if (do_p2p)     runBenchmark("Point-To-Point", &Benchmark::P2P);
  if (do_check_wilson)
  {
    grid_big_sep();
    std::cout << GridLogMessage << " Check Wilson" << std::endl;
    grid_big_sep();
    Benchmark::checkWilson();
  }

  int Ls = 12;
  if (do_flops)
  {
    auto runDeo = [&L_list, &pattern](const std::string& msg, int Ls, std::vector<double>& results, std::function<double(int,int,Grid::Coordinate)> fn)
    {
      grid_big_sep();
      std::cout << GridLogMessage << " " << msg << std::endl;
      for (int l = 0; l < L_list.size(); l++)
      {
        results.push_back(fn(Ls, L_list[l], pattern));
      }
    };

    // FP32 Deo FLOPS Benchmark
    // runDeo("fp32 SU(3) fundamental Wilson dslash 4d vectorised", 1, wilson_fp32, &Benchmark::DeoFlops<DomainWallFermionF>);
    runDeo("fp32 SU(3) fundamental Domain wall dslash 4d vectorised", Ls, dwf4_fp32, &Benchmark::DeoFlops<DomainWallFermionF>);
    if (do_flops_su4)
    {
      runDeo("fp32 SU(4) fundamental Domain wall dslash 4d vectorised", Ls, dwf4_su4_fp32, &Benchmark::DeoFlops<DomainWallFermionSU4F>);
    }
    if (do_flops_sp4_f)
    {
        // runDeo("fp32 Sp(4) fundamental Wilson dslash 4d vectorised", 1, wilson_sp4_fp32, &Benchmark::DeoFlops<DomainWallFermionSp4F>);
        runDeo("fp32 Sp(4) fundamental Domain wall dslash 4d vectorised", Ls, dwf4_sp4_fp32, &Benchmark::DeoFlops<DomainWallFermionSp4F>);
    }
    if (do_flops_sp4_2as)
    {
        // runDeo("fp32 Sp(4) two-index antisymmetric Wilson dslash 4d vectorised", 1, wilson_sp4_2as_fp32, &Benchmark::DeoFlops<DomainWallFermionSp4TwoIndexAntiSymmetricF>);
        runDeo("fp32 Sp(4) two-index antisymmetric Domain wall dslash 4d vectorised", Ls, dwf4_sp4_2as_fp32, &Benchmark::DeoFlops<DomainWallFermionSp4TwoIndexAntiSymmetricF>);
    }
    // runDeo("fp32 SU(3) fundamental Improved Staggered dslash 4d vectorised", 0, staggered_fp32, &Benchmark::DeoFlops<ImprovedStaggeredFermionF>);

    // FP64 Deo FLOPS Benchmark
    if (do_flops_fp64)
    {
      // runDeo("fp64 SU(3) fundamental Wilson dslash 4d vectorised", 1, wilson_fp64, &Benchmark::DeoFlops<DomainWallFermionD>);
      runDeo("fp64 SU(3) fundamental Domain wall dslash 4d vectorised", Ls, dwf4_fp64, &Benchmark::DeoFlops<DomainWallFermionD>);
      if (do_flops_su4)
      {
        runDeo("fp64 SU(4) fundamental Domain wall dslash 4d vectorised", Ls, dwf4_su4_fp64, &Benchmark::DeoFlops<DomainWallFermionSU4D>);
      }
      if (do_flops_sp4_f)
      {
        // runDeo("fp64 Sp(4) fundamental Wilson dslash 4d vectorised", 1, wilson_sp4_fp64, &Benchmark::DeoFlops<DomainWallFermionSp4D>);
        runDeo("fp64 Sp(4) fundamental Domain wall dslash 4d vectorised", Ls, dwf4_sp4_fp64, &Benchmark::DeoFlops<DomainWallFermionSp4D>);
      }
      if (do_flops_sp4_2as)
      {
        // runDeo("fp64 Sp(4) two-index antisymmetric Wilson dslash 4d vectorised", 1, wilson_sp4_2as_fp64, &Benchmark::DeoFlops<DomainWallFermionSp4TwoIndexAntiSymmetricD>);
        runDeo("fp64 Sp(4) two-index antisymmetric Domain wall dslash 4d vectorised", Ls, dwf4_sp4_2as_fp64, &Benchmark::DeoFlops<DomainWallFermionSp4TwoIndexAntiSymmetricD>);
      }
      // runDeo("fp64 SU(3) fundamental Improved Staggered dslash 4d vectorised", 0, staggered_fp64, &Benchmark::DeoFlops<ImprovedStaggeredFermionD>);
    }

    nlohmann::json tmp_flops;
    struct benchmarkResult
    {
      const char* header;
      const char* jsonName;
      const std::vector<double>& results;

      benchmarkResult(const char* header, const char* jsonName, const std::vector<double>& results)
        : header(header)
        , jsonName(jsonName)
        , results(results)
        {}
    };

    auto OutputDeoResults = [&NN, &tmp_flops, &L_list, &Ls]
      (const char* const precision,
       const std::vector<benchmarkResult>& results)
    {
      grid_big_sep();
      std::cout << GridLogMessage << "Gflop/s/node Summary table Ls=" << Ls << std::endl;
      std::cout << GridLogMessage << " * PRECISION: " << precision << std::endl;
      grid_big_sep();

      std::stringstream header;
      header << std::setw(5) << "L";
      for (auto & result : results)
      {
        if (!result.results.empty())
	        header << " " << std::setw(17) << std::string(result.header);
      }
      header << std::endl;

      grid_printf( "%s", header.str().c_str());
      for (int l = 0; l < L_list.size(); l++)
      {
        nlohmann::json tmp;
        tmp["L"] = L_list[l];
        tmp["Precision"] = precision;
        std::stringstream stream;
        stream.precision(2);
        stream << std::setw(5) << L_list[l];
        for (auto & result : results)
        {
          if (!result.results.empty())
          {
            stream << " " << std::setw(17) << std::fixed << result.results[l] / NN;
            tmp[std::string(result.jsonName)] = result.results[l] / NN;
          }
        }
	      stream << std::endl;
        tmp_flops["results"].push_back(tmp);
        grid_printf("%s", stream.str().c_str());
      }
    };

    OutputDeoResults("FP32",
    {
      {"Wilson", "Gflops_wilson", wilson_fp32},
      {"DWF", "Gflops_dwf4", dwf4_fp32},
      {"SU(4) DWF", "Gflops_dwf4_su4", dwf4_su4_fp32},
      {"Sp(4) Wilson fun", "Gflops_wilson_sp4_fund", wilson_sp4_fp32},
      {"Sp(4) DWF fun", "Gflops_dwf4_sp4_fund", dwf4_sp4_fp32},
      {"Sp(4) Wilson 2as", "Gflops_wilson_sp4_2as", wilson_sp4_2as_fp32},
      {"Sp(4) DWF 2as", "Gflops_dwf4_sp4_2as", dwf4_sp4_2as_fp32},
      {"Staggered", "Gflops_staggered", staggered_fp32}
    });

    if (do_flops_fp64)
    {
      OutputDeoResults("FP64",
      {
        {"Wilson", "Gflops_wilson", wilson_fp64},
        {"DWF", "Gflops_dwf4", dwf4_fp64},
        {"SU(4) DWF", "Gflops_dwf4_su4", dwf4_su4_fp64},
        {"Sp(4) Wilson fun", "Gflops_wilson_sp4_fund", wilson_sp4_fp64},
        {"Sp(4) DWF fun", "Gflops_dwf4_sp4_fund", dwf4_sp4_fp64},
        {"Sp(4) Wilson 2as", "Gflops_wilson_sp4_2as", wilson_sp4_2as_fp64},
        {"Sp(4) DWF 2as", "Gflops_dwf4_sp4_2as", dwf4_sp4_2as_fp64},
        {"Staggered", "Gflops_staggered", staggered_fp64},
      });
    }

    grid_big_sep();
    std::cout << GridLogMessage
              << " Comparison point     result: " << 0.5 * (dwf4_fp32[sel] + dwf4_fp32[selm1]) / NN
              << " Gflop/s per node" << std::endl;
    std::cout << GridLogMessage << " Comparison point is 0.5*(" << dwf4_fp32[sel] / NN << "+"
              << dwf4_fp32[selm1] / NN << ") " << std::endl;
    std::cout << std::setprecision(3);
    grid_big_sep();
    tmp_flops["comparison_point_Gflops"] = 0.5 * (dwf4_fp32[sel] + dwf4_fp32[selm1]) / NN;
    json_results["flops"] = tmp_flops;
  }

  // json_results["hostnames"] = get_mpi_hostnames();
  json_results["hostnames"] = "local\0";

  if (!json_filename.empty())
  {
    std::cout << GridLogMessage << "writing benchmark results to " << json_filename
              << std::endl;

    int me = 0;
    // MPI_Comm_rank(MPI_COMM_WORLD, &me);
    if (me == 0)
    {
      std::ofstream json_file(json_filename);
      json_file << std::setw(2) << json_results;
    }
  }

  Grid_finalize();
}
