#include "Benchmarks.hpp"

#include <Grid/Grid.h>
#include "Common.hpp"
#include "Benchmark-Grid/Utils.hpp"
#include "Benchmark-Grid/Roctx.hpp"

using namespace Grid;


void Benchmark::SU4(nlohmann::json& json_results)
{
  GRID_ROCTX_RANGE("SU4");
  const int Nwarmup = 50;
  const int Nc4 = 4;
  typedef Lattice<iMatrix<vComplexF, 4>> LatticeSU4;

  Coordinate simd_layout = GridDefaultSimd(Nd, vComplexF::Nsimd());
  Coordinate mpi_layout = GridDefaultMpi();

  std::cout << GridLogMessage << "Benchmarking z = y*x SU(4) bandwidth" << std::endl;
  grid_small_sep();
  grid_printf("%5s %15s %15s %15s %15s\n", "L", "size (MB/node)", "time (usec)",
              "GB/s/node", "Gflop/s/node");

  uint64_t NN;

  uint64_t lmax = 48;

  GridSerialRNG sRNG;
  sRNG.SeedFixedIntegers(std::vector<int>({45, 12, 81, 9}));
  for (int lat = 8; lat <= lmax; lat += 8)
  {
    GRID_ROCTX_RANGE((std::string("lat=") + std::to_string(lat)).c_str());

    Coordinate latt_size({lat * mpi_layout[0], lat * mpi_layout[1], lat * mpi_layout[2],
                          lat * mpi_layout[3]});
    double vol =
        static_cast<double>(latt_size[0]) * latt_size[1] * latt_size[2] * latt_size[3];

    GridCartesian Grid(latt_size, simd_layout, mpi_layout);

    NN = Grid.NodeCount();

    LatticeSU4 z(&Grid);
    z = Zero();
    LatticeSU4 x(&Grid);
    x = Zero();
    LatticeSU4 y(&Grid);
    y = Zero();

    const uint64_t Nloop = (200 * lmax * lmax * lmax / lat / lat / lat);

    {
      GRID_ROCTX_RANGE("Warmup");
      for (int i = 0; i < Nwarmup; i++)
      {
        z = x * y;
      }
    }
    double start = usecond();
    {
      GRID_ROCTX_RANGE("SU4 mul timed loop");
      for (int i = 0; i < Nloop; i++)
      {
        z = x * y;
      }
    }
    double stop = usecond();
    double time = (stop - start) / Nloop / 1.e6;

    double flops = vol * Nc4 * Nc4 * (6 + (Nc4 - 1) * 8) / 1.e9; // mul,add
    double bytes = 3.0 * vol * Nc4 * Nc4 * 2 * sizeof(RealF) / 1024. / 1024.;
    grid_printf("%5d %15.2f %15.2f %15.2f %15.2f\n", lat, bytes / NN, time * 1.e6,
                bytes / time / NN / 1024., flops / time / NN);

    nlohmann::json tmp;
    tmp["L"] = lat;
    tmp["size_MB"] = bytes / NN;
    tmp["GBps"] = bytes / time / NN / 1024.;
    tmp["GFlops"] = flops / time / NN;
    json_results["SU4"].push_back(tmp);
  }
};
