#include "Benchmarks.hpp"

#include <Grid/Grid.h>
#include "Common.hpp"
#include "Benchmark-Grid/Utils.hpp"
#include "Benchmark-Grid/Roctx.hpp"

using namespace Grid;


void Benchmark::Memory(nlohmann::json& json_results)
{
  GRID_ROCTX_RANGE("Memory");
  const int Nwarmup = 50;
  const int Nvec = 8;
  typedef Lattice<iVector<vReal, Nvec>> LatticeVec;
  typedef iVector<vReal, Nvec> Vec;

  Coordinate simd_layout = GridDefaultSimd(Nd, vReal::Nsimd());
  Coordinate mpi_layout = GridDefaultMpi();

  std::cout << GridLogMessage << "Benchmarking a*x + y bandwidth" << std::endl;
  grid_small_sep();
  grid_printf("%5s %15s %15s %15s %15s\n", "L", "size (MB/node)", "time (usec)",
              "GB/s/node", "Gflop/s/node");

  uint64_t NN;
  uint64_t lmax = 64;

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

    Vec rn;
    random(sRNG, rn);

    LatticeVec z(&Grid);
    z = Zero();
    LatticeVec x(&Grid);
    x = Zero();
    LatticeVec y(&Grid);
    y = Zero();
    double a = 2.0;

    const uint64_t Nloop = (200 * lmax * lmax * lmax / lat / lat / lat);

    {
      GRID_ROCTX_RANGE("Warmup");
      for (int i = 0; i < Nwarmup; i++)
      {
        z = a * x - y;
      }
    }
    double start = usecond();
    {
      GRID_ROCTX_RANGE("axpy timed loop");
      for (int i = 0; i < Nloop; i++)
      {
        z = a * x - y;
      }
    }
    double stop = usecond();
    double time = (stop - start) / Nloop / 1.e6;

    double flops = vol * Nvec * 2 / 1.e9; // mul,add
    double bytes = 3.0 * vol * Nvec * sizeof(Real) / 1024. / 1024.;

    grid_printf("%5d %15.2f %15.2f %15.2f %15.2f\n", lat, bytes / NN, time * 1.e6,
                bytes / time / NN / 1024., flops / time / NN);

    nlohmann::json tmp;
    tmp["L"] = lat;
    tmp["size_MB"] = bytes / NN;
    tmp["GBps"] = bytes / time / NN / 1024.;
    tmp["GFlops"] = flops / time / NN;
    json_results["axpy"].push_back(tmp);
  }
};
