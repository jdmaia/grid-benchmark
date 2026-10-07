#include "Benchmarks.hpp"

#include <Grid/Grid.h>
#include "Common.hpp"
#include "Benchmark-Grid/Utils.hpp"
#include "Benchmark-Grid/Roctx.hpp"

using namespace Grid;


void Benchmark::P2P(nlohmann::json& json_results)
{
  GRID_ROCTX_RANGE("P2P");
  // IMPORTANT: The P2P benchmark uses "MPI_COMM_WORLD" communicator, which is
  // not the quite the same as Grid.communicator. Practically speaking, the
  // latter one contains the same MPI-ranks but in a different order. Grid
  // does this make sure it can exploit ranks with shared memory (i.e.
  // multiple ranks on the same node) as best as possible.

  // buffer-size to benchmark. This number is the same as the largest one used
  // in the "Comms()" benchmark. ( L=48, Ls=12, double-prec-complex,
  // half-color-spin-vector. ). Mostly an arbitrary choice, but nice to match
  // it here
  size_t bytes = 127401984;

  int Nwarmup = 20;
  int Nloop = 100;

  std::cout << GridLogMessage << "Benchmarking point-to-point bandwidth" << std::endl;
  grid_small_sep();
  grid_printf("from to      mean(usec)           err           min           "
              "bytes    rate (GiB/s)\n");

  int ranks;
  int me;
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  void *buf_from = acceleratorAllocDevice(bytes);
  void *buf_to = acceleratorAllocDevice(bytes);
  nlohmann::json json_p2p;
  for (int from = 0; from < ranks; ++from)
    for (int to = 0; to < ranks; ++to)
    {
      if (from == to)
        continue;

      GRID_ROCTX_RANGE((std::string("ping-pong ") + std::to_string(from) + "->" +
                         std::to_string(to)).c_str());

      std::vector<double> t_time(Nloop);
      time_statistics timestat;
      MPI_Status status;

      for (int i = -Nwarmup; i < Nloop; ++i)
      {
        double start = usecond_precise();
        if (from == me)
        {
          auto err = MPI_Send(buf_from, bytes, MPI_CHAR, to, 0, MPI_COMM_WORLD);
          assert(err == MPI_SUCCESS);
        }
        if (to == me)
        {
          auto err =
              MPI_Recv(buf_to, bytes, MPI_CHAR, from, 0, MPI_COMM_WORLD, &status);
          assert(err == MPI_SUCCESS);
        }
        double stop = usecond_precise();
        if (i >= 0)
          t_time[i] = stop - start;
      }
      // important: only 'from' and 'to' have meaningful timings. we use
      // 'from's.
      MPI_Bcast(t_time.data(), Nloop, MPI_DOUBLE, from, MPI_COMM_WORLD);

      timestat.statistics(t_time);
      double rate = bytes / (timestat.mean / 1.e6) / 1024. / 1024. / 1024.;
      double rate_err = rate * timestat.err / timestat.mean;
      double rate_max = rate * timestat.mean / timestat.min;
      double rate_min = rate * timestat.mean / timestat.max;

      grid_printf("%2d %2d %15.4f %15.3f %15.4f %15zu %15.2f\n", from, to,
                  timestat.mean, timestat.err, timestat.min, bytes, rate);

      nlohmann::json tmp;
      tmp["from"] = from;
      tmp["to"] = to;
      tmp["bytes"] = bytes;
      tmp["time_usec"] = timestat.mean;
      tmp["time_usec_error"] = timestat.err;
      tmp["time_usec_min"] = timestat.min;
      tmp["time_usec_max"] = timestat.max;
      tmp["time_usec_full"] = t_time;
      nlohmann::json tmp_rate;
      tmp_rate["mean"] = rate;
      tmp_rate["error"] = rate_err;
      tmp_rate["max"] = rate_max;
      tmp_rate["min"] = rate_min;
      tmp["rate_GBps"] = tmp_rate;
     json_p2p.push_back(tmp);
    }
  json_results["p2p"] = json_p2p;

  acceleratorFreeDevice(buf_from);
  acceleratorFreeDevice(buf_to);
}
