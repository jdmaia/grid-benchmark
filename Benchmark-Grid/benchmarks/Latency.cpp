
#include "Benchmarks.hpp"

#include <Grid/Grid.h>
#include "Common.hpp"
#include "Benchmark-Grid/Utils.hpp"
#include "Benchmark-Grid/Roctx.hpp"

using namespace Grid;


void Benchmark::Latency(nlohmann::json& json_results)
{
  GRID_ROCTX_RANGE("Latency");
  int Nwarmup = 100;
  int Nloop = 300;

  std::cout << GridLogMessage << "Benchmarking point-to-point latency" << std::endl;
  grid_small_sep();
  grid_printf("from to      mean(usec)           err           max\n");

  int ranks;
  int me;
  MPI_Comm_size(MPI_COMM_WORLD, &ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  int bytes = 8;
  void *buf_from = acceleratorAllocDevice(bytes);
  void *buf_to = acceleratorAllocDevice(bytes);
  nlohmann::json json_latency;
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
      grid_printf("%2d %2d %15.4f %15.3f %15.4f\n", from, to, timestat.mean,
                  timestat.err, timestat.max);
      nlohmann::json tmp;
      tmp["from"] = from;
      tmp["to"] = to;
      tmp["time_usec"] = timestat.mean;
      tmp["time_usec_error"] = timestat.err;
      tmp["time_usec_min"] = timestat.min;
      tmp["time_usec_max"] = timestat.max;
      tmp["time_usec_full"] = t_time;
      json_latency.push_back(tmp);
    }
  json_results["latency"] = json_latency;

  acceleratorFreeDevice(buf_from);
  acceleratorFreeDevice(buf_to);
}
