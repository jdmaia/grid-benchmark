# Grid benchmarks

This repository contains benchmarks for the [Grid](https://github.com/paboyle/Grid) library.
The benchmarks can be summarised as follows.

- `Benchmark_Grid`: This benchmark measures floating-point performance for various fermion
  matrices, as well as bandwidth measurements for different operations. Measurements are
  performed for a fixed range of problem sizes.
- `Benchmark_IO`: Parallel I/O benchmark.

## TL;DR
Build and install Grid, all dependencies, and the benchmark with:
```bash
./bootstrap-env.sh <env_dir> <system>           # create benchmark environment
./build-grid.sh <env_dir> <config> <njobs>      # build Grid
./build-benchmark.sh <env_dir> <config> <njobs> # build benchmarks
```
where `<env_dir>` is an arbitrary directory where every product will be stored, `<system>`
is a sub-directory of `systems` containing system-specific configuration files 
(an existing preset or your own), and finally `<config>` is the name of a build config
in `systems/<system>/grid-config.json`. After a successful execution, the benchmark binaries
will be in `<env_dir>/prefix/gridbench_<config>`. Build tasks are executed using `<njobs>`
parallel processes.

## Developing configurations for additional systems

### System-specific directory
You can create a configuration for a new system by creating a new subdirectory in the
`systems` directory that should contain at least:
- a `grid-config.json` file according to the specification described below;
- a `files` directory containing any files that need to be copied to the environment directory.

### Configuration file
The system directory must contain a `grid-config.json` file specifying compilation flag
configurations for Grid. This file must contain a single `"configs"` JSON array, with all
elements of the form
```json
{
  "name": "foo",          // name of the configuration
  "env-script": "bar.sh", // additional script to source before building 
                          // (e.g. to load modules, path absolute or relative to the 
                          // environment directory, ignored if empty)
  "commit": "...",        // Grid commit to use 
                          // (anything that can be an argument of git checkout)
  "config-options": "..." // options to pass to the configure script,
  "env" : {               // environment variables
    "VAR": "value"        // export VAR="value" before building
  }
  "pixi-env": "..."       // Pixi environment to use for this configuration
}
```
Grid's dependencies are managed with [Pixi](https://pixi.sh/latest/) environments defined
in the [`pixi.toml`](pixi.toml) file. The following environments are available:
- `gpu-nvidia`: NVIDIA GPU Linux build
- `gpu-amd`: AMD GPU Linux build
- `cpu-linux`: CPU Linux build with LLVM
- `cpu-apple-silicon`: Apple Silicon build on macOS
and one of these strings must be used as a value for `"pixi-env"` above.

Please refer to [Grid's repository](https://github.com/paboyle/Grid) 
for documentation on the use of Grid's `configure` script.

### Environment setup
Once a complete system folder has been created as above, the associated environment can be
deployed with:
```bash
./bootstrap-env.sh ./env <system>
```
where `<system>` is the name of the system directory in `systems`. Here, `./env` is an
example of a deployment location, but any writable path can be used. This script will
install Pixi and deploy the relevant environments in `./env`, as well as copy all files
present in the system `files` directory. After successful completion, `./env` will contain
an `env.sh` file that can be sourced to activate a given environment:
```bash
source ./env/env.sh <config>
```
where `<config>` must match a `"name"` field from the `grid-config.json` file.

This script will:
1. make the embedded Pixi path available;
2. activate the Pixi environment specified in `"pixi-env"`; and
3. source the (optional) additional script specified in `"env-script"`.

The procedure above is used by the scripts `build-grid.sh` and `build-benchmark.sh`,
which can be used to build Grid and the benchmark as described above.

## Running the benchmarks
After building the benchmarks as described above, you can find the binaries in 
`<env_dir>/prefix/gridbench_<config>`. Depending on the system selected, the environment
directory might also contain example batch scripts. Each HPC system tends to have its own
runtime characteristics, and it is not possible to automate determining the best runtime
environment for the Grid benchmark. Examples of known supercomputing environments can be found:
- in the [`systems` directory](systems) of this repository; and
- in the [`systems` directory](https://github.com/paboyle/Grid/tree/develop/systems) of the Grid repository.

More information about the benchmark results is provided below.

### `Benchmark_Grid`
This benchmark performs flop/s measurements for typical lattice QCD sparse matrices, as
well as memory and inter-process bandwidth measurements using Grid routines. The benchmark
command accepts any Grid flag (see the complete list with `--help`), as well as any benchmark-specific flags. The Grid `--help` string for commit `6165931afaa53a9885b6183ff762fc2477f30b51` is reproduced below for convenience.
```
  --help : this message
Geometry:
  --mpi n.n.n.n   : default MPI decomposition
  --threads n     : default number of OMP threads
  --grid n.n.n.n  : default Grid size
  --shm  M        : allocate M megabytes of shared memory for comms
  --shm-mpi 0|1   : Force MPI usage under multi-rank per node 
  --shm-hugepages : use explicit huge pages in mmap call 
  --device-mem M  : Size of device software cache for lattice fields (MB) 
Verbose:
  --log list      : comma separated list from Error,Warning,Message,Performance,Iterative,Integrator,Debug,Colours
  --notimestamp   : suppress millisecond resolution stamps
  --decomposition : report on default omp,mpi and simd decomposition
Debug:
  --dylib-map     : print dynamic library map, useful for interpreting signal backtraces 
  --heartbeat     : periodic itimer wakeup (interrupts stuck system calls!) 
  --signal-delay n : pause for n seconds after signal handling (useful to get ALL nodes in stuck state) 
  --debug-stdout  : print stdout from EVERY node to file Grid.stdout/err.rank 
  --debug-signals : catch sigsegv and print a blame report, handle SIGHUP with a backtrace to stderr
  --debug-heartbeat : periodically report backtrace 
  --debug-mem     : print Grid allocator activity
Performance:
  --comms-overlap    : Overlap comms with compute
  --dslash-generic: Wilson kernel for generic Nc
  --dslash-unroll : Wilson kernel for Nc=3
  --dslash-asm    : Wilson kernel for AVX512
```
In addition, there are flags that are undocumented by the help string:
Flag | Result
--|--
--accelerator-threads | Multiplies the number of threads in a threadblock on the device.



The benchmark flags are as follows:

Flag | Result
--|--
--help                       | Display all benchmark CLI flags and Grid CLI flags.
--json-out \<file\>          | Save the measurement results in JSON format to `<file>`.
--benchmark-memory           | Enable axpy memory benchmark (default=on).
--no-benchmark-memory        | Disable axpy memory benchmark.
--benchmark-su4              | Enable SU(4) memory benchmark (default=on).
--no-benchmark-su4           | Disable SU(4) memory benchmark.
--benchmark-comms            | Enable communications benchmark (default=on).
--no-benchmark-comms         | Disable communications benchmark.
--benchmark-flops            | Enable all Dirac Matrix Flops benchmarks (default=on).
--no-benchmark-flops         | Disable all Dirac Matrix Flops benchmarks.
--benchmark-flops-su4        | Enable SU(4) Dirac Matrix Flops benchmark (default=on).
--no-benchmark-flops-su4     | Disable SU(4) Dirac Matrix Flops benchmark.
--benchmark-flops-sp4-f      | Enable Sp(4) Dirac Matrix Flops benchmark (default=on).
--no-benchmark-flops-sp4-f   | Disable Sp(4) Dirac Matrix Flops benchmark.
--benchmark-flops-sp4-2as    | Enable Sp(4) Two-Index AntiSymmetric Dirac Matrix Flops benchmark (default=on).
--no-benchmark-flops-sp4-2as | Disable Sp(4) Two-Index AntiSymmetric Dirac Matrix Flops benchmark.
--benchmark-flops-fp64       | Enable FP64 Dirac Matrix Flops benchmarks (default=on).
--no-benchmark-flops-fp64    | Disable FP64 Dirac Matrix Flops benchmarks.
--benchmark-latency          | Enable point-to-point communications latency benchmark (default=off).
--no-benchmark-latency       | Disable point-to-point communications latency benchmark.
--benchmark-p2p              | Enable point-to-point communications bandwidth benchmark (default=off).
--no-benchmark-p2p           | Disable point-to-point communications bandwidth benchmark.
--check-wilson               | Enable Wilson Fermion correctness check (default=on).
--no-check-wilson            | Disable Wilson Fermion correctness check.
--dwf4-only                  | Only run the SU(3) fundamental Domain wall dslash Flops benchmark, for L=24, 32 and 48, and only for the overlapped comms/compute scheme with the concurrent communicator policy (`G/O/C`) instead of the usual sweep over all four comms schemes. Disables every other benchmark and correctness check, and ignores `--max-L`. FP64 can still be turned off with `--no-benchmark-flops-fp64`.
--pattern \<x.y.z.t\>        | Scales the local lattice dimensions by the factors in the string x.y.z.t. This allows the volume-per-node to be matched between systems with dfferent GPUs-per-node if desired.
--max-L \<size\>             | Sets the maximum lattice size for the flops benchmarks. This must either be a power of 2, or 3 times a power of 2.

It is intended that only the `--json-out` flag should be necessary for standard use cases of the benchmark. The 
benchmarks are performed on a fixed set of problem sizes, and the Grid flag `--grid` will
be ignored. The *Floating-point performace* benchmark can run on additional larger problems sizes by specifying a valid `--max-L` flag.

The resulting metrics are as follows. All data-size units are in base 2 
(i.e. 1 kB = 1024 B) and flops are given in base 10 (i.e. 1 kB = 1000 B).

*Memory bandwidth*

One sub-benchmark measures the memory bandwidth using a lattice version of the `axpy` BLAS
routine, similar to the STREAM benchmark. This benchmark can be disabled with the `--no-benchmark-memory` CLI argument. The JSON entries under `"axpy"` have the form:
```json
{
  "GBps": 215.80653375861607,   // bandwidth in GB/s/node
  "GFlops": 19.310041765757834, // FP performance (double precision)
  "L": 8,                       // local lattice volume
  "size_MB": 3.0                // memory size in MB/node
}
```

A second benchmark performs site-wise SU(4) matrix multiplication and has a higher
arithmetic intensity than the `axpy` test (although it is still memory-bound). This benchmark can be disabled with the `--no-benchmark-su4` CLI argument. 
The JSON entries under `"SU4"` have the form:
```json
{
  "GBps": 394.76639187026865,  // bandwidth in GB/s/node
  "GFlops": 529.8464820758512, // FP performance (single precision)
  "L": 8,                      // local lattice size
  "size_MB": 6.0               // memory size in MB/node
}
```

*Inter-process bandwidth*

This sub-benchmark measures the achieved bidirectional bandwidth in threaded halo exchange
using routines in Grid. The exchange is performed in each direction on the MPI Cartesian
grid, which is parallelised across at least two processes. The resulting bandwidth is
related to node-local transfers (inter-CPU, NVLink, ...) or network transfers depending
on the MPI decomposition. This benchmark can be disabled with the `--no-benchmark-comms` CLI argument. The JSON entries under `"comms"` have the form:
```json
{
  "L": 40,                       // local lattice size
  "bytes": 73728000,             // payload size in B/rank
  "dir": 2,                      // direction of the exchange, 8 possible directions
                                 // (0: +x, 1: +y, ..., 5: -x, 6: -y, ...)
  "rate_GBps": {
    "error": 6.474271894240327,  // standard deviation across measurements (GB/s/node)
    "max": 183.10546875,         // maximum measured bandwidth (GB/s/node)
    "mean": 175.21747026766676   // average measured bandwidth (GB/s/node)
  },
  "time_usec": 3135.055          // average transfer time (microseconds)
}
```

*Floating-point performances*

This sub-benchmark measures the achieved floating-point performance for the Wilson-, 
domain-wall- and staggered-fermion sparse matrices provided by Grid. For Wilson and domain-wall, multiple matrix sizes per lattice site are tested. By default, both single- and double-precision versions of the benchmark are performed. The matrix size at each site is as follows:
name | element count
-- | --
Wilson | 3*4
Wilson Sp(4) | 4*4
Wilson Sp(4) Two-Index Antisymmetric | 5*4
domain-wall | 3*4
domain-wall SU(4) | 4*4
domain-wall Sp(4) | 4*4
domain-wall Sp(4) Two-Index Antisymmetric | 5*4
staggered | 3*4

The entire benchmark can be disabled with the `--no-benchmark-flops` flag. Individual sub-benchmarks can be disabled as follows:
- `--no-benchmark-flops-su4`: domain-wall SU(4)
- `--no-benchmark-flops-sp4-f`: Wilson Sp(4), domain-wall Sp(4)
- `--no-benchmark-flops-sp4-2as`: Wilson Sp(4) Two-Index Antisymmetric, domain-wall Sp(4) Two-Index Antisymmetric

Wilson, domain-wall, and staggered are not able to be disabled.

The double-precision versions of the benchmarks can be disabled with `--no-benchmark-fp64`.

For devices with large amounts of dedicated memory, it is possible to run the benchmark for additional problem sizes with the `--max-L <maximum size>` flag. Note that the maximum size *must* either be 2<sup>n</sup> or 3*2<sup>n</sup>. Since the default maximum size is 32 (the **minimum** acceptable `--max-L`) the next possible value is 48. The next highest and lowest valid values will be reported if an invalid value is entered.

In the `"flops"`
and `"results"` sections of the JSON output the best performances are recorded, e.g.:
```json
{
  "Gflops_dwf4": 930.5314401622717,             // domain-wall in Gflop/s/node
  "Gflops_dwf4_sp4_2as": 2029.311501559971,     // domain-wall Sp(4) Two-Index Antisymmetric in Gflop/s/node
  "Gflops_dwf4_sp4_fund": 1449.8298297273077,   // domain-wall Sp(4) in Gflop/s/node
  "Gflops_dwf4_su4": 1451.4278683482005,        // domain-wall SU(4) in Gflop/s/node
  "Gflops_staggered": 8.179053400638082,        // staggered in Gflop/s/node
  "Gflops_wilson": 55.21422625196834,           // Wilson in Gflop/s/node
  "Gflops_wilson_sp4_2as": 235.7911263021898,   // Wilson Sp(4) Two-Index Antisymmetric in Gflop/s/node
  "Gflops_wilson_sp4_fund": 156.44161527750148, // Wilson Sp(4) in Gflop/s/node
  "L": 8,                                       // local lattice size
  "Precision": "FP32"                           // data type (FP32 or FP64)
}
```
Here "best" means the best result across the different implementations of the routines.
Please see the benchmark log for a detailed breakdown. Finally, the JSON output contains
a "comparison point", which is the average of the L=24 and L=32 best domain-wall (`"Gflops_dwf4"`) single-precision performances.

*Peer-to-Peer latency*

This sub-benchmark is disabled by default and can be enabled with the `--benchmark-latency` CLI argument. A very small payload is transferred between every possible pair of MPI ranks to measure the communications latency between all ranks. The JSON entries under `"latency"` have the form:
```json
{
  "from": 0,               // The MPI rank sending data
  "to":   1,               // The MPI rank receiving data
  "time_usec": ...,        // The mean latency (microseconds)
  "time_usec_error": ...,  // The standard deviation of the latency (microseconds)
  "time_usec_min": ...,    // The smallest measured latency (microseconds)
  "time_usec_max": ...,    // The largest measured latency (microseconds)
  "time_usec_full": [...]  // All latency timings (microseconds)
}
```

*Peer-to-Peer bandwidth*

This sub-benchmark is disabled by default and can be enabled with the `--benchmark-p2p` CLI argument. A very large payload is transferred between every possible pair of MPI ranks to measure the communications bandwidth between all ranks. The JSON entries under `"p2p"` have the form:
```json
{
  "from": 0,               // The MPI rank sending data
  "to":   1,               // The MPI rank receiving data
  "bytes": 127401984,      // The payload size (B)
  "time_usec": ...,        // The mean transfer time (microseconds)
  "time_usec_error": ...,  // The standard deviation of the transfer time (microseconds)
  "time_usec_max": ...,    // The largest measured transfer time (microseconds)
  "time_usec_min": ...,    // The smallest measured transfer time (microseconds)
  "time_usec_full": [...], // All transfer time timings (microseconds)
  "rate_GBps": ...
}
```
The information about the transfer rate is separated to a sub-entry under `"rate_GBps"` with the following structure:
```json
{
  "mean": ...,   // The mean transfer rate (GB/s)
  "error": ...,  // The standard deviation of the transfer rate (GB/s)
  "max": ...,    // The largest measured transfer rate (GB/s)
  "min": ...     // The smallest measured transfer rate (GB/s)
}
```
As is consistent with the data units in the other benchmarks, these data rates are in base 2 (1 kB = 1024 B).

### `Benchmark_IO`

This benchmark tests the parallel I/O performance of Grid for both reading and writing, and compares this against the I/O performance of the C++ standard library. This benchmark is non-configurable and runs on a fixed set of problem sizes. Measurement results can be saved to a file with the `--json-out <file>` flag. All data-size units are in base 2 (i.e. kB = 1024 B). The benchmark flags are as follows:

Flag | Result
--|--
--help                       | Display all benchmark CLI flags and Grid CLI flags.
--json-out \<file\>          | Save the measurement results in JSON format to `<file>`.

The benchmark repeats 10 times. On each pass, the benchmark will:
1) Write a Grid structure to disk using `std::ofstream` using one file per process,
2) Read the same Grid structure back into memory using `std::ifstream` using one file per process,
3) Write a Grid structure to a single file using Grid's parallel I/O,
4) Read the same Grid structure back into memory using Grid's parallel I/O from the single file.

The results for each pass are written into the `"passes"` section of the output JSON. This consists of performance averages over all problem sizes 24 to 32 inclusive, as well as the performance on each problem size individually. Each pass has the following structure:
```json
{
  "Pass": 1,                           // pass index, 1-10
  "av_grid_read": 6426.877979887839,   // average Grid parallel I/O read speed in GB/s
  "av_grid_write": 1283.9060495508716, // average Grid parallel I/O write speed in GB/s
  "av_max_L": 32,                      // largest problem size used in the average
  "av_min_L": 24,                      // smallest problem size used in the average
  "av_std_read": 51489.06630810099,    // std::ifstream I/O read speed in GB/s
  "av_std_write": 5371.805416924437,   // std::ofstream I/O write speed in GB/s
  "volumes": [...]                     // individual results on each problem size
}
```
The subsection `"volumes"` of the `"passes"` section records results for each problem size on each pass in the following form:
```json
{
  "L": 8,                           // local lattice size
  "grid_read": 631.8449873631004,   // Grid parallel I/O read speed in GB/s
  "grid_write": 193.74838543012143, // Grid parallel I/O write speed in GB/s
  "std_read": 77922.07792207792,    // std::ifstream I/O read speed in GB/s
  "std_write": 4930.156121610517    // std::ofstream I/O write speed in GB/s
},
```
Summary statistics for the passes and the volume-averages are recorded in the `"average"` section of the JSON, which has the following structure:
```json
{
  "av_grid_read_mean": 6771.762525331016,     // mean of 'av_grid_read' over all 10 passes in GB/s
  "av_grid_read_rob": 96.36351333776568,      // robustness of 'av_grid_read' over all 10 passes in %
  "av_grid_read_stddev": 246.25424103184403,  // standard deviation of 'av_grid_read' over all 10 passes in GB/s
  "av_grid_write_mean": 1301.360669990908,    // mean of 'av_grid_write' over all 10 passes in GB/s
  "av_grid_write_rob": 98.05881576571363,     // robustness of 'av_grid_write' over all 10 passes in %
  "av_grid_write_stddev": 25.261808157067016, // standard deviation of 'av_grid_write' over all 10 passes in GB/s
  "av_max_L": 32,                             // largest problem size used in the average
  "av_min_L": 24,                             // smallest problem size used in the average
  "av_std_read_mean": 49775.25973673006,      // mean of 'av_std_read' over all 10 passes in GB/s
  "av_std_read_rob": 97.51189971101847,       // robustness of 'av_grid_read' over all 10 passes in %
  "av_std_read_stddev": 1238.4583813508862,   // standard deviation of 'av_grid_read' over all 10 passes in GB/s
  "av_std_write_mean": 5411.033884034247,     // mean of 'av_grid_write' over all 10 passes in GB/s
  "av_std_write_rob": 97.10663406573798,      // robustness of 'av_grid_write' over all 10 passes in %
  "av_std_write_stddev": 156.5610110920218,   // standard deviation of 'av_grid_write' over all 10 passes in GB/s,
  "volumes": [...]
}
```

'Robustness' is defined as `1 - std.dev. / |mean|`, where `|mean|` is the absolute value of the mean. Once again the pass-averages per-volume are individually recorded in the `"volumes"` subsection with the following structure:

```json
{
  "L": 8,                                  // local lattice size
  "grid_read_mean": 638.2219086851617,     // mean of 'grid_read' over all 10 passes in GB/s
  "grid_read_rob": 59.7194074010008,       // robustness of 'grid_read' over all 10 passes in GB/s
  "grid_read_stddev": 257.07956691502665,  // standard deviation of 'grid_read' over all 10 passes in %
  "grid_write_mean": 397.54878301223994,   // mean of 'grid_write' over all 10 passes in GB/s
  "grid_write_rob": 70.384645339958,       // robustness of 'grid_write' over all 10 passes in %
  "grid_write_stddev": 117.73548203575567, // standard deviation of 'grid_write' over all 10 passes in GB/s
  "std_read_mean": 63641.793997669476,     // mean of 'std_read' over all 10 passes in GB/s
  "std_read_rob": 84.87482271891673,       // robustness of 'std_read' over all 10 passes in %
  "std_read_stddev": 9625.934167009314,    // standard deviation of 'std_read' over all 10 passes in GB/s
  "std_write_mean": 5539.262098862089,     // mean of 'std_write' over all 10 passes in GB/s
  "std_write_rob": 93.89967914021757,      // robustness of 'std_write' over all 10 passes in %
  "std_write_stddev": 337.9127612949058    // standard deviation of 'std_write' over all 10 passes in GB/s
}
```
