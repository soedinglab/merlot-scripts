## benchmark scripts

A collection of scripts that was used to generate and benchmark the simulations provided with the tool MERLoT. Download the expression matrices and parameter files via the [GWDG FTP server](wwwuser.gwdg.de/~compbiol/merlot/). These will be used as examples for the various scripts.

Let `$scripts` be the directory where `merlot-scripts` is downloaded and `$dir` the directory where the simulations are downloaded and unpacked.

More details are available as comments in each script.

### benchmarkN_serial.sh
The entry point script. Will generate simulations, calculate predictions, and evaluate them. This command (provided line 31 is commented out) will run `destiny`, `SLICER`, `Monocle2` and `MERLoT` on the 100 simulated single bifurcation datasets provided for MERLoT:

`$scripts/benchmarkN_serial.sh $scripts test $dir/benchmark1 20 1`

This script calls a runner script for each simulation via the [`parallel`](https://www.gnu.org/software/parallel/) utility.

### run_as_one.sh
This script will calculate the embeddings - run diffusion maps for `destiny` and DDRTree embeddings for `Monocle2` in various configurations. This is done because the diffusion maps and DDRTree embeddings are used by multiple method flavors further in the pipeline and so it is more efficient to calculate them once. The script will then call the master evaluator script for each dataset.

(includes `run_destiny.R`, `run_monocle.R`)

### timed_benchmarksN.sh
This script evaluates the runs of each method. By default, it lets each method run for 60 minutes and terminates them if they need longer.

(includes `benchmark_MERLoT_dest.R`, `benchmark_MERLoT_mon.R`, `benchmark_SLICER.R`, `benchmark_destiny.R`, `benchmark_monocl2.R`, `benchmark_stopped.R`)

### simulating data
The repository also includes the scripts needed to run simulations (`run_simN.sh`) in the framework of `benchmarkN_serial.sh`. In particular, the parameters for the actual simulations are set in `generate_simN.py`.

### benchmark evaluation & plotting
The repository also includes two collections of scripts that were used to read and evaluate the output of the benchmark. These continue the evaluation indices used in the paper (F1 Measure, Jaccard Index, Matthews Correlation Coefficient, adjusted Mutual Information, Goodman-Kruskall Index) as well as others (`evaluate_method.R`, `various.R`).

## Dependencies

The following R packages are required:

From CRAN:

- `infotheo`
- `igraph`
- `entropy`
- `optparse`
- `FNN`
- `viridis`
- `ggplot2`
- `ggnetwork`
- `ggrepel`
- `network`
- `intergraph`
- `mclust`
- `pdist`
- `hashmap`
- `reshape2`

From Bioconductor:

- `destiny`
- `monocle`, version 2.7
- `AnnotationDbi`
- `org.Mm.eg.db`
- `org.Hs.eg.db`
- `topGO`

Other:
- [`rprosstt`](https://github.com/soedinglab/prosstt-r)
- [`merlot`](https://github.com/soedinglab/merlot)
- `slingshot`
- `ElPiGraph.R`
- `TSCAN`
- `SLICER`

Python dependencies:

- `numpy`
- `scipy`
- `pandas`
- `sklearn`
- 

- `prosstt`