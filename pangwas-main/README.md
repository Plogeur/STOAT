# pangwas


## About

`pangwas` finds associations in pangenomes using [snarls](https://github.com/vgteam/vg/wiki/Snarls-and-chains) to identify variants. 

## Dependencies

- [`htslib`](https://www.htslib.org/) (required by libvgio)
- Boost: Boost can be installed by running `sudo apt-get install libboost-all-dev`

`pangwas` includes some scripts that use [`vg`](https://github.com/vgteam/vg) to build and modify graphs.

## Installation

`pangwas` can be built from source using CMake.

First, clone the repo and its dependencies:

 ```
git clone --recursive https://forgemia.inra.fr/xian-hui.chang/pangwas.git
cd pangwas
```

Next, build using CMake in a directory `build`. 

```
mkdir build
cd build
cmake ..
make
```

The binary file will be placed in the directory `bin` and can be run with:

```
./bin/pangwas
```

The `bin` directory can be added to your `PATH` variable to allow `pangwas` to be run from any directory.
From the `pangwas` directory, run:

```
echo 'export PATH="${PATH}:'"$(pwd)"'/bin"' >>~/.bashrc
```

Then close your terminal and open it again, or run 

```
source ~/.bashrc
```

## Running in docker
`pangwas` is also available in a [docker container](https://hub.docker.com/repository/docker/xhchang/pangwas/general).
The docker container can be run with:

```
docker run -it -v [local_path_to_data]:/work/data -w /work/ -u `id -u $USER` --rm --entrypoint /bin/bash xhchang/pangwas
```

## Usage

### Preparing the graph

`pangwas` takes as input a variation graph and a [distance index](https://github.com/vgteam/vg/wiki/Index-Types#distance-index).
The variation graph can be in `.hg`, `.pg`, or `.gbz` format (from [libbdsg](https://github.com/vgteam/libbdsg) and [gbwtgraph](https://github.com/jltsiren/gbwtgraph)).
A graph in `.gfa` format can be converted to `.gbz` format and indexed using the script `convert_and_index_gfa.sh`. 

```
./scripts/convert_and_index_gfa.sh [input_graph_base] [output_graph_base] [reference_sample]
```

The input `graph_base` should not include the extension `.gfa`.
`reference_sample` is used to specify which paths to use as reference paths, if they are not already present as reference-sense paths in the graph.
It must match the sample/locus of one or more paths in the graph (See the [vg wiki](https://github.com/vgteam/vg/wiki/Path-Metadata-Model)). 


### Running `pangwas`

The base command for finding associated snarls is:
```
pangwas --graph [graph.hg] --distance-index [graph.dist] --sample-of-interest [sample name]  
```

`--sample-of-interest` is used to specify the samples with the trait of interest.
It may be repeated multiple times.
The sample name must match the sample of a path in the graph.
These can be found with `vg paths -M -x [graph].hg`.
Note that for generic-sense paths, the sample name may be in the `LOCUS` field. 

### Options

- `--method`

  `pangwas` currently has only one method 
  
  - `paths`: The `paths` method uses the paths in the graph itself to find associations. 
    Within each snarl, paths are partitioned based on the walk they take through the snarl.
    If one partition matches exactly the set of samples of interest, then the snarl is considered to be associated.

- `--output-format`
  
  The output format can be specified with `--output-format (tsv / fasta)`.
  
  - The tsv output contains one locus per line. It is formatted:
  
    ```
    #reference_path\tstart_offset\tend_offset\tvariant_length
    ```
    
    The start and end coordinates are given for the reference path. 
    If there isn't a reference-sense path in the graph, the reference specified by `--reference-sample` will be used. 
    If there is no reference-sense path and `--reference-sample` was not given, or if the reference did not traverse the variant, then any sample that traversed the variant may be used as `reference_path`. 
  
  - The fasta output produces a sequence for each path in each snarl.
    The description line is formatted:
  
    ```
    >snarl_name|sample_coordinates|reference_coordinates
    ```
    
    The snarl name is formatted `snarl:start_id-end_id` and the sample and reference coordinates are formatted `path_name:start_offset-end_offset`.
    As with the tsv output, the reference coordinates use the real reference if possible, and any path if not.
  
