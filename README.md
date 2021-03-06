# MovingFrame.jl

Combining moving frame and shape analysis.

## Dependency
The MovingFrame need working installation of python package igl and scipy. It can be installed using conda command:

```
conda install -c conda-forge igl
conda install -c anaconda scipy
```

or calling conda from julia:

```julia
julia> using Conda
julia> Conda.add("igl";channel="conda-forge")
julia> Conda.add("scipy")
```


## Installation
MovingFrame.jl requires Julia version 1.3.0 or above. To install MovingFrame.jl run the following command inside a Julia session:

```julia
julia> using Pkg
julia> Pkg.add("https://github.com/yusri-dh/MovingFrame.jl")
```

## Example
Sample code is available from the Pluto notebooks in [example.jl](example.jl) or Jupyter notebooks in [example.ipynb](example.ipynb). The example can also run directly on the cloud on [Nextjournal](https://nextjournal.com/yusri/integrated-analysis-of-cell-shape-and-movement-in-moving-frame)
