# Maximum_Flexi_Clique

This is the implementation of Maximum Flexi-Clique branch and bound algorithm, which is described in the following papaer submitted in VLDB 2027:
- Efficient Computation of Maximum Flexi-Clique in Networks

C++ implementation of exact and heuristic algorithms for finding the **maximum τ-flexi-clique** in an undirected graph.

A **τ-flexi-clique** of size *k* is a connected subgraph in which every node has degree at least ⌊k^τ⌋ within the subgraph (τ ∈ (0, 1]).

## Algorithms

| Binary | Description |
|--------|-------------|
| `flexi` | **EBA** — Efficient Branch-and-Bound (exact) |
| `flexi` | **FPA** — Flexi-Prune Algorithm (heuristic) |
| `flexi` | **NPA** — Naive Peeling Algorithm (heuristic) |

All three algorithms are compiled into a single binary and selected via `--algo`.

## Requirements

- C++17
- CMake ≥ 3.14

## Build

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

This produces a single binary `flexi` in the `build/` directory.

## Usage

```
./build/flexi --algo <npa|fpa|eba> --file <path/to/graph.dat> --tau <0.0–1.0> [--out <result.txt>]
```

**Example:**
```bash
./build/flexi --algo eba --file data/dolphins.dat --tau 0.9
```

**Output:**
```
=== Result ===
Algorithm : eba
File      : data/dolphins.dat
Tau       : 0.9
Size      : 12
Time(s)   : 0.003
Valid     : yes
Branches  : 847
```

## Input Format

Plain text edge list (whitespace-separated). Comment lines starting with `#` are ignored. Node IDs can be any non-negative integers; they are re-indexed to `[0, n)` internally.

```
# optional comment
0 1
0 2
1 3
...
```

