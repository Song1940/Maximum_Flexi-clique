#pragma once
#include "graph.h"
#include "npa.h"

// FPA: Flexi-Prune Algorithm
// K-core seed + Modified Greedy++ initialization + connectivity-aware peeling.
AlgoResult runFPA(const Graph& G, double tau);

// Modified Greedy++ heuristic used internally by FPA and EBA.
std::vector<int> modifiedGreedyPlusPlus(const Graph& G, double tau);
