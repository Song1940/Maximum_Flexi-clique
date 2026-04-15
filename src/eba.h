#pragma once
#include "graph.h"
#include "npa.h"

// EBA: Efficient Branch-and-Bound Algorithm for Maximum Flexi-Clique.
//
// Pruning rules can be disabled at compile time via preprocessor defines:
//   -DNO_RULE1        disable Rule 1 (min adjusted degree check)
//   -DNO_RULE2        disable Rule 2 (upper-bound size check)
//   -DNO_RULE3        disable Rule 3 (diameter-based candidate pruning)
//   -DNO_RULE4        disable Rule 4 (degree-distance combined pruning)
//   -DNO_RULE5        disable Rule 5 (follower-based sibling pruning)
//   -DNO_RULE6        disable Rule 6 (global degree threshold pruning)
//   -DNO_DEGREE_ORDER disable degree-ascending branching order

AlgoResult runEBA(const Graph& G, double tau);
