#pragma once
// Fully Dynamic Connectivity — HDT (Holm-de Lichtenberg-Thorup)
// Translated from Python FPA.py reference implementation
//
// Complexity (amortized):
//   insertEdge : O(log n)
//   deleteEdge : O(log² n)
//   connected  : O(log n)
//
// comp_count tracks number of connected components.
// AP check in FPA:
//   delete all edges of u  →  diff = new_comp - old_comp
//   diff == 1  →  only u isolated, rest connected  →  non-AP
//   diff  > 1  →  additional splits               →  AP

#include <vector>
#include <set>
#include <unordered_map>
#include <cmath>
#include "lct.h"

struct FDC {
    int n;
    int L;           // ceil(log2(n))
    int comp_count;  // number of connected components

    std::vector<LinkCutForest>                forests;        // forests[0..L]
    std::vector<std::set<std::pair<int,int>>> non_tree_edges; // per-level non-tree edges

    std::unordered_map<long long, int>  edge_level;
    std::unordered_map<long long, bool> is_tree_edge;

    static long long key(int u, int v) {
        if (u > v) std::swap(u, v);
        return (long long)u << 32 | (unsigned)v;
    }

    explicit FDC(int n) : n(n), comp_count(n) {
        L = (n > 1) ? (int)std::ceil(std::log2((double)n)) : 1;
        forests.reserve(L + 1);
        for (int i = 0; i <= L; ++i) forests.emplace_back(n);
        non_tree_edges.resize(L + 1);
    }

    bool connected(int u, int v) {
        if (u == v) return true;
        return forests[0].connected(u, v);
    }

    int get_comp_count() const { return comp_count; }

    // ── insertEdge ────────────────────────────────────────────────────────────
    void insertEdge(int u, int v) {
        if (u == v) return;
        if (u > v) std::swap(u, v);
        long long k = key(u, v);
        if (edge_level.count(k)) return;

        edge_level[k] = 0;
        if (!forests[0].connected(u, v)) {
            for (int i = 0; i <= L; ++i) forests[i].link(u, v);
            is_tree_edge[k] = true;
            --comp_count;
        } else {
            non_tree_edges[0].insert({u, v});
            is_tree_edge[k] = false;
        }
    }

    // ── deleteEdge ────────────────────────────────────────────────────────────
    void deleteEdge(int u, int v) {
        if (u == v) return;
        if (u > v) std::swap(u, v);
        long long k = key(u, v);
        auto it = edge_level.find(k);
        if (it == edge_level.end()) return;

        int level = it->second;
        edge_level.erase(it);
        bool was_tree = is_tree_edge[k];
        is_tree_edge.erase(k);

        if (!was_tree) {
            non_tree_edges[level].erase({u, v});
            return;
        }

        // Cut from all forests at level..L
        for (int i = level; i <= L; ++i) {
            if (forests[i].connected(u, v))
                forests[i].cut(u, v);
        }

        // Search for replacement non-tree edge (level down to 0)
        std::pair<int,int> replacement = {-1, -1};
        int rep_level = -1;
        for (int l = level; l >= 0 && rep_level < 0; --l) {
            for (auto& e : non_tree_edges[l]) {
                if (!forests[l].connected(e.first, e.second)) {
                    replacement = e;
                    rep_level = l;
                    break;
                }
            }
        }

        if (rep_level >= 0) {
            non_tree_edges[rep_level].erase(replacement);
            long long rk = key(replacement.first, replacement.second);
            is_tree_edge[rk] = true;
            for (int i = rep_level; i <= L; ++i) {
                if (!forests[i].connected(replacement.first, replacement.second))
                    forests[i].link(replacement.first, replacement.second);
            }
        } else {
            ++comp_count;
        }

        // Promote non-tree edges now within the same component
        for (int l = 0; l <= level; ++l) {
            std::vector<std::pair<int,int>> to_promote;
            for (auto& e : non_tree_edges[l]) {
                if (forests[l].connected(e.first, e.second))
                    to_promote.push_back(e);
            }
            for (auto& e : to_promote) {
                non_tree_edges[l].erase(e);
                int nl = l + 1;
                if (nl <= L) {
                    non_tree_edges[nl].insert(e);
                    edge_level[key(e.first, e.second)] = nl;
                } else {
                    edge_level.erase(key(e.first, e.second));
                }
            }
        }
    }
};
