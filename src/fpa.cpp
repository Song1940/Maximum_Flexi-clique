#include "fpa.h"
#include <algorithm>
#include <chrono>
#include <queue>
#include <set>
#include <unordered_map>
#include <unordered_set>

// Modified Greedy++: iteratively removes the minimum-score node per component.
// Score(v) = accumulated load + current degree within the subgraph.
std::vector<int> modifiedGreedyPlusPlus(const Graph& G, double tau) {
    std::vector<int> best;
    int global_best = 0;

    std::vector<bool> visited(G.n, false);
    std::vector<std::vector<int>> init_components;
    for (int s = 0; s < G.n; ++s) {
        if (visited[s]) continue;
        std::vector<int> comp;
        std::queue<int> q;
        q.push(s); visited[s] = true;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            comp.push_back(u);
            for (int w : G.adj[u])
                if (!visited[w]) { visited[w] = true; q.push(w); }
        }
        init_components.push_back(std::move(comp));
    }

    std::sort(init_components.begin(), init_components.end(),
              [](const auto& a, const auto& b){ return a.size() > b.size(); });

    std::queue<std::vector<int>> comp_q;
    for (auto& c : init_components) comp_q.push(std::move(c));

    while (!comp_q.empty()) {
        std::vector<int> comp = std::move(comp_q.front());
        comp_q.pop();
        if ((int)comp.size() <= global_best) continue;

        std::unordered_set<int> comp_set(comp.begin(), comp.end());
        std::unordered_map<int,int>    cur_deg;
        std::unordered_map<int,double> load;
        for (int v : comp) {
            cur_deg[v] = 0; load[v] = 0.0;
            for (int u : G.adj[v])
                if (comp_set.count(u)) ++cur_deg[v];
        }

        int comp_min_deg = INT_MAX;
        for (int v : comp) comp_min_deg = std::min(comp_min_deg, cur_deg[v]);
        int max_possible = (comp_min_deg > 0)
            ? (int)std::floor(std::pow((double)comp_min_deg, 1.0 / tau)) : 0;
        if (max_possible <= global_best) continue;

        std::set<int> remaining(comp.begin(), comp.end());

        while (!remaining.empty()) {
            int sz        = (int)remaining.size();
            int threshold = floorPow(sz, tau);

            int    min_node  = *remaining.begin();
            double min_score = load[min_node] + cur_deg[min_node];
            for (int v : remaining) {
                double sc = load[v] + cur_deg[v];
                if (sc < min_score) { min_score = sc; min_node = v; }
            }

            if (cur_deg[min_node] >= threshold) {
                std::vector<bool> rem_act(G.n, false);
                for (int v : remaining) rem_act[v] = true;

                if (G.isConnected(rem_act)) {
                    if (sz > global_best) {
                        global_best = sz;
                        best.assign(remaining.begin(), remaining.end());
                    }
                } else {
                    std::vector<bool> sub_vis(G.n, false);
                    for (int s : remaining) {
                        if (sub_vis[s]) continue;
                        std::vector<int> sub;
                        std::queue<int> bq;
                        bq.push(s); sub_vis[s] = true;
                        while (!bq.empty()) {
                            int u = bq.front(); bq.pop();
                            sub.push_back(u);
                            for (int w : G.adj[u])
                                if (rem_act[w] && !sub_vis[w]) { sub_vis[w] = true; bq.push(w); }
                        }
                        if ((int)sub.size() > global_best) comp_q.push(std::move(sub));
                    }
                }
                break;
            }

            load[min_node] += cur_deg[min_node];
            for (int u : G.adj[min_node])
                if (remaining.count(u)) --cur_deg[u];
            remaining.erase(min_node);
        }
    }

    return best;
}

// FPA: Flexi-Prune Algorithm
//   1. K-core seed selection.
//   2. Modified Greedy++ on the seed scope.
//   3. Connectivity-aware peeling on the scope.
AlgoResult runFPA(const Graph& G, double tau) {
    auto t0 = std::chrono::high_resolution_clock::now();

    // Step 1: k-core seed selection.
    std::vector<int> core = G.coreNumbers();
    int k_star = *std::max_element(core.begin(), core.end());

    std::vector<bool> seed_mask(G.n, false);
    int best_k = k_star;

    for (int k = 2; k <= k_star; ++k) {
        std::vector<bool> mask(G.n, false);
        for (int v = 0; v < G.n; ++v)
            if (core[v] >= k) mask[v] = true;
        std::vector<bool> lcc = G.largestCC(mask);
        int sz = 0; for (bool b : lcc) sz += b;
        if (floorPow(sz, tau) <= k) { best_k = k; seed_mask = lcc; break; }
    }

    // Scope = (best_k - 1)-core component containing the seed.
    std::vector<bool> scope(G.n, false);
    {
        std::vector<bool> lower(G.n, false);
        for (int v = 0; v < G.n; ++v)
            if (core[v] >= best_k - 1) lower[v] = true;

        int seed_node = -1;
        for (int v = 0; v < G.n; ++v) if (seed_mask[v]) { seed_node = v; break; }

        if (seed_node >= 0) {
            std::vector<bool> vis(G.n, false);
            std::queue<int> q;
            q.push(seed_node); vis[seed_node] = true;
            while (!q.empty()) {
                int u = q.front(); q.pop();
                scope[u] = true;
                for (int w : G.adj[u])
                    if (lower[w] && !vis[w]) { vis[w] = true; q.push(w); }
            }
        } else {
            scope = lower;
        }
    }

    // Step 2: run Modified Greedy++ on a subgraph induced by scope nodes.
    std::vector<int> scope_nodes;
    std::vector<int> remap(G.n, -1);
    for (int v = 0; v < G.n; ++v)
        if (scope[v]) { remap[v] = (int)scope_nodes.size(); scope_nodes.push_back(v); }

    int sn = (int)scope_nodes.size();
    Graph subG(sn);
    for (int i = 0; i < sn; ++i) {
        int v = scope_nodes[i];
        for (int u : G.adj[v])
            if (scope[u] && remap[u] > i) subG.addEdge(i, remap[u]);
    }
    subG.finalize();

    std::vector<int> greedy_local = modifiedGreedyPlusPlus(subG, tau);
    std::vector<int> greedy_result;
    for (int i : greedy_local) greedy_result.push_back(scope_nodes[i]);

    // Step 3: connectivity-aware peeling on scope.
    std::vector<bool> active = scope;
    std::vector<int>  loc_deg(G.n, 0);
    int active_cnt = 0;
    for (int v = 0; v < G.n; ++v) {
        if (!active[v]) continue;
        ++active_cnt;
        for (int u : G.adj[v]) if (active[u]) ++loc_deg[v];
    }

    std::vector<int> best_result = greedy_result;
    int best_size = (int)greedy_result.size();

    auto tryUpdate = [&]() {
        if (active_cnt <= best_size) return;
        int thr = floorPow(active_cnt, tau);
        for (int v = 0; v < G.n; ++v)
            if (active[v] && loc_deg[v] < thr) return;
        if (!G.isConnected(active)) return;
        best_size = active_cnt;
        best_result.clear();
        for (int v = 0; v < G.n; ++v) if (active[v]) best_result.push_back(v);
    };

    tryUpdate();

    bool changed = true;
    while (changed && active_cnt > 0) {
        changed = false;
        int thr = floorPow(active_cnt, tau);

        std::vector<int> cands;
        for (int v = 0; v < G.n; ++v) if (active[v]) cands.push_back(v);
        std::sort(cands.begin(), cands.end(),
                  [&](int a, int b){ return loc_deg[a] < loc_deg[b]; });

        for (int u : cands) {
            if (loc_deg[u] >= thr) break;
            if (!G.isArticulationPoint(u, active)) {
                active[u] = false;
                --active_cnt;
                for (int w : G.adj[u]) if (active[w]) --loc_deg[w];
                loc_deg[u] = 0;
                tryUpdate();
                changed = true;
                break;
            }
        }
    }

    auto t1 = std::chrono::high_resolution_clock::now();

    AlgoResult res;
    res.nodes    = best_result;
    res.time_sec = std::chrono::duration<double>(t1 - t0).count();
    return res;
}
