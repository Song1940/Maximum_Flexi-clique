#pragma once
// Link-Cut Tree (Splay-based) for dynamic connectivity
// Supports: link(u,v), cut(u,v), connected(u,v) in O(log n) amortized

struct LCTNode {
    LCTNode *ch[2], *par;
    bool rev;
    LCTNode() : par(nullptr), rev(false) { ch[0] = ch[1] = nullptr; }
};

namespace LCT {

inline bool isRoot(LCTNode* x) {
    return !x->par || (x->par->ch[0] != x && x->par->ch[1] != x);
}

inline void pushRev(LCTNode* x) {
    if (!x) return;
    std::swap(x->ch[0], x->ch[1]);
    x->rev ^= true;
}

inline void push(LCTNode* x) {
    if (x->rev) {
        pushRev(x->ch[0]);
        pushRev(x->ch[1]);
        x->rev = false;
    }
}

inline void rotate(LCTNode* x) {
    LCTNode* p = x->par;
    LCTNode* g = p->par;
    int i = (p->ch[1] == x), j = i ^ 1;

    if (!isRoot(p)) {
        if (g->ch[0] == p) g->ch[0] = x;
        else if (g->ch[1] == p) g->ch[1] = x;
    }
    x->par = g;
    p->ch[i] = x->ch[j];
    if (x->ch[j]) x->ch[j]->par = p;
    x->ch[j] = p;
    p->par = x;
}

// splay x to root of its auxiliary tree
void splay(LCTNode* x) {
    // push ancestors first
    static LCTNode* stk[128];
    int top = 0;
    LCTNode* y = x;
    stk[top++] = y;
    while (!isRoot(y)) { y = y->par; stk[top++] = y; }
    while (top) push(stk[--top]);

    while (!isRoot(x)) {
        LCTNode* p = x->par;
        if (!isRoot(p)) {
            LCTNode* g = p->par;
            if ((g->ch[0] == p) == (p->ch[0] == x)) rotate(p);
            else rotate(x);
        }
        rotate(x);
    }
}

void access(LCTNode* x) {
    LCTNode* last = nullptr;
    for (LCTNode* y = x; y; y = y->par) {
        splay(y);
        y->ch[1] = last;
        last = y;
    }
    splay(x);
}

void makeRoot(LCTNode* x) {
    access(x);
    pushRev(x);
    push(x);
}

LCTNode* findRoot(LCTNode* x) {
    access(x);
    while (x->ch[0]) { push(x); x = x->ch[0]; }
    splay(x);
    return x;
}

bool connected(LCTNode* u, LCTNode* v) {
    return findRoot(u) == findRoot(v);
}

void link(LCTNode* u, LCTNode* v) {
    makeRoot(u);
    if (findRoot(v) != u) u->par = v;
}

void cut(LCTNode* u, LCTNode* v) {
    makeRoot(u);
    access(v);
    if (v->ch[0] == u && !u->ch[1]) {
        v->ch[0] = nullptr;
        u->par = nullptr;
    }
}

} // namespace LCT


// Per-tree wrapper
struct LinkCutForest {
    int n;
    std::vector<LCTNode> nodes;

    LinkCutForest() : n(0) {}
    explicit LinkCutForest(int n) : n(n), nodes(n) {}

    bool connected(int u, int v) { return LCT::connected(&nodes[u], &nodes[v]); }
    void link(int u, int v)      { LCT::link(&nodes[u], &nodes[v]); }
    void cut(int u, int v)       { LCT::cut(&nodes[u], &nodes[v]); }
};
