#pragma once

#include <random>
#include <cmath>
#include <stdexcept>
#include <shared_mutex>
#include <mutex>
#include <vector>

struct TValue {
    double val;
    int id;
};

struct Node {
    TValue v;
    int prior;
    int sz;
    Node* l;
    Node* r;

    Node(double val, int id, int p)
        : v{val, id}, prior(p), sz(1), l(nullptr), r(nullptr) {}
};

inline int sz(Node* t) { return t ? t->sz : 0; }
inline void pull(Node* t) {
    if (t) t->sz = 1 + sz(t->l) + sz(t->r);
}

inline bool lessThan(const TValue& a, const TValue& b) {
    if (a.val != b.val) return a.val < b.val;
    return a.id < b.id;
}

inline void split(Node* t, const TValue& key, Node*& a, Node*& b) {
    if (!t) { a = b = nullptr; return; }
    if (lessThan(t->v, key)) {
        split(t->r, key, t->r, b);
        a = t;
        pull(a);
    } else {
        split(t->l, key, a, t->l);
        b = t;
        pull(b);
    }
}

inline Node* merge(Node* a, Node* b) {
    if (!a) return b;
    if (!b) return a;
    if (a->prior > b->prior) {
        a->r = merge(a->r, b);
        pull(a);
        return a;
    } else {
        b->l = merge(a, b->l);
        pull(b);
        return b;
    }
}

inline TValue find_k(Node* t, int k) {
    if (!t || k < 0 || k >= sz(t)) {
        throw std::out_of_range("find_k: k out of range");
    }
    int left = sz(t->l);
    if (k == left) return t->v;
    if (k < left) return find_k(t->l, k);
    return find_k(t->r, k - left - 1);
}

class TreapSamplerTS {
private:
    Node* root = nullptr;
    mutable std::shared_mutex mtx;

    std::mt19937 rng;
    std::uniform_int_distribution<int> prand;
    std::uniform_real_distribution<double> urand;

public:
    TreapSamplerTS()
        : prand(1, 2000000000),
          urand(0.0, 1.0)
    {
        rng.seed(std::random_device{}());
    }

    ~TreapSamplerTS() {
        destroy(root);
    }

    void insert(double val, int id) {
        std::unique_lock<std::shared_mutex> lock(mtx);

        TValue tv{val, id};
        Node *a, *b;
        split(root, tv, a, b);
        Node* nw = new Node(val, id, prand(rng));
        root = merge(merge(a, nw), b);
    }

    int size() const {
        std::shared_lock<std::shared_mutex> lock(mtx);
        return sz(root);
    }

    TValue find_by_rank(int k) const {
        std::shared_lock<std::shared_mutex> lock(mtx);
        return find_k(root, k);
    }

    TValue sample_min_x(int x) const {
        std::shared_lock<std::shared_mutex> lock(mtx);

        int n = sz(root);
        if (n == 0) throw std::runtime_error("sample_min_x: empty");

        static thread_local std::mt19937 trng(std::random_device{}());
        static thread_local std::uniform_real_distribution<double> uni(0.0, 1.0);

        double U = uni(trng);
        // 1-based: ceil(m - m*(1-U)^(1/x))
        double rawK = n - n * std::pow(1.0 - U, 1.0 / x);
        int K = (int)std::ceil(rawK) - 1; // 转 0-based

        if (K < 0) K = 0;
        if (K >= n) K = n - 1;

        return find_k(root, K);
    }

    TValue sample_top_y_min(double y_frac, int x) const {
        std::shared_lock<std::shared_mutex> lock(mtx);

        int n = sz(root);
        if (n == 0) return {0, -1};

        if (y_frac <= 0) y_frac = 0;
        if (y_frac > 1)  y_frac = 1;

        int m = (int)std::floor(y_frac * n);
        if (m < 1) m = 1;

        thread_local std::mt19937 trng(std::random_device{}());
        thread_local std::uniform_real_distribution<double> uni(0.0, 1.0);

        double U = uni(trng);

        double rawK = m - m * std::pow(1.0 - U, 1.0 / x);
        int K = (int)std::ceil(rawK) - 1;

        if (K < 0) K = 0;
        if (K >= m) K = m - 1;

        return find_k(root, K);
    }

private:
    static void destroy(Node* t) {
        if (!t) return;
        destroy(t->l);
        destroy(t->r);
        delete t;
    }
};


class TournamentSampler {
private:
    int mod;
    TreapSamplerTS global;
    std::vector<TreapSamplerTS> buckets;

public:
    TournamentSampler(int mod_)
        : mod(mod_), buckets(mod_) {}

    void insert(double val, int id) {
        global.insert(val, id);

        int b = id % mod;
        if (b < 0) b += mod;
        buckets[b].insert(val, id);
    }

    TValue sample_top_y_min(double y_frac, int x) const {
        return global.sample_top_y_min(y_frac, x);
    }

    TValue sample_min_dpt(int x, int module) const {
        int b = module % mod;
        if (b < 0) b += mod;
        return buckets[b].sample_min_x(x);
    }
};
