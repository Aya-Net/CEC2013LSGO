#pragma once

#include <random>
#include <cmath>
#include <stdexcept>
#include <shared_mutex>
#include <mutex>

// ===============================
// 节点值：(val, id)
// ===============================
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

// 按 key split：左子树 < key，右子树 >= key
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

// merge 两棵 Treap（保证所有 a < 所有 b）
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

// 找第 k 小（1-based）
inline TValue find_k(Node* t, int k) {
    if (!t || k < 1 || k > sz(t)) {
        throw std::out_of_range("find_k: k out of range");
    }
    int left = sz(t->l);
    if (k == left + 1) return t->v;
    if (k <= left) return find_k(t->l, k);
    return find_k(t->r, k - left - 1);
}

// ===============================
// 线程安全 top-y% 抽样器（多读单写）
// ===============================
class TopYTreapSampler {
private:
    Node* root = nullptr;

    mutable std::shared_mutex mtx;   // 读写锁

    std::mt19937 rng_;
    std::uniform_int_distribution<int> prand;
    std::uniform_real_distribution<double> urand;

public:
    TopYTreapSampler()
        : prand(1, 2000000000),
          urand(0.0, 1.0)
    {
        rng_.seed(std::random_device{}());
    }

    ~TopYTreapSampler() {
        // 简单递归释放（单线程销毁时调用）
        destroy(root);
    }

    // 写：需要独占锁
    void insert(double val, int id) {
        std::unique_lock<std::shared_mutex> lock(mtx);

        TValue tv{val, id};
        Node *a, *b;
        split(root, tv, a, b);
        Node* nw = new Node(val, id, prand(rng_));
        root = merge(merge(a, nw), b);
    }

    // 读：共享锁，多线程并发读安全
    int size() const {
        std::shared_lock<std::shared_mutex> lock(mtx);
        return sz(root);
    }

    // 读：按 rank 查询第 k 小（1-based）
    TValue find_by_rank(int k) const {
        std::shared_lock<std::shared_mutex> lock(mtx);
        return find_k(root, k);
    }

    // 读：从前 y_frac 部分，模拟“抽 x 次再取最小”
    TValue sample_top_y_min(double y_frac, int x, std::mt19937 &rng) {
        std::shared_lock<std::shared_mutex> lock(mtx);

        int n = sz(root);
        if (n == 0) throw std::runtime_error("no elements");

        if (y_frac <= 0.0) y_frac = 0.0;
        if (y_frac > 1.0)  y_frac = 1.0;

        int m = (int)std::floor(y_frac * n);
        if (m < 1) m = 1;

        double U = urand(const_cast<std::mt19937&>(rng));

        // K = ceil(m - m * (1-U)^(1/x))
        double rawK = m - m * std::pow(1.0 - U, 1.0 / x);
        int K = (int)std::ceil(rawK);
        if (K < 1) K = 1;
        if (K > m) K = m;

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
