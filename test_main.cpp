#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

using namespace std;
using namespace __gnu_pbds;

// Order Statistic Tree
template <class Key>
using ost = tree<Key, null_type, less<Key>, rb_tree_tag, tree_order_statistics_node_update>;

struct RNG {
    mt19937_64 gen;
    uniform_real_distribution<double> dist;
    RNG(uint64_t seed=123456789ULL) : gen(seed), dist(0.0,1.0) {}
    inline double uni01() { return dist(gen); }
    inline int randint(int n) { return (int)(gen() % n); } // fast, small bias ok here
};

// ---------- Naive algorithm ----------
inline int sample_rank_naive(RNG &rng, int n, int x) {
    int best = n;
    for (int i = 0; i < x; ++i) {
        int r = rng.randint(n) + 1;
        if (r < best) best = r;
    }
    return best;
}

// ---------- Inverse transform ----------
inline int sample_rank_inverse(RNG &rng, int n, int x) {
    double U = rng.uni01();
    double K = n - n * pow(1.0 - U, 1.0 / x);
    int Ki = (int)ceil(K);
    if (Ki < 1) Ki = 1;
    if (Ki > n) Ki = n;
    return Ki;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // -------- parameters --------
    int n = 20000;          // |P|
    int x = 2000;           // number of draws
    int trials = 2000000;  // number of comparisons

    RNG rng1(2025), rng2(2025);

    // OST initialization
    ost<int> T;
    for (int i = 1; i <= n; ++i) T.insert(i);

    vector<long long> cnt_naive(n+1, 0);
    vector<long long> cnt_ost(n+1, 0);

    // -------- naive timing --------
    auto t1 = chrono::high_resolution_clock::now();
    for (int i = 0; i < trials; ++i) {
        int k = sample_rank_naive(rng1, n, x);
        cnt_naive[k]++;
    }
    auto t2 = chrono::high_resolution_clock::now();

    // -------- OST timing --------
    auto t3 = chrono::high_resolution_clock::now();
    for (int i = 0; i < trials; ++i) {
        int k = sample_rank_inverse(rng2, n, x);
        auto it = T.find_by_order(k - 1); // rank-select
        (void)it;                         // suppress unused warning
        cnt_ost[k]++;
    }
    auto t4 = chrono::high_resolution_clock::now();

    // -------- output --------
    cout << "n = " << n << ", x = " << x << ", trials = " << trials << "\n\n";
    cout << "rank\tnaive\tost\n";
    for (int k = 1; k <= n; ++k) {
        cout << k << "\t" << cnt_naive[k] << "\t" << cnt_ost[k] << "\n";
    }

    auto naive_time = chrono::duration_cast<chrono::milliseconds>(t2 - t1).count();
    auto ost_time   = chrono::duration_cast<chrono::milliseconds>(t4 - t3).count();

    cout << "\nTime (ms):\n";
    cout << "Naive: " << naive_time << " ms\n";
    cout << "OST  : " << ost_time << " ms\n";

    return 0;
}