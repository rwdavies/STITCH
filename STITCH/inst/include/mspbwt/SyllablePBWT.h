#ifndef SYLLABLEPBWT_H_
#define SYLLABLEPBWT_H_

#include "vcfpp/vcfpp.h"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <unistd.h>
#include <unordered_map>
#include <vector>

using namespace std;

using IntMapU = unordered_map<int, int>;

// C++11 compatible
template<typename T, typename = typename std::enable_if<std::is_unsigned<T>::value>::type>
class SyllablePBWT
{
  public:
    static const unsigned long long MOD{-1ull - 58}, BASE{(1ull << 48) - 59};
    const int B{sizeof(T) * 8};
    int M{0}, N{0}, n{0};
    int minSiteL{B * 2 - 1};
    vector<vector<int>> x, a;
    vector<vector<T>> r;
    vector<vector<unsigned long long>> h;
    vector<T> z_;
    vector<unsigned __int128> xp;
    vector<unsigned long long> hz;
    vector<int> z, up_end, dn_end;

  public:
    SyllablePBWT(){};
    virtual ~SyllablePBWT(){};

    int build(std::string vcfpanel, std::string samples, std::string region)
    {
        vcfpp::BcfReader vcf(vcfpanel, samples, region);
        vcfpp::BcfRecord var(vcf.header);
        // get M = # haplotypes, N = # sites, n = # syllables, IDs
        M = vcf.nsamples * 2;
        N = 0;
        vector<bool> gt;
        while(vcf.getNextVariant(var))
        {
            var.getGenotypes(gt);
            if(!var.isNoneMissing() || !var.allPhased())
                throw runtime_error("there is a site with non-phased or missing value!\n");
            N++;
        }
        n = (N + B - 1) / B;
        // resize data structures for syllabic panel
        x.resize(M, vector<int>(n));
        r.resize(n);
        a.resize(n + 1, vector<int>(M));
        iota(a[0].begin(), a[0].end(), 0);
        vcf.setRegion(region); // seek back region
        {
            vector<T> x_(M);

            // K = current site, k = current syllable, kp1 = k+1
            for(int K = 0, k = 0, kp1 = 1; K < N; K++)
            {
                vcf.getNextVariant(var);
                var.getGenotypes(gt);
                // update current syllable with current site
                for(int index = 0; index < M; index++)
                {
                    x_[index] = (x_[index] << 1) | (gt[index] != 0);
                }

                if(K == N - 1)
                { // if last site, pad last syllable with 0s
                    int pad = n * B - N;
                    for(int i = 0; i < M; i++) x_[i] <<= pad;
                    K += pad;
                }
                if((K + 1) % B > 0) continue; // skip if not end of syllable

                // coordinate-compress syllable values
                vector<T> x_copy = x_;
                sort(x_copy.begin(), x_copy.end());
                x_copy.resize(unique(x_copy.begin(), x_copy.end()) - x_copy.begin());
                r[k] = x_copy;
                for(int i = 0; i < M; i++)
                {
                    x[i][k] = lower_bound(r[k].begin(), r[k].end(), x_[i]) - r[k].begin();
                }
                memset(&x_[0], 0, M * sizeof(T));

                // counting sort to build a[k+1] with a[k] and x[][k]
                int r_sz = r[k].size();
                vector<vector<int>> counter(r_sz);
                for(int i = 0; i < M; i++)
                {
                    int a_ = a[k][i];
                    counter[x[a_][k]].push_back(a_);
                }
                for(int i = 0, p = 0; i < r_sz; i++)
                {
                    for(int v : counter[i])
                    {
                        a[kp1][p++] = v;
                    }
                }
                k++, kp1++;
            }
        }

        xp.resize(n), hz.resize(n + 1), up_end.resize(n + 1), dn_end.resize(n + 1), z_.resize(n), z.resize(n);
        // xp[i] = pow(BASE, i) % MOD
        xp[0] = 1;
        for(int k = 1; k < n; k++)
        {
            xp[k] = xp[k - 1] * BASE % MOD;
        }

        // build prefix hashes; h[i][k+1] = hash of x[i][0, k]
        h.resize(M, vector<unsigned long long>(n + 1));
        for(int i = 0; i < M; i++)
        {
            for(int k = 0; k < n; k++)
            {
                h[i][k + 1] = (h[i][k] + xp[k] * (x[i][k] + 1)) % MOD;
            }
        }

        return 0;
    }

    int save(std::string save_file)
    {
        ofstream out(save_file, ios::binary);
        if(out.fail()) return 1;

        out.write((char *)&M, sizeof(M));
        out.write((char *)&N, sizeof(N));
        out.write((char *)&B, sizeof(B));
        for(int i = 0; i < M; i++)
        {
            for(int j = 0; j < n; j++)
            {
                out.write((char *)&x[i][j], sizeof(int));
            }
        }
        for(int i = 0; i < n; i++)
        {
            size_t sz = r[i].size();
            out.write((char *)&sz, sizeof(sz));
            for(T v : r[i])
            {
                out.write((char *)&v, sizeof(T));
            }
        }
        for(int i = 0; i <= n; i++)
        {
            for(int j = 0; j < M; j++)
            {
                out.write((char *)&a[i][j], sizeof(int));
            }
        }
        return 0;
    }

    int load(std::string load_file)
    {
        ifstream in(load_file, ios::binary);
        if(in.fail()) return 1;

        if(!in.read((char *)&M, sizeof(M))) return 2;
        if(!in.read((char *)&N, sizeof(N))) return 2;
        if(M < 1 || N < 1) return 2;
        if(!in.read((char *)&B, sizeof(B))) return 2;
        if(B != sizeof(T) * 8) return 3;
        n = (N + B - 1) / B, minSiteL = B * 2 - 1;
        xp.resize(n), hz.resize(n + 1), up_end.resize(n + 1), dn_end.resize(n + 1), z_.resize(n), z.resize(n);
        xp[0] = 1;
        for(int i = 1; i < n; i++)
        {
            xp[i] = xp[i - 1] * BASE % MOD;
        }
        x.resize(M + 1, vector<int>(n));
        h.resize(M, vector<unsigned long long>(n + 1));
        for(int i = 0; i < M; i++)
        {
            for(int j = 0; j < n; j++)
            {
                if(!in.read((char *)&x[i][j], sizeof(int))) return 2;
                h[i][j + 1] = (h[i][j] + xp[j] * (x[i][j] + 1)) % MOD;
            }
        }
        r.resize(n);
        for(int i = 0; i < n; i++)
        {
            size_t sz;
            if(!in.read((char *)&sz, sizeof(sz)) || sz < 1) return 2;
            r[i].resize(sz);
            for(int j = 0; j < (int)sz; j++)
            {
                if(!in.read((char *)&r[i][j], sizeof(T))) return 2;
            }
        }
        a.resize(n + 1, vector<int>(M));
        for(int i = 0; i <= n; i++)
        {
            for(int j = 0; j < M; j++)
            {
                if(!in.read((char *)&a[i][j], sizeof(int)) || !(-1 < a[i][j] && a[i][j] < M)) return 2;
            }
        }

        return 0;
    }

    vector<int> randhapz()
    {
        std::random_device rd; // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> dist(0, 101);
        vector<int> rz(N);
        for(int i = 0; i < N; i++) rz[i] = dist(gen) > 95;
        return rz;
    }

    vector<T> encode_zg(const vector<int> & rz)
    {
        vector<T> zg(n);
        for(int K = 0; K < N; K++)
        {
            int k = K / B;
            zg[k] = (zg[k] << 1) | (rz[K] != 0);
        }
        zg[n - 1] <= n * B - N;
        return zg;
    }

    void query(std::string vcfquery,
               std::string samples,
               std::string region,
               int L,
               IntMapU & haplens,
               IntMapU & hapends)
    {
        vcfpp::BcfReader vcf(vcfquery, samples, region);
        vcfpp::BcfRecord var(vcf.header);
        vector<bool> gt;
        int Q = vcf.nsamples * 2;
        vector<vector<T>> Z(Q, vector<T>(n));
        for(int K = 0; K < N; K++)
        {
            vcf.getNextVariant(var);
            var.getGenotypes(gt);
            if(!var.isNoneMissing() || !var.allPhased())
                throw runtime_error("there is a site with non-phased or missing value!\n");
            int k = K / B;
            for(int i = 0; i < Q; i++) Z[i][k] = (Z[i][k] << 1) | (gt[i] != 0);
        }
        for(int i = 0; i < Q; i++) Z[i][n - 1] <= n * B - N; // pad last syllable with 0s
        {
            query_zg(Z[1], L, haplens, hapends);
        }
    }

    int query_zg(const vector<T> & zg, int L, IntMapU & haplens, IntMapU & hapends)
    { // 1-length sites
        if(L < minSiteL) L = minSiteL + 1;

        const int l = (L - B + 1) / B; // min # syllables that must be covered by an L-site match

        for(int k = 0; k < n; k++)
        {
            z[k] = lower_bound(r[k].begin(), r[k].end(), zg[k]) - r[k].begin();
            if(z[k] < (int)r[k].size() && zg[k] != r[k][z[k]]) z[k] = r[k].size();
            hz[k + 1] = (hz[k] + xp[k] * (z[k] + 1)) % MOD;
        }
        memset(&up_end[1], 0, n * sizeof(int));
        memset(&dn_end[1], 0, n * sizeof(int));

        for(int k = l - 1, kp1 = l, t = 0, up_on = 0, dn_on = 0, no_beg = 0; k < n; k = kp1 - 1)
        {
            if(z[k] == (int)r[k].size())
            {
                t = M;
                no_beg = k + 1;
            }
            else
            {
                int lo_i = 0, hi_i = M, tu_beg = k, td_beg = k;
                while(lo_i < hi_i)
                {
                    int g_i = (lo_i + hi_i) >> 1, xi = a[kp1][g_i];

                    if(z[k] < x[xi][k])
                    {
                        hi_i = g_i;
                        td_beg = k;
                    }
                    else if(z[k] > x[xi][k])
                    {
                        lo_i = g_i + 1;
                        tu_beg = k;
                    }
                    else
                    {
                        int lo_k;
                        if(k - no_beg < 10)
                        {
                            lo_k = k;
                            while(lo_k-- > no_beg && z[lo_k] == x[xi][lo_k])
                            {
                            }
                        }
                        else
                        {
                            lo_k = no_beg;
                            int hi_k = k;
                            unsigned long long rt =
                                hz[kp1] < h[xi][kp1] ? MOD - h[xi][kp1] + hz[kp1] : hz[kp1] - h[xi][kp1];
                            while(lo_k < hi_k)
                            {
                                int g_k = (lo_k + hi_k) >> 1;

                                if((hz[g_k] < h[xi][g_k] ? MOD - h[xi][g_k] + hz[g_k] : hz[g_k] - h[xi][g_k])
                                   == rt)
                                {
                                    hi_k = g_k;
                                }
                                else
                                {
                                    lo_k = g_k + 1;
                                }
                            }
                            lo_k--;
                        }

                        if(lo_k > -1 && z[lo_k] < x[xi][lo_k])
                        {
                            hi_i = g_i;
                            td_beg = lo_k;
                        }
                        else
                        {
                            lo_i = g_i + 1;
                            tu_beg = lo_k;
                        }
                    }
                }

                t = lo_i;
                no_beg = max(no_beg, min(tu_beg, td_beg) + 1);
            }

            up_on -= up_end[k];
            dn_on -= dn_end[k];
            if(kp1 < no_beg + l)
            {
                kp1 = no_beg + l;
                continue;
            }
            int s_idx = k - l;

            unsigned long long zh =
                hz[kp1] < hz[kp1 - l] ? MOD - hz[kp1 - l] + hz[kp1] : hz[kp1] - hz[kp1 - l];

            int p = t - up_on, lo_p = 0, hi_p = p;
            while(lo_p < hi_p)
            {
                int g_p = (lo_p + hi_p) >> 1, xi = a[kp1][g_p];
                if((h[xi][kp1] < h[xi][kp1 - l] ? MOD - h[xi][kp1 - l] + h[xi][kp1]
                                                : h[xi][kp1] - h[xi][kp1 - l])
                   == zh)
                {
                    hi_p = g_p;
                }
                else
                {
                    lo_p = g_p + 1;
                }
            }

            while(p > lo_p)
            {
                int xi = a[kp1][--p];

                int lo = kp1, pivot = min(n, lo + 10);
                while(lo < pivot && z[lo] == x[xi][lo]) lo++;
                if(lo == pivot)
                {
                    int hi = n;
                    unsigned long long lt = hz[k] < h[xi][k] ? MOD - h[xi][k] + hz[k] : hz[k] - h[xi][k];

                    while(lo < hi)
                    {
                        int g = (lo + hi + 1) >> 1;

                        if(lt == (hz[g] < h[xi][g] ? MOD - h[xi][g] + hz[g] : hz[g] - h[xi][g]))
                        {
                            lo = g;
                        }
                        else
                        {
                            hi = g - 1;
                        }
                    }
                }
                up_end[lo]++;

                int start, end;
                refine(s_idx, lo, xi, &start, &end);
                if(end - start >= L)
                {
                    // out << qID << '\t' << IDs[xi] << '\t' << start << '\t' << end << '\n';
                    if(haplens.count(n) == 0)
                    {
                        haplens[xi] = end - start;
                        hapends[xi] = end;
                    }
                    else if(end - start > haplens[n])
                    {
                        haplens[xi] = end - start;
                        hapends[xi] = end;
                    }
                }
            }
            up_on = t - p;

            p = t + dn_on, lo_p = p, hi_p = M;
            while(lo_p < hi_p)
            {
                int g_p = (lo_p + hi_p) >> 1, xi = a[kp1][g_p];
                if((h[xi][kp1] < h[xi][kp1 - l] ? MOD - h[xi][kp1 - l] + h[xi][kp1]
                                                : h[xi][kp1] - h[xi][kp1 - l])
                   == zh)
                {
                    lo_p = g_p + 1;
                }
                else
                {
                    hi_p = g_p;
                }
            }

            while(p < lo_p)
            {
                int xi = a[kp1][p++];

                int lo = kp1, pivot = min(n, lo + 10);
                while(lo < pivot && z[lo] == x[xi][lo]) lo++;
                if(lo == pivot)
                {
                    int hi = n;
                    unsigned long long lt = hz[k] < h[xi][k] ? MOD - h[xi][k] + hz[k] : hz[k] - h[xi][k];

                    while(lo < hi)
                    {
                        int g = (lo + hi + 1) >> 1;

                        if(lt == (hz[g] < h[xi][g] ? MOD - h[xi][g] + hz[g] : hz[g] - h[xi][g]))
                        {
                            lo = g;
                        }
                        else
                        {
                            hi = g - 1;
                        }
                    }
                }
                dn_end[lo]++;

                int start, end;
                refine(s_idx, lo, xi, &start, &end);
                if(end - start >= L)
                {
                    // out << qID << '\t' << IDs[xi] << '\t' << start << '\t' << end << '\n';
                    if(haplens.count(n) == 0)
                    {
                        haplens[xi] = end - start;
                        hapends[xi] = end;
                    }
                    else if(end - start > haplens[n])
                    {
                        haplens[xi] = end - start;
                        hapends[xi] = end;
                    }
                }
            }
            dn_on = p - t;

            kp1++;
        }

        return 0;
    }

    inline void refine(int s_idx, int e_idx, int xi, int * start, int * end);
};

template<>
inline void SyllablePBWT<unsigned long long>::refine(int s_idx, int e_idx, int xi, int * start, int * end)
{
    if(s_idx > -1)
    {
        *start = (s_idx + 1) * B - __builtin_ctzll(z_[s_idx] ^ r[s_idx][x[xi][s_idx]]);
    }
    else
        *start = 0;
    if(e_idx < n)
    {
        *end = e_idx * B + __builtin_clzll(z_[e_idx] ^ r[e_idx][x[xi][e_idx]]);
    }
    else
        *end = N;
}

template<>
inline void SyllablePBWT<unsigned __int128>::refine(int s_idx, int e_idx, int xi, int * start, int * end)
{
    if(s_idx > -1)
    {
        unsigned __int128 x_val = r[s_idx][x[xi][s_idx]];
        unsigned long long right = z_[s_idx] ^ x_val;
        if(right == 0)
        {
            *start = (s_idx + 1) * B - 64 - __builtin_ctzll((z_[s_idx] >> 64) ^ (x_val >> 64));
        }
        else
        {
            *start = (s_idx + 1) * B - __builtin_ctzll(right);
        }
    }
    else
        *start = 0;
    if(e_idx < n)
    {
        unsigned __int128 x_val = r[e_idx][x[xi][e_idx]];
        unsigned long long left = (z_[e_idx] >> 64) ^ (x_val >> 64);
        if(left == 0)
        {
            *end = e_idx * B + 64 + __builtin_clzll(z_[e_idx] ^ x_val);
        }
        else
        {
            *end = e_idx * B + __builtin_clzll(left);
        }
    }
    else
        *end = N;
}

#endif // SYLLABLEPBWT_H_
