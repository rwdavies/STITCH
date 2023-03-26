/*******************************************************************************
 * @file        mspbwt.h
 * @author      Zilong Li
 * @email       zilong.dk@gmail.com
 * @version     v0.1.0
 * @breif       Multiple Symbols Positional Burrows–Wheeler Transform (msPBWT)
 * Copyright (C) 2022. The use of this code is governed by the LICENSE file.
 ******************************************************************************/

#ifndef MSPBWT_H_
#define MSPBWT_H_

#include "vcfpp/vcfpp.h"

#include <algorithm>
#include <bitset>
#include <fstream>
#include <iterator>
#include <map>
#include <numeric>
#include <random>
#include <unordered_map>
#include <unordered_set>


using namespace std;
using namespace vcfpp;

using IntMapU = unordered_map<int, int>;
using IntVec1D = vector<int>;
using IntVec2D = vector<IntVec1D>;
using IntVec3D = vector<IntVec2D>;

template <typename T>
T reverseBits(T n, size_t B = sizeof(T) * 8)
{
    assert(B <= std::numeric_limits<T>::digits);
    T rv = 0;
    for (size_t i = 0; i < B; ++i, n >>= 1)
        rv = (rv << 1) | (n & 0x01);
    return rv;
}

// 0-based
inline vector<int> seq_by(int start, int end, int by)
{
    int n = (end - start + 1) % by == 0 ? (end - start + 1) / by : ((end - start + 1) / by + 1);
    vector<int> seq(n);
    int i{0}, x;
    for (x = start; x <= end; x = x + by)
        seq[i++] = x;
    return seq;
}

// C++11 compatible
template <typename T, typename = typename std::enable_if<std::is_unsigned<T>::value>::type>
class msPBWT
{
private:
    using grid_t = T;
    using GridVec = vector<grid_t>;
    using GridVec2D = vector<GridVec>;
    using GridVec3D = vector<GridVec2D>;
    using GridSetU = unordered_set<grid_t>;
    using GridMapU = unordered_map<grid_t, int>;                 // {symbol : index}
    using GridVecMapU = unordered_map<grid_t, std::vector<int>>; // {symbol : vec(index)}
    using SymbolIdxMap = map<int, int, less<int>>;               // {index: rank}
    using WgSymbolMap = map<grid_t, SymbolIdxMap, less<grid_t>>; // {symbol:{index:rank}}

    int B{sizeof(T) * 8}, N{0}, M{0}, G{0}, G1{0}, G2{0}, nindices{4};
    bool is_save_X{1}, is_save_D{1};
    vector<GridMapU> GMV; // for debug purpose
    GridVec2D X;          // Grids x Haps
    GridVec2D S;          // Grids x Sorted and Unique symbols
    IntVec3D W;
    IntVec2D C;
    IntVec2D A;          // nindices x Grids x Haps
    IntVec2D D;          // nindices x Grids x Haps
    vector<int> reorder; // (M) reorder SNPs or just subset SNPs

public:
    msPBWT(int nindices_ = 4) : nindices(nindices_)
    {
    }

    virtual ~msPBWT(){};

    bool verbose{0};
    bool debug{0};

    IntVec1D randhapz()
    {
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> dist(0, 101);
        IntVec1D z(M);
        for (int i = 0; i < M; i++)
        {
            z[i] = dist(gen) > 95;
        }
        return z;
    }

    GridVec encodezg(const vector<int>& z_)
    {
        assert(z_.size() == M);
        int k{0};
        GridVec zg(G);
        vector<bool> z(M);
        for (k = 0; k < M; k++)
            z[k] = z_[reorder[k]] != 0;
        k = 0;
        for (size_t m = 0; m < z.size(); m++)
        {
            zg[k] = (zg[k] << 1) | (z[m] != 0);
            if ((m + 1) % B == 0)
            {
                zg[k] = reverseBits(zg[k]);
                k++;
            }
        }
        if (G == k + 1)
        {
            zg[k] <<= G * B - M; // padding zeros
            zg[k] = reverseBits(zg[k]);
        }
        else if (G == k)
        {
            if (verbose)
                cerr << "no need padding\n";
        }
        else
        {
            throw std::runtime_error("something wrong\n");
        }
        return zg;
    }

    GridVec randzg()
    {
        auto z = randhapz();
        return encodezg(z);
    }

    GridMapU build_C(const GridVec& y, const GridVec& s)
    {
        int c{0};
        GridMapU C;
        for (const auto& si : s)
        {
            c = 0;
            for (const auto& xi : y)
            {
                if (xi < si)
                    c++;
            }
            C[si] = c;
        }
        return C;
    }

    WgSymbolMap save_W(const GridVec& y, const GridVec& s)
    {
        int k{0};
        WgSymbolMap W;
        for (const auto& si : s)
        {
            k = 0;
            for (size_t i = 0; i < y.size(); i++)
            {
                if (y[i] == si)
                    W[si][i] = k++; // k: how many occurrences before i
            }
        }

        return W;
    }

    IntVec1D save_C(const GridVec& y, const GridVec& s)
    {
        int c{0}, i{0};
        vector<int> C(s.size());
        for (const auto& si : s)
        {
            c = 0;
            for (const auto& yi : y)
            {
                if (yi < si)
                    c++;
            }
            C[i++] = c;
        }

        return C;
    }

    IntVec2D save_Occ(const GridVec& y, const GridVec& s)
    {
        vector<vector<int>> Occ(s.size());
        vector<int> idx;
        size_t i{0};
        int k{0};
        for (const auto& si : s)
        {
            for (i = 0; i < y.size(); i++)
            {
                if (y[i] == si)
                    idx.push_back(i);
            }
            Occ[k++] = idx;
            idx.clear();
        }

        return Occ;
    }

    GridVecMapU build_W(const GridVec& y, const GridVec& s)
    {
        int c{0};
        size_t i{0};
        auto C = build_C(y, s);
        GridVecMapU W;
        for (const auto& si : s)
        {
            W[si] = vector<int>(y.size());
            c = 0;
            for (i = 0; i < y.size(); i++)
            {
                if (y[i] == si)
                    c++;
                W[si][i] = c + C.at(si);
            }
        }

        return W;
    }

    void build(const std::string& vcfpanel, const std::string& samples, const std::string& region,
               double maf = 0)
    {
        int k{0}, m{0}, i{0};
        BcfReader vcf(vcfpanel, samples, region);
        BcfRecord var(vcf.header);
        N = vcf.nsamples * 2;
        M = 0;
        vector<bool> gt;
        vector<vector<bool>> allgts, gt_rares, gt_commons;
        vector<int> snp_rares, snp_commons;
        double af;
        while (vcf.getNextVariant(var))
        {
            var.getGenotypes(gt);
            if (!var.isNoneMissing() || !var.allPhased())
                continue;
            // keep track of snp index with AF < minaf
            af = 0;
            for (auto g : gt)
                af += g;
            af /= N;
            if (af < maf)
            {
                snp_rares.push_back(M);
                gt_rares.push_back(gt);
            }
            else
            {
                snp_commons.push_back(M);
                gt_commons.push_back(gt);
            }
            if ((M + 1) % B == 0)
            {
                reorder.insert(reorder.end(), snp_rares.begin(), snp_rares.end());
                snp_rares.clear();
                reorder.insert(reorder.end(), snp_commons.begin(), snp_commons.end());
                snp_commons.clear();
                allgts.insert(allgts.cend(), gt_rares.begin(), gt_rares.end());
                gt_rares.clear();
                allgts.insert(allgts.cend(), gt_commons.begin(), gt_commons.end());
                gt_commons.clear();
            }
            M++;
        }
        reorder.insert(reorder.end(), snp_rares.begin(), snp_rares.end());
        snp_rares.clear();
        reorder.insert(reorder.end(), snp_commons.begin(), snp_commons.end());
        snp_commons.clear();
        allgts.insert(allgts.cend(), gt_rares.begin(), gt_rares.end());
        gt_rares.clear();
        allgts.insert(allgts.cend(), gt_commons.begin(), gt_commons.end());
        gt_commons.clear();
        G = (M + B - 1) / B;
        if (verbose)
            cerr << "N:" << N << ",M:" << M << ",G:" << G << ",B:" << B << ",nindices:" << nindices << endl;
        X.resize(G, GridVec(N));
        k = 0;
        for (m = 0; m < M; m++)
        {
            for (i = 0; i < N; i++)
                X[k][i] = (X[k][i] << 1) | (allgts[m][i] != 0);
            if ((m + 1) % B == 0)
            {
                for (i = 0; i < N; i++)
                    X[k][i] = reverseBits(X[k][i]); // reverset bits
                k++;                                // update next grid
            }
        }
        if (G == k + 1)
        {
            for (i = 0; i < N; i++)
            {
                X[k][i] <<= G * B - M;
                X[k][i] = reverseBits(X[k][i]); // reverset bits
            }
        }
        else if (G == k)
        {
            if (verbose)
                cerr << "no need padding\n";
        }
        else
        {
            throw std::runtime_error("something wrong\n");
        }

        // building A and save indices now
        A.resize(G);
        S.resize(G);
        C.resize(G);
        W.resize(G);
        if (is_save_D)
            D.resize(G);
        // create some intermediate varibales
        vector<int> Occ;
        vector<vector<int>> kas; // for building A
        vector<vector<int>> kds; // for building D
        vector<int> sqp;         // for building D
        GridSetU symbols;        // keep track of unique symbols at k
        GridVec y1(N);
        vector<int> a0(N);     // initilize a[k]
        vector<int> d0(N + 1); // sentinels, d[0],d[N] = k + 1
        int Gi, ki, ni{-1};
        for (i = 0; i < nindices; i++)
        {
            auto Gv = seq_by(i, G - 1, nindices);
            Gi = Gv.size();
            // initial a and d at ki=-1
            std::iota(a0.begin(), a0.end(), 0);
            if (is_save_D)
            {
                std::fill(d0.begin(), d0.end(), 0);
                d0[0] = 1;
                d0[N] = 1;
            }
            for (ki = 0; ki < Gi; ki++)
            {
                ni++;
                k = Gv[ki];
                build_indices_k(ni, ki, X[k], y1, S[ni], C[ni], W[ni], Occ, kas, kds, sqp, a0, d0, symbols);
            }
        }
    }


    void build_indices_k(int ni, int ki, const GridVec& xk, GridVec& yk, GridVec& sk, vector<int>& ck,
                         vector<vector<int>>& wk, vector<int>& Occ, vector<vector<int>>& kas,
                         vector<vector<int>>& kds, vector<int>& sqp, vector<int>& a0, vector<int>& d0,
                         GridSetU& symbols)
    {
        int s, n, c, e, i;
        for (n = 0; n < N; n++)
        {
            yk[n] = xk[a0[n]];
            symbols.insert(yk[n]);
        }
        sk = GridVec(symbols.begin(), symbols.end());
        std::sort(sk.begin(), sk.end());
        symbols.clear();
        // C[k] = save_C(y1, S[k]);
        // W[k] = save_Occ(y1, S[k]);
        // auto Wg = build_W(y1, S[k]); // here Wg is S x N
        // for (n = 0; n < N; n++)
        //     A[k][Wg[y1[n]][n] - 1] = a0[n];
        // next run
        // a0 = A[k];

        int sn = sk.size();
        ck.resize(sn);
        wk.resize(sn);
        for (s = 0; s < sn; s++)
        {
            c = 0;
            for (n = 0; n < N; n++)
            {
                if (yk[n] < sk[s])
                    c++;
                if (yk[n] == sk[s])
                    Occ.push_back(n);
            }
            ck[s] = c;
            wk[s] = Occ;
            Occ.clear();
        }

        // save A and D
        kas.resize(sn);
        if (is_save_D)
        {
            kds.resize(sn);
            sqp.resize(sn, ki + 1);
        }
        // iterate over all haplotyps in reverese prefix sorted order
        for (n = 0; n < N; n++)
        {
            i = a0[n]; // index ihap at Y[k]
            if (is_save_D)
                e = d0[n];
            for (s = 0; s < sn; s++)
            {
                if (is_save_D)
                {
                    if (e > sqp[s])
                        sqp[s] = e;
                }
                if (xk[i] == sk[s])
                {
                    kas[s].push_back(i);
                    if (is_save_D)
                    {
                        kds[s].push_back(sqp[s]);
                        sqp[s] = 0;
                    }
                }
            }
        }
        // make sure A[k], D[k] are empty
        for (s = 0; s < sn; s++)
        {
            A[ni].insert(A[ni].end(), kas[s].begin(), kas[s].end());
            if (is_save_D)
                D[ni].insert(D[ni].end(), kds[s].begin(), kds[s].end());
        }
        kas.clear();
        a0 = A[ni];
        if (is_save_D)
        {
            kds.clear();
            sqp.clear();
            D[ni][N] = ki + 1; // add sentinel
            d0 = D[ni];
        }
    }

    int save(const std::string& filename)
    {
        ofstream out(filename, ios::out | ios::binary);
        if (out.fail())
            return 1;
        out.write((char*)&B, sizeof(B));
        out.write((char*)&M, sizeof(M));
        out.write((char*)&N, sizeof(N));
        out.write((char*)&G, sizeof(G));
        out.write((char*)&nindices, sizeof(nindices));
        // write reorder (M)
        int n, k, m;
        for (m = 0; m < M; m++)
            out.write((char*)&reorder[m], sizeof(int));
        if (is_save_X)
        {
            // write X
            for (k = 0; k < G; k++)
            {
                for (n = 0; n < N; n++)
                    out.write((char*)&X[k][n], sizeof(T));
            }
        }
        // write S
        for (k = 0; k < G; k++)
        {
            size_t sz = S[k].size();
            out.write((char*)&sz, sizeof(sz));
            for (size_t i = 0; i < sz; i++)
                out.write((char*)&S[k][i], sizeof(T));
        }
        // write A
        for (k = 0; k < G; k++)
        {
            for (n = 0; n < N; n++)
                out.write((char*)&A[k][n], sizeof(int));
        }
        if (is_save_D)
        {
            // write D
            for (k = 0; k < G; k++)
            {
                for (n = 0; n < N + 1; n++)
                    out.write((char*)&D[k][n], sizeof(int));
            }
        }
        // write C
        for (k = 0; k < G; k++)
        {
            size_t sz = C[k].size();
            out.write((char*)&sz, sizeof(sz));
            for (size_t i = 0; i < sz; i++)
                out.write((char*)&C[k][i], sizeof(int));
        }
        // write W
        for (k = 0; k < G; k++)
        {
            size_t sz = W[k].size();
            out.write((char*)&sz, sizeof(sz));
            for (size_t i = 0; i < sz; i++)
            {
                size_t sz1 = W[k][i].size();
                out.write((char*)&sz1, sizeof(sz1));
                for (size_t j = 0; j < sz1; j++)
                    out.write((char*)&W[k][i][j], sizeof(int));
            }
        }
        return 0;
    }

    int load(const std::string& filename)
    {
        ifstream in(filename, ios::in | ios::binary);
        if (in.fail())
            return 1;
        if (!in.read((char*)&B, sizeof(B)))
            return 2;
        if (B != sizeof(T) * 8)
            throw invalid_argument("the binary file may be created with different B!\n");
        if (!in.read((char*)&M, sizeof(M)))
            return 2;
        if (!in.read((char*)&N, sizeof(N)))
            return 2;
        if (!in.read((char*)&G, sizeof(G)))
            return 2;
        if (!in.read((char*)&nindices, sizeof(nindices)))
            return 2;
        cerr << "N: " << N << ",M: " << M << ",G: " << G << ",B: " << B << ",nindices " << nindices << endl;
        // read reorder (M)
        reorder.resize(M);
        int n, m, k;
        for (m = 0; m < M; m++)
            in.read((char*)&reorder[m], sizeof(int));
        if (is_save_X)
        {
            // read X
            X.resize(G, GridVec(N));
            for (k = 0; k < G; k++)
            {
                for (n = 0; n < N; n++)
                {
                    if (!in.read((char*)&X[k][n], sizeof(T)))
                        return 2;
                }
            }
            if (verbose)
                cerr << "load X done" << endl;
        }
        // read S
        S.resize(G);
        if (debug)
            GMV.resize(G);
        for (k = 0; k < G; k++)
        {
            size_t sz;
            if (!in.read((char*)&sz, sizeof(sz)) || sz < 1)
                return 2;
            S[k].resize(sz);
            for (size_t i = 0; i < sz; i++)
            {
                if (!in.read((char*)&S[k][i], sizeof(T)))
                    return 2;
                if (debug)
                    GMV[k][S[k][i]] = i;
            }
        }
        if (verbose)
            cerr << "load S done" << endl;
        // read A
        A.resize(G, vector<int>(N));
        for (k = 0; k < G; k++)
        {
            for (n = 0; n < N; n++)
            {
                if (!in.read((char*)&A[k][n], sizeof(int)))
                    return 2;
            }
        }
        if (verbose)
            cerr << "load A done" << endl;
        if (is_save_D)
        {
            // read D
            D.resize(G, vector<int>(N + 1));
            for (k = 0; k < G; k++)
            {
                for (n = 0; n < N + 1; n++)
                {
                    if (!in.read((char*)&D[k][n], sizeof(int)))
                        return 2;
                }
            }
            if (verbose)
                cerr << "load D done" << endl;
            // for (int i = 0; i < N + 1; i++)
            // {
            //     for (int j = 0; j < G; j++)
            //         cerr << D[j][i] << ",";
            //     cerr << endl;
            // }
        }
        // read C
        C.resize(G);
        for (k = 0; k < G; k++)
        {
            size_t sz;
            if (!in.read((char*)&sz, sizeof(sz)) || sz < 1)
                return 2;
            C[k].resize(sz);
            for (size_t i = 0; i < sz; i++)
            {
                if (!in.read((char*)&C[k][i], sizeof(int)))
                    return 2;
            }
        }
        if (verbose)
            cerr << "load C done" << endl;
        // read W
        W.resize(G);
        for (k = 0; k < G; k++)
        {
            size_t sz;
            if (!in.read((char*)&sz, sizeof(sz)) || sz < 1)
                return 2;
            W[k].resize(sz);
            for (size_t i = 0; i < sz; i++)
            {
                size_t sz1;
                if (!in.read((char*)&sz1, sizeof(sz1)) || sz1 < 1)
                    return 2;
                W[k][i].resize(sz1);
                for (size_t j = 0; j < sz1; j++)
                {
                    if (!in.read((char*)&W[k][i][j], sizeof(int)))
                        return 2;
                }
            }
        }
        if (verbose)
            cerr << "load W done" << endl;

        return 0;
    }

    void query(const std::string& vcfquery, const std::string& samples, const std::string& region,
               int viewk = 0)
    {
        BcfReader vcf(vcfquery, samples, region);
        BcfRecord var(vcf.header);
        vector<bool> gt;
        vector<GridVec> Z(vcf.nsamples * 2, GridVec(G));
        size_t i{0};
        int k{0}, m{0};
        while (vcf.getNextVariant(var))
        {
            var.getGenotypes(gt);
            if (!var.isNoneMissing() || !var.allPhased())
                continue;
            for (i = 0; i < gt.size(); i++)
                Z[i][k] = (Z[i][k] << 1) | (gt[i] != 0);
            m++;
            if (m % B == 0)
            {
                for (i = 0; i < gt.size(); i++)
                    Z[i][k] = reverseBits(Z[i][k]); // reverset bits
                k++;                                // update next grid
            }
        }
        assert(m == M);
        if (G == k + 1)
        {
            for (i = 0; i < gt.size(); i++)
            {
                Z[i][k] <<= G * B - M;
                Z[i][k] = reverseBits(Z[i][k]); // reverset bits
            }
        }
        else if (G == k)
        {
            if (verbose)
                cerr << "no need padding\n";
        }
        else
        {
            throw std::runtime_error("something wrong\n");
        }
        IntMapU haplens, hapends, hapnindicies;
        report_setmaximal(haplens, hapends, hapnindicies, Z[1], viewk);
    }

    void real_setmaximal(int iind, IntMapU& haplens, IntMapU& hapends, IntMapU& hapnindicies, GridVec& zg,
                         const IntVec1D& gv, int step, int viewk = -1)
    {
        int klen, ks, k, s, n, i, j, count, e, f, g, e1, f1, g1, valid_grid_start{0};
        bool matches_lower, matches_upper, re_count{0};
        IntVec1D ghost;
        // bool first_valid_grid_start = true;
        for (k = 0; k < gv.size(); k++)
        {
            ks = k + step;
            auto kzs = std::lower_bound(S[ks].begin(), S[ks].end(), zg[gv[k]]);
            s = std::fmin(std::distance(S[ks].begin(), kzs), S[ks].size() - 1);
            if (S[ks][s] != zg[gv[k]])
            {
                ghost.push_back(k);
                if (verbose)
                    cerr << "ghost symbol at k: " << k << endl;
                // // what to do if having ghost symbols
                // if (kzs == S[ks].end() || s == 0)
                //     zg[gv[k]] = S[ks][s]; // already lower bound!
                // else
                //     zg[gv[k]] = S[ks][s - 1]; // corce to be lower bound!
            }

            if (k == valid_grid_start)
            {
                f1 = C[ks][s];
                g1 = C[ks][s] + W[ks][s].size(); // could be N
                e1 = valid_grid_start;
            }
            else if (k > valid_grid_start)
            {
                // assume symbol exists
                // g1 >= f1 >= C[ks][s] if g >= f
                auto fzk = std::lower_bound(W[ks][s].begin(), W[ks][s].end(), f);
                auto gzk = std::lower_bound(W[ks][s].begin(), W[ks][s].end(), g);
                f1 = C[ks][s] + std::fmin(std::distance(W[ks][s].begin(), fzk), W[ks][s].size());
                g1 = C[ks][s] + std::fmin(std::distance(W[ks][s].begin(), gzk), W[ks][s].size());
                if (f1 == g1)
                {

                    if (debug && viewk == ks)
                    {
                        for (i = 0; i < N; i++)
                        {
                            if (i == f)
                            {
                                cout << "zg is below" << endl;
                                for (int j = 0; j < k; j++)
                                    if (GMV[j + step].count(zg[gv[j]]))
                                        cout << GMV[j + step][zg[gv[j]]] << " ";
                                    else
                                        cout << "-1 ";
                                cout << endl;
                                cout << "f is below - " << f << endl;
                            }
                            for (int j = 0; j < k; j++)
                                cout << GMV[j + step][X[j][A[ks - 1][i]]] << " ";
                            cout << endl;
                            if (i == g - 1)
                                cout << "g - 1 is above - " << g - 1 << endl;
                            if (i == g)
                            {
                                cout << "g is above - " << g << endl;
                                for (int j = 0; j < k; j++)
                                    if (GMV[j + step].count(zg[gv[j]]))
                                        cout << GMV[j + step][zg[gv[j]]] << " ";
                                    else
                                        cout << "-1 ";
                                cout << endl;
                                cout << "zg is above" << endl;
                            }
                        }
                    }

                    // report matches from e to k for [f, g)
                    for (i = f; i < g; i++)
                    {
                        n = A[ks - 1][i];
                        klen = k - e; // k-1-e+1
                        // okay if ghost symbols in [e, k)
                        re_count = false;
                        for (auto it = ghost.rbegin(); it != ghost.rend(); ++it)
                        {
                            if (*it < e)
                                break;
                            if (*it >= e && *it < k)
                            {
                                re_count = true;
                                break;
                            }
                        }
                        if (re_count)
                        {
                            j = 0, klen = 0;
                            while (k >= ++j && X[gv[k - j]][n] == zg[gv[k - j]])
                                klen++;
                        }
                        // check if all equals
                        if (debug)
                        {
                            j = 0, count = 0;
                            while (k >= ++j && X[gv[k - j]][n] == zg[gv[k - j]])
                                count++;
                            cerr << k << "," << i << "," << klen << "," << count << "," << n << endl;
                        }
                        if (haplens.count(n) == 0)
                        {
                            haplens[n] = klen;
                            hapends[n] = k - 1;
                            hapnindicies[n] = 1;
                        }
                        else if (klen > haplens[n])
                        {
                            haplens[n] = klen;
                            hapends[n] = k - 1;
                        }
                        if (hapnindicies.count(n))
                            hapnindicies[n] = iind >= hapnindicies[n] ? (iind + 1) : hapnindicies[n];
                    }
                    // finding new e1, f1, g1
                    e1 = D[ks][f1] - 1;     // y[f1] and y[f1-1] diverge here, so upper bound for e
                    if (f1 == N && e1 == k) // recall sentinels d0[N] = k + 1;
                        e1 = k - 1;

                    // cerr << X[gv[e1]][A[ks][f1]] << "," << X[gv[e1]][A[ks][f1 - 1]] << endl;
                    // cerr << std::bitset<sizeof(T) * 8>(reverseBits(X[gv[e1]][A[ks][f1 - 1]])) << ","
                    //      << std::bitset<sizeof(T) * 8>(reverseBits(X[gv[e1]][A[ks][f1]])) << endl;

                    matches_lower = false;
                    matches_upper = false;
                    // if there is a ghost symbol in zg, e1++ will continue forever so add e1 < k
                    while ((e1 < k) && (!matches_lower) && (!matches_upper))
                    {
                        if (f1 > 0)
                            matches_upper = (zg[gv[e1]] == X[gv[e1]][A[ks][f1 - 1]]);
                        else
                            matches_upper = false;
                        if (f1 < N)
                            matches_lower = (zg[gv[e1]] == X[gv[e1]][A[ks][f1]]);
                        else
                            matches_lower = false;
                        // if matches neither y[f1] or y[f1-1], eg. symbol missing or just happens, e1++
                        if ((!matches_lower) && (!matches_upper))
                            ++e1;
                    }

                    // if ghost symbols exists and loop --e1 stops as it is
                    // this will shorten the potential longest matches
                    if (matches_upper)
                    {
                        --f1;
                        // make sure e1 > 0
                        while (e1 > valid_grid_start && zg[gv[e1 - 1]] == X[gv[e1 - 1]][A[ks][f1]])
                            --e1;
                        // we can't skip missing symbols otherwise here is not true
                        while (f1 > 0 && D[ks][f1] <= e1)
                            --f1;
                    }
                    if (matches_lower)
                    {
                        ++g1;
                        // make sure e1 > 0
                        while (e1 > valid_grid_start && zg[gv[e1 - 1]] == X[gv[e1 - 1]][A[ks][f1]])
                            --e1;
                        while (g1 < N && D[ks][g1] <= e1)
                            ++g1;
                    }

                    if (debug && (matches_lower) && (matches_upper))
                    {
                        cerr << "matches both upper and lower: " << e1 << "," << f1 << "," << g1 << endl;
                        for (i = 0; i < N; i++)
                        {
                            if (i == f1)
                            {
                                cout << "zg is below" << endl;
                                for (int j = 0; j < k; j++)
                                    if (GMV[j + step].count(zg[gv[j]]))
                                        cout << GMV[j + step][zg[gv[j]]] << " ";
                                    else
                                        cout << "-1 ";
                                cout << endl;
                                cout << "f1 is below - " << f << endl;
                            }
                            for (int j = 0; j < k; j++)
                                cout << GMV[j + step][X[j][A[ks][i]]] << " ";
                            cout << endl;
                            if (i == g1 - 1)
                                cout << "g1 - 1 is above - " << g1 - 1 << endl;
                            if (i == g1)
                            {
                                cout << "g1 is above - " << g1 << endl;
                                for (int j = 0; j < k; j++)
                                    if (GMV[j + step].count(zg[gv[j]]))
                                        cout << GMV[j + step][zg[gv[j]]] << " ";
                                    else
                                        cout << "-1 ";
                                cout << endl;
                                cout << "zg is above" << endl;
                            }
                        }
                    }
                }
                else
                {
                    assert(f1 < g1);
                    // nothing to do
                    e1 = e;
                }
            }
            // update for ext run
            f = f1;
            g = g1;
            e = e1;
        }

        // final report matches from e to k for [f, g)
        assert(ks == gv.size() - 1 + step);
        for (i = f; i < g; i++)
        {
            n = A[ks][i];
            klen = k - e;
            // okay if ghost symbols in [e, k)
            re_count = false;
            for (auto it = ghost.rbegin(); it != ghost.rend(); ++it)
            {
                if (*it < e)
                    break;
                if (*it >= e && *it < k)
                {
                    re_count = true;
                    break;
                }
            }
            if (re_count)
            {
                j = 0, klen = 0;
                while (k >= ++j && X[gv[k - j]][n] == zg[gv[k - j]])
                    klen++;
            }
            if (debug)
            {
                j = 0, count = 0;
                while (k >= ++j && X[gv[k - j]][n] == zg[gv[k - j]])
                    count++;
                cerr << k << "," << i << "," << klen << "," << count << "," << n << endl;
            }
            if (haplens.count(n) == 0)
            {
                haplens[n] = klen;
                hapends[n] = k - 1;
                hapnindicies[n] = 1;
            }
            else if (klen > haplens[n])
            {
                haplens[n] = klen;
                hapends[n] = k - 1;
            }
            if (hapnindicies.count(n))
                hapnindicies[n] = iind >= hapnindicies[n] ? (iind + 1) : hapnindicies[n];
        }
    }

    void impute_setmaximal(int iind, IntMapU& haplens, IntMapU& hapends, IntMapU& hapnindicies, GridVec& zg,
                         const IntVec1D& gv, int step)
    {
        int klen, ks, k, s, n, i, j, e, f, g, e1, f1, g1, valid_grid_start{0};
        bool re_count;
        IntVec1D ghost;
        for (k = 0; k < gv.size(); k++)
        {
            ks = k + step;
            auto kzs = std::lower_bound(S[ks].begin(), S[ks].end(), zg[gv[k]]);
            s = std::fmin(std::distance(S[ks].begin(), kzs), S[ks].size() - 1);
            if (S[ks][s] != zg[gv[k]])
            {
                ghost.push_back(k);
                if (verbose)
                    cerr << "ghost symbol at k: " << k << endl;
            }

            if (k == valid_grid_start)
            {
                f1 = C[ks][s];
                g1 = C[ks][s] + W[ks][s].size(); // could be N
                e1 = valid_grid_start;
            }
            else if (k > valid_grid_start)
            {
                // assume symbol exists
                // g1 >= f1 >= C[ks][s] if g >= f
                auto fzk = std::lower_bound(W[ks][s].begin(), W[ks][s].end(), f);
                auto gzk = std::lower_bound(W[ks][s].begin(), W[ks][s].end(), g);
                f1 = C[ks][s] + std::fmin(std::distance(W[ks][s].begin(), fzk), W[ks][s].size());
                g1 = C[ks][s] + std::fmin(std::distance(W[ks][s].begin(), gzk), W[ks][s].size());
                if (f1 == g1)
                {
                    // report matches from e to k for [f, g)
                    for (i = f; i < g; i++)
                    {
                        n = A[ks - 1][i];
                        klen = k - e; // k-1-e+1
                        // okay if ghost symbols in [e, k)
                        re_count = false;
                        for (auto it = ghost.rbegin(); it != ghost.rend(); ++it)
                        {
                            if (*it < e)
                                break;
                            if (*it >= e && *it < k)
                            {
                                re_count = true;
                                break;
                            }
                        }
                        if (re_count)
                        {
                            j = 0, klen = 0;
                            while (k >= ++j && X[gv[k - j]][n] == zg[gv[k - j]])
                                klen++;
                        }
                        if (haplens.count(n) == 0)
                        {
                            haplens[n] = klen;
                            hapends[n] = k - 1;
                            hapnindicies[n] = 1;
                        }
                        else if (klen > haplens[n])
                        {
                            haplens[n] = klen;
                            hapends[n] = k - 1;
                        }
                        if (hapnindicies.count(n))
                            hapnindicies[n] = iind >= hapnindicies[n] ? (iind + 1) : hapnindicies[n];
                    }
                    // finding new e1, f1, g1
                    e1 = D[ks][f1] - 1;     // y[f1] and y[f1-1] diverge here, so upper bound for e
                    if (f1 == N && e1 == k) // recall sentinels d0[N] = k + 1;
                        e1 = k - 1;

                    // if ghost symbols exists and loop --e1 stops as it is
                    // this will shorten the potential longest matches
                    // follow standard process as Durbin
                    if (zg[gv[e1]] <= X[gv[e1]][A[ks][f1 - 1]] && f1 > 0)
                    {
                        --f1;
                        // make sure e1 > 0
                        while (e1 > valid_grid_start && zg[gv[e1 - 1]] == X[gv[e1 - 1]][A[ks][f1]])
                            --e1;
                        // we can't skip missing symbols otherwise here is not true
                        while (f1 > 0 && D[ks][f1] <= e1)
                            --f1;
                    }
                    else
                    {
                        ++g1;
                        // make sure e1 > 0
                        while (e1 > valid_grid_start && zg[gv[e1 - 1]] == X[gv[e1 - 1]][A[ks][f1]])
                            --e1;
                        while (g1 < N && D[ks][g1] <= e1)
                            ++g1;
                    }

                }
                else
                {
                    assert(f1 < g1);
                    // nothing to do
                    e1 = e;
                }
            }
            // update for ext run
            f = f1;
            g = g1;
            e = e1;
        }

        // final report matches from e to k for [f, g)
        assert(ks == gv.size() - 1 + step);
        for (i = f; i < g; i++)
        {
            n = A[ks][i];
            klen = k - e;
            // okay if ghost symbols in [e, k)
            re_count = false;
            for (auto it = ghost.rbegin(); it != ghost.rend(); ++it)
            {
                if (*it < e)
                    break;
                if (*it >= e && *it < k)
                {
                    re_count = true;
                    break;
                }
            }
            if (re_count)
            {
                j = 0, klen = 0;
                while (k >= ++j && X[gv[k - j]][n] == zg[gv[k - j]])
                    klen++;
            }
            if (haplens.count(n) == 0)
            {
                haplens[n] = klen;
                hapends[n] = k - 1;
                hapnindicies[n] = 1;
            }
            else if (klen > haplens[n])
            {
                haplens[n] = klen;
                hapends[n] = k - 1;
            }
            if (hapnindicies.count(n))
                hapnindicies[n] = iind >= hapnindicies[n] ? (iind + 1) : hapnindicies[n];
        }
    }

    void report_setmaximal(IntMapU& haplens, IntMapU& hapends, IntMapU& hapnindicies, GridVec& zg,
                           int viewk = -1)
    {
        int iind, step{0};
        for (iind = 0; iind < nindices; iind++)
        {
            auto gv = seq_by(iind, G - 1, nindices);
            real_setmaximal(iind, haplens, hapends, hapnindicies, zg, gv, step);
            step += gv.size();
        }
    }

    void report_neighourings(IntMapU& haplens, IntMapU& hapends, IntMapU& hapnindicies, const GridVec& zg,
                             int L = 32)
    {
        int k, s, klen, j, n, l, Gi, ki, iind, ni{-1};
        for (iind = 0; iind < nindices; iind++)
        {
            auto Gv = seq_by(iind, G - 1, nindices);
            Gi = Gv.size();
            int zak_prev, zak_curr, valid_grid_start = 0;
            bool first_valid_grid_start = true;
            for (ki = 0; ki < Gi; ki++)
            {
                k = Gv[ki];
                ni++;
                auto kzs = std::upper_bound(S[ni].begin(), S[ni].end(), zg[k]);
                kzs = kzs == S[ni].begin() ? S[ni].begin() : std::prev(kzs);
                s = std::distance(S[ni].begin(), kzs);
                zak_curr = C[ni][s];

                // if (zg[k] == S[ni][s])
                // {
                //     if (first_valid_grid_start)
                //         valid_grid_start = ki;
                //     first_valid_grid_start = false;
                // }
                // else
                // {
                //     if (verbose)
                //         cerr << "skip: " << ki << endl;
                //     // if zg[k] symbol not exists, skip this grid and start over.
                //     first_valid_grid_start = true;
                //     continue;
                // }

                if (ki > valid_grid_start)
                {
                    if (zak_prev >= zak_curr)
                    {
                        auto kzi = std::upper_bound(W[ni][s].begin(), W[ni][s].end(), zak_prev);
                        kzi = kzi == W[ni][s].begin() ? W[ni][s].begin() : std::prev(kzi);
                        zak_curr += std::distance(W[ni][s].begin(), kzi);
                    }

                    // if ( zg[k] != S[ni][s] ) {
                    //     // re-confirm new position with longest matches if symbol mis-match
                    //     // should work with mis-match symbol at certain grid backwards?
                    //     int i = 2;
                    //     while (ki >= i && X[Gv[ki - i + 1]][A[ni][zak_curr]] == zg[Gv[ki - i + 1]] &&
                    //            X[Gv[ki - i]][A[ni][zak_curr]] < zg[Gv[ki - i]])
                    //     {
                    //         zak_curr++;
                    //         if (X[Gv[ki - i]][A[ni][zak_curr]] == zg[Gv[ki - i]])
                    //             i++;
                    //         if (zak_curr == N)
                    //             break;
                    //     }
                    // }


                    for (l = 0; l < L; l++)
                    {
                        n = A[ni][std::fmax(zak_curr - l - 1, 0)];
                        klen = 0;
                        for (j = ki; j >= 0; j--)
                        {
                            if (X[Gv[j]][n] == zg[Gv[j]])
                            {
                                klen++;
                            }
                            else
                            {
                                // TODO : update this with combining all nindicies
                                if (haplens.count(n) == 0)
                                {
                                    haplens[n] = klen;
                                    hapends[n] = ki;
                                    hapnindicies[n] = 1;
                                }
                                else if (klen >= haplens[n])
                                {
                                    haplens[n] = klen;
                                    hapends[n] = ki;
                                }
                                if (hapnindicies.count(n))
                                    hapnindicies[n] = iind >= hapnindicies[n] ? (iind + 1) : hapnindicies[n];
                                break;
                            }
                        }
                        if (zak_curr - l - 1 == 0)
                            break;
                    }

                    for (l = 0; l < L; l++)
                    {
                        n = A[ni][std::fmin(zak_curr + l, N - 1)];
                        klen = 0;
                        for (j = ki; j >= 0; j--)
                        {
                            if (X[Gv[j]][n] == zg[Gv[j]])
                            {
                                klen++;
                            }
                            else
                            {
                                if (haplens.count(n) == 0)
                                {
                                    haplens[n] = klen;
                                    hapends[n] = ki;
                                    hapnindicies[n] = 1;
                                }
                                else if (klen >= haplens[n])
                                {
                                    haplens[n] = klen;
                                    hapends[n] = ki;
                                }
                                if (hapnindicies.count(n))
                                    hapnindicies[n] = iind >= hapnindicies[n] ? (iind + 1) : hapnindicies[n];
                                break;
                            }
                        }
                        if (zak_curr + l == N - 1)
                            break;
                    }
                }
                zak_prev = zak_curr;
            }
        }
    }

    IntVec1D insert_and_match(const GridVec& zg)
    {
        vector<int> za(G);
        size_t i, s;
        int k;
        for (k = 0; k < G; k++)
        {
            auto kzs = std::lower_bound(S[k].begin(), S[k].end(), zg[k]);
            s = std::fmin(std::distance(S[k].begin(), kzs), S[k].size() - 1);
            za[k] = C[k][s];
            if (k > 0)
            {
                if (za[k - 1] >= za[k])
                {
                    auto kzi = std::lower_bound(W[k][s].begin(), W[k][s].end(), za[k - 1]);
                    za[k] += std::fmin(std::distance(W[k][s].begin(), kzi), W[k][s].size() - 1);
                }
                i = 1;
                while (k >= i && X[k - i + 1][A[k][za[k]]] == zg[k - i + 1] &&
                       X[k - i][A[k][za[k]]] < zg[k - i])
                {
                    za[k]++;
                    if (X[k - i][A[k][za[k]]] == zg[k - i])
                        i++;
                }
                cerr << k << "," << i << endl;
            }
        }
        return za;
    }

    IntVec1D insert_at_last(const GridVec& zg, bool start_with_last = false)
    {
        IntVec1D za(G);
        int k, s;
        for (k = 0; k < G; k++)
        {
            auto kzs = std::lower_bound(S[k].begin(), S[k].end(), zg[k]);
            s = std::fmin(std::distance(S[k].begin(), kzs), S[k].size() - 1);
            if (start_with_last && k == 0)
                za[k] = C[k][s] + W[k][s].size() - 1;
            else
                za[k] = C[k][s] + 1;
            if (k > 0)
            {
                if (za[k - 1] >= za[k])
                {
                    auto kzi = std::lower_bound(W[k][s].begin(), W[k][s].end(), za[k - 1]);
                    za[k] += std::fmin(std::distance(W[k][s].begin(), kzi), W[k][s].size() - 1);
                }
            }
        }
        return za;
    }

    IntVec1D insert(const GridVec& zg)
    {
        IntVec1D za(G);
        int k, s;
        for (k = 0; k < G; k++)
        {
            auto kzs = std::lower_bound(S[k].begin(), S[k].end(), zg[k]);
            s = std::fmin(std::distance(S[k].begin(), kzs), S[k].size() - 1);
            za[k] = C[k][s];
            if (k > 0)
            {
                if (za[k - 1] >= za[k])
                {
                    auto kzi = std::lower_bound(W[k][s].begin(), W[k][s].end(), za[k - 1]);
                    za[k] += std::fmin(std::distance(W[k][s].begin(), kzi), W[k][s].size() - 1);
                }
            }
        }
        return za;
    }

    void view_panel(int k, bool bit = true)
    {
        for (size_t i = 0; i < X[0].size(); i++)
        {
            for (int j = 0; j <= k + 1; j++)
            {
                if (bit)
                {
                    auto rb = reverseBits(X[j][A[k][i]]);
                    cout << std::bitset<sizeof(T) * 8>(rb) << " ";
                }
                else
                {
                    cout << X[j][A[k][i]] << " ";
                }
            }
            cout << endl;
        }
    }

    void view_zg(const GridVec& zg, int k, bool bit = true, int L = 0)
    {
        auto za1 = insert(zg);
        auto za2 = insert_and_match(zg);
        for (size_t i = 0; i < X[0].size(); i++)
        {
            if (((i < za1[k] - L) || (i > za1[k] + L)) && L != 0)
                continue;
            for (int j = 0; j <= k + 2; j++)
            {
                if (bit)
                {
                    auto rb = reverseBits(X[j][A[k][i]]);
                    cout << std::bitset<sizeof(T) * 8>(rb) << " ";
                }
                else
                {
                    cout << X[j][A[k][i]] << " ";
                }
            }
            cout << endl;
            if (i == za1[k])
            {
                // print out original Z bits
                cout << "========= fisrt zg is inserting here ========  k=" << k << ", za1[k]=" << za1[k]
                     << endl;
                for (int j = 0; j <= k + 2; j++)
                {
                    if (bit)
                    {
                        auto rb = reverseBits(zg[j]);
                        cout << std::bitset<sizeof(T) * 8>(rb) << " ";
                    }
                    else
                    {
                        cout << zg[j] << " ";
                    }
                }
                cout << endl;
                cout << "========= fisrt zg is inserting here ========  k=" << k << ", za1[k]=" << za1[k]
                     << endl;
            }
            if (i == za2[k])
            {
                // print out original Z bits
                cout << "========= last zg is inserting here ========  k=" << k << ", za2[k]=" << za2[k]
                     << endl;
                for (int j = 0; j <= k + 2; j++)
                {
                    if (bit)
                    {
                        auto rb = reverseBits(zg[j]);
                        cout << std::bitset<sizeof(T) * 8>(rb) << " ";
                    }
                    else
                    {
                        cout << zg[j] << " ";
                    }
                }
                cout << endl;
                cout << "========= last zg is inserting here ========  k=" << k << ", za2[k]=" << za2[k]
                     << endl;
            }
        }
    }
};


#endif // MSPBWT_H_
