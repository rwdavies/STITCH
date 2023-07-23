/*******************************************************************************
 * @file        mspbwt.h
 * @author      Zilong Li
 * @email       zilong.dk@gmail.com
 * @version     v0.1.0
 * @breif       Multiple Symbols Positional Burrows–Wheeler Transform (MSPBWT)
 * Copyright (C) 2022. The use of this code is governed by the LICENSE file.
 ******************************************************************************/

#ifndef MSPBWT_H_
#define MSPBWT_H_

#include "vcfpp/vcfpp.h"
#include <algorithm>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

using IntMapU = std::unordered_map<int, int>;
using Int1D = std::vector<int>;
using Int2D = std::vector<Int1D>;
using Int3D = std::vector<Int2D>;
using Bool1D = std::vector<bool>;
using Bool2D = std::vector<Bool1D>;
using String1D = std::vector<std::string>;

const int MAXSYMBOLS = 1024;

template<typename T>
T reverseBits(T n, size_t B = sizeof(T) * 8)
{
    assert(B <= std::numeric_limits<T>::digits);
    T rv = 0;
    for(size_t i = 0; i < B; ++i, n >>= 1) rv = (rv << 1) | (n & 0x01);
    return rv;
}

// 0-based
inline Int1D seq_by(int start, int end, int by)
{
    int n = (end - start + 1) % by == 0 ? (end - start + 1) / by : ((end - start + 1) / by + 1);
    Int1D seq(n);
    int i{0}, x;
    for(x = start; x <= end; x = x + by) seq[i++] = x;
    return seq;
}

// 0-based
inline Int1D seq_by2(int ichunk, int end, int nchunks)
{
    int end_one = end + 1;
    int isize = end_one / nchunks; // size of chunks except the last chunk
    int lsize = end_one - (nchunks - 1) * isize; // size of last chunk
    Int1D seq;
    if(ichunk == nchunks - 1)
    {
        seq.resize(lsize);
        int j{0}, i;
        for(i = isize * ichunk; i < end_one; i++) seq[j++] = i;
        assert(j == lsize);
    }
    else
    {
        seq.resize(isize);
        int j{0}, i;
        for(i = isize * ichunk; i < isize * (1 + ichunk); i++) seq[j++] = i;
        assert(j == isize);
    }
    return seq;
}

struct QUILT_RHB
{
    Int2D rhb_t; // nhaps x ngrids for quilt
    Int2D rare_per_hap_info; // store rare SNP index for each hap, 1-based
    Int1D pos, ac;
    String1D ref, alt;
    Bool1D snp_is_common;
    int n_skipped = 0, nhaps, nsnps, nGrids;
};

// C++11 compatible
template<typename T, typename = typename std::enable_if<std::is_unsigned<T>::value>::type>
class MSPBWT
{
  private:
    using grid_t = T;
    using GridVec = std::vector<grid_t>;
    using GridVec2D = std::vector<GridVec>;
    using GridVec3D = std::vector<GridVec2D>;
    using GridSetU = std::unordered_set<grid_t>;
    using GridMapU = std::unordered_map<grid_t, int>; // {symbol : index}
    using GridVecMapU = std::unordered_map<grid_t, std::vector<int>>; // {symbol : vec(index)}
    using SymbolIdxMap = std::map<int, int, std::less<int>>; // {index: rank}
    using WgSymbolMap = std::map<grid_t, SymbolIdxMap, std::less<grid_t>>; // {symbol:{index:rank}}

    int B{sizeof(T) * 8}, N{0}, M{0}, Mtotal{0}, G{0}, nindices{4};
    bool is_save_X{1}, is_save_D{0};
    GridVec2D X; // Grids x Haps
    GridVec2D S; // Grids x Sorted and Unique symbols
    Int3D W;
    Int2D C;
    Int2D A; // nindices x Grids x Haps
    Int2D D; // nindices x Grids x Haps
    Int1D keep; // keep only index of common variants
    Int1D pos;

  public:
    MSPBWT(int nindices_ = 4) : nindices(nindices_) {}

    virtual ~MSPBWT(){};

    bool verbose{0};
    bool debug{0};
    bool is_quilt_rhb{0};

    QUILT_RHB quilt;

    Int1D randhapz()
    {
        std::random_device rd; // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> dist(0, 101);
        Int1D z(M);
        for(int i = 0; i < M; i++) z[i] = dist(gen) > 95;
        return z;
    }

    int count1bits(T n)
    {
        return std::bitset<sizeof(T)>(n).count();
    }

    GridVec encodezg(const Int1D & z)
    {
        bool is_z_commom;
        if(z.size() == M)
            is_z_commom = true;
        else if(z.size() > M)
            is_z_commom = false;
        else
            throw std::runtime_error("something wrong! can't encode z as a grid!");

        int k{0}, m{0};
        GridVec zg(G);
        for(m = 0; m < M; m++)
        {
            if(is_z_commom)
                zg[k] = (zg[k] << 1) | (z[m] != 0);
            else
                zg[k] = (zg[k] << 1) | (z[keep[m]] != 0);
            if((m + 1) % B == 0)
            {
                zg[k] = reverseBits(zg[k]);
                k++;
            }
        }
        if(G == k + 1)
        {
            zg[k] <<= G * B - M; // padding zeros
            zg[k] = reverseBits(zg[k]);
        }
        else if(G == k)
        {
            if(verbose) std::cerr << "no need padding\n";
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

    GridMapU build_C(const GridVec & y, const GridVec & s)
    {
        int c{0};
        GridMapU C;
        for(const auto & si : s)
        {
            c = 0;
            for(const auto & xi : y)
            {
                if(xi < si) c++;
            }
            C[si] = c;
        }
        return C;
    }

    WgSymbolMap save_W(const GridVec & y, const GridVec & s)
    {
        int k{0};
        WgSymbolMap W;
        for(const auto & si : s)
        {
            k = 0;
            for(size_t i = 0; i < y.size(); i++)
            {
                if(y[i] == si) W[si][i] = k++; // k: how many occurrences before i
            }
        }

        return W;
    }

    Int1D save_C(const GridVec & y, const GridVec & s)
    {
        int c{0}, i{0};
        Int1D C(s.size());
        for(const auto & si : s)
        {
            c = 0;
            for(const auto & yi : y)
            {
                if(yi < si) c++;
            }
            C[i++] = c;
        }

        return C;
    }

    Int2D save_Occ(const GridVec & y, const GridVec & s)
    {
        Int2D Occ(s.size());
        Int1D idx;
        size_t i{0};
        int k{0};
        for(const auto & si : s)
        {
            for(i = 0; i < y.size(); i++)
            {
                if(y[i] == si) idx.push_back(i);
            }
            Occ[k++] = idx;
            idx.clear();
        }

        return Occ;
    }

    GridVecMapU build_W(const GridVec & y, const GridVec & s)
    {
        int c{0};
        size_t i{0};
        auto C = build_C(y, s);
        GridVecMapU W;
        for(const auto & si : s)
        {
            W[si] = Int1D(y.size());
            c = 0;
            for(i = 0; i < y.size(); i++)
            {
                if(y[i] == si) c++;
                W[si][i] = c + C.at(si);
            }
        }

        return W;
    }

    void build(const std::string & vcfpanel, const std::string & samples, const std::string & region, double maf = 0)
    {
        int k{0}, m{0}, i{0};
        vcfpp::BcfReader vcf(vcfpanel, samples, region);
        vcfpp::BcfRecord var(vcf.header);
        N = vcf.nsamples * 2, Mtotal = 0;
        if(is_quilt_rhb) quilt.rare_per_hap_info.resize(N);
        {
            Bool1D gt;
            Bool2D allgts;
            double af;
            int n_skipped = 0, prev_pos = -1;
            while(vcf.getNextVariant(var))
            {
                var.getGenotypes(gt);
                // only keep if meets conditions
                //  - bi-allelic
                //  - snp
                //  - position increased from previous site
                if(!var.isSNP() || !var.isNoneMissing() || !var.allPhased() || !(var.POS() > prev_pos))
                {
                    n_skipped += 1;
                    continue;
                }
                // keep track of snp index with AF < minaf
                int ac = 0, i = 0;
                for(bool g : gt)
                {
                    ac += g;
                    if(g && is_quilt_rhb) quilt.rare_per_hap_info[i].push_back(Mtotal + 1); // 1-based
                    i++;
                }
                af = (double)ac / N;
                if(af >= maf)
                {
                    keep.push_back(Mtotal);
                    allgts.push_back(gt);
                }
                Mtotal++;
                prev_pos = var.POS();
                if(is_quilt_rhb)
                {
                    quilt.ac.push_back(ac);
                    quilt.pos.push_back(var.POS());
                    quilt.ref.push_back(var.REF());
                    quilt.alt.push_back(var.ALT());
                    if(af >= maf)
                        quilt.snp_is_common.push_back(true);
                    else
                        quilt.snp_is_common.push_back(false);
                }
            }

            if(is_quilt_rhb)
            {
                quilt.nhaps = N;
                quilt.nsnps = Mtotal;
                quilt.n_skipped = n_skipped;
                const int QB = 32;
                const int n_common_snps = keep.size();
                const int nGrids = (n_common_snps + QB - 1) / QB; // keep grids of common SNPs only
                quilt.nGrids = nGrids;
                quilt.rhb_t.resize(N, Int1D(nGrids));
                int d32_times_bs, imax, ihap, k;
                // check rcpp_int_contract
                for(int bs = 0; bs < nGrids; bs++)
                {
                    for(ihap = 0; ihap < N; ihap++)
                    {
                        d32_times_bs = 32 * bs;
                        if(bs < (nGrids - 1))
                        {
                            imax = 31;
                        }
                        else
                        {
                            imax = n_common_snps - d32_times_bs - 1; // final one!
                        }
                        std::uint32_t itmp = 0;
                        for(k = imax; k >= 0; k--)
                        {
                            itmp <<= 1;
                            int j = allgts[d32_times_bs + k][ihap];
                            itmp |= j & 0x1;
                        }
                        quilt.rhb_t[ihap][bs] = itmp;
                    }
                }
            }

            M = keep.size(); // common SNPs only
            G = (M + B - 1) / B;
            if(G < nindices)
            {
                std::cerr << "resize mspbwt_nindices from " << nindices << " to 1" << std::endl;
                nindices = 1;
            }

            if(verbose)
                std::cerr << "N:" << N << ", M_total:" << Mtotal << ", M_common:" << M << ", G:" << G << ", B:" << B
                          << ", nindices:" << nindices << std::endl;
            X.resize(G, GridVec(N));
            k = 0;
            for(m = 0; m < M; m++)
            {
                for(i = 0; i < N; i++) X[k][i] = (X[k][i] << 1) | (allgts[m][i] != 0);
                if((m + 1) % B == 0)
                {
                    for(i = 0; i < N; i++) X[k][i] = reverseBits(X[k][i]); // reverse bits
                    k++; // update next grid
                }
            }
            if(G == k + 1)
            {
                for(i = 0; i < N; i++)
                {
                    X[k][i] <<= G * B - M;
                    X[k][i] = reverseBits(X[k][i]); // reverset bits
                }
            }
            else if(G == k)
            {
                if(verbose) std::cerr << "no need padding\n";
            }
            else
            {
                throw std::runtime_error("something wrong\n");
            }
        }

        // building A and save indices now
        A.resize(G);
        S.resize(G);
        C.resize(G);
        W.resize(G);
        if(is_save_D) D.resize(G);
        // create some intermediate varibales
        Int1D Occ;
        Int1D sqp; // for building D
        Int2D kas; // for building A
        Int2D kds; // for building D
        GridSetU symbols; // keep track of unique symbols at k
        GridVec y1(N);
        Int1D a0(N); // initilize a[k]
        Int1D d0(N + 1); // sentinels, d[0],d[N] = k + 1
        int Gi, ki, ni{-1};
        for(i = 0; i < nindices; i++)
        {
            auto Gv = seq_by(i, G - 1, nindices);
            Gi = Gv.size();
            // initial a and d at ki=-1
            std::iota(a0.begin(), a0.end(), 0);
            if(is_save_D)
            {
                std::fill(d0.begin(), d0.end(), 0);
                d0[0] = 1;
                d0[N] = 1;
            }
            if(verbose) std::cerr << "build_indices_k\n";
            for(ki = 0; ki < Gi; ki++)
            {
                ni++;
                k = Gv[ki];
                build_indices_k(ni, ki, X[k], y1, S[ni], C[ni], W[ni], Occ, kas, kds, sqp, a0, d0, symbols);
            }
        }
    }

    void build_indices_k(int ni,
                         int ki,
                         GridVec & xk,
                         GridVec & yk,
                         GridVec & sk,
                         Int1D & ck,
                         Int2D & wk,
                         Int1D & Occ,
                         Int2D & kas,
                         Int2D & kds,
                         Int1D & sqp,
                         Int1D & a0,
                         Int1D & d0,
                         GridSetU & symbols)
    {
        int s, n, c, e, i;
        for(n = 0; n < N; n++)
        {
            yk[n] = xk[a0[n]];
            symbols.insert(yk[n]);
        }
        sk = GridVec(symbols.begin(), symbols.end());
        symbols.clear();
        std::sort(sk.begin(), sk.end());
        int sn = sk.size();

        if(sn > MAXSYMBOLS) // keep only maximum 256 symbols
        {
            sk.erase(sk.end() - (sk.size() - MAXSYMBOLS), sk.end());
            sn = sk.size();
            for(n = 0; n < N; n++)
            {
                if(xk[a0[n]] > sk[sn - 1]) yk[n] = sk[sn - 1];
            }
        }

        // C[k] = save_C(y1, S[k]);
        // W[k] = save_Occ(y1, S[k]);
        // auto Wg = build_W(y1, S[k]); // here Wg is S x N
        // for (n = 0; n < N; n++)
        //     A[k][Wg[y1[n]][n] - 1] = a0[n];
        // next run
        // a0 = A[k];

        ck.resize(sn);
        wk.resize(sn);
        for(s = 0; s < sn; s++)
        {
            c = 0;
            for(n = 0; n < N; n++)
            {
                if(yk[n] < sk[s]) c++;
                if(yk[n] == sk[s]) Occ.push_back(n);
            }
            ck[s] = c;
            wk[s] = Occ;
            Occ.clear();
        }

        // save A and D
        kas.resize(sn);
        if(is_save_D)
        {
            kds.resize(sn);
            sqp.resize(sn, ki + 1);
        }
        // iterate over all haplotyps in reverese prefix sorted order
        for(n = 0; n < N; n++)
        {
            i = a0[n]; // index ihap at Y[k]
            if(is_save_D) e = d0[n];
            for(s = 0; s < sn; s++)
            {
                if(is_save_D)
                {
                    if(e > sqp[s]) sqp[s] = e;
                }
                if(xk[i] == sk[s])
                {
                    kas[s].push_back(i);
                    if(is_save_D)
                    {
                        kds[s].push_back(sqp[s]);
                        sqp[s] = 0;
                    }
                }
            }
            // allow xk > sk[sn-1] is useful when we want the number of symbols to be capped
            s = sn - 1;
            if(xk[i] > sk[s])
            {
                kas[s].push_back(i);
                if(is_save_D)
                {
                    kds[s].push_back(sqp[s]);
                    sqp[s] = 0;
                }
            }
        }
        // make sure A[k], D[k] are empty
        for(s = 0; s < sn; s++)
        {
            A[ni].insert(A[ni].end(), kas[s].begin(), kas[s].end());
            if(is_save_D) D[ni].insert(D[ni].end(), kds[s].begin(), kds[s].end());
        }
        kas.clear();
        a0 = A[ni];
        if(is_save_D)
        {
            kds.clear();
            sqp.clear();
            D[ni][N] = ki + 1; // add sentinel
            d0 = D[ni];
        }
    }

    void report_neighourings(Int1D & haps, Int1D & ends, Int1D & lens, GridVec & zg, int mspbwtL = 32, int mspbwtM = 3)
    {
        int k, s, klen, j, n, l, Gi, ki, iind, ks{-1};
        for(iind = 0; iind < nindices; iind++)
        {
            auto Gv = seq_by(iind, G - 1, nindices);
            Gi = Gv.size();
            int f, g, e, f1, g1, e1, fg, cur, prev, valid_grid_start{0};
            Int1D hap_up_prev(mspbwtL, -1), hap_up_cur(mspbwtL, -1), len_up_prev(mspbwtL, -1), len_up_cur(mspbwtL, -1);
            Int1D hap_down_prev(mspbwtL, -1), hap_down_cur(mspbwtL, -1), len_down_prev(mspbwtL, -1),
                len_down_cur(mspbwtL, -1);
            for(ki = 0; ki < Gi; ki++)
            {
                ks++;
                k = Gv[ki];
                auto kzs = std::lower_bound(S[ks].begin(), S[ks].end(), zg[k]);
                s = std::fmin(std::distance(S[ks].begin(), kzs), S[ks].size() - 1);
                if(S[ks][s] != zg[k])
                {
                    if(verbose) std::cerr << "ghost symbol at k: " << k << std::endl;
                    // coerce zg[k] to the closest symbol in S[ks]
                    Int1D tmp(S[ks].size());
                    for(size_t i = 0; i < S[ks].size(); i++) tmp[i] = count1bits(zg[k] ^ S[ks][i]);
                    l = std::distance(tmp.begin(), std::min_element(tmp.begin(), tmp.end()));
                    zg[k] = S[ks][l];
                    auto kzs = std::lower_bound(S[ks].begin(), S[ks].end(), zg[k]);
                    s = std::fmin(std::distance(S[ks].begin(), kzs), S[ks].size() - 1);
                }

                if(ki == valid_grid_start)
                {
                    f1 = C[ks][s];
                    g1 = C[ks][s] + W[ks][s].size(); // could be N
                    e1 = valid_grid_start;
                    /* fg = (f1 + g1 - 1) / 2; */
                    fg = std::fmin(f1 + mspbwtL, N - 1);
                    // initilize prev haps/lens first
                    for(l = 0; l < mspbwtL; l++)
                    {
                        n = A[ks][std::fmax(fg - l, 0)]; // go up
                        hap_up_prev[l] = n;
                        len_up_prev[l] = 0;
                        if(fg - l == 0) break;
                    }
                    for(l = 0; l < mspbwtL; l++)
                    {
                        n = A[ks][std::fmin(fg + l + 1, N - 1)]; // go down
                        hap_down_prev[l] = n;
                        len_down_prev[l] = 0;
                        if(fg + l + 1 == N - 1) break;
                    }
                }
                else if(ki > valid_grid_start)
                {
                    // assume symbol exists, then g1 >= f1 >= C[ks][s] if g >= f
                    auto fzk = std::lower_bound(W[ks][s].begin(), W[ks][s].end(), f);
                    auto gzk = std::lower_bound(W[ks][s].begin(), W[ks][s].end(), g);
                    f1 = C[ks][s] + std::fmin(std::distance(W[ks][s].begin(), fzk), W[ks][s].size());
                    g1 = C[ks][s] + std::fmin(std::distance(W[ks][s].begin(), gzk), W[ks][s].size());
                    fg = f1;

                    for(prev = 0, cur = 0; prev < mspbwtL; prev++)
                    {
                        n = A[ks][std::fmax(fg - cur, 0)]; // scan up
                        if(n == hap_up_prev[prev])
                        {
                            // it is a match, save and increment
                            hap_up_cur[cur] = n;
                            len_up_cur[cur] = len_up_prev[prev] + 1;
                            cur++;
                        }
                        else
                        {
                            // this is not a match, report it
                            klen = len_up_prev[prev];
                            if(klen >= mspbwtM)
                            {
                                haps.push_back(hap_up_prev[prev]);
                                lens.push_back(klen);
                                ends.push_back(ki - 1);
                            }
                        }
                        if(fg - cur == 0) break;
                    }

                    for(cur = 0; cur < mspbwtL; cur++)
                    {
                        n = A[ks][std::fmax(fg - cur, 0)]; // scan up
                        hap_up_cur[cur] = n;
                        // go backwards to find the lens
                        for(klen = 0, j = ki; j >= 0; j--)
                        {
                            if(X[Gv[j]][n] == zg[Gv[j]])
                                klen++;
                            else
                                break;
                        }
                        len_up_cur[cur] = klen;
                    }
                    hap_up_prev = hap_up_cur; // next
                    len_up_prev = len_up_cur; // next

                    for(prev = 0, cur = 0; prev < mspbwtL; prev++)
                    {
                        n = A[ks][std::fmax(fg + cur + 1, N - 1)]; // scan down
                        if(n == hap_down_prev[prev])
                        {
                            // it is a match, save and increment
                            hap_down_cur[cur] = n;
                            len_down_cur[cur] = len_down_prev[prev] + 1;
                            cur++;
                        }
                        else
                        {
                            // this is not a match, report it
                            klen = len_down_prev[prev];
                            if(klen >= mspbwtM)
                            {
                                haps.push_back(hap_down_prev[prev]);
                                lens.push_back(klen);
                                ends.push_back(ki - 1);
                            }
                        }
                        if(fg + cur + 1 == N - 1) break;
                    }

                    for(cur = 0; cur < mspbwtL; cur++)
                    {
                        n = A[ks][std::fmax(fg + cur + 1, N - 1)]; // scan down
                        hap_down_cur[cur] = n;
                        // go backwards to find the lens
                        for(klen = 0, j = ki; j >= 0; j--)
                        {
                            if(X[Gv[j]][n] == zg[Gv[j]])
                                klen++;
                            else
                                break;
                        }
                        len_down_cur[cur] = klen;
                    }
                    hap_down_prev = hap_down_cur; // next
                    len_down_prev = len_down_cur; // next

                    if(!(g1 > f1)) // g1==f1
                    {
                        // restart to find new f1 and g1
                        bool matches_lower = false, matches_upper = false;
                        if(is_save_D)
                            e1 = D[ks][f1] - 1; // y[f1] and y[f1-1] diverge here, so upper bound for e
                        else
                            e1 = ki - 1; // without D, go backwards find the divergence
                        // if there is a ghost symbol in zg, e1++ will continue forever so add e1 < k
                        while((e1 < ki) && (!matches_lower) && (!matches_upper))
                        {
                            if(f1 > 0)
                                matches_upper = (zg[Gv[e1]] == X[Gv[e1]][A[ks][f1 - 1]]);
                            else
                                matches_upper = false;
                            if(f1 < N)
                                matches_lower = (zg[Gv[e1]] == X[Gv[e1]][A[ks][f1]]);
                            else
                                matches_lower = false;
                            // if matches neither y[f1] or y[f1-1], eg. symbol missing or just happens, e1++
                            if((!matches_lower) && (!matches_upper)) ++e1;
                        }
                        // if ghost symbols exists and loop --e1 stops as it is
                        // this will shorten the potential longest matches
                        if(matches_upper)
                        {
                            --f1;
                            // make sure e1 > 0
                            while(e1 > valid_grid_start && zg[Gv[e1 - 1]] == X[Gv[e1 - 1]][A[ks][f1]]) --e1;
                            if(is_save_D)
                            {
                                while(f1 > 0 && D[ks][f1] <= e1) --f1;
                            }
                            else
                            {
                                int ei = e1;
                                while(f1 > 0 && ei <= e1)
                                {
                                    ei = ki;
                                    while(ei > valid_grid_start && X[Gv[ei]][A[ks][f1 - 1]] == X[Gv[ei]][A[ks][f1]])
                                        --ei;
                                    --f1;
                                }
                            }
                        }
                        if(matches_lower)
                        {
                            ++g1;
                            // make sure e1 > 0
                            while(e1 > valid_grid_start && zg[Gv[e1 - 1]] == X[Gv[e1 - 1]][A[ks][f1]]) --e1;
                            if(is_save_D)
                            {
                                while(g1 < N && D[ks][g1] <= e1) ++g1;
                            }
                            else
                            {
                                int ei = e1;
                                while(g1 > N && ei <= e1)
                                {
                                    ei = ki;
                                    while(ei > valid_grid_start && X[Gv[ei]][A[ks][g1]] == X[Gv[ei]][A[ks][g1 + 1]])
                                        --ei;
                                    ++g1;
                                }
                            }
                        }
                    }
                }
                // update for ext run
                f = f1;
                g = g1;
                e = e1;
            }
            // final report everything
            for(l = 0; l < mspbwtL; l++)
            {
                haps.push_back(hap_up_cur[l]);
                lens.push_back(len_up_cur[l]);
                ends.push_back(ki - 1);
                haps.push_back(hap_down_cur[l]);
                lens.push_back(len_down_cur[l]);
                ends.push_back(ki - 1);
            }
        }
    }

    int save(const std::string & filename)
    {
        std::ofstream out(filename, std::ios::out | std::ios::binary);
        if(out.fail()) return 1;
        out.write((char *)&B, sizeof(B));
        out.write((char *)&M, sizeof(M));
        out.write((char *)&N, sizeof(N));
        out.write((char *)&G, sizeof(G));
        out.write((char *)&nindices, sizeof(nindices));
        // write reorder (M)
        int n, k, m;
        for(m = 0; m < M; m++) out.write((char *)&keep[m], sizeof(int));
        if(is_save_X)
        {
            // write X
            for(k = 0; k < G; k++)
            {
                for(n = 0; n < N; n++) out.write((char *)&X[k][n], sizeof(T));
            }
        }
        // write S
        for(k = 0; k < G; k++)
        {
            size_t sz = S[k].size();
            out.write((char *)&sz, sizeof(sz));
            for(size_t i = 0; i < sz; i++) out.write((char *)&S[k][i], sizeof(T));
        }
        // write A
        for(k = 0; k < G; k++)
        {
            for(n = 0; n < N; n++) out.write((char *)&A[k][n], sizeof(int));
        }
        if(is_save_D)
        {
            // write D
            for(k = 0; k < G; k++)
            {
                for(n = 0; n < N + 1; n++) out.write((char *)&D[k][n], sizeof(int));
            }
        }
        // write C
        for(k = 0; k < G; k++)
        {
            size_t sz = C[k].size();
            out.write((char *)&sz, sizeof(sz));
            for(size_t i = 0; i < sz; i++) out.write((char *)&C[k][i], sizeof(int));
        }
        // write W
        for(k = 0; k < G; k++)
        {
            size_t sz = W[k].size();
            out.write((char *)&sz, sizeof(sz));
            for(size_t i = 0; i < sz; i++)
            {
                size_t sz1 = W[k][i].size();
                out.write((char *)&sz1, sizeof(sz1));
                for(size_t j = 0; j < sz1; j++) out.write((char *)&W[k][i][j], sizeof(int));
            }
        }
        return 0;
    }

    int load(const std::string & filename)
    {
        std::ifstream in(filename, std::ios::in | std::ios::binary);
        if(in.fail()) return 1;
        if(!in.read((char *)&B, sizeof(B))) return 2;
        if(B != sizeof(T) * 8) throw std::invalid_argument("the binary file may be created with different B!\n");
        if(!in.read((char *)&M, sizeof(M))) return 2;
        if(!in.read((char *)&N, sizeof(N))) return 2;
        if(!in.read((char *)&G, sizeof(G))) return 2;
        if(!in.read((char *)&nindices, sizeof(nindices))) return 2;
        std::cerr << "N: " << N << ",M: " << M << ",G: " << G << ",B: " << B << ",nindices " << nindices << std::endl;
        keep.resize(M);
        int n, m, k;
        for(m = 0; m < M; m++) in.read((char *)&keep[m], sizeof(int));
        if(is_save_X)
        {
            // read X
            X.resize(G, GridVec(N));
            for(k = 0; k < G; k++)
            {
                for(n = 0; n < N; n++)
                {
                    if(!in.read((char *)&X[k][n], sizeof(T))) return 2;
                }
            }
            if(verbose) std::cerr << "load X done" << std::endl;
        }
        // read S
        S.resize(G);
        for(k = 0; k < G; k++)
        {
            size_t sz;
            if(!in.read((char *)&sz, sizeof(sz)) || sz < 1) return 2;
            S[k].resize(sz);
            for(size_t i = 0; i < sz; i++)
            {
                if(!in.read((char *)&S[k][i], sizeof(T))) return 2;
            }
        }
        if(verbose) std::cerr << "load S done" << std::endl;
        // read A
        A.resize(G, Int1D(N));
        for(k = 0; k < G; k++)
        {
            for(n = 0; n < N; n++)
            {
                if(!in.read((char *)&A[k][n], sizeof(int))) return 2;
            }
        }
        if(verbose) std::cerr << "load A done" << std::endl;
        if(is_save_D)
        {
            // read D
            D.resize(G, Int1D(N + 1));
            for(k = 0; k < G; k++)
            {
                for(n = 0; n < N + 1; n++)
                {
                    if(!in.read((char *)&D[k][n], sizeof(int))) return 2;
                }
            }
            if(verbose) std::cerr << "load D done" << std::endl;
            // for (int i = 0; i < N + 1; i++)
            // {
            //     for (int j = 0; j < G; j++)
            //         cerr << D[j][i] << ",";
            //     cerr << endl;
            // }
        }
        // read C
        C.resize(G);
        for(k = 0; k < G; k++)
        {
            size_t sz;
            if(!in.read((char *)&sz, sizeof(sz)) || sz < 1) return 2;
            C[k].resize(sz);
            for(size_t i = 0; i < sz; i++)
            {
                if(!in.read((char *)&C[k][i], sizeof(int))) return 2;
            }
        }
        if(verbose) std::cerr << "load C done" << std::endl;
        // read W
        W.resize(G);
        for(k = 0; k < G; k++)
        {
            size_t sz;
            if(!in.read((char *)&sz, sizeof(sz)) || sz < 1) return 2;
            W[k].resize(sz);
            for(size_t i = 0; i < sz; i++)
            {
                size_t sz1;
                if(!in.read((char *)&sz1, sizeof(sz1)) || sz1 < 1) return 2;
                W[k][i].resize(sz1);
                for(size_t j = 0; j < sz1; j++)
                {
                    if(!in.read((char *)&W[k][i][j], sizeof(int))) return 2;
                }
            }
        }
        if(verbose) std::cerr << "load W done" << std::endl;

        return 0;
    }
};

#endif // MSPBWT_H_
