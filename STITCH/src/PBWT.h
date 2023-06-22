#ifndef PBWT_H_
#define PBWT_H_

#include "vcfpp/vcfpp.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <fstream>
#include <iterator>
#include <random>
#include <stdexcept>

using namespace std;

using Bool1D = vector<bool>;
using Bool2D = vector<Bool1D>;
using Int1D = vector<int>;
using Int2D = vector<Int1D>;
using Int3D = vector<Int2D>;

class PBWT
{
  public:
    Bool2D X;
    Int2D A; // M x N
    Int2D D; // M x N
    Int2D U; // M x (N+1)
    int M{0}, N{0};
    bool verbose{0};

  public:
    PBWT() {}

    virtual ~PBWT() {}

    int save(const std::string & filename)
    {
        ofstream out(filename, ios::out | ios::binary);
        if(out.fail()) return 1;
        out.write((char *)&M, sizeof(M));
        out.write((char *)&N, sizeof(N));
        int k, n;
        // write X
        for(k = 0; k < M; k++)
        {
            for(bool b : X[k]) out.write(reinterpret_cast<const char *>(&b), sizeof(bool));
        }
        // write A
        for(k = 0; k < M; k++)
        {
            for(n = 0; n < N; n++) out.write((char *)&A[k][n], sizeof(int));
        }
        // write D
        for(k = 0; k < M; k++)
        {
            for(n = 0; n < N; n++) out.write((char *)&D[k][n], sizeof(int));
        }
        // write U
        for(k = 0; k < M; k++)
        {
            for(n = 0; n < N + 1; n++) out.write((char *)&U[k][n], sizeof(int));
        }
        return 0;
    }

    int load(const std::string & filename)
    {
        ifstream in(filename, ios::in | ios::binary);
        if(in.fail()) return 1;
        in.read((char *)&M, sizeof(M));
        in.read((char *)&N, sizeof(N));
        int k, n;
        // write X
        X.resize(M, Bool1D(N));
        bool b;
        for(k = 0; k < M; k++)
        {
            for(n = 0; n < N; n++)
            {
                in.read((char *)&b, sizeof(bool));
                X[k][n] = b;
            }
        }
        // write A
        A.resize(M, Int1D(N));
        for(k = 0; k < M; k++)
        {
            for(n = 0; n < N; n++) in.read((char *)&A[k][n], sizeof(int));
        }
        // write D
        D.resize(M, Int1D(N));
        for(k = 0; k < M; k++)
        {
            for(n = 0; n < N; n++) in.read((char *)&D[k][n], sizeof(int));
        }
        // write A
        U.resize(M, Int1D(N + 1));
        for(k = 0; k < M; k++)
        {
            for(n = 0; n < N + 1; n++) in.read((char *)&U[k][n], sizeof(int));
        }
        return 0;
    }

    void build(std::string vcfpanel, std::string samples, std::string region)
    {
        vcfpp::BcfReader vcf(vcfpanel, samples, region);
        vcfpp::BcfRecord var(vcf.header);
        // get N = # haplotypes, M = # sites
        N = vcf.nsamples * 2, M = 0;
        Bool1D gt;
        while(vcf.getNextVariant(var))
        {
            var.getGenotypes(gt);
            if(!var.isNoneMissing() || !var.allPhased())
                throw runtime_error("there is a site with non-phased or missing value!\n");
            M++;
            X.push_back(gt);
        }
        A.resize(M, Int1D(N));
        D.resize(M, Int1D(N));
        U.resize(M, Int1D(N + 1));
        // V.resize(M, Int1D(N + 1));
        Int1D a1(N), d1(N), a0(N), d0(N);
        int i, k, u_, v_, p, q;
        for(k = 0; k < M; k++)
        {
            u_ = 0, v_ = 0, p = k + 1, q = k + 1;
            for(i = 0; i < N; i++)
            {
                int d_ = k > 0 ? D[k - 1][i] : 0;
                int a_ = k > 0 ? A[k - 1][i] : i;
                p = std::max(p, d_);
                q = std::max(q, d_);
                U[k][i] = u_;
                // V[k][i] = v_;
                if(X[k][a_])
                {
                    a1[v_] = a_;
                    d1[v_] = q;
                    v_++;
                    q = 0;
                }
                else
                {
                    a0[u_] = a_;
                    d0[u_] = p;
                    u_++;
                    p = 0;
                }
            }
            U[k][N] = u_;
            // V[k][N] = N;
            for(i = 0; i < N; i++)
            {
                // V[k][i] += u_;
                if(i < u_)
                {
                    A[k][i] = a0[i];
                    D[k][i] = d0[i];
                }
                else
                {
                    A[k][i] = a1[i - u_];
                    D[k][i] = d1[i - u_];
                }
            }
        }
        cerr << "build pbwt indices done\n";
    }

    Bool1D encodezg(const Int1D & z_)
    {
        Bool1D z;
        for(int i : z_) z.push_back(i);
        return z;
    }

    void report_setmaximal(Bool1D z, Int1D & haps, Int1D & ends, Int1D & lens)
    {
        int j, k, e, f, g, f1, g1, e1;
        e = 0, f = 0, g = N;
        for(k = 0; k < M; k++)
        {
            f1 = z[k] ? U[k][N] + (f - U[k][f]) : U[k][f];
            g1 = z[k] ? U[k][N] + (g - U[k][g]) : U[k][g];
            if(g1 > f1) // no change to e
            {
                f = f1;
                g = g1;
            }
            else // we have reached a maximum -- need to report and update e, f, g
            {
                // first report matches
                for(j = f; j < g; ++j)
                {
                    // klen = (k - 1) - e + 1;
                    // hap_idx, end, len
                    // cout << A[k - 1][j] << "\t" << e << "\t" << k - 1 << endl;
                    haps.push_back(A[k - 1][j]);
                    lens.push_back(k - e);
                    ends.push_back(k - 1);
                    // check matches
                    if(verbose)
                    {
                        for(int i = e; i < k; i++) assert(z[i] == X[i][A[k - 1][j]]);
                        if(k < M && z[k] == X[k][A[k - 1][j]])
                            throw runtime_error("match not maximal - can extend forwards\n");
                        if(e > 1 && z[e - 1] == X[e - 1][A[k - 1][j]])
                            throw runtime_error("match not maximal - can extend backwards\n");
                    }
                }
                // then update e,f,g
                e1 = D[k][f1] - 1; /* y[f1] and y[f1-1] diverge here, so upper bound for e */
                if((z[e1] == 0 && f1 > 0) || f1 == N)
                { // match upper
                    f1 = g1 - 1;
                    while(e1 > 0 && z[e1 - 1] == X[e1 - 1][A[k][f1]]) --e1;
                    while(D[k][f1] <= e1) --f1;
                }
                else if(f1 < N)
                {
                    g1 = f1 + 1;
                    while(e1 > 0 && z[e1 - 1] == X[e1 - 1][A[k][f1]]) --e1;
                    while(g1 < M && D[k][g1] <= e1) ++g1;
                }
                e = e1, f = f1, g = g1;
            }
        }
        // report matches to the end
        for(j = f; j < g; ++j)
        {
            // hap_idx, end, len
            // cout << A[k - 1][j] << "\t" << e << "\t" << k - 1 << endl;
            haps.push_back(A[k - 1][j]);
            lens.push_back(k - e);
            ends.push_back(k - 1);
            if(verbose)
            {
                for(int i = e; i < k; i++) assert(z[i] == X[i][A[k - 1][j]]);
                if(k < M && z[k] == X[k][A[k - 1][j]]) throw runtime_error("match not maximal - can extend forwards\n");
                if(e > 1 && z[e - 1] == X[e - 1][A[k - 1][j]])
                    throw runtime_error("match not maximal - can extend backwards\n");
            }
        }
    }

    void query(std::string vcfquery, std::string samples, std::string region, std::string outfile)
    {
        vcfpp::BcfReader vcf(vcfquery, samples, region);
        vcfpp::BcfRecord var(vcf.header);
        Bool1D gt;
        Bool2D G, GT;
        while(vcf.getNextVariant(var))
        {
            var.getGenotypes(gt);
            if(!var.isNoneMissing() || !var.allPhased())
                throw runtime_error("there is a site with non-phased or missing value!\n");
            G.push_back(gt);
        }
        GT.resize(G[0].size(), Bool1D(G.size()));
        for(size_t i = 0; i < G.size(); i++)
        {
            for(size_t j = 0; j < G[i].size(); j++) GT[j][i] = G[i][j];
        }
        std::ofstream ofs(outfile);
        for(size_t j = 0; j < GT.size(); j++)
        {
            Int1D haps, lens, ends;
            report_setmaximal(GT[j], haps, ends, lens);
            for(size_t i = 0; i < haps.size(); i++)
            {
                ofs << j << "\t" << haps[i] << "\t" << ends[i] - lens[i] + 1 << "\t" << ends[i] << "\t" << lens[i]
                    << std::endl;
            }
        }
    }
};

#endif // PBWT_H_
