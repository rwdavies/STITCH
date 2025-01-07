/*******************************************************************************
 * @file        https://github.com/Zilong-Li/vcfpp/vcfpp.h
 * @author      Zilong Li
 * @email       zilong.dk@gmail.com
 * @version     v0.6.1
 * @breif       a single C++ file for manipulating VCF
 * Copyright (C) 2022-2023.The use of this code is governed by the LICENSE file.
 ******************************************************************************/

/*! \mainpage The documentation of the single C++ file *vcfpp.h* for manipulating VCF/BCF
 *
 * \section intro_sec Introduction
 *
 * This project https://github.com/Zilong-Li/vcfpp introduces a single C++ file as interface to the basic
 * htslib, which can be easily included in a C++ program for scripting high-performance genomic analyses.
 *
 * - vcfpp.BcfHeader keeps track of the header information in  VCF/BCF
 * - vcfpp.BcfRecord keeps track of the variants information in VCF/BCF
 * - vcfpp.BcfReader streams in variants from VCF/BCF file or stdin
 * - vcfpp.BcfWriter streams out variants to VCF/BCF file or stdout
 *
 * \section install_sec Installation
 *
 * - <EM> include "vcfpp.h" </EM> to your program and compile it by <EM> g++ my.cpp -std=c++11 -Wall -I. -lhts
 * - </EM>
 * - make sure you have https://github.com/samtools/htslib installed on your system and the it is in your
 * environment.
 *
 *
 * \copyright Copyright (C) 2022 Zilong Li . This project is governed by the LICENSE file in
 * https://github.com/Zilong-Li/vcfpp.
 *
 *
 */

#pragma once

#include <cstddef>
#include <stdexcept>
#ifndef VCFPP_H_
#    define VCFPP_H_

#    include <iostream>
#    include <memory>
#    include <string>
#    include <type_traits>
#    include <vector>

// make sure you have htslib installed
extern "C"
{
#    include <htslib/hts.h>
#    include <htslib/kstring.h>
#    include <htslib/tbx.h>
#    include <htslib/vcf.h>
#    include <htslib/vcfutils.h>
}

namespace vcfpp
{
template<typename T>
using isValidFMT =
    typename std::enable_if<std::is_same<T, std::string>::value || std::is_same<T, std::vector<char>>::value
                                || std::is_same<T, std::vector<int>>::value
                                || std::is_same<T, std::vector<float>>::value,
                            bool>::type;

template<typename T>
using isValidInfo =
    typename std::enable_if<std::is_same<T, std::string>::value || std::is_same<T, std::vector<int>>::value
                                || std::is_same<T, std::vector<float>>::value,
                            bool>::type;

template<typename T>
using isInfoVector = typename std::enable_if<std::is_same<T, std::vector<int>>::value
                                                 || std::is_same<T, std::vector<float>>::value,
                                             bool>::type;

template<typename T>
using isScalar = typename std::enable_if<std::is_same<T, int>::value || std::is_same<T, float>::value
                                             || std::is_same<T, double>::value,
                                         bool>::type;

template<typename T>
using isString = typename std::enable_if<std::is_same<T, std::string>::value, bool>::type;

template<typename T>
using isValidGT = typename std::enable_if<std::is_same<T, std::vector<bool>>::value
                                              || std::is_same<T, std::vector<char>>::value,
                                          bool>::type;

template<typename T>
using isFormatVector = typename std::enable_if<std::is_same<T, std::vector<float>>::value
                                                   || std::is_same<T, std::vector<char>>::value
                                                   || std::is_same<T, std::vector<int>>::value,
                                               bool>::type;

namespace details
{

template<typename T>
isScalar<T> isMissing(T const & v)
{
    if(std::is_same<T, float>::value)
    {
        return bcf_float_is_missing(v);
    }
    else if(std::is_same<T, int>::value)
    {
        return bcf_int32_missing(v);
    }
}

// test if a string is end with another string
inline bool isEndWith(std::string const & s, std::string const & e)
{
    if(s.length() >= e.length())
    {
        return (0 == s.compare(s.length() - e.length(), e.length(), e));
    }
    else
    {
        return false;
    }
}

// determinate the mode for read/write the compressed/uncompressed VCF/BCF
inline std::string getMode(std::string const & fname, std::string mode = "r")
{
    if(isEndWith(fname, "bcf.gz")) mode += "b";
    if(isEndWith(fname, "bcf")) mode += "bu";
    if(isEndWith(fname, "vcf.gz")) mode += "z";
    return mode;
}

// string split by separator
inline std::vector<std::string> split_string(const std::string & s, const std::string & separators)
{
    std::vector<std::string> ret;
    bool is_seperator[256] = {false};
    for(auto & ch : separators)
    {
        is_seperator[(unsigned int)ch] = true;
    }
    int begin = 0;
    for(int i = 0; i <= (int)s.size(); i++)
    {
        if(is_seperator[(uint8_t)s[i]] || i == (int)s.size())
        {
            ret.push_back(std::string(s.begin() + begin, s.begin() + i));
            begin = i + 1;
        }
    }
    return ret;
}

// deleter for htsFile
struct hts_file_close
{
    void operator()(htsFile * x)
    {
        if(x) hts_close(x);
    }
};

// deleter for hts iterator
struct hts_iter_close
{
    void operator()(hts_itr_t * x)
    {
        if(x) hts_itr_destroy(x);
    }
};

// deleter for hts idx
struct hts_idx_close
{
    void operator()(hts_idx_t * x)
    {
        if(x) hts_idx_destroy(x);
    }
};

// deleter for tabix idx
struct tabix_idx_close
{
    void operator()(tbx_t * x)
    {
        if(x) tbx_destroy(x);
    }
};

// deleter for variant
struct bcf_line_close
{
    void operator()(bcf1_t * x)
    {
        if(x) bcf_destroy(x);
    }
};

// deleter for bcf header
struct bcf_hdr_close
{
    void operator()(bcf_hdr_t * x)
    {
        if(x) bcf_hdr_destroy(x);
    }
};

} // namespace details

/**
 * @class BcfHeader
 * @brief Object represents header of the VCF, offering methods to access and modify the tags
 * @note  BcfHeader has 3 friends, BcfReader, BcfWriter and BcfRecord.
 **/
class BcfHeader
{
    friend class BcfRecord;
    friend class BcfReader;
    friend class BcfWriter;

  private:
    bcf_hdr_t * hdr = NULL; // bcf header
    bcf_hrec_t * hrec = NULL; // populate header

  public:
    BcfHeader() {}

    ~BcfHeader()
    {
        if(hrec) bcf_hrec_destroy(hrec);
        if(hdr) bcf_hdr_destroy(hdr);
    }

    /** @brief stream out the header */
    friend std::ostream & operator<<(std::ostream & out, const BcfHeader & h)
    {
        out << h.asString();
        return out;
    }

    // TODO: check if the value is valid for vcf specification

    /** @brief add INFO field to header
     *  @param id tag name in INFO field
     *  @param number Number of elements
     *  @param type Integer, Floast or String
     *  @param description Description of the tag
     *  */
    inline void addINFO(const std::string & id,
                        const std::string & number,
                        const std::string & type,
                        const std::string & description)
    {
        addLine("##INFO=<ID=" + id + ",Number=" + number + ",Type=" + type + ",Description=\"" + description
                + "\">");
    }

    /** @brief add FORMAT field to header
     *  @param id tag name in FORMAT field
     *  @param number Number of elements
     *  @param type Integer, Floast or String
     *  @param description Description of the tag
     *  */
    inline void addFORMAT(const std::string & id,
                          const std::string & number,
                          const std::string & type,
                          const std::string & description)
    {
        addLine("##FORMAT=<ID=" + id + ",Number=" + number + ",Type=" + type + ",Description=\"" + description
                + "\">");
    }

    /**
     * @brief add one FILTER line to header
     * @param id FILTER name
     * @param description Description of the FILTER
     * */
    inline void addFILTER(const std::string & id, const std::string & description)
    {
        addLine("##FILTER=<ID=" + id + ",Description=\"" + description + "\">");
    }

    /** @brief add contig to header
     *  @param id contig or chromosome name
     *  */
    inline void addContig(const std::string & id)
    {
        addLine("##contig=<ID=" + id + ">");
    }

    /**
     * @brief add one line to header
     * */
    inline void addLine(const std::string & str)
    {
        int ret = 0;
        ret = bcf_hdr_append(hdr, str.c_str());
        if(ret != 0) throw std::runtime_error("could not add " + str + " to header\n");
        ret = bcf_hdr_sync(hdr);
        if(ret != 0) throw std::runtime_error("could not add " + str + " to header\n");
    }

    /** @brief add one sample name to header */
    inline void addSample(const std::string & sample) const
    {
        bcf_hdr_add_sample(hdr, sample.c_str());
        if(bcf_hdr_sync(hdr) != 0)
        {
            throw std::runtime_error("couldn't add the sample.\n");
        }
    }

    /** @brief return header as a string */
    inline std::string asString() const
    {
        kstring_t s = {0, 0, NULL}; // kstring
        if(bcf_hdr_format(hdr, 0, &s) == 0) // append header string to s.s! append!
        {
            std::string out(s.s, s.l);
            free(s.s);
            return out;
        }
        else
        {
            throw std::runtime_error("failed to convert formatted header to string");
        }
    }

    /** @brief return all samples in a string vector */
    std::vector<std::string> getSamples() const
    {
        std::vector<std::string> vec;
        for(int i = 0; i < bcf_hdr_nsamples(hdr); i++)
        {
            vec.push_back(std::string(hdr->samples[i]));
        }
        return vec;
    }

    /** @brief return all contig/chromosome names in a string vector */
    std::vector<std::string> getSeqnames() const
    {
        int ret = 0;
        const char ** seqs = bcf_hdr_seqnames(hdr, &ret);
        if(ret == 0) printf("there is no contig id in the header!\n");
        std::vector<std::string> vec;
        for(int i = 0; i < ret; i++)
        {
            vec.push_back(std::string(seqs[i]));
        }
        free(seqs);
        return vec;
    }

    /**
     *  @brief get the type of a given tag
     *  @param tag in the FORMAT
     *  @return 1: int; 2: float; 3: string; 0: error;
     */
    inline int getFormatType(std::string tag) const
    {
        int tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, tag.c_str());
        if(tag_id < 0) return 0;
        if(bcf_hdr_id2type(hdr, BCF_HL_FMT, tag_id) == (BCF_HT_INT & 0xff))
        {
            return 1;
        }
        else if(bcf_hdr_id2type(hdr, BCF_HL_FMT, tag_id) == (BCF_HT_REAL & 0xff))
        {
            return 2;
        }
        else if(bcf_hdr_id2type(hdr, BCF_HL_FMT, tag_id) == (BCF_HT_STR & 0xff))
        {
            return 3;
        }
        else
        {
            return 0;
        }
    }

    /** @brief remove a contig tag from header */
    inline void removeContig(std::string tag) const
    {
        bcf_hdr_remove(hdr, BCF_HL_CTG, tag.c_str());
    }

    /** @brief remove a INFO tag from header */
    inline void removeInfo(std::string tag) const
    {
        bcf_hdr_remove(hdr, BCF_HL_INFO, tag.c_str());
    }

    /** @brief remove a FORMAT tag from header */
    inline void removeFormat(std::string tag) const
    {
        bcf_hdr_remove(hdr, BCF_HL_FMT, tag.c_str());
    }

    /** @brief remove a FILTER tag from header */
    inline void removeFilter(std::string tag) const
    {
        bcf_hdr_remove(hdr, BCF_HL_FLT, tag.c_str());
    }

    /**
     * @brief update the sample id in the header
     * @param samples a comma-separated string for multiple new samples
     * @note this only update the samples name in a given sequential order
     * */
    void updateSamples(const std::string & samples)
    {
        auto ss = details::split_string(samples, ",");
        const int nsamples = nSamples();
        if(nsamples != (int)ss.size())
            throw std::runtime_error("You provide either too few or too many samples");
        kstring_t htxt = {0, 0, 0};
        bcf_hdr_format(hdr, 1, &htxt);
        // Find the beginning of the #CHROM line
        int i = htxt.l - 2, ncols = 0;
        while(i >= 0 && htxt.s[i] != '\n')
        {
            if(htxt.s[i] == '\t') ncols++;
            i--;
        }
        if(i < 0 || strncmp(htxt.s + i + 1, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", 45))
        {
            if(i > 0 && !strncmp(htxt.s + i + 1, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", 38))
                throw std::runtime_error("Error: missing FORMAT fields, cowardly refusing to add samples\n");
            throw std::runtime_error("Could not parse the header: " + std::string(htxt.s));
        }

        ncols = 0;
        while(ncols != 9)
        {
            i++;
            if(htxt.s[i] == '\t') ncols++;
        }
        htxt.l = i;

        // Replace all samples
        for(i = 0; i < nsamples; i++)
        {
            kputc('\t', &htxt);
            kputs(ss[i].c_str(), &htxt);
        }
        kputc('\n', &htxt);

        // destroy the old and make new header
        bcf_hdr_destroy(hdr);
        hdr = bcf_hdr_init("r");
        if(bcf_hdr_parse(hdr, htxt.s) < 0)
            throw std::runtime_error("An error occurred while parsing the header\n");
        // free everything
        free(htxt.s);
    }

    /**
     * @brief explicitly set samples to be extracted
     * @param samples samples to include or exclude as a comma-separated string
     * */
    inline void setSamples(const std::string & samples) const
    {
        int ret = 0;
        ret = bcf_hdr_set_samples(hdr, samples.c_str(), 0);
        if(ret != 0)
        {
            throw std::runtime_error("the " + std::to_string(ret)
                                     + "-th sample are not in the VCF.\nparameter samples:" + samples);
        }
    }

    /** @brief set the VCF version */
    inline void setVersion(const std::string & version) const
    {
        bcf_hdr_set_version(hdr, version.c_str());
    }

    /** @brief return the number of samples in the VCF */
    inline int nSamples() const
    {
        return bcf_hdr_nsamples(hdr);
    }
};

/**
 * @class BcfRecord
 * @brief Object represents a variant record in the VCF, offering methods to access and modify fields.
 * @note  BcfRecord has to be associated with a BcfHeader object and needs to be filled in by calling
 *BcfReader.getNextVariant function.
 **/
class BcfRecord
{
    friend class BcfReader;
    friend class BcfWriter;

  private:
    BcfHeader * header;
    std::shared_ptr<bcf1_t> line = std::shared_ptr<bcf1_t>(bcf_init(), details::bcf_line_close()); // variant
    bcf_hdr_t * hdr_d = NULL; // a dup header by bcf_hdr_dup(header->hdr)
    bcf_fmt_t * fmt = NULL;
    bcf_info_t * info = NULL;
    int32_t * gts = NULL;
    int ndst, ret, nsamples;
    bool noneMissing = true; // whenever parsing a tag have to reset this variable
    bool isAllPhased = false;
    int nploidy = 0; // the number of ploidy
    int nvalues = 0; // the number of values for a tag in FORMAT

  public:
    /// if there is "." in GT for the sample, then it's coded as missing (TRUE)
    std::vector<char> isGenoMissing;

  public:
    /// empty constructor. call init() afterwards
    BcfRecord() {}

    /// constructor with a given BcfHeader object
    BcfRecord(BcfHeader & h)
    {
        initHeader(h);
    }

    ~BcfRecord()
    {
        if(gts) free(gts);
        if(hdr_d) bcf_hdr_destroy(hdr_d);
    }

    /// initilize the header associated with BcfRecord object by pointing to another BcfHeader object
    void initHeader(BcfHeader & h)
    {
        header = &h;
        if(!header->hdr) throw std::runtime_error("please initial header first\n");
        nsamples = header->nSamples();
        typeOfGT.resize(nsamples);
        gtPhase.resize(nsamples, 0);
    }

    /// reset the header associated with BcfRecord object by pointing to another BcfHeader object
    void resetHeader(BcfHeader & h)
    {
        header = &h;
    }

    /** @brief stream out the variant */
    friend std::ostream & operator<<(std::ostream & out, const BcfRecord & v)
    {
        out << v.asString();
        return out;
    }

    /** @brief return current variant as raw string */
    inline std::string asString() const
    {
        kstring_t s = {0, 0, NULL}; // kstring
        if(vcf_format(header->hdr, line.get(), &s) == 0)
        {
            std::string out(s.s, s.l);
            free(s.s);
            return out;
        }
        else
        {
            throw std::runtime_error("couldn't format current record into a string.\n");
        }
    }

    /**
     * @brief fill in the input vector with genotypes of 0 and 1. only works for ploidy<=2. Genotypes with
     * missing allele is coded as heterozygous
     * @param v valid input includes vector<bool> and vector<char> type
     * @return bool
     * @note  use isNoneMissing() to check if all genotypes are with no missingness. Alternatively, one can
     * use vector<int> as the input type as noted in the other overloading function getGenotypes().
     * */
    template<typename T>
    isValidGT<T> getGenotypes(T & v)
    {
        ndst = 0;
        ret = bcf_get_genotypes(header->hdr, line.get(), &gts, &ndst);
        if(ret <= 0)
        {
#    if defined(VERBOSE)
            std::cerr << "GT not present for current site. did you initilize the variant object?\n";
#    endif
            return false;
        }
        // if nploidy is not set manually. find the max nploidy using the first variant (eg. 2) resize v as
        // max(nploidy)
        // * nsamples (ret) NOTE: if ret == nsamples and only one sample then haploid
        if(ret != nploidy * nsamples && nploidy > 0)
        {
            // rare case if nploidy is set manually. eg. only one sample. the first variant is '1' but the
            // second is '1/0'. ret = 1 but nploidy should be 2
            v.resize(nploidy * nsamples);
        }
        else
        {
            v.resize(ret);
            nploidy = ret / nsamples;
        }
        // work with nploidy == 1, haploid?
        isGenoMissing.assign(nsamples, 0);
        int i, j, nphased = 0;
        noneMissing = true;
        fmt = bcf_get_fmt(header->hdr, line.get(), "GT");
        int nploidy_cur = ret / nsamples; // requires nploidy_cur <= nploidy
        for(i = 0; i < nsamples; i++)
        {
            // check and fill in typeOfGT; only supports SNPs now. check vcfstats.c for inspiration
            typeOfGT[i] = bcf_gt_type(fmt, i, NULL, NULL);
            if(typeOfGT[i] == GT_UNKN)
            {
                noneMissing = false; // set missing as het, user should check if isNoneMissing();
                isGenoMissing[i] = 1;
                v[i * nploidy + 0] = 1;
                for(j = 1; j < nploidy_cur; j++) v[i * nploidy + j] = 0;
                continue;
            }

            for(j = 0; j < nploidy_cur; j++)
            {
                // TODO: right now only parse 0 and 1, ie max(nploidy)=2; other values 2,3... will be set to 1
                v[i * nploidy + j] = bcf_gt_allele(gts[j + i * nploidy_cur]) != 0;
            }
            if(nploidy == 2)
            {
                gtPhase[i] = (gts[1 + i * nploidy_cur] & 1) == 1;
                nphased += gtPhase[i];
            }
        }
        if(nphased == nsamples)
        {
            isAllPhased = true;
        }
        else
        {
            isAllPhased = false;
        }
        return true;
    }

    /**
     * @brief fill in the input vector with genotype values, 0, 1 or -9 (missing).
     * @param v valid input is vector<int> type
     * @return bool
     * @note this function provides full capability to handle all kinds of genotypes
     * in multi-ploidy data costing more spae than the other function. missing allele is set as -9.
     * */
    bool getGenotypes(std::vector<int> & v)
    {
        ndst = 0;
        ret = bcf_get_genotypes(header->hdr, line.get(), &gts, &ndst);
        if(ret <= 0)
        {
#    if defined(VERBOSE)
            std::cerr << "GT not present for current site. did you initilize the variant object?\n";
#    endif
            return false;
        }
        v.resize(ret);
        isGenoMissing.assign(nsamples, 0);
        nploidy = ret / nsamples;
        int i, j, nphased = 0;
        noneMissing = true;
        for(i = 0; i < nsamples; i++)
        {
            int32_t * ptr = gts + i * nploidy;
            int is_phased = 0;
            for(j = 0; j < nploidy; j++)
            {
                // if true, the sample has smaller ploidy
                if(ptr[j] == bcf_int32_vector_end) break;
                // missing allele
                if(bcf_gt_is_missing(ptr[j]))
                {
                    noneMissing = false;
                    isGenoMissing[i] = 1;
                    v[i * nploidy + j] = -9;
                    continue;
                }
                v[i * nploidy + j] = bcf_gt_allele(ptr[j]);
                is_phased += bcf_gt_is_phased(ptr[j]);
            }
            if(nploidy == is_phased)
            {
                gtPhase[i] = true;
                nphased += gtPhase[i];
            }
        }
        if(nphased == nsamples)
        {
            isAllPhased = true;
        }
        else
        {
            isAllPhased = false;
        }
        return true;
    }

    /**
     * @brief get tag value in FORMAT
     * @param tag valid tag name in FORMAT column declared in the VCF header
     * @param v valid input include vector of float, char, int type
     * @return bool
     * */
    template<typename T, typename S = typename T::value_type>
    isFormatVector<T> getFORMAT(std::string tag, T & v)
    {
        fmt = bcf_get_fmt(header->hdr, line.get(), tag.c_str());
        if(!fmt)
        {
            throw std::invalid_argument("no FORMAT=" + tag + " in the VCF header.\n");
        }
        nvalues = fmt->n;
        ndst = 0;
        S * dst = NULL;
        int tagid = header->getFormatType(tag);
        if(tagid == 1)
        {
            ret = bcf_get_format_int32(header->hdr, line.get(), tag.c_str(), &dst, &ndst);
        }
        else if(tagid == 2)
        {
            ret = bcf_get_format_float(header->hdr, line.get(), tag.c_str(), &dst, &ndst);
        }
        else if(tagid == 3)
        {
            ret = bcf_get_format_char(header->hdr, line.get(), tag.c_str(), &dst, &ndst);
        }
        else
        {
            ret = -1;
        }

        if(ret >= 0)
        {
            // user have to check if there is missing in the return v;
            v = std::vector<S>(dst, dst + ret);
            free(dst);
            return true;
        }
        else
        {
            free(dst);
            return false;
        }
    }

    /**
     * @brief get tag value in FORMAT
     * @param tag valid tag name in FORMAT column declared in the VCF header
     * @param v vector of string
     * @return bool
     * */
    bool getFORMAT(std::string tag, std::vector<std::string> & v)
    {
        fmt = bcf_get_fmt(header->hdr, line.get(), tag.c_str());
        if(!fmt)
        {
            throw std::invalid_argument("no FORMAT=" + tag + " in the VCF header.\n");
        }
        nvalues = fmt->n;
        // if ndst < (fmt->n+1)*nsmpl; then realloc is involved
        ret = -1, ndst = 0;
        char ** dst = NULL;
        if(header->getFormatType(tag) == 3)
            ret = bcf_get_format_string(header->hdr, line.get(), tag.c_str(), &dst, &ndst);
        if(ret > 0)
        {
            v.clear();
            for(int i = 0; i < nsamples; i++)
            {
                // Ugly: GT field is considered to be a string by the VCF header but BCF represents it as INT.
                v.emplace_back(dst[i]);
            };
            free(dst[0]);
            free(dst);
            return true;
        }
        else
        {
            free(dst[0]);
            free(dst);
            return false;
        }
    }

    /**
     * @brief get tag value in INFO
     * @param tag valid tag name in INFO column declared in the VCF header
     * @param v valid input include vector of float, int type
     * @return bool
     * */
    template<typename T, typename S = typename T::value_type>
    isInfoVector<T> getINFO(std::string tag, T & v)
    {
        info = bcf_get_info(header->hdr, line.get(), tag.c_str());
        if(!info)
        {
            throw std::invalid_argument("no INFO=" + tag + " in the VCF header.\n");
        }
        S * dst = NULL;
        ndst = 0;
        ret = -1;
        if(info->type == BCF_BT_INT8 || info->type == BCF_BT_INT16 || info->type == BCF_BT_INT32)
        {
            ret = bcf_get_info_int32(header->hdr, line.get(), tag.c_str(), &dst, &ndst);
        }
        else if(info->type == BCF_BT_FLOAT)
        {
            ret = bcf_get_info_float(header->hdr, line.get(), tag.c_str(), &dst, &ndst);
        }
        if(ret >= 0)
        {
            v = std::vector<S>(dst, dst + ret); // user have to check if there is missing in the return v;
            free(dst);
            return true;
        }
        else
        {
            free(dst);
            return false;
        }
    }

    /**
     * @brief get tag value in INFO
     * @param tag valid tag name in INFO column declared in the VCF header
     * @param v valid input include scalar value of float or int type
     * @return bool
     * */
    template<typename T>
    isScalar<T> getINFO(std::string tag, T & v)
    {
        info = bcf_get_info(header->hdr, line.get(), tag.c_str());
        if(!info)
        {
            throw std::invalid_argument("no INFO=" + tag + " in the VCF header.\n");
        }
        // scalar value
        if(info->len == 1)
        {
            if(info->type == BCF_BT_INT8 || info->type == BCF_BT_INT16 || info->type == BCF_BT_INT32)
            {
                v = info->v1.i;
            }
            else if(info->type == BCF_BT_FLOAT)
            {
                v = info->v1.f;
            }
            return true;
        }
        else
        {
#    if defined(VERBOSE)
            std::cerr << "there are multiple values for " + tag
                             + " in INFO for current site. please use vector instead\n";
#    endif
            return false;
        }
    }

    /**
     * @brief get tag value in INFO
     * @param tag valid tag name in INFO column declared in the VCF header
     * @param v valid input is std::string
     * @return bool
     * */
    template<typename T>
    isString<T> getINFO(std::string tag, T & v)
    {
        info = bcf_get_info(header->hdr, line.get(), tag.c_str());
        if(!info)
        {
            throw std::invalid_argument("no INFO=" + tag + " in the VCF header.\n");
        }
        if(info->type == BCF_BT_CHAR)
        {
            v = std::string(reinterpret_cast<char *>(info->vptr), info->vptr_len);
            return true;
        }
        else
        {
            return false;
        }
    }

    /**
     * @brief set tag value for INFO
     * @param tag valid tag name in INFO column declared in the VCF header
     * @param v valid input include scalar value of float or int type
     * @return bool
     * */
    template<typename T>
    isScalar<T> setINFO(std::string tag, const T & v)
    {
        // bcf_hrec_set_val
        // bcf_update_info_flag(header.hdr, line, tag.c_str(), NULL, 1);
        int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
        if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_INT & 0xff))
        {
            ret = bcf_update_info_int32(header->hdr, line.get(), tag.c_str(), &v, 1);
        }
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_REAL & 0xff))
        {
            float v2 = static_cast<float>(v);
            ret = bcf_update_info_float(header->hdr, line.get(), tag.c_str(), &v2, 1);
        }
        else
        {
            ret = -1;
        }
        if(ret < 0)
        {
#    if defined(VERBOSE)
            std::cerr << "couldn't set " + tag + " for this variant.\nplease add the tag in headerfirst.\n";
#    endif
            return false;
        }
        return true;
    }

    /**
     * @brief set tag value for INFO
     * @param tag valid tag name in INFO column declared in the VCF header
     * @param v valid input include vector<int> vector<float> std::string
     * @return bool
     * */
    template<typename T>
    isValidInfo<T> setINFO(std::string tag, const T & v)
    {
        // bcf_update_info_flag(header.hdr, line, tag.c_str(), NULL, 1);
        int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
        if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_INT & 0xff))
        {
            ret = bcf_update_info_int32(header->hdr, line.get(), tag.c_str(), v.data(), v.size());
        }
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_REAL & 0xff))
        {
            ret = bcf_update_info_float(header->hdr, line.get(), tag.c_str(), v.data(), v.size());
        }
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_STR & 0xff))
        {
            ret = bcf_update_info_string(header->hdr, line.get(), tag.c_str(), v.data());
        }
        else
        {
            ret = -1;
        }

        if(ret < 0)
        {
#    if defined(VERBOSE)
            std::cerr << "couldn't set " + tag + " for this variant.\nplease add the tag in headerfirst.\n";
#    endif
            return false;
        }
        return true;
    }

    /// remove the given tag from INFO of the variant
    void removeINFO(std::string tag)
    {
        int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
        if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_INT & 0xff))
        {
            ret = bcf_update_info_int32(header->hdr, line.get(), tag.c_str(), NULL, 0);
        }
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_REAL & 0xff))
        {
            ret = bcf_update_info_float(header->hdr, line.get(), tag.c_str(), NULL, 0);
        }
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_INFO, tag_id) == (BCF_HT_STR & 0xff))
        {
            ret = bcf_update_info_string(header->hdr, line.get(), tag.c_str(), NULL);
        }
        else
        {
            ret = -1;
        }

        if(ret < 0)
        {
            throw std::runtime_error("couldn't remove " + tag + " for this variant.\n");
        }
    }

    /**
     * @brief set genotypes from scratch even if genotypes not present
     * @param v  the genotypes of vector<int> type
     * @return bool
     * */
    bool setGenotypes(const std::vector<int> & v)
    {
        // bcf_gt_type
        int i, j, k;
        nploidy = v.size() / nsamples;
        int32_t * gt = (int32_t *)malloc(v.size() * sizeof(int32_t));
        for(i = 0; i < nsamples; i++)
        {
            for(j = 0; j < nploidy; j++)
            {
                k = i * nploidy + j;
                if(v[k] == -9 || v[k] == bcf_int32_missing)
                {
                    gt[k] = bcf_gt_missing;
                }
                else if(gtPhase[i])
                {
                    gt[k] = bcf_gt_phased(v[k]);
                }
                else
                {
                    gt[k] = bcf_gt_unphased(v[k]);
                }
            }
        }
        if(bcf_update_genotypes(header->hdr, line.get(), gt, v.size()) < 0)
        {
            free(gt);
#    if defined(VERBOSE)
            std::cerr << "couldn't set genotypes correctly.\n";
#    endif
            return false;
        }
        free(gt);
        return true;
    }

    /**
     * @brief set phasing status for all diploid samples using given vector
     * @param v valid input includes vector<char>
     * */
    void setPhasing(const std::vector<char> & v)
    {
        if((int)v.size() != nsamples)
            throw std::runtime_error("the size of input vector is not matching the size of genotypes");
        gtPhase = v;
    }

    /// remove the given tag from FORMAT of the variant
    void removeFORMAT(std::string tag)
    {
        ret = -1;
        int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
        if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_INT & 0xff))
        {
            ret = bcf_update_format_int32(header->hdr, line.get(), tag.c_str(), NULL, 0);
        }
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_STR & 0xff))
        {
            ret = bcf_update_format_char(header->hdr, line.get(), tag.c_str(), NULL, 0);
        }
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_REAL & 0xff))
        {
            ret = bcf_update_format_float(header->hdr, line.get(), tag.c_str(), NULL, 0);
        }

        if(ret < 0) throw std::runtime_error("couldn't remove " + tag + " correctly.\n");
    }

    /**
     * @brief set tag values for all samples in FORMAT using given vector
     * @param tag valid tag name in FORMAT column declared in the VCF header
     * @param v valid input include vector<int>, vector<float>, vector<char>, std::string
     * @return bool
     * */
    template<typename T>
    isValidFMT<T> setFORMAT(std::string tag, const T & v)
    {
        int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
        if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_INT & 0xff))
        {
            ret = bcf_update_format_int32(header->hdr, line.get(), tag.c_str(), v.data(), v.size());
        }
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_STR & 0xff))
        {
            ret = bcf_update_format_char(header->hdr, line.get(), tag.c_str(), v.data(), v.size());
        }
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_REAL & 0xff))
        {
            ret = bcf_update_format_float(header->hdr, line.get(), tag.c_str(), v.data(), v.size());
        }
        else
        {
            ret = -1;
        }

        if(ret < 0)
        {
#    if defined(VERBOSE)
            std::cerr << "couldn't set format " + tag + " corectly.\n";
#    endif
            return false;
        }
        return true;
    }

    /**
     * @brief set tag for a single sample in FORMAT using given singular value. this works only when there is
     * one sample in the vcf
     * @param tag valid tag name in FORMAT column declared in the VCF header
     * @param v valid input include int, float or double
     * @return void
     * */
    template<typename T>
    isScalar<T> setFORMAT(std::string tag, const T & v)
    {
        float v2 = v;
        int tag_id = bcf_hdr_id2int(header->hdr, BCF_DT_ID, tag.c_str());
        if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_INT & 0xff))
        {
            ret = bcf_update_format_int32(header->hdr, line.get(), tag.c_str(), &v, 1);
        }
        else if(bcf_hdr_id2type(header->hdr, BCF_HL_FMT, tag_id) == (BCF_HT_REAL & 0xff))
        {
            ret = bcf_update_format_float(header->hdr, line.get(), tag.c_str(), &v2, 1);
        }
        else
        {
            ret = -1;
        }
        if(ret < 0)
        {
#    if defined(VERBOSE)
            std::cerr << "couldn't set format " + tag + " corectly.\n";
#    endif
            return false;
        }
        return true;
    }

    /** @brief add one variant record from given string*/
    void addLineFromString(const std::string & vcfline)
    {
        kstring_t s = {0, 0, NULL};
        kputsn(vcfline.c_str(), vcfline.length(), &s);
        ret = vcf_parse1(&s, header->hdr, line.get());
        free(s.s);
        if(ret > 0) throw std::runtime_error("error parsing: " + vcfline + "\n");
        if(line->errcode == BCF_ERR_CTG_UNDEF)
        {
            std::string contig(bcf_hdr_id2name(header->hdr, line->rid));
            hdr_d = bcf_hdr_dup(header->hdr);
            header->hrec = bcf_hdr_id2hrec(hdr_d, BCF_DT_CTG, 0, line->rid);
            if(header->hrec == NULL)
                throw std::runtime_error("contig" + contig + " unknow and not found in the header.\n");
            ret = bcf_hdr_add_hrec(header->hdr, header->hrec);
            if(ret < 0) throw std::runtime_error("error adding contig " + contig + " to header.\n");
            ret = bcf_hdr_sync(header->hdr);
        }
    }

    /** @brief if all samples have non missing values for the tag in FORMAT */
    inline bool isNoneMissing() const
    {
        return noneMissing;
    }

    /** @brief return boolean value indicates if current variant is Structual Variant or not */
    inline bool isSV() const
    {
        if(bcf_get_info(header->hdr, line.get(), "SVTYPE") == NULL)
        {
            return false;
        }
        else
        {
            return true;
        }
    }

    /** @brief return boolean value indicates if current variant is exclusively INDEL */
    inline bool isIndel() const
    {
        // REF has multiple allels
        if(REF().length() > 1 && !isSV()) return true;
        for(int i = 1; i < line->n_allele; i++)
        {
            std::string alt(line->d.allele[i]);
            if(alt == ".") return true;
            if(alt.length() != REF().length() && !isSV()) return true;
        }
        return false;
    }

    /** @brief return boolean value indicates if current variant is multiallelic sites */
    inline bool isMultiAllelics() const
    {
        if(line->n_allele <= 2) return false;
        return true;
    }

    /** @brief return boolean value indicates if current variant is exclusively multiallelic SNP sites */
    inline bool isMultiAllelicSNP() const
    {
        // skip REF with length > 1, i.e. INDEL
        if(REF().length() > 1 || line->n_allele <= 2) return false;
        for(int i = 1; i < line->n_allele; i++)
        {
            std::string snp(line->d.allele[i]);
            if(snp.length() != 1)
            {
                return false;
            }
        }
        return true;
    }

    /** @brief return boolean value indicates if current variant is exclusively biallelic SNP. Note ALT=* are
     * skipped */
    inline bool isSNP() const
    {
        // REF and ALT have multiple allels
        if(REF().length() > 1 || line->n_allele > 2) return false;
        std::string snp(line->d.allele[1]);
        if(!(snp == "A" || snp == "C" || snp == "G" || snp == "T"))
        {
            return false;
        }
        return true;
    }

    /** @brief return boolean value indicates if current variant has SNP type defined in vcf.h (htslib>=1.16)
     */
    inline bool hasSNP() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_SNP, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /// return boolean value indicates if current variant has INDEL type defined in vcf.h (htslib>=1.16)
    inline bool hasINDEL() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_INDEL, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /// return boolean value indicates if current variant has INS type defined in vcf.h (htslib>=1.16)
    inline bool hasINS() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_INS, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /// return boolean value indicates if current variant has DEL type defined in vcf.h (htslib>=1.16)
    inline bool hasDEL() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_DEL, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /// return boolean value indicates if current variant has MNP type defined in vcf.h (htslib>=1.16)
    inline bool hasMNP() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_MNP, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /// return boolean value indicates if current variant has BND type defined in vcf.h (htslib>=1.16)
    inline bool hasBND() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_BND, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /// return boolean value indicates if current variant has OTHER type defined in vcf.h (htslib>=1.16)
    inline bool hasOTHER() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_OTHER, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /// return boolean value indicates if current variant has OVERLAP type defined in vcf.h (htslib>=1.16)
    inline bool hasOVERLAP() const
    {
        int type = bcf_has_variant_types(line.get(), VCF_OVERLAP, bcf_match_overlap);
        if(type < 0) throw std::runtime_error("something wrong with variant type\n");
        if(type == 0) return false;
        return true;
    }

    /** @brief return CHROM name */
    inline std::string CHROM() const
    {
        return std::string(bcf_hdr_id2name(header->hdr, line->rid));
    }

    /** @brief return ID field */
    inline std::string ID() const
    {
        return std::string(line->d.id);
    }

    /** @brief return 1-base position */
    inline int64_t POS() const
    {
        return line->pos + 1;
    }

    /** @brief modify CHROM value */
    inline void setCHR(const std::string & s)
    {
        line->rid = bcf_hdr_name2id(header->hdr, s.c_str());
    }

    /** @brief modify position given 1-based value */
    inline void setPOS(int64_t p)
    {
        line->pos = p - 1;
    }

    /** @brief update ID */
    inline void setID(const std::string & s)
    {
        bcf_update_id(header->hdr, line.get(), s.c_str());
    }

    /** @brief set REF and ALT alleles given a string seperated by comma */
    inline void setRefAlt(const std::string & s)
    {
        bcf_update_alleles_str(header->hdr, line.get(), s.c_str());
    }

    /** @brief modify the QUAL value */
    inline void setQUAL(float q)
    {
        line->qual = q;
    }

    /** @brief modify the QUAL value */
    inline void setQUAL(char q)
    {
        bcf_float_set_missing(line->qual);
    }

    /** @brief modify the FILTER value */
    inline void setFILTER(const std::string & s)
    {
        int32_t tmpi = bcf_hdr_id2int(header->hdr, BCF_DT_ID, s.c_str());
        bcf_add_filter(header->hdr, line.get(), tmpi);
    }

    /** @brief return 0-base start of the variant (can be any type) */
    inline int64_t Start() const
    {
        return line->pos;
    }

    /** @brief return 0-base end of the variant (can be any type) */
    inline int64_t End() const
    {
        return line->pos + line->rlen;
    }

    /** @brief return raw REF alleles as string */
    inline std::string REF() const
    {
        return std::string(line->d.allele[0]);
    }

    /** @brief swap REF and ALT for biallelic SNP */
    inline void swap_REF_ALT()
    {
        if(!isMultiAllelicSNP()) std::swap(line->d.allele[0], line->d.allele[1]);
    }

    /** @brief return raw ALT alleles as string */
    inline std::string ALT() const
    {
        std::string s;
        for(int i = 1; i < line->n_allele; i++)
        {
            s += std::string(line->d.allele[i]) + ",";
        }
        if(s.length() > 1) s.pop_back();
        return s;
    }

    /** @brief return QUAL value */
    inline float QUAL()
    {
        if(bcf_float_is_missing(line->qual))
        {
            noneMissing = false;
            return bcf_float_missing;
        }
        else
        {
            return line->qual;
        }
    }

    /** @brief return raw FILTER column as string */
    inline std::string FILTER()
    {
        if(line->d.n_flt == 0)
        {
            return ".";
        }
        else if(line->d.n_flt == 1)
        {
            return std::string(bcf_hdr_int2id(header->hdr, BCF_DT_ID, line->d.flt[0]));
        }
        else
        {
            std::string s;
            for(int i = 1; i < line->d.n_flt; i++)
            {
                s += std::string(bcf_hdr_int2id(header->hdr, BCF_DT_ID, line->d.flt[i])) + ",";
            }
            s.pop_back();
            return s;
        }
    }

    /** @brief return raw INFO column as string. recommend to use getINFO for specific tag. */
    inline std::string allINFO()
    {
        int32_t max_dt_id = header->hdr->n[BCF_DT_ID];
        kstring_t * s = (kstring_t *)calloc(1, sizeof(kstring_t));
        if(line->n_info)
        {
            int first = 1;
            for(int i = 0; i < line->n_info; ++i)
            {
                bcf_info_t * z = &line->d.info[i];
                if(!z->vptr) continue;
                if(!first) kputc(';', s);
                first = 0;
                if(z->key < 0 || z->key >= max_dt_id || header->hdr->id[BCF_DT_ID][z->key].key == NULL)
                    throw std::runtime_error("Invalid BCF and wrong INFO tag");
                kputs(header->hdr->id[BCF_DT_ID][z->key].key, s);
                if(z->len <= 0) continue;
                kputc('=', s);
                if(z->len == 1)
                {
                    switch(z->type)
                    {
                        case BCF_BT_INT8:
                            if(z->v1.i == bcf_int8_missing)
                                kputc('.', s);
                            else
                                kputw(z->v1.i, s);
                            break;
                        case BCF_BT_INT16:
                            if(z->v1.i == bcf_int16_missing)
                                kputc('.', s);
                            else
                                kputw(z->v1.i, s);
                            break;
                        case BCF_BT_INT32:
                            if(z->v1.i == bcf_int32_missing)
                                kputc('.', s);
                            else
                                kputw(z->v1.i, s);
                            break;
                        case BCF_BT_INT64:
                            if(z->v1.i == bcf_int64_missing)
                                kputc('.', s);
                            else
                                kputll(z->v1.i, s);
                            break;
                        case BCF_BT_FLOAT:
                            if(bcf_float_is_missing(z->v1.f))
                                kputc('.', s);
                            else
                                kputd(z->v1.f, s);
                            break;
                        case BCF_BT_CHAR:
                            kputc(z->v1.i, s);
                            break;
                        default:
                            throw std::runtime_error("Unexpected type in INFO");
                    }
                }
                else
                    bcf_fmt_array(s, z->len, z->type, z->vptr);
            }
            if(first) kputc('.', s);
        }
        else
            kputc('.', s);
        std::string out = std::string(s->s, s->l);
        free(s->s);
        free(s);
        return out;
    }

    /** @brief return boolean value indicates if genotypes of all samples are phased */
    inline bool allPhased() const
    {
        return isAllPhased;
    }

    /** @brief return the number of ploidy of current variant */
    inline int ploidy() const
    {
        return nploidy;
    }

    /** @brief in a rare case, one may want to set the number of ploidy manually */
    inline void setPloidy(int v)
    {
        nploidy = v;
    }

    /**
     * @brief vector of nsamples length. keep track of the type of genotype (one of GT_HOM_RR, GT_HET_RA,
     *        GT_HOM_AA, GT_HET_AA, GT_HAPL_R, GT_HAPL_A or GT_UNKN).
     * @note GT_HOM_RR 0 \n
     *       GT_HOM_AA 1 \n
     *       GT_HET_RA 2 \n
     *       GT_HET_AA 3 \n
     *       GT_HAPL_R 4 \n
     *       GT_HAPL_A 5 \n
     *       GT_UNKN   6 \n
     * */
    std::vector<char> typeOfGT;

    /** @brief vector of nsamples length. keep track of the phasing status of each sample */
    std::vector<char> gtPhase;
};

/**
 * @class BcfReader
 * @brief Stream in variants from compressed/uncompressed VCF/BCF file or stdin
 **/
class BcfReader
{
  private:
    std::shared_ptr<htsFile> fp; // hts file
    std::shared_ptr<hts_idx_t> hidx; // hts index file
    std::shared_ptr<tbx_t> tidx; // .tbi .csi index file for vcf files
    std::shared_ptr<hts_itr_t> itr; // hts iterator
    kstring_t s = {0, 0, NULL}; // kstring
    std::string fname;
    bool isBcf = false; // if the input file is bcf or vcf;

  public:
    /// a BcfHeader object
    BcfHeader header;
    /// number of samples in the VCF
    int nsamples;
    /// a vector of samples name in the VCF
    std::vector<std::string> SamplesName;

    /// Construct an empty BcfReader
    BcfReader() {}

    /**
     *  @brief construct a vcf/bcf reader from file.
     *  @param file   the input vcf/bcf with suffix vcf(.gz)/bcf(.gz) or stdin "-"
     */
    BcfReader(const std::string & file) : fname(file)
    {
        open(file);
    }

    /**
     *  @brief construct a vcf/bcf reader with subset samples
     *  @param file   the input vcf/bcf with suffix vcf(.gz)/bcf(.gz) or stdin "-"
     *  @param region samtools-like region "chr:start-end", skip if empty
     */
    BcfReader(const std::string & file, const std::string & region) : fname(file)
    {
        open(file);
        if(!region.empty()) setRegion(region);
        SamplesName = header.getSamples();
    }

    /**
     *  @brief construct a vcf/bcf reader with subset samples in target region
     *  @param file   the input vcf/bcf with suffix vcf(.gz) or bcf(.gz)
     *  @param region samtools-like region "chr:start-end", skip if empty
     *  @param samples  LIST samples to include or exclude as a comma-separated string. \n
     *                  LIST : select samples in list \n
     *                  ^LIST : exclude samples from list \n
     *                  "-" : include all samples \n
     *                  "" : exclude all samples
     */
    BcfReader(const std::string & file, const std::string & region, const std::string & samples) : fname(file)
    {
        open(file);
        if(!region.empty()) setRegion(region);
        if(!samples.empty()) setSamples(samples);
    }

    ~BcfReader()
    {
        if(s.s) free(s.s);
    }

    /// Close the VCF file and its associated files
    void close()
    {
        if(fp) fp.reset();
        if(hidx) hidx.reset();
        if(itr) itr.reset();
        if(tidx) tidx.reset();
    }

    /// Open a VCF/BCF/STDIN file for streaming in
    void open(const std::string & file)
    {
        if(!fname.empty() && fname != file)
            throw std::runtime_error("does not support re-open a file yet. " + fname);
        fname = file;
        fp = std::shared_ptr<htsFile>(hts_open(fname.c_str(), "r"), details::hts_file_close());
        if(!fp) throw std::invalid_argument("I/O error: input file is invalid");
        enum htsExactFormat hts_format = hts_get_format(fp.get())->format;
        if(hts_format == bcf) isBcf = true;
        header.hdr = bcf_hdr_read(fp.get());
        nsamples = bcf_hdr_nsamples(header.hdr);
        SamplesName = header.getSamples();
    }

    /** @brief set the number of threads to use */
    inline int setThreads(int n)
    {
        return hts_set_threads(fp.get(), n);
    }

    /// return a BcfHeader object
    const BcfHeader & getHeader() const
    {
        return header;
    }

    /**
     * @brief query the status of a given region in the VCF
     * @return -2: the region is not a valid bcftools-like format,
     *             or it is not presenting in the VCF even though it's bcftols-like format. \n
     *         -1: there is no index file found. \n
     *          0: the region is valid but empty. \n
     *          1: vaild and not empty. \n
     */
    int getStatus(const std::string & region)
    {
        try
        {
            setRegion(region);
            BcfRecord v(header);
            if(!getNextVariant(v)) return 0;
        }
        catch(const std::invalid_argument & e)
        {
            return -1;
        }
        catch(const std::runtime_error & e)
        {
            return -2;
        }
        return 1;
    }

    /**
     * @brief count the number of variants by iterating through a given region.
     * @note If you want to continue work on that region, remember to reset the region by setRegion()! \n
     *        Also, check the status of the region first to handle the different cases!
     */
    int getVariantsCount(const std::string & region)
    {
        int c{0};
        setRegion(region);
        BcfRecord v(header);
        while(getNextVariant(v)) c++;
        return c;
    }

    /**
     * @brief explicitly stream to specific samples
     * @param samples the string is bcftools-like format, which is comma separated list of samples to include
     * (or exclude with "^" prefix).
     * */
    void setSamples(const std::string & samples)
    {
        header.setSamples(samples);
        nsamples = bcf_hdr_nsamples(header.hdr);
        SamplesName = header.getSamples();
    }

    /**
     * @brief explicitly stream to specific region. throw invalid_argument error if index file not found.
     * throw runtime_error if the region was not a valid bcftools-like format or was not presenting in the
     * VCF.
     * @param region the string for region is samtools-like format, which can be 'chr', 'chr:start' and
     * 'chr:start-end'
     * */
    void setRegion(const std::string & region)
    {
        // 1. check and load index first
        // 2. query iterval region
        // 3. if region is empty, use "."
        if(isBcf)
        {
            hidx = std::shared_ptr<hts_idx_t>(bcf_index_load(fname.c_str()), details::hts_idx_close());
            if(itr) itr.reset(); // reset current region.
            if(region.empty())
                itr = std::shared_ptr<hts_itr_t>(bcf_itr_querys(hidx.get(), header.hdr, "."),
                                                 details::hts_iter_close());
            else
                itr = std::shared_ptr<hts_itr_t>(bcf_itr_querys(hidx.get(), header.hdr, region.c_str()),
                                                 details::hts_iter_close());
        }
        else
        {
            tidx = std::shared_ptr<tbx_t>(tbx_index_load(fname.c_str()), details::tabix_idx_close());
            if(tidx.get() == NULL) throw std::invalid_argument(" no tabix index found!");
            if(itr) itr.reset(); // reset
            if(region.empty())
                itr = std::shared_ptr<hts_itr_t>(tbx_itr_querys(tidx.get(), "."), details::hts_iter_close());
            else
                itr = std::shared_ptr<hts_itr_t>(tbx_itr_querys(tidx.get(), region.c_str()),
                                                 details::hts_iter_close());
        }
        if(itr.get() == NULL)
            throw std::runtime_error("region was not found! make sure the region format is correct");
    }

    /** @brief read in the next variant
     *  @param r the BcfRecord object to be filled in. */
    bool getNextVariant(BcfRecord & r)
    {
        int ret = -1;
        if(itr.get() != NULL)
        {
            if(isBcf)
            {
                ret = bcf_itr_next(fp.get(), itr.get(), r.line.get());
                bcf_subset_format(r.header->hdr, r.line.get()); // has to be called explicitly for bcf
                bcf_unpack(r.line.get(), BCF_UN_ALL);
                return (ret >= 0);
            }
            else
            {
                int slen = tbx_itr_next(fp.get(), tidx.get(), itr.get(), &s);
                if(slen > 0)
                {
                    ret = vcf_parse1(&s, r.header->hdr, r.line.get()); // ret > 0, error
                    bcf_unpack(r.line.get(), BCF_UN_ALL);
                }
                return (ret <= 0) && (slen > 0);
            }
        }
        else
        {
            ret = bcf_read(fp.get(), r.header->hdr, r.line.get());
            // unpack record immediately. not lazy
            bcf_unpack(r.line.get(), BCF_UN_ALL);
            return (ret == 0);
        }
    }
};

/**
 * @class BcfWriter
 * @brief Stream out variants to compressed/uncompressed VCF/BCF file or stdout
 **/
class BcfWriter
{
  private:
    std::shared_ptr<htsFile> fp; // hts file
    std::shared_ptr<bcf1_t> b = std::shared_ptr<bcf1_t>(bcf_init(), details::bcf_line_close()); // variant
    int ret;
    bool isHeaderWritten = false;
    const BcfHeader * hp;

  public:
    /// header object initialized by initalHeader
    BcfHeader header;

    /// Construct an empty BcfWriter
    BcfWriter() {}

    /**
     * @brief          Open VCF/BCF file for writing. The format is infered from file's suffix
     * @param fname    The file name or "-" for stdin/stdout. For indexed files
     * @param version  The output header is constructed with the internal template given a specific version
     */
    BcfWriter(const std::string & fname, std::string version = "VCF4.1")
    {
        open(fname);
        initalHeader(version);
        initalHeader(header);
    }

    /**
     * @brief          Open VCF/BCF file for writing. The format is infered from file's suffix
     * @param fname    The file name or "-" for stdin/stdout. For indexed files
     * @param h        The output header is pointing to the given BcfHeader object
     */
    BcfWriter(const std::string & fname, const BcfHeader & h)
    {
        open(fname);
        initalHeader(h);
    }

    /**
     * @brief          Open VCF/BCF file for writing using given mode
     * @param fname    The file name or "-" for stdin/stdout. For indexed files
     * @param version  The output header is constructed with the internal template given a specific version
     * @param mode     Mode matching \n
     *                 [w]b  .. compressed BCF \n
     *                 [w]bu .. uncompressed BCF \n
     *                 [w]z  .. compressed VCF \n
     *                 [w]   .. uncompressed VCF
     */
    BcfWriter(const std::string & fname, const std::string & version, const std::string & mode)
    {
        open(fname, mode);
        initalHeader(version);
        initalHeader(header);
    }

    /**
     * @brief          Open VCF/BCF file for writing using given mode
     * @param fname    The file name or "-" for stdin/stdout. For indexed files
     * @param h        The output header is pointing to the given BcfHeader object
     * @param mode     Mode matching \n
     *                 [w]b  .. compressed BCF \n
     *                 [w]bu .. uncompressed BCF \n
     *                 [w]z  .. compressed VCF \n
     *                 [w]   .. uncompressed VCF
     */
    BcfWriter(const std::string & fname, const BcfHeader & h, const std::string & mode)
    {
        open(fname, mode);
        initalHeader(h);
    }

    ~BcfWriter() {}

    /**
     * @brief          Open VCF/BCF file for writing. The format is infered from file's suffix
     * @param fname    The file name or "-" for stdin/stdout. For indexed files
     */
    void open(const std::string & fname)
    {
        auto mode = details::getMode(fname, "w");
        fp = std::shared_ptr<htsFile>(hts_open(fname.c_str(), mode.c_str()), details::hts_file_close());
        if(!fp) throw std::invalid_argument("I/O error: input file is invalid");
    }

    /**
     * @brief          Open VCF/BCF file for writing using given mode
     * @param fname    The file name or "-" for stdin/stdout. For indexed files
     * @param mode     Mode matching \n
     *                 [w]b  .. compressed BCF \n
     *                 [w]bu .. uncompressed BCF \n
     *                 [w]z  .. compressed VCF \n
     *                 [w]   .. uncompressed VCF
     */
    void open(const std::string & fname, const std::string & mode)
    {
        fp = std::shared_ptr<htsFile>(hts_open(fname.c_str(), mode.c_str()), details::hts_file_close());
        if(!fp) throw std::invalid_argument("I/O error: input file is invalid");
    }

    /// close the BcfWriter object.
    void close()
    {
        if(!isHeaderWritten) writeHeader();
        if(b) b.reset();
        if(fp) fp.reset();
    }

    /// initial a VCF header using the internal template given a specific version. VCF4.1 is the default
    void initalHeader(std::string version = "VCF4.1")
    {
        header.hdr = bcf_hdr_init("w");
        header.setVersion(version);
    }

    /// initial a VCF header by pointing to header of another VCF
    void initalHeader(const BcfHeader & h)
    {
        hp = &h;
    }

    /// copy header of given VCF and restrict on samples. if samples=="", FORMAT removed and only sites left
    void copyHeader(const std::string & vcffile, std::string samples = "-")
    {
        htsFile * fp2 = hts_open(vcffile.c_str(), "r");
        if(!fp2) throw std::invalid_argument("I/O error: input file is invalid");
        if(samples == "")
        { // site-only
            bcf_hdr_t * hfull = bcf_hdr_read(fp2);
            header.hdr = bcf_hdr_subset(hfull, 0, 0, 0);
            bcf_hdr_remove(header.hdr, BCF_HL_FMT, NULL);
            bcf_hdr_destroy(hfull);
        }
        else
        {
            header.hdr = bcf_hdr_read(fp2);
            header.setSamples(samples);
        }
        hts_close(fp2);
        initalHeader(header);
    }

    /// copy a string to a vcf line
    void writeLine(const std::string & vcfline)
    {
        if(!isHeaderWritten && !writeHeader()) throw std::runtime_error("could not write header\n");
        kstring_t s = {0, 0, NULL};
        kputsn(vcfline.c_str(), vcfline.length(), &s);
        ret = vcf_parse1(&s, hp->hdr, b.get());
        free(s.s);
        if(ret > 0) throw std::runtime_error("error parsing: " + vcfline + "\n");
        if(b->errcode == BCF_ERR_CTG_UNDEF)
        {
            throw std::runtime_error("contig id " + (std::string)bcf_hdr_id2name(hp->hdr, b->rid)
                                     + " not found in the header. please run header->AddContig() first.\n");
        }
        ret = bcf_write(fp.get(), hp->hdr, b.get());
        if(ret != 0) throw std::runtime_error("error writing: " + vcfline + "\n");
    }

    /// streams out the header
    bool writeHeader()
    {
        ret = bcf_hdr_write(fp.get(), hp->hdr);
        if(ret == 0) return isHeaderWritten = true;
        return false;
    }

    /// streams out the given variant of BcfRecord type
    bool writeRecord(BcfRecord & v)
    {
        if(!isHeaderWritten) writeHeader();
        if(bcf_write(fp.get(), v.header->hdr, v.line.get()) < 0) return false;
        return true;
    }
};

} // namespace vcfpp

#endif // VCFPP_H_
