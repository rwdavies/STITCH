#ifndef SEQLIB_COMMON_H
#define SEQLIB_COMMON_H

/*! \mainpage SeqLib 1.0
 *
 * \section intro_sec Introduction
 *
 * SeqLib is a C++ package for querying BAM/SAM/CRAM files with HTSlib, performing
 * BWA-MEM operations in memory, and peforming sequence assembly with FermiKit.
 * See https://github.com/walaj/SeqLib for
 * full description.
 */

#include <string>
#include <vector>

/** HTSlib/BWA-MEM/BLAT/Fermi operations */
namespace SeqLib
{

    static const char RCOMPLEMENT_TABLE[128] = {' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
                                                ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
                                                ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 'T',
                                                ' ', 'G', ' ', ' ', ' ', 'C', ' ', ' ', ' ', ' ', ' ', ' ', 'N', ' ', ' ', ' ', ' ', ' ', 'A', ' ', ' ', ' ',
                                                ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', 't', ' ', 'g', ' ', ' ', ' ', 'c', ' ', ' ', ' ', ' ', ' ', ' ',
                                                'n', ' ', ' ', ' ', ' ', ' ', 'a', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};

    /*
    static const std::vector<std::string> CHR_NAME {"1", "2", "3", "4", "5", "6", "7", "8", "9",
        "10", "11", "12", "13", "14", "15", "16", "17",
        "18", "19", "20", "21", "22", "X", "Y", "M"};
    static const std::vector<std::string> CHR_NAME_NUM {"1", "2", "3", "4", "5", "6", "7", "8", "9",
        "10", "11", "12", "13", "14", "15", "16", "17",
        "18", "19", "20", "21", "22", "23", "24"};

    static const std::vector<int> CHR_LEN_VEC = {249250621, 243199373, 198022430, 191154276, //1-4
                                     180915260, 171115067, //5-6
                                     159138663, 146364022, 141213431, 135534747, 135006516, 133851895, //7-12
                                     115169878, 107349540, 102531392, 90354753,  81195210,  78077248, //13-18
                                     59128983,  63025520,  48129895,  51304566,  155270560, 59373566, //19-24
                                     16571}; //25

    static const std::vector<double> CHR_CUMSUM_WEIGHT_X = {0.08209014, 0.16218732, 0.22740558, 0.29036182, 0.34994586, 0.40630223, 0.45871420,
                                                            0.50691887, 0.55342720, 0.59806527, 0.64252937, 0.68661320, 0.72454415, 0.75989948,
                                                            0.79366797, 0.82342611, 0.85016757, 0.87588214, 0.89535614, 0.91611346, 0.93196494,
                                                            0.94886198, 1.00000000};

    static const std::vector<double> CHR_WEIGHT_X = {0.08209014, 0.08009718, 0.06521825, 0.06295624, //1-4
                                     0.05958404, 0.05635637, //5-6
                                     0.05241197, 0.04820467, 0.04650833, 0.04463807, 0.04446410, 0.04408383, //7-12
                                     0.03793095, 0.03535534, 0.03376849, 0.02975814, 0.02674146, 0.02571457,
                                     0.01947400, 0.02075732, 0.01585148, 0.01689705, 0.05113802};

    static const int CHR_LEN [25] = {249250621, 243199373, 198022430, 191154276, //1-4
                                     180915260, 171115067, //5-6
                                     159138663, 146364022, 141213431, 135534747, 135006516, 133851895, //7-12
                                     115169878, 107349540, 102531392, 90354753,  81195210,  78077248, //13-18
                                     59128983,  63025520,  48129895,  51304566,  155270560, 59373566, //19-24
                                     16571}; //25

    static const uint32_t CHR_CLEN [25] = {0, 249250621, 492449994,  690472424, 881626700, 1062541960, 1233657027,
                                           1392795690,1539159712,1680373143,1815907890,1950914406,2084766301,
                                           2199936179, 2307285719, 2409817111, 2500171864, 2581367074, 2659444322,
                                           2718573305, 2781598825, 2829728720, 2881033286, 3036303846, 3095677412};

    static const uint32_t genome_size_XY = 3095677411;

    static std::string const REFHG19 = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";

    static const int NONCENT_CHR [44] = {1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,
                                         11,11,12,12,13,14,15,16,16,17,17,18,18,19,19,20,20,21,21,22,23,23,24,24};

    */
} // namespace SeqLib

#endif
