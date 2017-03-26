#include <Rcpp.h>

// the following "workaround" is needed as I cannot include both
// seqLib and Rcpp in the same .cpp for reasons I do not understand

std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<int>, std::vector<int>, std::vector<std::string>, std::vector<int>, std::vector<std::string>, std::vector<std::string>> test_get_reads(std::string region, std::string file_name, std::string ref_fa);

//' @export
// [[Rcpp::export]]
Rcpp::List get_sample_data_from_SeqLib(std::string region, std::string file_name, std::string reference = "") {
    std::vector<std::string> qname;
    std::vector<std::string> strand;
    std::vector<int> pos;
    std::vector<int> mapq;
    std::vector<std::string> cigar;
    std::vector<int> isize;
    std::vector<std::string> seq;
    std::vector<std::string> qual;
    std::tie(qname, strand, pos, mapq, cigar, isize, seq, qual) = test_get_reads(region, file_name, reference) ;
    return Rcpp::List::create(
        Rcpp::Named("qname") = qname,
	Rcpp::Named("strand") = strand,
	Rcpp::Named("pos") = pos,
	Rcpp::Named("mapq") = mapq,
	Rcpp::Named("cigar") = cigar,
	Rcpp::Named("isize") = isize,
	Rcpp::Named("seq") = seq,
	Rcpp::Named("qual") = qual
    );
}

