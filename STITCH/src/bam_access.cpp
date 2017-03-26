#include <SeqLib/BamReader.h>

//' @export
// [[Rcpp::export]]
std::string get_header_using_SeqLib(std::string file_name) {
    SeqLib::BamReader reader;
    reader.Open(file_name);
    std::string header_text = reader.HeaderConcat();
    return header_text;
}



std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<int>, std::vector<int>, std::vector<std::string>, std::vector<int>, std::vector<std::string>, std::vector<std::string>> test_get_reads(std::string region, std::string file_name, std::string reference = "") {
    SeqLib::BamReader reader;
    reader.Open(file_name);
    if (reference != "") {
        reader.SetCramReference(reference);
    }
    SeqLib::GenomicRegion gr(region, reader.Header());
    reader.SetRegion(gr);
    SeqLib::BamRecord record;
    std::vector<std::string> qname;
    std::vector<std::string> strand;
    std::vector<int> pos;
    std::vector<int> mapq;
    std::vector<std::string> cigar;
    std::vector<int> isize;
    std::vector<std::string> seq;
    std::vector<std::string> qual;
    while(reader.GetNextRecord(record)) {
        if (record.DuplicateFlag())
            continue;
        qname.push_back(record.Qname());
	// convert boolean to +, - used by Rsamtools
	bool s = record.ReverseFlag();
	if (s) {
	  strand.push_back("-");
	} else {
	  strand.push_back("+");	  
	}
	pos.push_back(record.Position() + 1); // make 1-based
	mapq.push_back(record.MapQuality());
	cigar.push_back(record.CigarString());
	isize.push_back(record.InsertSize());
        seq.push_back(record.Sequence());
	qual.push_back(record.Qualities());
    }
    return std::make_tuple(qname, strand, pos, mapq, cigar, isize, seq, qual);
}
