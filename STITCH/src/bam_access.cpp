#include <vector>
#include <string>
#include <SeqLib/BamReader.h>
#include <SeqLib/RefGenome.h>


//' @export
// [[Rcpp::export]]
std::string query_region(const std::string& file_name, const std::string& chrom, int32_t p1, int32_t p2) {
    SeqLib::RefGenome rg;
    rg.LoadIndex(file_name);
    if (rg.IsEmpty())
        return "";
    return rg.QueryRegion(chrom, p1, p2);
}


//' @export
// [[Rcpp::export]]
std::string get_header_using_SeqLib(std::string file_name) {
    SeqLib::BamReader reader;
    reader.Open(file_name);
    std::string header_text = reader.HeaderConcat();
    return header_text;
}


//' @export
// [[Rcpp::export]]
int get_read_span(std::vector<int> cigarLengthVec, std::vector<std::string> cigarTypeVec) {
    int readLength = 0;
    int cigarLength;
    std::string cigarType;
    for(std::size_t iM=0; iM < cigarLengthVec.size(); iM++) {
        cigarLength = cigarLengthVec[iM];
        cigarType = cigarTypeVec[iM];
	    // keep if it is an M or D. basically,
        // keep X and = as well
        if((cigarType == "M") || (cigarType == "D") || (cigarType == "=") || (cigarType == "X")) {
	        readLength = readLength + cigarLength;
	    }
    }
    return readLength;
}



std::tuple<std::vector<int>, std::vector<std::string>> split_cigar(std::string cigarRead) {
    int p = 0; // now 0-based
    std::vector<int> cigarLengthVec;
    std::vector<std::string> cigarTypeVec;
    for(std::size_t j=0; j < cigarRead.length(); j++) {
        if (! isdigit( cigarRead.substr(j, 1)[0] )) {
	  cigarLengthVec.push_back( std::stoi(cigarRead.substr(p, j-p)) );
	  cigarTypeVec.push_back( cigarRead.substr(j, 1));
	  p = j + 1;
	}
    }
    return std::make_tuple(cigarLengthVec, cigarTypeVec);
}


// now, take cigarRead, and if first or last cigar is an S, operate on it
// we will want to return the original components, plus pos and seq,
//    which we may operate on
std::tuple<std::vector<int>, std::vector<std::string>, std::string, std::string, int> deal_with_soft_clipped_bases(std::vector<int> cigarLengthVec, std::vector<std::string> cigarTypeVec, std::string seq, std::string qual, int pos, bool useSoftClippedBases) {
    if (useSoftClippedBases == true) {
        if (cigarTypeVec[0] == "S") {
  	    // change pos to match soft clipped bases	  
	    pos = pos - cigarLengthVec[0];
	    cigarTypeVec[0] = "M";
	}
	if (cigarTypeVec[cigarTypeVec.size() - 1] == "S") {
	    // change to be a useable cigar. pos is fine
	    cigarTypeVec[cigarTypeVec.size() - 1] = "M";
	}
    } else {
        // remove bases specified in cigar
        if (cigarTypeVec[0] == "S") {
	  // string.erase: first entry is 0-based position, second is 1-based n to remove
	  seq.erase(0, cigarLengthVec[0]);
	  qual.erase(0, cigarLengthVec[0]);
	  cigarLengthVec.erase(cigarLengthVec.begin()); // remove first entry in vector
	  cigarTypeVec.erase(cigarTypeVec.begin());	  	  
	}
	if (cigarTypeVec[cigarTypeVec.size() - 1] == "S") {
	  cigarLengthVec.erase(cigarLengthVec.begin() + cigarLengthVec.size() - 1); // remove last entry in vector
	  cigarTypeVec.erase(cigarTypeVec.begin() + cigarTypeVec.size() - 1);
	}
    }
    return std::make_tuple(cigarLengthVec, cigarTypeVec, seq, qual, pos);
}



std::tuple<std::vector<std::string>, std::vector<std::string>, std::vector<int>, std::vector<int>, std::vector<std::string>, std::vector<int>, std::vector<std::string>, std::vector<std::string>> get_reads_from_seqLib(std::string region, std::string file_name, std::string reference = "") {
    SeqLib::BamReader reader;
    // have to set cram reference first before opening bam!
    if (reference != "") {
       reader.SetCramReference(reference);
    }
    reader.Open(file_name);
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


std::tuple<std::vector<int>, std::vector<int>, std::vector<int>, std::vector<int>, std::vector<std::string>, std::vector<std::string>, std::vector<std::string>, std::vector<int>, std::vector<int>, std::vector<std::string>, std::vector<int>, std::vector<int>> get_sampleReadsRaw_using_SeqLib(const bool useSoftClippedBases, const int bqFilter, const int iSizeUpperLimit, std::vector<std::string> ref, std::vector<std::string> alt, const int nSNPs, std::vector<int> L, std::string region, std::string file_name, std::string reference, const bool save_sampleReadsInfo, const bool use_bx_tag) {
    //
    // initialize SeqLib stuff
    //
    SeqLib::BamReader reader;
    // have to set cram reference first before opening bam!
    if (reference != "") {
        reader.SetCramReference(reference);
    }
    reader.Open(file_name);
    SeqLib::GenomicRegion gr(region, reader.Header());
    reader.SetRegion(gr);
    SeqLib::BamRecord record;
    int isize, mapq;
    std::string seq, seq_temp, qual, qual_temp;
    bool got_seq_already;
    //
    // declare output variables
    //
    std::vector<int> out_num_SNPs; // 0-based
    std::vector<int> out_BQs; // pred-scaled base qualities, scaled with <0 -> ref, >0 -> alt
    std::vector<int> out_SNP_pos; // 0-based position of SNPs
    std::vector<int> out_iRead; // 0-based read SNPs came from
    std::vector<int> out_readStart; // 1-based start of read
    std::vector<int> out_readEnd; // 1-based end of read
    // slightly separate
    std::vector<std::string> qname; // read name for reads with SNPs in them
    std::vector<std::string> bxtag; // for recording bx tag
    std::vector<std::string> strand; // strand for reads with SNPs in them
    //
    // if save_sampleReadsInfo
    //
    std::vector<std::string> qname_all; // read name for reads with SNPs in them
    std::vector<int> out_readStart_all; // 1-based start of read
    std::vector<int> out_readEnd_all; // 1-based end of read
    //
    // new variables
    //
    int x1, x2, y, iM, t, tMin, tMax, readStart, readEnd, readLength;
    int nSNPInRead = -1;
    std::vector<int> cigarLengthVec;
    std::vector<int> cigarLengthVec_temp;    
    std::vector<std::string> cigarTypeVec;
    std::vector<std::string> cigarTypeVec_temp;
    int maxnSNPInRead = 1000;
    std::vector<int> qualLocal(maxnSNPInRead);
    std::vector<int> posLocal(maxnSNPInRead); // there shouldnt be this many SNPs
    int refPosition, refOffset, strandOffset, refPosition_temp;
    int iNumOfMs, cigarLength;
    char s;
    std::string cigarType;
    int whileVar, localbq;
    //
    // for bx tag
    //
    const std::string bx_string_tag = "BX";
    std::string for_bx_tag_s;
    bool string_z_tag_result;
    //
    // begin
    //
    // want - tMin - first SNP after last SNP before read
    // also - tMax - first SNP before last SNP after read
    // if the distance is >=0 between them, it is in that read
    //
    tMin = 0; 
    tMax = -1; // left and right boundaries of what SNPs to look at
    //
    // replaces for(iRead=0; iRead<=numberOfReads-1; iRead++)
    //
    while(reader.GetNextRecord(record)) {
        if (record.DuplicateFlag())
            continue;
	// convert boolean to +, - used by Rsamtools
	//
	// now - skip if isize or bq is not met
	//
	isize = record.InsertSize();
	mapq = record.MapQuality();
	if ((mapq < bqFilter) | (abs(isize) > iSizeUpperLimit))
	  continue;
	//
	// okay, if we're here, we're considering this read
	//
	refPosition = record.Position() + 1; // 1-based
	std::string cigarRead = record.CigarString();
	// if unmapped, skip
	if (cigarRead == "") // this is how "*" cigarReads are returned by seqLib
	  continue;
	std::tie(cigarLengthVec, cigarTypeVec) = split_cigar(cigarRead);
	// okay, probably want to clean this up at some point, might be wrong way of doing this
	got_seq_already = false;
	if ((cigarTypeVec[0] == "S") || (cigarTypeVec[cigarTypeVec.size() - 1] == "S")) {
	  // get seq. will get again later, am assuming most reads don't have clips!
	  //std::cout << "inside clipping" << std::endl;
	  //std::cout << "seq=" << seq << std::endl;
	  //std::cout << "qual=" << qual << std::endl;	  
	  seq = record.Sequence();
	  qual = record.Qualities();	  
	  got_seq_already = true;
	  std::tie(cigarLengthVec_temp, cigarTypeVec_temp, seq_temp, qual_temp, refPosition_temp) = deal_with_soft_clipped_bases(cigarLengthVec, cigarTypeVec, seq, qual, refPosition, useSoftClippedBases);
	  // reset!
	  cigarLengthVec = cigarLengthVec_temp;
	  cigarTypeVec = cigarTypeVec_temp;
	  seq = seq_temp;
	  qual = qual_temp;	  
	  refPosition = refPosition_temp;
	  //std::cout << "seq=" << seq << std::endl; 	  	  
	  //std::cout << "qual=" << qual << std::endl; 	  	  
	}
	readLength = get_read_span(cigarLengthVec, cigarTypeVec);
	readStart = refPosition; // simple this way
	readEnd = readStart + readLength - 1; // includes final base	
	// get current read position
	whileVar=0;
	// for t-Min
	while(whileVar==0) {
	  // dont continue if too far
	  if(tMin<(nSNPs-1)) {
	    // continue 
	    if(L[tMin] < readStart) {
	      tMin++;
	    } else {
	      whileVar=1; // break loop - done!
	    }
	  } else {
	    whileVar=1; // break loop
	  }
	}
	// now - go right
	whileVar=0;
	while(whileVar==0) {
	  // dont continue if too far
	  if(tMax < (nSNPs - 1)) {
	    if(L[tMax + 1] <= readEnd) { 
	      // number is only to deal with weird cigars	  
	      tMax++;
	    } else {
	      whileVar=1; // break loop - done!
	    }
	  } else {
	    tMax = nSNPs - 1; // done, cannot be larger
	    whileVar=1; // break loop
	  }
	}
	//
	// no matter what, save some minimal info on read
	//
	if (save_sampleReadsInfo) {
	  //
  	    qname_all.push_back(record.Qname());
	    out_readStart_all.push_back(readStart);
	    out_readEnd_all.push_back(readEnd);
	}
	//
	// now, for this read, calculate whether there are SNPs - only bother if tMin>=tMax
	//
	//
	nSNPInRead = -1;
	if(tMin <= tMax) {
  	    // OK, we've got to do stuff
  	    // might have already gotten (and modified) sequence above
  	    // if that is the case, don't re-get
	    if (got_seq_already == false) {
	        seq = record.Sequence();
  	        qual = record.Qualities();
	    }
	    //
	    // do cigar stuff here
	    //
	    iNumOfMs = cigarLengthVec.size();
	    // set some things
	    refOffset = 0; // offset against the reference sequence
	    strandOffset = 0; // offset in the strand
	    //
	    // now, loop over each part of the read (M, D=del, I=ins)
	    //
	    for(iM=0; iM < iNumOfMs; iM++) {
	      cigarLength = cigarLengthVec[iM];
	      cigarType = cigarTypeVec[iM];
	      // if its an M - scan
	    if(cigarType == "M" || cigarType == "=" || cigarType == "X") {
		x1 = refPosition + refOffset; // left part of M
		x2 = refPosition + refOffset + cigarLength-1; // right part of M
		// determine whether that snps is spanned by the read
		for(t=tMin; t<=tMax; t++) {
		    y = L[t];
		    // if this is true - have a SNP!
		    if(x1 <= y && y <= x2) {
		      s = seq[y-refPosition-refOffset+strandOffset];
		      // check if ref or ALT - only keep if true
		      // also only use if BQ at least bqFilter (17) (as in 17 or greater)
		      localbq=int(qual[y-refPosition-refOffset+strandOffset])-33;
		      // also bound BQ above by MQ
		      if(localbq > mapq) // if greater, than reduce
			localbq = mapq;
		      if ((s==ref[t][0] || s==alt[t][0]) && (localbq>=bqFilter)) {
			// is this the reference or alternate?
			nSNPInRead = nSNPInRead+1;
                        if (nSNPInRead >= maxnSNPInRead) {
                          qualLocal.resize(maxnSNPInRead * 2);
                          posLocal.resize(maxnSNPInRead * 2);
                          maxnSNPInRead *= 2;
                        }
			if(s==ref[t][0]) {
			  qualLocal[nSNPInRead] = - localbq;
			}
			if(s==alt[t][0]) {
			  qualLocal[nSNPInRead] = localbq;
			}
			posLocal[nSNPInRead] = t;
		      } // end of check if ref or alt
		    } // end of whether this SNP intersects read
		} // end of loop on SNP
		// now, bump up ref and pos offset by x1
		refOffset = refOffset + cigarLength;
		strandOffset = strandOffset + cigarLength;
	      } // end of if statement on whether cigar type is M
	      // if it is an insertion - bump the strand offset
	      if (cigarType == "I")
		strandOffset = strandOffset + cigarLength;
	      // if it is a deletion - bump the reference position
	      if (cigarType == "D")
		refOffset = refOffset + cigarLength;
	    } // close loop on M
	} // end of check on whether there can be results to run
	//
	// save result if it is worth saving!
	//
      if(nSNPInRead > -1) {
	  qname.push_back(record.Qname());
	  // from seqlib cpp
	  //   /** Get a string (Z) tag 
	  //   * @param tag Name of the tag. eg "XP"
	  //   * @param s The string to be filled in with the tag information
	  //   * @return Returns true if the tag is present, even if empty. Return false if no tag or not a Z tag.
	  //  */
          if (use_bx_tag) {
	      string_z_tag_result = record.GetZTag(bx_string_tag, for_bx_tag_s);
	      if (string_z_tag_result) {
                  bxtag.push_back(for_bx_tag_s);
	      } else {
		  bxtag.push_back("");
	      }
	  }
	  //std::cout << "qname = " << record.Qname() << std::endl;	  
	  //std::cout << "bx_string_tag = " << bx_string_tag << std::endl;
	  //std::cout << "for_bx_tag_s = " << for_bx_tag_s << std::endl;	  
	  //
	  bool ss = record.ReverseFlag();
	  if (ss) {
	    strand.push_back("-");
	  } else {
	    strand.push_back("+");	  
	  }
	  out_readStart.push_back(readStart);
	  out_readEnd.push_back(readEnd);
	  // get start, end of fragment here
	  for(int i=0; i <= nSNPInRead; i++) {
	    out_num_SNPs.push_back(nSNPInRead);
	    out_BQs.push_back(qualLocal[i]); // base qualities	    
	    out_SNP_pos.push_back(posLocal[i]); // positions
	    out_iRead.push_back(qname.size() - 1); // 0-based
	  }
	}
    } // end of loop on read
    return std::make_tuple(out_num_SNPs, out_BQs, out_SNP_pos, out_iRead, qname, bxtag, strand, out_readStart, out_readEnd, qname_all, out_readStart_all, out_readEnd_all);
}
