#ifndef OPTIMOTU_PIPELINE_FASTQ_H
#define OPTIMOTU_PIPELINE_FASTQ_H

#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <deque>

// disable warnings for use of auto_ptr and unnecessary parentheses in boost
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#pragma GCC diagnostic ignored "-Wparentheses"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string/predicate.hpp>

#pragma GCC diagnostic pop

typedef boost::iostreams::filtering_stream<boost::iostreams::input>
  filter_stream_in;
typedef boost::iostreams::filtering_stream<boost::iostreams::output>
  filter_stream_out;
typedef std::string string;

#define NO_LIMIT (std::numeric_limits<std::streamsize>::max())

void open_fastx_in(filter_stream_in & instream, const std::string & infile);

void open_fastx_out(filter_stream_out & outstream, const std::string & outfile, bool compress);

// [[Rcpp::depends(BH)]]

#endif
