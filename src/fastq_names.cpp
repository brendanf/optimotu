#include "fastq.h"

//' Get the names of reads in a FASTQ file
//' @param x (`character`) FASTQ file path, optionally gzipped.
//' @return (`character`) Read names.
// [[Rcpp::export]]
Rcpp::CharacterVector fastq_names(std::string x) {
  std::deque<string> names;
  filter_stream_in instream;
  open_fastx_in(instream, x);
  string line;
  while (instream >> line) {
    if (line.size() == 0) break;
    line = line.substr(1);
    names.push_back(line);
    instream.ignore(NO_LIMIT, '\n');
    instream.ignore(NO_LIMIT, '\n');
    instream.ignore(NO_LIMIT, '\n');
    instream.ignore(NO_LIMIT, '\n');
  }
  return Rcpp::wrap(names);
}
