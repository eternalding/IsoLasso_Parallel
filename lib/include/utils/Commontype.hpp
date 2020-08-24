#ifndef COMMONTYPE_HPP
#define COMMONTYPE_HPP
#include <vector>
#include <string>

constexpr std::string_view version = "1.0.0";

template<typename T>
using TwoDimVec = std::vector<std::vector<T>>;


typedef std::pair<u_int32_t,int32_t> range_type;


/*
 * Parameters used to remove reads that spans too long and spans too many reads.
 * if an exon spans >=MAX_EXON_SPAN and covers >=MAX_JUNCTION_COVERAGE reads, this read is considered as an error read.
 */
constexpr uint32_t    MAX_EXON_SPAN         {50000}  ;
constexpr uint32_t    MAX_JUNCTION_COVERAGE {100000} ;

constexpr uint32_t    MIN_GAP_SPAN          {1000}   ; // The minimum length of the gap between two reads to be considered as separate genes.

constexpr uint32_t    MIN_RG_SIZE           {4}      ; // The minimum number of clustered reads to output.
constexpr uint32_t    MAX_PE_SPAN           {700000} ; // The minimum span between two ends of paired-end reads


#endif