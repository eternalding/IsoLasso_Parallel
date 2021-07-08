#ifndef COMMONTYPE_HPP
#define COMMONTYPE_HPP
#include <vector>
#include <string>
#include <fstream>

//#define DEBUG

constexpr std::string_view version = "1.0.0";

template<typename T>
using TwoDimVec = std::vector<std::vector<T>>;


typedef std::pair<uint32_t,int32_t> range_type;


/*
 * Parameters used to remove reads that spans too long and spans too many reads.
 * if an exon spans >=MAX_EXON_SPAN and covers >=MAX_JUNCTION_COVERAGE reads, this read is considered as an error read.
 */
constexpr uint32_t    MAX_EXON_SPAN         {50000}  ;
constexpr uint32_t    MAX_JUNCTION_COVERAGE {100000} ;
constexpr uint32_t    MIN_GAP_SPAN          {100}    ; // The minimum length of the gap between two reads to be considered as separate genes.
constexpr uint32_t    MIN_RG_SIZE           {4}      ; // The minimum number of clustered reads to output.
constexpr uint32_t    MAX_PE_SPAN           {700000} ; // The minimum span between two ends of paired-end reads
constexpr uint32_t    MIN_EXON_COV          {1}      ; // The minimum coverage requirement for exon, setting this value small will increase the sensitivity of exon boundary
/*
 * The coverage of a boundary shall be larger than MIN_JUNC_COV 
 * to be considered as a junction.
 */ 
constexpr uint32_t    MIN_JUNC_COV          {1}      ; 
constexpr uint32_t    MIN_OVERLAP           {10}     ;
constexpr double      MIN_CVG_FRAC          {0.25}   ;

constexpr uint32_t    PE_READ_DISTANCE      {200}    ;
constexpr uint32_t    PE_READ_STD           {20}     ;

/* EM */
// GenerateCandidateIsoform
constexpr uint32_t    MIN_JUNC_READ          {1}      ;
constexpr uint32_t    MAX_ISOFORM_NUM        {1024}   ; //TODO : Check the growth of memory i.t.f ISONUM
constexpr double      EXON_MIN_FRAC          {0.003}  ;
constexpr double      JUNC_EXP_FRAC          {0.002}  ;
constexpr double      INTRON_RETENTION_FRAC  {0.1}    ;
constexpr double      EXON_READLEN_FRAC      {2.0}    ;
constexpr double      THERSHOLD_GROWTH       {1.2}    ;
constexpr uint32_t    MAXPATH_PERNODE        {128}    ;



namespace IsoLasso::utils
{
    extern  std::ofstream RG_STATS_FS,GTF_FS          ;
}


#endif