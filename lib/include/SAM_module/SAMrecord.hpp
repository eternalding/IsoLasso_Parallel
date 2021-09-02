#ifndef SAM_RECORD_HPP
#define SAM_RECORD_HPP

#include<iostream>
#include"utils/Commontype.hpp"
#include<string>
#include<vector>
#include<variant>
#include"utils/ArgParser.hpp"
#include<iostream>
#include<ranges>
#include<sstream>
#include<iterator>
#include<numeric>

namespace IsoLasso::format
{
    struct MetaData
    {
        std::string                        tag;
        std::variant<uint32_t,std::string> value;
    };

    class Sam_record
    {
        public:
        // Required fields
        std::string         QName;  //Query template NAME
        uint32_t            Flag;   //bitwise FLAG
        std::string         RName;  //Reference sequence NAME
        uint32_t            Pos;    //1-based leftmost mapping POSition
        uint32_t            MapQ;   //MAPping Quality
        std::string         CIGAR;  //CIGAR string
        std::string         RNext;  //Reference name of the mate/next read
        uint32_t            PNext;  //Position of the mate/next read
        int32_t             TLen;   // observed Template LENgth
        std::string         Seq;    //segment SEQuence
        std::string         Qual;   //ASCII of Phred-scaled base QUALity+33
        
        // Optional fields
        int32_t             SpliceDir   {0};

        // Indicator fields 
        bool                ValidBit    {false};
        bool                isPairedEnd {false};
        // Containers
        std::vector<uint32_t>   SegmentStart,SegmentEnd;


        /*
         * I/O functions
         */
        friend bool
        operator>>(std::istream& fin,IsoLasso::format::Sam_record& Record);

        friend std::ostream&
        operator<<(std::ostream& os,const IsoLasso::format::Sam_record& Record);

        /*Utility functions*/
          
        inline void
        reset() 
        {
            SegmentStart.resize(0);
            SegmentStart.shrink_to_fit();
            SegmentEnd.resize(0);
            SegmentEnd.shrink_to_fit();
            ValidBit=false;
            isPairedEnd=false;
            return;
        }
        range_type
        GetRange() const;
    };

    class Header_record
    {
        std::string                             Header_tag; //@HD,@SQ,@RG,@PG,@CO
        std::vector<IsoLasso::format::MetaData> HRecords;

        friend bool
        operator>>(std::istream& fin,IsoLasso::format::Header_record& HRecord);

        friend std::ostream&
        operator<<(std::ostream& os,const IsoLasso::format::Header_record& HRecord);

    };
};//end of namespace IsoLasso::format

namespace IsoLasso::utils
{
    void
    ParseCIGAR(IsoLasso::format::Sam_record&);

    inline uint32_t
    GetEfficientLen(const IsoLasso::format::Sam_record& Record)
    {
        uint32_t EfficientLen {0};
        for(auto i =0;i<Record.SegmentStart.size();i++) 
            EfficientLen+=(Record.SegmentEnd[i]-Record.SegmentStart[i]+1);
        return EfficientLen;
    }

    void
    Setfields(format::Sam_record&);

    inline uint32_t 
    FileRead( std::istream & is, std::vector<char> & buff )
    {
        is.read( &buff[0], buff.size() );
        return is.gcount();
    }

    inline uint64_t 
    CountLines( const std::vector <char> & buff, uint32_t sz )
    {
        uint64_t newlines = 0;
        const char * p = &buff[0];
        for ( uint32_t i = 0; i < sz; i++ ) 
        {
            if ( p[i] == '\n' ) {
                newlines++;
            }
        }
        return newlines;
    }

    inline const uint64_t 
    TotalReadCnt(const std::string FileName)
    {
        constexpr uint32_t SZ = 1024 * 1024;
        std::vector <char> buff( SZ );
        std::ifstream ifs(FileName);
        uint64_t LineCnt {0};
        while(auto cc = FileRead( ifs, buff ) ) 
            LineCnt += CountLines( buff, cc );
        return LineCnt;
    }



}//end of namespace IsoLasso::utils


#endif