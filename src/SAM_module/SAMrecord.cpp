#include<SAM_module/SAMrecord.hpp>
#include<utils/ArgParser.hpp>
#include<iostream>
#include<ranges>
#include<sstream>
#include<iterator>
#include<cmath>

namespace IsoLasso::format
{
    ///START OF HEADRECORD
    bool
    operator>>(std::istream& fin,Header_record& HRecord)
    {
        if(fin.eof())
            return false;
        std::string oneline;
        if (fin.peek()=='@')
        {
            std::getline(fin,oneline);
            return true;
        }
        else
            return false;    
    }
    ///END OF HEADRECORD
    ///START OF RECORD
    std::istream&
    operator>>(std::istream& fin,Sam_record& Record)
    {
        IsoLasso::format::Args arguments;
        Record.reset();

        std::string oneline;
        std::getline(fin,oneline);

        if(oneline[0]=='@')
        {
            Record.ValidBit = false;
            return fin;
        }
        std::istringstream iss(oneline);
        if(iss.fail())
            throw std::invalid_argument("istringstream failure!");

        iss>>Record.QName>>Record.Flag >>Record.RName>>Record.Pos
        >>Record.MapQ >>Record.CIGAR>>Record.RNext>>Record.PNext
        >>Record.TLen >>Record.Seq  >>Record.Qual;
        
        IsoLasso::utils::Setfields(Record);
        IsoLasso::utils::ParseCIGAR(Record);

        //Splice Direction
        Record.SpliceDir = 0;
        if (iss.rdbuf()->in_avail())
        {
            std::string Opt_field;
            while(iss>>Opt_field)
            {
                if (auto pos=Opt_field.find("XS:A:") ; pos!=std::string::npos)
                {
                    if(Opt_field[pos+5]=='+')
                        Record.SpliceDir = 1;
                    else if(Opt_field[pos+5]=='-')
                        Record.SpliceDir = -1;
                }
            }
        }
        return fin;

    }

    std::ostream&
    operator<<(std::ostream& os,const Sam_record& Record)
    {
        os<<Record.QName<<'\t'<<Record.Flag<<'\t'<<Record.RName<<'\t'<<Record.Pos<<'\t'
          <<Record.MapQ <<'\t'<<Record.CIGAR<<'\t'<<Record.RNext<<'\t'<<Record.PNext<<'\t'
          <<Record.TLen <<'\t'<<Record.Seq<<'\t'<<Record.Qual<<'\t'<<Record.SpliceDir;
        return os;
    }
    
    range_type
    Sam_record::GetRange() const
    {
        if(SegmentStart.size())
        {
            if(isPairedEnd)
            {
                if(Pos<PNext) // Left-end
                    return range_type(SegmentStart.front(),SegmentStart.front()+TLen-1);
                else 
                    return range_type(PNext,SegmentEnd.back());
            }
            else 
                return range_type(SegmentStart.front(),SegmentEnd.back()); 
        }
        else
            return range_type(0,0);
    }

}//end of namespace IsoLasso::format

/*
 * Parese CIGAR strings
 * notice that the range should be interpreted as [exonstartpos, exonendpos], not [exonstartpos,exonendpos).
 * Handle I and D Cigar characters
 * exonstart / exonend shall record exonstartpos/exonendpos for all segments, which means
 * each read shall have a exonstart/exonend
 *
 * op    Description
 * M     Alignment match (can be a sequence match or mismatch
 * I     Insertion to the reference
 * D     Deletion from the reference
 * N     Skipped region from the reference
 * S     Soft clip on the read (clipped sequence present in <seq>)
 * H     Hard clip on the read (clipped sequence NOT present in <seq>)
 * P     Padding (silent deletion from the padded reference sequence) s
 *
 * RefPos:     1  2  3  4  5  6  7     8  9 10 11 12 13 14 15 16 17 18 19
 * Reference:  C  C  A  T  A  C  T     G  A  A  C  T  G  A  C  T  A  A  C
 * Read:                   A  C  T  A  G  A  A     T  G  G  C  T
 * CIGAR: 3M1I3M1D5M
 */
namespace IsoLasso::utils
{
  void
  ParseCIGAR(IsoLasso::format::Sam_record& Record)
  {
    std::istringstream CIGARStream(Record.CIGAR);

    long start_pos = Record.Pos;
    int number; 
    char operation;
    bool junction_flag = true;
   
    while(CIGARStream.rdbuf()->in_avail())
    {
        CIGARStream>>number;
        CIGARStream>>operation;

        switch(operation)
        {
            case 'M':
                if(junction_flag)//New segment
                {
                    Record.SegmentStart.push_back(start_pos);
                    Record.SegmentEnd.push_back(start_pos+number-1);
                }
                else 
                    Record.SegmentEnd.back()+=number;
                junction_flag = false;
                start_pos+=number;
                break;                
            case 'N':
                junction_flag=true;
                start_pos+=number;
                break;
            case 'D':
                if(junction_flag)//New segment
                {
                    Record.SegmentStart.push_back(start_pos);
                    Record.SegmentEnd.push_back(start_pos+number-1);
                }
                else 
                    Record.SegmentEnd.back()+=number;
                Record.SegmentEnd.back()+=number;
                start_pos+=number; 
                break;
            case 'I':
            case 'H':
            case 'S':
            case 'P':
                break;
            default:
                std::cerr<<"Invalid CIGAR character "<<operation<<std::endl;
                Record.ValidBit=false;
                break;
        }
    }
    return;
  }  

  void
  Setfields(IsoLasso::format::Sam_record& Record)
  {
    //Valid record
    if((Record.Flag & 0x4))
    {
        std::cout<<"Ignore unmapped record."<<std::endl;
        Record.ValidBit=false;
        return;
    }

    std::string_view QNameSuffix(Record.QName.substr(Record.QName.length()-2));
    if((QNameSuffix=="/1")||(QNameSuffix=="/2"))
    {
        Record.QName       = Record.QName.substr(0,Record.QName.length()-2);
        Record.isPairedEnd = true;
    }
    else
    {
        //Concordant paired-end pairs : (99,147) & (83,163)
        Record.isPairedEnd = (Record.Flag==99)||(Record.Flag==147)||(Record.Flag==83)||(Record.Flag==163)||
                             (Record.Flag==419)||(Record.Flag==339)||(Record.Flag==355)||(Record.Flag==403);

        /*
        * If :
        * 1. the paired-end distance is too large
        * 2. RNEXT is identical to RNAME (=) or rnext is unavailable (*)
        * then force to use single-end version.
        */ 
        if(std::abs(Record.TLen)>MAX_PE_SPAN||(Record.RNext!="*" && Record.RNext!="="))
            Record.isPairedEnd = false;
    }
    Record.ValidBit=true;

    return;
  }
} 

