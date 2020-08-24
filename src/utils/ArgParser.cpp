#include<utils/ArgParser.hpp>

namespace IsoLasso::utils
{
    void 
    printHelpMsg()
    {
        std::cerr<<"IsoLassoParallel "<<version<<std::endl;
        std::cerr<<"Usage: processsam {options} <in.sam|->\n\n";
        std::cerr<<"Input:\n";
        std::cerr<<"A SAM format file containing the read mapping information, or '-' to read from STDIN. See NOTE for further information.\n";
        std::cerr<<"\nOptions:\n\n";
        std::cerr<<"PARAMETERS:\n\n";
        std::cerr<<"  -i/\tThe path of Input SAM file.\n\n";
        std::cerr<<"  -g/--min-gap-length <int>\tThe minimum length of the gap between two reads to be considered as separate genes. Default 100.\n\n";
        std::cerr<<"  -c/--min-read-num <int>\tThe minimum number of clustered reads to output. Default 4.\n\n";
        std::cerr<<"  -k/--max-pe-span <int> \tThe maximum pair-end spanning. Paired-end reads whose spanning exceeds this number will be discarded. Default 700000.\n\n";
        std::cerr<<"  -s/--max-num-instance \tThe maximum number of instances be written to the file. Default -1 (no limit)\n\n";
        std::cerr<<"  -u/--min-cvg-cut <0.0-1.0>\tThe fraction for coverage cutoff, should be between 0-1. A higher value will be more sensitive to coverage discrepancies in one gene. Default 0.05.\n\n";
        std::cerr<<"  -b/--single-only \t\tTreat reads as single-end reads, even if they are paired-end reads.\n\n";
        std::cerr<<"  -j/--min-junc-count <int>\tMinimum junction count. Only junctions  with no less than this number of supporting reads are considered. Default 1.\n\n";
        std::cerr<<"  -d/--direction <+|->\tIf this is the stranded RNA-Seq data, specify the direction of the data. Default '.' (non-stranded).\n\n";
        std::cerr<<"IO OPTIONS:\n\n";
        std::cerr<<"  -n/--isoinfer\t\t\tGenerate IsoInfer input files (.readinfo, .bound and .generange). Default off.\n\n";
        std::cerr<<"  -x/--annotation <string>\tProvide existing gene annotation file (in BED format). Adding this parameter will automatically incorporate existing gene annotation information into instance file. The bed file should be sorted according to the chromosome name and starting position of isoforms. This option is mutually exclusive to the -r/--range option.\n\n";
        std::cerr<<"  -r/--range <string>\t\tUse the provided gene ranges specified by the file (in BED format). This option is mutually exclusive to the -x/--annotation option.\n\n";
        std::cerr<<"  -e/--segment-bound <string>\tProvide the exon-intron boundary information specified by the filename. See NOTE for more information about the file format.\n\n";
        std::cerr<<"  -a/--annotation\t\tOutput annoation files, including read coverage (.real.wig), read coverage considering junctions and paired-end read spans (.wig), instance range and boundary (.bound.bed), junctions (.bed) and  junction summary (.junction.bed).\n\n";
        std::cerr<<"  -o/--prefix <string>\t\tSpecify the prefix of all generated files. The default value is the provided file name. \n\n";
        std::cerr<<"  -v/--no-coverage\t\tDon't output coverage information to the instance file.\n\n";
        std::cerr<<"NOTE\n\n\t1. processsam acceptes STDIN input of sam file by using '-' as filename. This is especially useful if you have the .bam file (e.g., from Tophat output), or you want to do some read filtering before running IsoLasso.  For example, if Samtools is installed, then use the following command to run processsam on only chromosome 1 reads:\n";
        std::cerr<<"\tsamtools view accepted_hits.bam chr1 | processsam -a -o accepted_hits -\n\n";
        std::cerr<<"\t2. The sam/bam file must be sorted according to the chromosome name and starting position. The bam file format can be sorted using 'samtools sort' command, while for the sam file, you can use the sort command. In Unix or Mac systems, use the following command: \n";
        std::cerr<<"\tsort -k 3,3 -k 4,4n in.sam > in.sorted.sam\n";
        std::cerr<<"to sort in.sam  into in.sorted.sam, or use the pipe:\n";
        std::cerr<<"\tsort -k 3,3 -k 4,4n in.sam | processsam -a -o accepted_hits -\n";
        std::cerr<<"\t3. The exon-intron boundary file (specified by -e/--segment-bound option) records the exon-intron boundary used by IsoLasso. Each line in the file represents one boundary information, and should include chromosome name, start position, end position (equal to start position) and direction (+/-). These fields should be tab-separated, and only the first 4 fields are used. For example,\n";
        std::cerr<<"\tchr1\t15796\t15796\t+\n";
        return;
    }

    bool
    ParseArg(const int& argc ,char* argv[],IsoLasso::format::Args& Args)
    {
        std::string arg_type; 
        for(auto i=1;i<argc;i++)
        {
            if (i%2==1)
                arg_type=argv[i];
            else
            {
                if (arg_type=="-i") 
                    Args.SAMFILE=argv[i];
                else if (arg_type=="-o") 
                    Args.OUTPUT_PREFIX=argv[i];
                else
                    std::__throw_invalid_argument("Invalid argument!");
                                                  
            }
        }
        return true;
    }

}//end of namespace IsoLasso::utils