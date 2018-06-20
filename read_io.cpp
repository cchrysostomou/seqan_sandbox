#include <zlib.h>
#include <seqan/seq_io.h>
// #include <seqan/stream.h>

using namespace std;

int main(int argc, char *argv[])
{

    seqan::CharString seqFileName = argv[1]; //"/home/costa/Documents/projects/seqan/test.fq.gz";
    //seqan::CharString seqFileName = "/home/costa/Documents/files/Undetermined_S0_L001_R1_01.fastq.gz";
    std::cout<<seqFileName<<std::endl;
    seqan::CharString id;
    seqan::Dna5String seq;

    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;
    
    seqan::SeqFileIn seqFileIn; //(seqan::toCString(seqFileName));
    
    if (!seqan::open(seqFileIn,seqan::toCString(seqFileName)))
    {
        std::cerr << "ERROR: could not open input file.\n";
        return 1;
    }

    while (!seqan::atEnd(seqFileIn)){
        seqan::readRecords(ids, seqs, seqFileIn, 100);
        for (int i = 0; i<seqan::length(ids); i++){
            cout<<i<<endl;
            cout << ids[i]<<endl;
            cout << seqs[i]<<endl;
        }
        break;
    }
    //seqan::readRecord(id, seq, seqFileIn);
    //std::cout << id << '\t' << seq << '\n';
    
    return 0;
}
