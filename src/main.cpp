#include <iostream>
#include "seqiter.hpp"
#include "edlign.hpp"

int main(int argc, char** argv) {
    assert(argc == 4);
    std::string queries(argv[1]);
    std::string targets(argv[2]);
    uint64_t segment_length = std::stoi(argv[3]);
    //std::cout << "query.name" << "\t" << "target.name" << "\t" << "x" << "\t" << "y" << "\t" << "edit.dist" << "\t" << "identity" << "\n";
    seqiter::for_each_seq_in_file(
        queries,
        [&](const std::string& qname,
            const std::string& qseq) {
            seqiter::for_each_seq_in_file(
                targets,
                [&](const std::string& tname,
                    const std::string& tseq) {
                    edlign::edlign_full(
                        qname, qseq,
                        tname, tseq,
                        segment_length);                    
                });
        });
    

    return 0;
}
