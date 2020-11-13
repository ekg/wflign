#include <iostream>
#include "seqiter.hpp"
#include "wflign.hpp"

int main(int argc, char** argv) {
    assert(argc == 6);
    std::string targets(argv[1]);
    std::string queries(argv[2]);
    uint64_t segment_length = std::stoi(argv[3]);
    uint64_t min_wavefront_length = std::stoi(argv[4]);
    uint64_t max_distance_threshold = std::stoi(argv[5]);
    seqiter::for_each_seq_in_file(
        queries,
        [&](const std::string& qname,
            const std::string& qseq) {
            seqiter::for_each_seq_in_file(
                targets,
                [&](const std::string& tname,
                    const std::string& tseq) {
                    wflign::wflign_affine_wavefront_reduced(
                        qname, qseq,
                        tname, tseq,
                        segment_length,
                        min_wavefront_length,
                        max_distance_threshold);
                });
        });
    

    return 0;
}
