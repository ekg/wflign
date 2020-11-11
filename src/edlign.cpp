#include "edlign.hpp"

namespace edlign {

void edlign_full(
        const std::string& query_name,
        const std::string& query,
        const std::string& target_name,
        const std::string& target,
        const uint64_t& segment_length) {
    uint64_t target_step = segment_length / 20;
    uint64_t query_step = segment_length / 20;
    for (uint64_t i = 0; i < target.size() - segment_length - 1; i += target_step) {
        for (uint64_t j = 0; j < query.size() - segment_length - 1; j += query_step) {
            do_alignment(query_name, query, j, target_name, target, i, segment_length, query_step, std::cout);
        }
        // do the last alignment in the row
        do_alignment(query_name, query, query.size()-segment_length, target_name, target, i, segment_length, query_step, std::cout);
    }
    // do the last alignment in the final column
    do_alignment(query_name, query, query.size()-segment_length, target_name, target, target.size()-segment_length, segment_length, query_step, std::cout);
}

void edlign_wavefront(
    const std::string& query_name,
    const std::string& query,
    const std::string& target_name,
    const std::string& target,
    const uint64_t& segment_length) {

    // set up our implicit matrix
    uint64_t steps_per_segment = 2;
    uint64_t step_size = segment_length / steps_per_segment;

    // Pattern & Text
    const int pattern_length = query.size() / step_size - steps_per_segment;
    const int text_length = target.size() / step_size - steps_per_segment;

    //std::cerr << "running WFA on " << pattern_length << " x " << text_length << std::endl;
    // Init Wavefronts
    edit_wavefronts_t wavefronts;
    edit_wavefronts_init(&wavefronts, pattern_length, text_length);
    edit_wavefronts_clean(&wavefronts);
    auto extend_match =
        [&](int& v,
            int& h,
            int& offset,
            const int& distance) {
            bool aligned = do_alignment(
                query_name,
                query,
                v * step_size,
                target_name,
                target,
                h * step_size,
                segment_length,
                step_size,
                std::cout);
            std::cerr << "v/h" << "\t"
                      << v << "\t"
                      << h << "\t"
                      << offset << "\t"
                      << distance << "\t"
                      << aligned << std::endl;
            /*
            if (aligned) {
                v += steps_per_segment / 2;
                h += steps_per_segment / 2;
                offset += steps_per_segment / 2;
            } else {
                ++v;
                ++h;
            }
            v = std::min(v, pattern_length);
            h = std::min(h, text_length);
            offset = std::min(offset, text_length);
            */
            ++v;
            ++h;
            if (aligned) ++offset;
            return aligned;
        };
    edit_wavefronts_align(&wavefronts,
                          extend_match,
                          pattern_length,
                          text_length);
    // todo write the edit cigar, use this to bound the alignment
}

bool do_alignment(
    const std::string& query_name,
    const std::string& query,
    const uint64_t& j,
    const std::string& target_name,
    const std::string& target,
    const uint64_t& i,
    const uint64_t& segment_length,
    const uint64_t& step_size,
    std::ostream& output) {

    auto edlib_config = edlibNewAlignConfig(step_size, // * 1.618034,
                                            EDLIB_MODE_NW,
                                            EDLIB_TASK_PATH,
                                            NULL, 0);

    EdlibAlignResult result = edlibAlign(query.c_str()+j, segment_length,
                                         target.c_str()+i, segment_length,
                                         edlib_config);

    bool aligned = false;
    //Output to file
    if (result.status == EDLIB_STATUS_OK
        && result.alignmentLength != 0
        && result.editDistance >= 0) {

        uint64_t matches = 0;
        uint64_t mismatches = 0;
        uint64_t insertions = 0;
        uint64_t deletions = 0;
        uint64_t softclips = 0;
        uint64_t refAlignedLength = 0;
        uint64_t qAlignedLength = 0;

        char* cigar = alignmentToCigar(result.alignment,
                                       result.alignmentLength,
                                       refAlignedLength,
                                       qAlignedLength,
                                       matches,
                                       mismatches,
                                       insertions,
                                       deletions,
                                       softclips);

        size_t alignmentRefPos = i + result.startLocations[0];
        double total = refAlignedLength + (qAlignedLength - softclips);
        double identity = (double)(total - mismatches * 2 - insertions - deletions) / total;

        output << query_name
               << "\t" << query.size()
               << "\t" << j
               << "\t" << j + qAlignedLength
               << "\t" << "+" // todo (currentRecord.strand == skch::strnd::FWD ? "+" : "-")
               << "\t" << target_name
               << "\t" << target.size()
               << "\t" << alignmentRefPos
               << "\t" << alignmentRefPos + refAlignedLength
               << "\t" << matches
               << "\t" << std::max(refAlignedLength, qAlignedLength)
               << "\t" << std::round(float2phred(1.0-identity))
               << "\t" << "id:f:" << identity
               << "\t" << "ma:i:" << matches
               << "\t" << "mm:i:" << mismatches
               << "\t" << "ni:i:" << insertions
               << "\t" << "nd:i:" << deletions
               << "\t" << "ns:i:" << softclips
               << "\t" << "ed:i:" << result.editDistance
               << "\t" << "al:i:" << result.alignmentLength
               << "\t" << "se:f:" << result.editDistance / (double)result.alignmentLength
               << "\t" << "cg:Z:" << cigar
               << "\n";

        free(cigar);
        aligned = true;
    }

    edlibFreeAlignResult(result);

    return aligned;
}

char* alignmentToCigar(
    const unsigned char* const alignment,
    const int alignmentLength,
    uint64_t& refAlignedLength,
    uint64_t& qAlignedLength,
    uint64_t& matches,
    uint64_t& mismatches,
    uint64_t& insertions,
    uint64_t& deletions,
    uint64_t& softclips) {

    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char moveCodeToChar[] = {'=', 'I', 'D', 'X'};

    std::vector<char>* cigar = new std::vector<char>();
    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    for (int i = 0; i <= alignmentLength; i++) {
        // if new sequence of same moves started
        if (i == alignmentLength || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0)) {
            // calculate matches, mismatches, insertions, deletions
            switch (lastMove) {
            case 'I':
                // assume that starting and ending insertions are softclips
                if (i == alignmentLength || cigar->empty()) {
                    softclips += numOfSameMoves;
                } else {
                    insertions += numOfSameMoves;
                }
                qAlignedLength += numOfSameMoves;
                break;
            case '=':
                matches += numOfSameMoves;
                qAlignedLength += numOfSameMoves;
                refAlignedLength += numOfSameMoves;
                break;
            case 'X':
                mismatches += numOfSameMoves;
                qAlignedLength += numOfSameMoves;
                refAlignedLength += numOfSameMoves;
                break;
            case 'D':
                deletions += numOfSameMoves;
                refAlignedLength += numOfSameMoves;
                break;
            default:
                break;
            }
                  
            // Write number of moves to cigar string.
            int numDigits = 0;
            for (; numOfSameMoves; numOfSameMoves /= 10) {
                cigar->push_back('0' + numOfSameMoves % 10);
                numDigits++;
            }
            std::reverse(cigar->end() - numDigits, cigar->end());
            // Write code of move to cigar string.
            cigar->push_back(lastMove);
            // If not at the end, start new sequence of moves.
            if (i < alignmentLength) {
                // Check if alignment has valid values.
                if (alignment[i] > 3) {
                    delete cigar;
                    return 0;
                }
                numOfSameMoves = 0;
            }
        }
        if (i < alignmentLength) {
            lastMove = moveCodeToChar[alignment[i]];
            numOfSameMoves++;
        }
    }
    cigar->push_back(0);  // Null character termination.
    char* cigar_ = (char*) malloc(cigar->size() * sizeof(char));
    std::memcpy(cigar_, &(*cigar)[0], cigar->size() * sizeof(char));
    delete cigar;

    return cigar_;
}

double float2phred(const double& prob) {
    if (prob == 1)
        return 255;  // guards against "-0"
    double p = -10 * (double) log10(prob);
    if (p < 0 || p > 255) // int overflow guard
        return 255;
    else
        return p;
}

}
