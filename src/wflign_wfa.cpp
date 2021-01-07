#include "wflign_wfa.hpp"

namespace wflign {

namespace wavefront {

void wflign_affine_wavefront(
    std::ostream& out,
    const std::string& query_name,
    const char* query,
    const uint64_t& query_total_length,
    const uint64_t& query_offset, // todo this is broken for reverse alignment reporting!!!!
    const uint64_t& query_length,
    const bool& query_is_rev,
    const std::string& target_name,
    const char* target,
    const uint64_t& target_total_length,
    const uint64_t& target_offset,
    const uint64_t& target_length,
    const uint64_t& segment_length,
    const float& min_identity,
    const int& wflambda_min_wavefront_length, // with these set at 0 we do exact WFA for wflambda
    const int& wflambda_max_distance_threshold) {
    //const int& wfa_min_wavefront_length, // with these set at 0 we do exact WFA for WFA itself
    //const int& wfa_max_distance_threshold) {

    // set up our implicit matrix
    uint64_t steps_per_segment = 2;
    uint64_t step_size = segment_length / steps_per_segment;

    // Pattern & Text
    const int pattern_length = query_length / step_size;
    const int text_length = target_length / step_size;

    const int wfa_min_wavefront_length = segment_length / 10;
    const int wfa_max_distance_threshold = segment_length / 2;

    // Allocate MM
    wflambda::mm_allocator_t* const wflambda_mm_allocator = wflambda::mm_allocator_new(BUFFER_SIZE_8M);
    // Set penalties
    wflambda::affine_penalties_t wflambda_affine_penalties = {
        .match = 0,
        .mismatch = 4,
        .gap_opening = 6,
        .gap_extension = 2,
    };
    // Init Affine wflambda
    wflambda::affine_wavefronts_t* affine_wavefronts;
    if (wflambda_min_wavefront_length || wflambda_max_distance_threshold) {
        affine_wavefronts = wflambda::affine_wavefronts_new_reduced(
            pattern_length, text_length, &wflambda_affine_penalties,
            wflambda_min_wavefront_length, wflambda_max_distance_threshold,
            NULL, wflambda_mm_allocator);
    } else {
        affine_wavefronts = wflambda::affine_wavefronts_new_complete(
            pattern_length, text_length, &wflambda_affine_penalties, NULL, wflambda_mm_allocator);
    }

    whash::patchmap<uint64_t,alignment_t*> alignments;
    // save this in a pair-indexed patchmap

    // setup affine WFA
    wfa::mm_allocator_t* const wfa_mm_allocator = wfa::mm_allocator_new(BUFFER_SIZE_8M);
    wfa::affine_penalties_t wfa_affine_penalties = {
        .match = 0,
        .mismatch = 4,
        .gap_opening = 6,
        .gap_extension = 2,
    };

    auto extend_match =
        [&](const int& v,
            const int& h) {
            bool aligned = false;
            uint64_t k = encode_pair(v, h);
            if (v >= 0 && h >= 0
                && v < pattern_length
                && h < text_length) {
                auto f = alignments.find(k);
                if (f != alignments.end()) {
                    aligned = true;
                } else  {
                    uint64_t query_begin = (v < pattern_length-1 ? v * step_size
                                            : query_length - segment_length);
                    uint64_t target_begin = (h < text_length-1 ? h * step_size
                                             : target_length - segment_length);
                    alignment_t* aln = new alignment_t(
                        segment_length,
                        segment_length,
                        wfa_min_wavefront_length,
                        wfa_max_distance_threshold,
                        wfa_mm_allocator,
                        wfa_affine_penalties);
                    aligned =
                        do_alignment(
                            query_name,
                            query,
                            query_length,
                            query_begin,
                            target_name,
                            target,
                            target_length,
                            target_begin,
                            segment_length,
                            step_size,
                            *aln);
                    if (aligned) {
                        alignments[k] = aln;
                    } else {
                        delete aln;
                    }
                }
            }
            return aligned;
        };

    // accumulate runs of matches in reverse order
    // then trim the cigars of successive mappings
    //
    std::vector<alignment_t*> trace;
    
    auto trace_match =
        [&](const int& v, const int& h) {
            uint64_t k = encode_pair(v, h);
            auto f = alignments.find(k);
            if (f != alignments.end()) {
                trace.push_back(f->second);
                return true;
            } else {
                return false;
            }
        };

    // Align
    wflambda::affine_wavefronts_align(
        affine_wavefronts,
        extend_match,
        trace_match,
        pattern_length,
        text_length);

    // get alignment score
    const int score = wflambda::edit_cigar_score_gap_affine(
        &affine_wavefronts->edit_cigar, &wflambda_affine_penalties);
#ifdef WFLIGN_DEBUG
    std::cerr << "[wflign::wflign_affine_wavefront] alignment score " << score << " for query: " << query_name << " target: " << target_name << std::endl;
#endif

    // todo: implement alignment identifier based on hash of the input, params, and commit
    // annotate each PAF record with it and the full alignment score

    // Trim alignments that overlap in the query
    if (trace.size()) {
        int last_i = 0;
        int last_j = 0;
        for (auto x = trace.begin(); x != trace.end(); ++x) {
            auto& curr = **x;
            curr.keep_query_length = segment_length;
        }
        for (auto x = trace.rbegin()+1; x != trace.rend(); ++x) {
            auto& curr = **x;
            auto& last = **(x-1);
            int ovlp_j = last.j + segment_length - curr.j;
            if (ovlp_j > 0) {
                int trim_last = ovlp_j / 2;
                int trim_curr = ovlp_j - trim_last;
                last.keep_query_length -= trim_last;
                curr.keep_query_length -= trim_curr;
                curr.skip_query_start += trim_curr;
            }
        }
        for (auto x = trace.rbegin(); x != trace.rend(); ++x) {
            write_alignment(out, **x,
                            query_name, query_total_length, query_offset, query_length,
                            query_is_rev,
                            target_name, target_total_length, target_offset, target_length,
                            min_identity);
        }
    }

    // Free
    wflambda::affine_wavefronts_delete(affine_wavefronts);
    wflambda::mm_allocator_delete(wflambda_mm_allocator);
    for (const auto& p : alignments) {
        delete p.second;
    }
}

// accumulate alignment objects
// run the traceback determine which are part of the main chain
// order them and write them out
// needed--- 0-cost deduplication of alignment regions (how????)
//     --- trim the alignment back to the first 1/2 of the query

bool do_alignment(
    const std::string& query_name,
    const char* query,
    const uint64_t& query_length,
    const uint64_t& j,
    const std::string& target_name,
    const char* target,
    const uint64_t& target_length,
    const uint64_t& i,
    const uint64_t& segment_length,
    const uint64_t& step_size,
    alignment_t& aln) {

    aln.score = wfa::affine_wavefronts_align(
        aln.affine_wavefronts, query+j, segment_length, target+i, segment_length);

    aln.j = j;
    aln.i = i;

    return aln.score < segment_length;
}


void write_alignment(
    std::ostream& out,
    const alignment_t& aln,
    const std::string& query_name,
    const uint64_t& query_total_length,
    const uint64_t& query_offset, // query offset on the forward strand
    const uint64_t& query_length, // used to compute the coordinates for reversed alignments
    const bool& query_is_rev,
    const std::string& target_name,
    const uint64_t& target_total_length,
    const uint64_t& target_offset,
    const uint64_t& target_length, // unused
    const float& min_identity,
    const bool& with_endline) {
//    bool aligned = false;
    //Output to file
    //auto& result = aln.result;
    //if (aln.score < 10) {
    {

        /*
        if (result.status == EDLIB_STATUS_OK
        && result.alignmentLength != 0
        && result.editDistance >= 0) {
        */

        uint64_t matches = 0;
        uint64_t mismatches = 0;
        uint64_t insertions = 0;
        uint64_t deletions = 0;
        uint64_t softclips = 0;
        uint64_t refAlignedLength = 0;
        uint64_t qAlignedLength = 0;
        int skipped_target_start = 0;
        int kept_target_length = 0;

        char* cigar = alignmentToCigar(&aln.affine_wavefronts->edit_cigar,
                                       aln.skip_query_start,
                                       aln.keep_query_length,
                                       skipped_target_start,
                                       kept_target_length,
                                       refAlignedLength,
                                       qAlignedLength,
                                       matches,
                                       mismatches,
                                       insertions,
                                       deletions,
                                       softclips);

        size_t alignmentRefPos = aln.i;
        double total = refAlignedLength + (qAlignedLength - softclips);
        double identity = (double)(total - mismatches * 2 - insertions - deletions) / total;
        // convert our coordinates to be relative to forward strand (matching PAF standard)
        uint64_t q_start;
        if (query_is_rev) {
            q_start = query_offset + (query_length - (aln.j + aln.skip_query_start + qAlignedLength));
        } else {
            q_start = query_offset + aln.j + aln.skip_query_start;
        }
        if (identity >= min_identity) {
            out << query_name
                << "\t" << query_total_length
                << "\t" << q_start
                << "\t" << q_start + qAlignedLength
                << "\t" << (query_is_rev ? "-" : "+")
                << "\t" << target_name
                << "\t" << target_total_length
                << "\t" << target_offset + alignmentRefPos + skipped_target_start
                << "\t" << target_offset + alignmentRefPos + skipped_target_start + refAlignedLength
                << "\t" << matches
                << "\t" << std::max(refAlignedLength, qAlignedLength)
                << "\t" << std::round(float2phred(1.0-identity))
                << "\t" << "as:i:" << aln.score
                << "\t" << "id:f:" << identity
                << "\t" << "ma:i:" << matches
                << "\t" << "mm:i:" << mismatches
                << "\t" << "ni:i:" << insertions
                << "\t" << "nd:i:" << deletions
                << "\t" << "ns:i:" << softclips
                << "\t" << "cg:Z:" << cigar;
            if (with_endline) {
                out << std::endl;
            }
        }
        free(cigar);
    }
}

char* alignmentToCigar(
    const wfa::edit_cigar_t* const edit_cigar,
    const int skip_query_start,
    const int keep_query_length,
    int& skipped_target_start,
    int& kept_target_length,
    uint64_t& refAlignedLength,
    uint64_t& qAlignedLength,
    uint64_t& matches,
    uint64_t& mismatches,
    uint64_t& insertions,
    uint64_t& deletions,
    uint64_t& softclips) {

    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    //char moveCodeToChar[] = {'=', 'I', 'D', 'X'};

    std::vector<char>* cigar = new std::vector<char>();
    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    //bool in_keep = false;
    int seen_query = 0;
    int start_idx = edit_cigar->begin_offset;

    //for (int i = edit_cigar->begin_offset; i <= edit_cigar->end_offset; ++i) {
    //    if (i == edit_cigar->end_offset || (edit_cigar->operations[i] != lastMove && lastMove != 0)) {
    
    while (start_idx < edit_cigar->end_offset
           && seen_query < skip_query_start) {
        switch (edit_cigar->operations[start_idx++]) {
        case 'M':
        case 'X':
        case '=':
            ++skipped_target_start;
            ++seen_query;
            break;
        case 'I':
            ++seen_query;
            break;
        case 'D':
            ++skipped_target_start;
            break;
        default:
            break;
        }
    }
    int end_idx = start_idx;
    seen_query = 0;
    while (end_idx < edit_cigar->end_offset
        && seen_query < keep_query_length) {
        switch (edit_cigar->operations[end_idx++]) {
        case 'M':
        case 'X':
        case '=':
            ++kept_target_length;
            ++seen_query;
            break;
        case 'I':
            ++seen_query;
            break;
        case 'D':
            ++kept_target_length;
            break;
        default:
            break;
        }
    }
    if (end_idx == start_idx) {
        end_idx = edit_cigar->end_offset;
    }
    for (int i = start_idx; i <= end_idx; i++) {
        // if new sequence of same moves started
        if (i == end_idx || (edit_cigar->operations[i] != lastMove && lastMove != 0)) {
            // calculate matches, mismatches, insertions, deletions
            switch (lastMove) {
            case 'I':
                // assume that starting and ending insertions are softclips
                if (i == end_idx || cigar->empty()) {
                    softclips += numOfSameMoves;
                } else {
                    insertions += numOfSameMoves;
                }
                qAlignedLength += numOfSameMoves;
                break;
            case 'M':
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
            if (i < end_idx) {
                numOfSameMoves = 0;
            }
        }
        if (i < end_idx) {
            lastMove = edit_cigar->operations[i];
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

}
