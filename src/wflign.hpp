#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include "edlib.h"
#include "wfa_edit_callback.hpp"
#include "gap_affine/affine_wavefront_align.hpp"

namespace wflign {

struct alignment_t {
    int v = 0;
    int h = 0;
    double identity = 0;
    const std::string* query_name;
    uint64_t query_size;
    const std::string* target_name;
    uint64_t target_size;
    EdlibAlignResult result;
    ~alignment_t(void) {
        edlibFreeAlignResult(result);
    }
};

void wflign_full(
    const std::string& query_name,
    const std::string& query,
    const std::string& target_name,
    const std::string& target,
    const uint64_t& segment_length);

void wflign_wavefront(
    const std::string& query_name,
    const std::string& query,
    const std::string& target_name,
    const std::string& target,
    const uint64_t& segment_length);

void wflign_affine_wavefront(
    const std::string& query_name,
    const std::string& query,
    const std::string& target_name,
    const std::string& target,
    const uint64_t& segment_length);

void wflign_affine_wavefront_reduced(
    const std::string& query_name,
    const std::string& query,
    const std::string& target_name,
    const std::string& target,
    const uint64_t& segment_length,
    const int& min_wavefront_length,
    const int& max_distance_threshold);

bool do_alignment(
    const std::string& query_name,
    const std::string& query,
    const uint64_t& j,
    const std::string& target_name,
    const std::string& target,
    const uint64_t& i,
    const uint64_t& segment_length,
    const uint64_t& step_size,
    alignment_t& alignment);

std::ostream& operator<<(
    std::ostream& out,
    const alignment_t& aln);

char* alignmentToCigar(
    const unsigned char* const alignment,
    const int alignmentLength,
    uint64_t& refAlignedLength,
    uint64_t& qAlignedLength,
    uint64_t& matches,
    uint64_t& mismatches,
    uint64_t& insertions,
    uint64_t& deletions,
    uint64_t& softclips);

double float2phred(const double& prob);

}
