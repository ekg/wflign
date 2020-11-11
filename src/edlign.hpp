#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include "edlib.h"
#include "wfa_edit_callback.hpp"

namespace edlign {

void edlign_full(
    const std::string& query_name,
    const std::string& query,
    const std::string& target_name,
    const std::string& target,
    const uint64_t& segment_length);

void edlign_wavefront(
    const std::string& query_name,
    const std::string& query,
    const std::string& target_name,
    const std::string& target,
    const uint64_t& segment_length);

void do_alignment(
    const std::string& query_name,
    const std::string& query,
    const uint64_t& j,
    const std::string& target,
    const std::string& target_name,
    const uint64_t& i,
    const uint64_t& segment_length,
    const uint64_t& step_size);

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
