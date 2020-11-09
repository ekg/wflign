#include "edlign.hpp"

namespace edlign {

void edlign(
        const std::string& query_name,
        const std::string& query,
        const std::string& target_name,
        const std::string& target,
        const uint64_t& segment_length) {
    for (uint64_t i = 0; i < target.size(); i += segment_length) {
        for (uint64_t j = 0; j < query.size(); j += segment_length) {
            uint64_t query_size = std::min(segment_length, query.size()-j);
            uint64_t target_size = std::min(segment_length, target.size()-i);
            EdlibAlignResult result = edlibAlign(query.c_str()+j, query_size,
                                                 target.c_str()+i, target_size,
                                                 edlibDefaultAlignConfig());
            if (result.status == EDLIB_STATUS_OK) {
                uint64_t max_size = std::max(query_size, target_size);
                std::cout
                    << query_name << "\t" << target_name << "\t"
                    << i << "\t" << j << "\t" << result.editDistance << "\t"
                    << (double)(max_size - result.editDistance) / (double) max_size << "\n";
            }
            edlibFreeAlignResult(result);
        }
    }
}

}
