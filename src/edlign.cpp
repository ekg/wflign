#include "edlign.hpp"

namespace edlign {

void edlign(
        const std::string& query_name,
        const std::string& query,
        const std::string& target_name,
        const std::string& target,
        const uint64_t& segment_length) {
    uint64_t step_size = segment_length / 4;
    auto edlib_config = edlibNewAlignConfig(segment_length / 2, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
    for (uint64_t i = 0; i < target.size(); i += step_size) {
        for (uint64_t j = 0; j < query.size(); j += step_size) {
            uint64_t query_size = std::min(segment_length, query.size()-j);
            uint64_t target_size = std::min(segment_length, target.size()-i);\
            EdlibAlignResult result = edlibAlign(query.c_str()+j, query_size,
                                                 target.c_str()+i, target_size,
                                                 edlib_config);
            if (result.status == EDLIB_STATUS_OK && result.editDistance >= 0) {
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
