#ifndef KRAKEN2_ADDITIONAL_HASH_H_
#define KRAKEN2_ADDITIONAL_HASH_H_

#include "kraken2_headers.h"
#include "kraken2_data.h"
#include "taxonomy.h"

namespace kraken2
{

    class AdditionalMap
    {

    public:
        void clearTemp();
        void addk_v2Temp(taxid_t taxon, uint64_t minimizer);
        void saveTemp(vector<taxid_t> &taxa);
        uint16_t GetWeight(uint64_t minimizer);
        size_t GetSize();
        bool IsEmpty();

    private:
        void AddPair(uint64_t minimizer);
        std::unordered_map<uint64_t, uint16_t> ump;
        // 一次read的临时kmer存储
        std::unordered_map<taxid_t, vector<uint64_t>> temp;
        void Add(uint64_t minimizer, uint16_t weight);
    };

}

#endif