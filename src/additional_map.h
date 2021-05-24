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
        // 对于冲突来说的
        void ReadConflictFile(const char *filename);
        void WriteConflictMap(const char *filename);
        void saveTemp(unordered_map<taxid_t, vector<uint64_t>> &temp, set<taxid_t> &ancestor);
        double GetWeight(uint64_t minimizer);
        size_t GetConflictSize();
        unordered_map<uint64_t, uint64_t> conflict_ump;
        // 对于未识别的来说的.
        void ReadAdditionalFile(const char *filename);
        void WriteAdditionalMap(const char *filename);
        void AddAdditionalPair(uint64_t minimizer, taxid_t tax_id, Taxonomy &taxonomy);
        taxid_t GetTax(uint64_t minimizer);
        size_t GetAdditionaSize();

    private:
        // 对于冲突来说的
        void AddMinimizer(uint64_t minimizer);

        // 对于未识别的来说的.
        unordered_map<uint64_t, taxid_t> additional_ump;
        void AddTaxon(uint64_t minimizer, taxid_t tax_id);
    };

}

#endif