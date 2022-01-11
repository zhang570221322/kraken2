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
        // Test
        unordered_map<string, uint64_t> MyCounter1;
        unordered_map<string, uint64_t> MyCounter2;
        //金标准数据
        unordered_map<string, uint64_t> _seqID_taxID;
        void loadGod_data(const char *filename);
        // 对于冲突来说的
        unordered_map<uint64_t, uint64_t> conflict_ump;
        vector<taxid_t> conflict_ancestor;
        unordered_map<taxid_t, vector<uint64_t>> conflict_temp;
        unordered_map<taxid_t, bool> leaf;
        unordered_map<taxid_t, int> leaf_count;
        unordered_map<double, vector<taxid_t>> score_taxid;
        void clearConflictObj();
        void ReadConflictFile(const char *filename);
        void WriteConflictMap(const char *filename);
        void saveTemp(taxon_counts_t &hit_counts, double max_score, uint64_t &conflicts, Taxonomy &taxonomy);
        double GetWeight(uint64_t minimizer, uint64_t child_count);
        size_t GetConflictSize();
        // 对于未识别的来说的.
        void ReadAdditionalFile(const char *filename);
        void WriteAdditionalMap(const char *filename);
        void AddAdditionalPair(uint64_t minimizer, taxid_t tax_id, Taxonomy &taxonomy);
        taxid_t GetTax(uint64_t minimizer);
        size_t GetAdditionaSize();

    private:
        // 对于冲突来说的
        void AddMinimizer(uint64_t minimizer);
        // 计算需要降低权重的tax. 存入conflict_ancestor
        void findConflictLTR(taxon_counts_t &hit_counts, double max_score, uint64_t &conflicts, Taxonomy &taxonomy);
        // 对于未识别的来说的.
        unordered_map<uint64_t, taxid_t> additional_ump;
        void AddTaxon(uint64_t minimizer, taxid_t tax_id);
    };

}

#endif