#include "additional_map.h"

using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::unordered_map;
using std::vector;

namespace kraken2
{

    void AdditionalMap::AddPair(uint64_t minimizer)
    {
        uint16_t defalut_weight = 2;
        //如果已经存在次数,再继续提升次数
        if (ump.find(minimizer) != ump.end())
        {
            ump[minimizer] = 1 + ump[minimizer];
        }
        else
        {
            Add(minimizer, defalut_weight);
        }
    }

    void AdditionalMap::Add(uint64_t minimizer, uint16_t weight)
    {
        ump.emplace(minimizer, weight);
    }
    void AdditionalMap::clearTemp()
    {

        temp.clear();
    }

    void AdditionalMap::addk_v2Temp(taxid_t taxon, uint64_t minimizer)
    {
        temp[taxon].push_back(minimizer);
    }
    void AdditionalMap::saveTemp(vector<taxid_t> &taxa)
    {
        for (size_t i = 0; i < taxa.size(); i++)
        {
            auto code = taxa[i];
            if (temp.find(code) != temp.end())
            {
                vector<uint64_t> kmers = temp[code];
                for (size_t i = 0; i < kmers.size(); i++)
                {
                    AddPair(kmers[i]);
                }
            }
        }
        clearTemp();
    }

    uint16_t AdditionalMap::GetWeight(uint64_t minimizer)
    {
        uint16_t weight = 1;

        if (ump.find(minimizer) != ump.end())
        {
            weight = ump[minimizer];
        }

        return weight;
    }

    size_t AdditionalMap::GetSize()
    {
        return ump.size();
    }

    bool AdditionalMap::IsEmpty()
    {
        return ump.size() == 0;
    }

}