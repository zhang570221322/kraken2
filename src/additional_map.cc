#include "additional_map.h"

using std::cout;
using std::endl;
using std::getline;
using std::ifstream;
using std::istringstream;
using std::ofstream;
using std::string;
using std::unordered_map;
using std::vector;

namespace kraken2
{

    bool isFileEmpty(std::ifstream &pFile)
    {
        return pFile.peek() == ifstream::traits_type::eof();
    }

    void AdditionalMap::loadGod_data(const char *filename)
    {
        ifstream ifile(filename);

        string line;
        string first;
        taxid_t second;
        taxid_t third;
        string forth;
        if (ifile.is_open() && !isFileEmpty(ifile))
        {
            while (!ifile.eof())
            {

                getline(ifile, line);
                istringstream linestream(line);
                linestream >> first >> second >> third >> forth;
                _seqID_taxID.emplace(first, second);
            }
        }
        else if (!ifile.is_open())
        {
            ofstream outfile(filename);
            outfile << "";
            outfile.close();
        }
        ifile.close();
    }

    void AdditionalMap::clearConflictObj()
    {
        conflict_ump.clear();
        conflict_temp.clear();
        conflict_ancestor.clear();
        leaf.clear();
        leaf_count.clear();
        score_taxid.clear();
    }
    void AdditionalMap::ReadConflictFile(const char *filename)
    {
        ifstream ifile(filename);

        if (ifile.is_open() && !isFileEmpty(ifile))
        {
            while (!ifile.eof())
            {
                string line;
                uint64_t key;
                taxid_t value;
                getline(ifile, line);
                istringstream linestream(line);
                linestream >> key >> value;
                conflict_ump.emplace(key, value);
                ;
            }
        }
        else if (!ifile.is_open())
        {
            ofstream outfile(filename);
            outfile << "";
            outfile.close();
        }
        ifile.close();
    }
    void AdditionalMap::WriteConflictMap(const char *filename)
    {
        if (GetConflictSize())
        {
            ofstream ofile(filename);
            if (ofile.is_open())
            {
                for (auto x : conflict_ump)
                {
                    ofile << x.first << "\t" << x.second;
                    ofile << "\n";
                }
            }
            ofile.close();
        }
    }
    double AdditionalMap::GetWeight(uint64_t minimizer, uint64_t child_count)
    //score = 1.05 - (0.05*x)-log(y),  X为这个kmer冲突次数,y为此taxo的子节点数,
    {
        uint16_t weight = 1;

        if (conflict_ump.find(minimizer) != conflict_ump.end())
        {
            weight = conflict_ump[minimizer];
            return 13.0 - 0.05 * weight - log(child_count + 1);
        }
        else
        {
            return 1.0;
        }
    }

    void AdditionalMap::AddMinimizer(uint64_t minimizer)
    {

        uint16_t defalut_weight = 2;
        //如果已经存在次数,再继续提升次数
        if (conflict_ump.find(minimizer) != conflict_ump.end())
        {
            conflict_ump[minimizer] = 1 + conflict_ump[minimizer];
        }
        else
        {
            conflict_ump.emplace(minimizer, defalut_weight);
        }
    }
    // Func 寻找到冲突路径,将冲突路径上的taxid_t填充到thread_conflict_ancestor中
    // 1. leaf中存储当前read的叶子节点,判断当前叶子节点分数是否符合条件，如果符合，设置为ture
    // 2. 对于刚好等于max_score的,直接设置为ture
    void AdditionalMap::findConflictLTR(taxon_counts_t &hit_counts, double max_score, uint64_t &conflicts, Taxonomy &taxonomy)
    {
        vector<taxid_t> &max_taxons = score_taxid[max_score];
        for (auto &kv_pair1 : score_taxid)
        {
            vector<taxid_t> &taxons = kv_pair1.second;
            // 处理相近的叶子结点
            for (auto taxon : taxons)
            {

                if (labs(kv_pair1.first - max_score) <= max_score / 20.0 && kv_pair1.first != max_score)
                {

                    for (auto &leaf_taxon_pair : leaf)
                    {
                        if (taxon == leaf_taxon_pair.first)
                        {
                            // 设置为冲突路径
                            conflicts += 1;
                            leaf[leaf_taxon_pair.first] = true;
                        }
                    }
                }
            }
        }

        if (max_taxons.size() > 1)
        {

            for (auto taxon : max_taxons)
            {
                leaf[taxon] = true;
                conflicts += 1;
            }
        }
        else
        {
            // 如果有相近的路径,则原来路径也是冲突路径
            if (conflicts)
            {
                leaf[max_taxons.front()] = true;
                conflicts += 1;
            }
        }

        for (auto &kv_pair2 : hit_counts)
        {
            taxid_t taxon2 = kv_pair2.first;
            for (auto &leaf_taxon : leaf)
            {
                if (taxon2 != leaf_taxon.first && leaf_taxon.second && taxonomy.IsAAncestorOfB(taxon2, leaf_taxon.first))
                {
                    conflict_ancestor.push_back(taxon2);
                }
            }
        }
    }

    void AdditionalMap::saveTemp(taxon_counts_t &hit_counts, double max_score, uint64_t &conflicts, Taxonomy &taxonomy)
    {
        // 先计算出需要降低的tax,存入conflict_ancestor
        findConflictLTR(hit_counts, max_score, conflicts, taxonomy);
        for (auto &code : conflict_ancestor)
        {

            if (conflict_temp.find(code) != conflict_temp.end())
            {
                vector<uint64_t> &kmers = conflict_temp[code];
                for (size_t i = 0; i < kmers.size(); i++)
                {
                    // 持久化,将每个read的冲突kmer存入conflict_ump,但是对于多线程必须是串行进行的.
                    AddMinimizer(kmers[i]);
                }
            }
        }
    }

    size_t AdditionalMap::GetConflictSize()
    {
        return conflict_ump.size();
    }

    //对于未识别来说

    void AdditionalMap::ReadAdditionalFile(const char *filename)
    {
        ifstream ifile(filename);

        if (ifile.is_open() && !isFileEmpty(ifile))
        {

            while (!ifile.eof())
            {
                string line;
                uint64_t key;
                taxid_t value;

                getline(ifile, line);
                istringstream linestream(line);
                linestream >> key >> value;

                AddTaxon(key, value);
            }
        }
        else if (!ifile.is_open())
        {
            ofstream outfile(filename);
            outfile << "";
            outfile.close();
        }

        ifile.close();
    }
    void AdditionalMap::WriteAdditionalMap(const char *filename)
    {
        if (GetAdditionaSize())
        {
            ofstream ofile(filename);
            if (ofile.is_open())
            {
                for (auto x : additional_ump)
                {
                    ofile << x.first << "\t" << x.second;
                    ofile << "\n";
                }
            }

            ofile.close();
        }
    }
    void AdditionalMap::AddAdditionalPair(uint64_t minimizer, taxid_t tax_id, Taxonomy &taxonomy)
    {
        taxid_t new_tax = tax_id;
        taxid_t old_tax = 0;

        if (additional_ump.find(minimizer) != additional_ump.end())
        {
            old_tax = additional_ump[minimizer];
            new_tax = taxonomy.LowestCommonAncestor(new_tax, old_tax);
            additional_ump[minimizer] = new_tax;
        }
        else
        {
            AddTaxon(minimizer, new_tax);
        }
    }
    taxid_t AdditionalMap::GetTax(uint64_t minimizer)
    {
        taxid_t tax = 0;

        if (additional_ump.find(minimizer) != additional_ump.end())
        {
            tax = additional_ump[minimizer];
        }

        return tax;
    }
    size_t AdditionalMap::GetAdditionaSize()
    {
        return additional_ump.size();
    }

    void AdditionalMap::AddTaxon(uint64_t minimizer, taxid_t tax_id)
    {
        additional_ump.emplace(minimizer, tax_id);
    }
}