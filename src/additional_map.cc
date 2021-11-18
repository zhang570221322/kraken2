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

        if (ifile.is_open() && !isFileEmpty(ifile))
        {
            while (!ifile.eof())
            {
                string line;
                string seqID;
                taxid_t empty;
                taxid_t taxID;
                getline(ifile, line);
                istringstream linestream(line);
                linestream >> seqID >> empty >> taxID;
                _seqID_taxID.emplace(seqID, empty);
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
        }

        // return 1.0 / (0.81 + 0.05 * weight * weight);

        return 13.0 - 0.5 * weight - log(child_count + 1);
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
    // 1.找到score_taxid的map中的较为相近的score对应的taxid_t
    // 2.从hit_counts中找到taxid_t的祖先,填充到thread_conflict_ancestor中
    void AdditionalMap::findConflictLTR(taxon_counts_t &hit_counts, double max_score, uint64_t &conflicts, Taxonomy &taxonomy)
    {

        for (auto &kv_pair1 : score_taxid)
        {

            if (kv_pair1.first == max_score)
            {
                vector<taxid_t> &taxons = kv_pair1.second;

                if (taxons.size() > 1)
                {
                    conflicts += taxons.size();
                    for (auto taxon : taxons)
                    {
                        leaf.push_back(taxon);
                    }
                }
            }
            else if (labs(kv_pair1.first - max_score) <= max_score / 20.0)
            {
                vector<taxid_t> &taxons = kv_pair1.second;
                conflicts += taxons.size();
                for (auto taxon : taxons)
                {
                    leaf.push_back(taxon);
                }
            }
        }

        for (auto &kv_pair2 : hit_counts)
        {
            taxid_t taxon2 = kv_pair2.first;
            for (auto taxon : leaf)
            {
                if (taxonomy.IsAAncestorOfB(taxon2, taxon))
                {
                    conflict_ancestor.push_back(taxon2);
                }
            }
        }
    }
    // 持久化,将临时的temp,但是对于多线程必须是串行进行的.
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