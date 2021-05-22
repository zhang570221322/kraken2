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

                AddWeight(key, value);
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
    double AdditionalMap::GetWeight(uint64_t minimizer)
    //y =  1.0 / (0.81 + 0.1 * weight * weight);  X为这个kmer冲突次数
    {
        uint16_t weight = 1;

        if (conflict_ump.find(minimizer) != conflict_ump.end())
        {
            weight = conflict_ump[minimizer];
        }
        return 1.0 / (0.81 + 0.05 * weight * weight);
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
            AddWeight(minimizer, defalut_weight);
        }
    }

    void AdditionalMap::AddWeight(uint64_t minimizer, uint16_t weight)
    {
        conflict_ump.emplace(minimizer, weight);
    }
    void AdditionalMap::clearTemp()
    {

        temp.clear();
    }

    void AdditionalMap::addk_v2Temp(taxid_t taxon, uint64_t minimizer)
    {
        temp[taxon].push_back(minimizer);
    }
    void AdditionalMap::saveTemp()
    {
        for (size_t i = 0; i < ancestor.size(); i++)
        {
            auto code = ancestor[i];
            if (temp.find(code) != temp.end())
            {
                vector<uint64_t> kmers = temp[code];
                for (size_t i = 0; i < kmers.size(); i++)
                {
                    AddMinimizer(kmers[i]);
                }
            }
        }
        clearTemp();
    }

    size_t AdditionalMap::GetConflictSize()
    {
        return conflict_ump.size();
    }

    void AdditionalMap::AddAncestor(taxid_t taxa)
    {
        ancestor.push_back(taxa);
    }
    void AdditionalMap::clearAncestor()
    {
        ancestor.clear();
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