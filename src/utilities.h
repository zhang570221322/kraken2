/*
 * Copyright 2013-2020, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */

#ifndef KRAKEN2_UTILITIES_H_
#define KRAKEN2_UTILITIES_H_

#include "kraken2_headers.h"

#include "taxonomy.h"
// Functions used by 2+ programs that I couldn't think of a better place for.
const std::string species_types[] = {
    "species",
    "subspecies",
    "variety",
    "subvariety",
    "form",
    "subform",
    "strain",
    "no rank",
    "species subgroup",
    "species group",
    "forma specialis",
    "solate"};
const std::string genus_types[] = {
    "genus",
    "subgenus",
    "section",
    "subsection",
    "series",
    "subseries"};
const std::string others_types[] = {
    "root",
    "superkingdom",
    "kingdom",
    "subkingdom ",
    "phylum",
    "subphylum",
    "class",
    "subclass",
    "superorder",
    "order",
    "suborder",
    "family",
    "subfamily",
    "tribe",
    "subtribe",
    "biotype"};
namespace kraken2
{

    // Turns a simple bitstring of 0s and 1s into one where the 1s and 0s are expanded,
    // e.g. 010110 expanded by a factor of 2 would become 001100111100
    // Allows specification of the spaced seed to represent positions in the sequence
    // rather than actual bits in the internal representation
    void ExpandSpacedSeedMask(uint64_t &spaced_seed_mask, const int bit_expansion_factor);

    std::vector<std::string> SplitString(const std::string &str,
                                         const std::string &delim = "\t", const size_t max_fields = (size_t)-1);

    bool IsSpecies(Taxonomy &tax, TaxonomyNode &node);

    bool IsGenus(Taxonomy &tax, TaxonomyNode &node);

    bool IsOther(Taxonomy &tax, TaxonomyNode &node);
}

#endif
