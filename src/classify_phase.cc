/*
 * Copyright 2013-2020, Derrick Wood <dwood@cs.jhu.edu>
 *
 * This file is part of the Kraken 2 taxonomic sequence classification system.
 */
#include <math.h>
#include "kraken2_headers.h"
#include "kv_store.h"
#include "taxonomy.h"
#include "seqreader.h"
#include "mmscanner.h"
#include "compact_hash.h"
#include "kraken2_data.h"
#include "aa_translate.h"
#include "reports.h"
#include "utilities.h"
#include "readcounts.h"
//新加的类.
#include "additional_map.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::map;
using std::ofstream;
using std::ostringstream;
using std::set;
using std::string;
using std::vector;
using namespace kraken2;

static const size_t NUM_FRAGMENTS_PER_THREAD = 10000;
static const taxid_t MATE_PAIR_BORDER_TAXON = TAXID_MAX;
static const taxid_t READING_FRAME_BORDER_TAXON = TAXID_MAX - 1;
static const taxid_t AMBIGUOUS_SPAN_TAXON = TAXID_MAX - 2;

struct Options
{
  string index_filename;
  string taxonomy_filename;
  string options_filename;
  string report_filename;
  string classified_output_filename;
  string unclassified_output_filename;
  string kraken_output_filename;
  string conflict_filename;
  bool mpa_style_report;
  bool report_kmer_data;
  bool quick_mode;
  bool report_zero_counts;
  bool use_translated_search;
  bool print_scientific_name;
  double confidence_threshold;
  int num_threads;
  bool paired_end_processing;
  bool single_file_pairs;
  int minimum_quality_score;
  int minimum_hit_groups;
  bool use_memory_mapping;
  bool match_input_order;
};

struct ClassificationStats
{
  uint64_t total_sequences;
  uint64_t total_bases;
  uint64_t total_classified;
  uint64_t total_conflict;
};

struct OutputStreamData
{
  bool initialized;
  bool printing_sequences;
  std::ostream *classified_output1;
  std::ostream *classified_output2;
  std::ostream *unclassified_output1;
  std::ostream *unclassified_output2;
  std::ostream *kraken_output;
};

struct OutputData
{
  uint64_t block_id;
  string kraken_str;
  string classified_out1_str;
  string classified_out2_str;
  string unclassified_out1_str;
  string unclassified_out2_str;
};

void ParseCommandLine(int argc, char **argv, Options &opts);
void usage(int exit_code = EX_USAGE);
void ProcessFiles(const char *filename1, const char *filename2,
                  KeyValueStore *hash, Taxonomy &tax,
                  IndexOptions &idx_opts, Options &opts, ClassificationStats &stats,
                  OutputStreamData &outputs, taxon_counters_t &total_taxon_counters, AdditionalMap &add_map);
taxid_t ClassifySequence(Sequence &dna, Sequence &dna2, ostringstream &koss,
                         KeyValueStore *hash, Taxonomy &tax, IndexOptions &idx_opts,
                         Options &opts, ClassificationStats &stats, MinimizerScanner &scanner,
                         vector<taxid_t> &taxa, taxon_counts_t &hit_counts,
                         vector<string> &tx_frames, taxon_counters_t &my_taxon_counts, AdditionalMap &add_map);
taxid_t ResolveTree(taxon_counts_t &hit_counts, Taxonomy &tax, size_t total_minimizers, Options &opts,
                    AdditionalMap &add_map, ostringstream &koss, ClassificationStats &stats);
void ReportStats(struct timeval time1, struct timeval time2, ClassificationStats &stats);

void MaskLowQualityBases(Sequence &dna, int minimum_quality_score);

int main(int argc, char **argv)
{
  Options opts;
  opts.quick_mode = false;
  opts.confidence_threshold = 0;
  opts.paired_end_processing = false;
  opts.single_file_pairs = false;
  opts.num_threads = 1;
  opts.mpa_style_report = false;
  opts.report_kmer_data = false;
  opts.report_zero_counts = false;
  opts.use_translated_search = false;
  opts.print_scientific_name = false;
  opts.minimum_quality_score = 0;
  opts.minimum_hit_groups = 0;
  opts.use_memory_mapping = false;

  taxon_counters_t taxon_counters; // stats per taxon
  ParseCommandLine(argc, argv, opts);

  omp_set_num_threads(opts.num_threads);
  cerr << "Create/update conflict hash map." << endl;
  IndexOptions idx_opts = {0};
  ifstream idx_opt_fs(opts.options_filename);
  struct stat sb;
  if (stat(opts.options_filename.c_str(), &sb) < 0)
    errx(EX_OSERR, "unable to get filesize of %s", opts.options_filename.c_str());

  cerr << "Loading database information...";
  auto opts_filesize = sb.st_size;
  idx_opt_fs.read((char *)&idx_opts, opts_filesize);
  opts.use_translated_search = !idx_opts.dna_db;
  Taxonomy taxonomy(opts.taxonomy_filename, opts.use_memory_mapping);
  KeyValueStore *hash_ptr = new CompactHashTable(opts.index_filename, opts.use_memory_mapping);
  cerr << "done." << endl;

  cerr << "Loading conflict hashmap...";
  AdditionalMap add_map;
  add_map.ReadConflictFile(opts.conflict_filename.c_str());
  cerr << "done." << endl;
  cerr << "Exec classification conflict analysis." << endl;
  ClassificationStats stats = {0, 0, 0, 0};

  OutputStreamData outputs = {false, false, nullptr, nullptr, nullptr, nullptr, &std::cout};

  struct timeval tv1, tv2;
  gettimeofday(&tv1, nullptr);
  if (optind == argc)
  {
    if (opts.paired_end_processing && !opts.single_file_pairs)
      errx(EX_USAGE, "paired end processing used with no files specified");
    ProcessFiles(nullptr, nullptr, hash_ptr, taxonomy, idx_opts, opts, stats, outputs, taxon_counters, add_map);
  }
  else
  {
    for (int i = optind; i < argc; i++)
    {
      if (opts.paired_end_processing && !opts.single_file_pairs)
      {
        if (i + 1 == argc)
        {
          errx(EX_USAGE, "paired end processing used with unpaired file");
        }
        ProcessFiles(argv[i], argv[i + 1], hash_ptr, taxonomy, idx_opts, opts, stats, outputs, taxon_counters, add_map);
        i += 1;
      }
      else
      {
        ProcessFiles(argv[i], nullptr, hash_ptr, taxonomy, idx_opts, opts, stats, outputs, taxon_counters, add_map);
      }
    }
  }
  gettimeofday(&tv2, nullptr);
  cerr << "\nWrite conflict kmer into " << opts.conflict_filename.c_str() << " (kmers:" << add_map.GetConflictSize() << ")" << endl;
  add_map.WriteConflictMap(opts.conflict_filename.c_str());

  delete hash_ptr;
  ReportStats(tv1, tv2, stats);
  return 0;
}

void ReportStats(struct timeval time1, struct timeval time2,
                 ClassificationStats &stats)
{
  time2.tv_usec -= time1.tv_usec;
  time2.tv_sec -= time1.tv_sec;
  if (time2.tv_usec < 0)
  {
    time2.tv_sec--;
    time2.tv_usec += 1000000;
  }
  double seconds = time2.tv_usec;
  seconds /= 1e6;
  seconds += time2.tv_sec;

  if (isatty(fileno(stderr)))
    cerr << "\r";
  fprintf(stderr,
          "%llu sequences (%.2f Mbp) (Find conflicts %lu times)  processed in %.3fs (%.1f Kseq/m, %.2f Mbp/m).  \n",
          (unsigned long long)stats.total_sequences,
          stats.total_bases / 1.0e6,
          stats.total_conflict,
          seconds,
          stats.total_sequences / 1.0e3 / (seconds / 60),
          stats.total_bases / 1.0e6 / (seconds / 60));
}
void ProcessFiles(const char *filename1, const char *filename2,
                  KeyValueStore *hash, Taxonomy &tax,
                  IndexOptions &idx_opts, Options &opts, ClassificationStats &stats,
                  OutputStreamData &outputs,
                  taxon_counters_t &total_taxon_counters, AdditionalMap &add_map)
{
  std::istream *fptr1 = nullptr, *fptr2 = nullptr;

  if (filename1 == nullptr)
    fptr1 = &std::cin;
  else
  {
    fptr1 = new std::ifstream(filename1);
  }
  if (opts.paired_end_processing && !opts.single_file_pairs)
  {
    fptr2 = new std::ifstream(filename2);
  }

  // The priority queue for output is designed to ensure fragment data
  // is output in the same order it was input
  auto comparator = [](const OutputData &a, const OutputData &b)
  {
    return a.block_id > b.block_id;
  };
  std::priority_queue<OutputData, vector<OutputData>, decltype(comparator)>
      output_queue(comparator);
  uint64_t next_input_block_id = 0;
  uint64_t next_output_block_id = 0;
  omp_lock_t output_lock;
  omp_init_lock(&output_lock);

#pragma omp parallel
  {
    MinimizerScanner scanner(idx_opts.k, idx_opts.l, idx_opts.spaced_seed_mask,
                             idx_opts.dna_db, idx_opts.toggle_mask,
                             idx_opts.revcom_version);
    vector<taxid_t> taxa;
    taxon_counts_t hit_counts;
    ostringstream kraken_oss, c1_oss, c2_oss, u1_oss, u2_oss;
    ClassificationStats thread_stats = {0, 0, 0};
    vector<string> translated_frames(6);
    BatchSequenceReader reader1, reader2;
    Sequence seq1, seq2;
    uint64_t block_id;
    OutputData out_data;
    taxon_counters_t thread_taxon_counters;
    //多线程添加冲突
    AdditionalMap thread_add_map;

    while (true)
    {
      thread_stats.total_sequences = 0;
      thread_stats.total_bases = 0;
      thread_stats.total_conflict = 0;
      auto ok_read = false;
      thread_add_map.clearConflictObj();
#pragma omp critical(seqread)
      { // Input processing block
        if (!opts.paired_end_processing)
        {
          // Unpaired data?  Just read in a sized block
          ok_read = reader1.LoadBlock(*fptr1, (size_t)(3 * 1024 * 1024));
        }
        else if (!opts.single_file_pairs)
        {
          // Paired data in 2 files?  Read a line-counted batch from each file.
          ok_read = reader1.LoadBatch(*fptr1, NUM_FRAGMENTS_PER_THREAD);
          if (ok_read && opts.paired_end_processing)
            ok_read = reader2.LoadBatch(*fptr2, NUM_FRAGMENTS_PER_THREAD);
        }
        else
        {
          auto frags = NUM_FRAGMENTS_PER_THREAD * 2;
          // Ensure frag count is even - just in case above line is changed
          if (frags % 2 == 1)
            frags++;
          ok_read = reader1.LoadBatch(*fptr1, frags);
        }
        block_id = next_input_block_id++;
      }

      if (!ok_read)
        break;

      // Reset all dynamically-growing things
      kraken_oss.str("");
      c1_oss.str("");
      c2_oss.str("");
      u1_oss.str("");
      u2_oss.str("");
      thread_taxon_counters.clear();

      while (true)
      {
        auto valid_fragment = reader1.NextSequence(seq1);
        if (opts.paired_end_processing && valid_fragment)
        {
          if (opts.single_file_pairs)
            valid_fragment = reader1.NextSequence(seq2);
          else
            valid_fragment = reader2.NextSequence(seq2);
        }
        if (!valid_fragment)
          break;
        thread_stats.total_sequences++;
        if (opts.minimum_quality_score > 0)
        {
          MaskLowQualityBases(seq1, opts.minimum_quality_score);
          if (opts.paired_end_processing)
            MaskLowQualityBases(seq2, opts.minimum_quality_score);
        }

        uint64_t read_conflicts = ClassifySequence(seq1, seq2,
                                                   kraken_oss, hash, tax, idx_opts, opts, thread_stats, scanner,
                                                   taxa, hit_counts, translated_frames, thread_taxon_counters, thread_add_map);
        // if (read_conflicts)
        // {
        //   // 打印seq_id
        //   cout << "seq_id: " << seq1.id << ", read_conflicts_times: " << read_conflicts << endl;
        // }
        thread_stats.total_bases += seq1.seq.size();
        thread_stats.total_conflict += read_conflicts;
        if (opts.paired_end_processing)
          thread_stats.total_bases += seq2.seq.size();
      }

#pragma omp atomic
      stats.total_sequences += thread_stats.total_sequences;
#pragma omp atomic
      stats.total_bases += thread_stats.total_bases;
#pragma omp atomic
      stats.total_conflict += thread_stats.total_conflict;
#pragma omp critical(output_stats)
      {
        if (isatty(fileno(stderr)))
          cerr << "\rProcessed " << stats.total_sequences
               << " sequences (" << stats.total_bases << " bp)"
               << " conflicts (" << stats.total_conflict << " times) ... ";
      }

      out_data.block_id = block_id;
      out_data.kraken_str.assign(kraken_oss.str());
      out_data.classified_out1_str.assign(c1_oss.str());
      out_data.classified_out2_str.assign(c2_oss.str());
      out_data.unclassified_out1_str.assign(u1_oss.str());
      out_data.unclassified_out2_str.assign(u2_oss.str());

#pragma omp critical(output_queue)
      {
        output_queue.push(out_data);
      }
      // 合并所有的冲突的kmer
#pragma omp critical(merge_conflict_map)
      {
        unordered_map<uint64_t, uint64_t> &shared_conflict_ump = add_map.conflict_ump;

        for (auto &kv_pair : thread_add_map.conflict_ump)
        {
          shared_conflict_ump[kv_pair.first] += kv_pair.second;
        }
      }
#pragma omp critical(update_taxon_counters)
      {
        for (auto &kv_pair : thread_taxon_counters)
        {
          total_taxon_counters[kv_pair.first] += std::move(kv_pair.second);
        }
      }

      bool output_loop = true;
      while (output_loop)
      {
#pragma omp critical(output_queue)
        {
          output_loop = !output_queue.empty();
          if (output_loop)
          {
            out_data = output_queue.top();
            if (out_data.block_id == next_output_block_id)
            {
              output_queue.pop();
              // Acquiring output lock obligates thread to print out
              // next output data block, contained in out_data
              omp_set_lock(&output_lock);
              next_output_block_id++;
            }
            else
              output_loop = false;
          }
        }
        if (!output_loop)
          break;
        if (outputs.kraken_output != nullptr)
          (*outputs.kraken_output) << out_data.kraken_str;
        if (outputs.classified_output1 != nullptr)
          (*outputs.classified_output1) << out_data.classified_out1_str;
        if (outputs.classified_output2 != nullptr)
          (*outputs.classified_output2) << out_data.classified_out2_str;
        if (outputs.unclassified_output1 != nullptr)
          (*outputs.unclassified_output1) << out_data.unclassified_out1_str;
        if (outputs.unclassified_output2 != nullptr)
          (*outputs.unclassified_output2) << out_data.unclassified_out2_str;
        omp_unset_lock(&output_lock);
      } // end while output loop
    }   // end while
  }     // end parallel block
  omp_destroy_lock(&output_lock);
  if (fptr1 != nullptr)
    delete fptr1;
  if (fptr2 != nullptr)
    delete fptr2;
  if (outputs.kraken_output != nullptr)
    (*outputs.kraken_output) << std::flush;
  if (outputs.classified_output1 != nullptr)
    (*outputs.classified_output1) << std::flush;
  if (outputs.classified_output2 != nullptr)
    (*outputs.classified_output2) << std::flush;
  if (outputs.unclassified_output1 != nullptr)
    (*outputs.unclassified_output1) << std::flush;
  if (outputs.unclassified_output2 != nullptr)
    (*outputs.unclassified_output2) << std::flush;
}

taxid_t ResolveTree(taxon_counts_t &hit_counts, Taxonomy &taxonomy, size_t total_minimizers, Options &opts, AdditionalMap &thread_add_map,
                    ostringstream &koss, ClassificationStats &stats)
{

  double max_score = 0;
  // 取到引用,存储最终的剪枝树分数.
  unordered_map<double, vector<taxid_t>> &score_taxid = thread_add_map.score_taxid;
  uint64_t conflicts = 0;
  // Sum each taxon's LTR path, find taxon with highest LTR score
  for (auto &kv_pair : hit_counts)
  {
    taxid_t taxon = kv_pair.first;
    double score = 0;
    bool isLeaf = true;
    for (auto &kv_pair2 : hit_counts)
    {
      taxid_t taxon2 = kv_pair2.first;

      if (taxonomy.IsAAncestorOfB(taxon2, taxon))
      {

        score += kv_pair2.second;
      }
      // 如果taxon是其中一个的祖先,则不是叶子节点
      if (isLeaf && taxon != taxon2 && taxonomy.IsAAncestorOfB(taxon, taxon2))
      {
        isLeaf = false;
      }
    }
    // 如果是叶子节点,存入叶子vector
    if (isLeaf)
    {
      // 加入叶子节点,但是不被选中
      thread_add_map.leaf[taxon] = false;
    }
    if (score > max_score)
    {
      max_score = score;
    }

    score_taxid[score].push_back(taxon);
  }
  // 打印分数和对应的taxon
  // for (auto &kv_pair : score_taxid)
  // {
  //   cout << "分数:" << kv_pair.first << "路径=";
  //   for (auto temp_tax : kv_pair.second)
  //   {
  //     cout << taxonomy.nodes()[temp_tax].external_id << ",";
  //   }
  //   cout << endl;
  // }
  // 打印所有命中
  // cout << "所有命中";
  // for (auto &ht_pair : hit_counts)
  // {
  //   cout << taxonomy.nodes()[ht_pair.first].external_id << ":" << ht_pair.second << " ";
  // }
  // cout << endl;
  // 打印叶子节点
  // cout << "叶子节点";
  // for (auto &leaf_pair : thread_add_map.leaf)
  // {
  //   cout << taxonomy.nodes()[leaf_pair.first].external_id << ":" << leaf_pair.second << " ";
  // }
  // cout << endl;
  // 达到阈值,计算并写入所有的冲突kmer
  thread_add_map.saveTemp(hit_counts, max_score, conflicts, taxonomy);
  return conflicts;
}

taxid_t ClassifySequence(Sequence &dna, Sequence &dna2, ostringstream &koss,
                         KeyValueStore *hash, Taxonomy &taxonomy, IndexOptions &idx_opts,
                         Options &opts, ClassificationStats &stats, MinimizerScanner &scanner,
                         vector<taxid_t> &taxa, taxon_counts_t &hit_counts, vector<string> &tx_frames,
                         taxon_counters_t &curr_taxon_counts, AdditionalMap &thread_add_map)
{
  uint64_t *minimizer_ptr;
  uint64_t conflict = 0;
  taxa.clear();
  hit_counts.clear();
  //  在这里做初始化清空操作
  thread_add_map.conflict_temp.clear();
  thread_add_map.conflict_ancestor.clear();
  thread_add_map.leaf.clear();
  thread_add_map.score_taxid.clear();

  auto frame_ct = opts.use_translated_search ? 6 : 1;
  int64_t minimizer_hit_groups = 0;

  for (int mate_num = 0; mate_num < 2; mate_num++)
  {
    if (mate_num == 1 && !opts.paired_end_processing)
      break;

    if (opts.use_translated_search)
    {
      TranslateToAllFrames(mate_num == 0 ? dna.seq : dna2.seq, tx_frames);
    }
    // index of frame is 0 - 5 w/ tx search (or 0 if no tx search)
    for (int frame_idx = 0; frame_idx < frame_ct; frame_idx++)
    {
      if (opts.use_translated_search)
      {
        scanner.LoadSequence(tx_frames[frame_idx]);
      }
      else
      {
        scanner.LoadSequence(mate_num == 0 ? dna.seq : dna2.seq);
      }
      uint64_t last_minimizer = UINT64_MAX;
      taxid_t last_taxon = TAXID_MAX;

      while ((minimizer_ptr = scanner.NextMinimizer()) != nullptr)
      {

        taxid_t taxon;
        if (scanner.is_ambiguous())
        {
          taxon = AMBIGUOUS_SPAN_TAXON;
        }
        else
        {
          if (*minimizer_ptr != last_minimizer)
          {
            bool skip_lookup = false;
            if (idx_opts.minimum_acceptable_hash_value)
            {
              if (MurmurHash3(*minimizer_ptr) < idx_opts.minimum_acceptable_hash_value)
                skip_lookup = true;
            }
            taxon = 0;

            if (!skip_lookup)
            {
              taxon = hash->Get(*minimizer_ptr);
            }
            last_taxon = taxon;
            last_minimizer = *minimizer_ptr;

            // Increment this only if (a) we have DB hit and
            // (b) minimizer != last minimizer
            if (taxon)
            {
              minimizer_hit_groups++;
              // New minimizer should trigger registering minimizer in RC/HLL
              auto last = scanner.last_minimizer();
              // 命中为非leaf的kmer,将其添加到temp.LCA
              if (taxonomy.nodes()[taxon].child_count)
                thread_add_map.conflict_temp[taxon].push_back(last);
              curr_taxon_counts[taxon].add_kmer(last);
            }
          }
          else
          {
            taxon = last_taxon;
          }
          if (taxon)
          {
            if (opts.quick_mode && minimizer_hit_groups >= opts.minimum_hit_groups)
            {

              goto finished_searching; // need to break 3 loops here
            }
            hit_counts[taxon] += 1;
          }
        }
        taxa.push_back(taxon);
      }
      if (opts.use_translated_search && frame_idx != 5)
        taxa.push_back(READING_FRAME_BORDER_TAXON);
    }
    if (opts.paired_end_processing && mate_num == 0)
      taxa.push_back(MATE_PAIR_BORDER_TAXON);
  }

finished_searching:

  auto total_kmers = taxa.size();
  if (opts.paired_end_processing)
    total_kmers--;                // account for the mate pair marker
  if (opts.use_translated_search) // account for reading frame markers
    total_kmers -= opts.paired_end_processing ? 4 : 2;
  conflict = ResolveTree(hit_counts, taxonomy, total_kmers, opts, thread_add_map, koss, stats);

  return conflict;
}

void MaskLowQualityBases(Sequence &dna, int minimum_quality_score)
{
  if (dna.format != FORMAT_FASTQ)
    return;
  if (dna.seq.size() != dna.quals.size())
    errx(EX_DATAERR, "%s: Sequence length (%d) != Quality string length (%d)",
         dna.id.c_str(), (int)dna.seq.size(), (int)dna.quals.size());
  for (size_t i = 0; i < dna.seq.size(); i++)
  {
    if ((dna.quals[i] - '!') < minimum_quality_score)
      dna.seq[i] = 'x';
  }
}

void ParseCommandLine(int argc, char **argv, Options &opts)
{
  int opt;

  while ((opt = getopt(argc, argv, "h?H:A:t:o:T:p:R:C:U:O:Q:g:nmzqPSMK")) != -1)
  {
    switch (opt)
    {
    case 'h':
    case '?':
      usage(0);
      break;
    case 'H':
      opts.index_filename = optarg;
      break;
    case 'A':
      opts.conflict_filename = optarg;
      break;
    case 't':
      opts.taxonomy_filename = optarg;
      break;
    case 'T':
      opts.confidence_threshold = std::stod(optarg);
      if (opts.confidence_threshold < 0 || opts.confidence_threshold > 1)
      {
        errx(EX_USAGE, "confidence threshold must be in [0, 1]");
      }
      break;
    case 'o':
      opts.options_filename = optarg;
      break;
    case 'q':
      opts.quick_mode = true;
      break;
    case 'p':
      opts.num_threads = atoi(optarg);
      if (opts.num_threads < 1)
        errx(EX_USAGE, "number of threads can't be less than 1");
      break;
    case 'g':
      opts.minimum_hit_groups = atoi(optarg);
      break;
    case 'P':
      opts.paired_end_processing = true;
      break;
    case 'S':
      opts.paired_end_processing = true;
      opts.single_file_pairs = true;
      break;
    case 'm':
      opts.mpa_style_report = true;
      break;
    case 'K':
      opts.report_kmer_data = true;
      break;
    case 'R':
      opts.report_filename = optarg;
      break;
    case 'z':
      opts.report_zero_counts = true;
      break;
    case 'C':
      opts.classified_output_filename = optarg;
      break;
    case 'U':
      opts.unclassified_output_filename = optarg;
      break;
    case 'O':
      opts.kraken_output_filename = optarg;
      break;
    case 'n':
      opts.print_scientific_name = true;
      break;
    case 'Q':
      opts.minimum_quality_score = atoi(optarg);
      break;
    case 'M':
      opts.use_memory_mapping = true;
      break;
    }
  }

  if (opts.index_filename.empty() ||
      opts.taxonomy_filename.empty() ||
      opts.options_filename.empty())
  {
    warnx("mandatory filename missing");
    usage();
  }

  if (opts.mpa_style_report && opts.report_filename.empty())
  {
    warnx("-m requires -R be used");
    usage();
  }
}

void usage(int exit_code)
{
  cerr << "Usage: classify [options] <fasta/fastq file(s)>" << endl
       << endl
       << "Options: (*mandatory)" << endl
       << "* -H filename      Kraken 2 index filename" << endl
       << "* -t filename      Kraken 2 taxonomy filename" << endl
       << "* -o filename      Kraken 2 options filename" << endl
       << "  -q               Quick mode" << endl
       << "  -M               Use memory mapping to access hash & taxonomy" << endl
       << "  -T NUM           Confidence score threshold (def. 0)" << endl
       << "  -p NUM           Number of threads (def. 1)" << endl
       << "  -Q NUM           Minimum quality score (FASTQ only, def. 0)" << endl
       << "  -P               Process pairs of reads" << endl
       << "  -S               Process pairs with mates in same file" << endl
       << "  -R filename      Print report to filename" << endl
       << "  -m               In comb. w/ -R, use mpa-style report" << endl
       << "  -z               In comb. w/ -R, report taxa w/ 0 count" << endl
       << "  -n               Print scientific name instead of taxid in Kraken output" << endl
       << "  -g NUM           Minimum number of hit groups needed for call" << endl
       << "  -C filename      Filename/format to have classified sequences" << endl
       << "  -U filename      Filename/format to have unclassified sequences" << endl
       << "  -O filename      Output file for normal Kraken output" << endl
       << "  -K               In comb. w/ -R, provide minimizer information in report" << endl;
  std::exit(exit_code);
}
