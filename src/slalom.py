import sys, argparse
from slalom_structures import GlobalState, EnrichmentCountType
from slalom_auxiliar import CustomHelpFormatter, ArgumentProcessor, CSVParser, InputFileProcessor, DataProcessor

if __name__ != '__main__':
    sys.exit()

#Parsing input arguments
usage = '%(prog)s [options] [-s SEQ_LEN_DB_FILE] [-m GROUP_MAP_FILE] -a1 ANNO_1_FILE -a2 ANNO_2_FILE -o OUTPUT_FILE'
version = '%(prog)s SLALOM version 2.1.1b'
arg_parser = argparse.ArgumentParser(usage = usage, allow_abbrev = False, formatter_class=CustomHelpFormatter)
arg_parser._optionals.title = None
arg_parser.description = 'Welcome to SLALOM (StatisticaL Analysis of Locus Overlap Method)! Abbreviations: SID = sequence identifier; GID = group identifier'
arg_parser.add_argument('--version', action = 'version', version = version)
main_files = arg_parser.add_argument_group('----- Basic options -----\n\nMain input/output files')
simplified_mode = arg_parser.add_argument_group('Simplified modes').add_mutually_exclusive_group()
operating_mode = arg_parser.add_argument_group('Operating mode setup')
core_controls = arg_parser.add_argument_group('Core algorithm controls')
input_format = arg_parser.add_argument_group('----- Advanced options -----\n\nInput file format')
input_alternatives = arg_parser.add_argument_group('Input alternatives (reduce the need in files/columns for simplified cases)')
input_controls = arg_parser.add_argument_group('Input controls')
output_files_extra = arg_parser.add_argument_group('Additional output files')
output_controls = arg_parser.add_argument_group('Output controls')
other_options = arg_parser.add_argument_group('Other options')
main_files.add_argument('-s', '--seqlenfile', metavar = 'SEQ_LEN_DB_FILE', dest = 'len_db', type = str, default = '', help = 'Input file with the database of sequence lengths')
main_files.add_argument('-m', '--mapfile', metavar = 'GROUP_MAP_FILE', dest = 'group_map', type = str, default = '', help = 'Input file with the sequence group mapping')
main_files.add_argument('-a1', '--anno1file', metavar = 'ANNO_1_FILE', dest = 'anno1', type = str, required = True, help = 'Input file with the first annotation')
main_files.add_argument('-a2', '--anno2file', metavar = 'ANNO_2_FILE', dest = 'anno2', type = str, required = True, help = 'Input file with the second annotation')
main_files.add_argument('-o', '--outfile', dest = 'output_file', type = str, required = True, help = 'Output TSV file with calculated similarity/performance measures')
simplified_mode.add_argument('--genbank', dest = 'genbank', action = 'store_true', help = 'Compare a pair of genomes in GenBank format')
simplified_mode.add_argument('--bed', dest = 'bed', action = 'store_true', help = 'Compare a pair of genomes in BED format')
operating_mode.add_argument('-b', '--benchmarking', dest = 'benchmark', action = 'store_true', help = 'Treat the first annotation as benchmark (default: the annotations are equal)')
operating_mode.add_argument('-E', '--enrichment_count', dest = 'enrichment_count', type = EnrichmentCountType, default = 0,
                            help = "If 0: consider the input sites separately - symbol-resolved mode; if positive int: minimal number of sites with occurrence to consider a position enriched - enrichment mode; " +\
                            "if 'gross': count all the occurrencies - gross mode (default: 0)")
core_controls.add_argument('-Os', '--overlap_symbols', dest = 'overlap_symbols', type = int, default = 1, help = 'Minimal overlaping symbols of a site to be considered as match, [1,inf.)')
core_controls.add_argument('-Op', '--overlap_part', dest = 'overlap_part', type = float, default = 0.0, help = 'Minimal overlaping part of a site to be considered as match, [0,1]')
core_controls.add_argument('-Oa', '--overlap_apply', dest = 'overlap_apply', default = 'shortest', choices = ['shortest', 'longest', 'current', 'patched'],
                           help = "Apply the overlapping criteia to the {shortest} of two site, the {longest}, the {current} or the current allowing {patched} matches (default: 'shortest')")
core_controls.add_argument('-On', '--overlap_nature', dest = 'predictor_nature', default = 'neutral', choices = ['lagging', 'any', 'leading'],
                           help = "Required overlap nature (benchmark mode only): the predictor is {lagging} (predicted sites start not earlier than the benchmark ones), {any} or {leading} (default: 'any')")
core_controls.add_argument('-a', '--averaging', dest = 'averaging', default = 'group', choices = ['sequence', 'group', 'dataset'],
                           help = "Averaging of basic measures: {sequence}-wise (macro-macro), {group}-wise (micro-macro) or {dataset}-wise (micro-micro) (default: 'group')")
core_controls.add_argument('-A', '--adjust_for_seqlen', dest = 'len_adjust', action = 'store_true', help = "Adjust residue counts for the sequence length (default: plain-sum the counts)")
input_format.add_argument('-sd', '--seqlenfile_delim', dest = 'len_db_delimiter', type = str, default = '\t', help = 'Delimiter in the length database file (default: <tab>)')
input_format.add_argument('-sh', '--seqlenfile_headers', dest = 'len_db_headers', type = int, default = 0, help = 'Number of header lines to discard in the length database file')
input_format.add_argument('-sc', '--seqlenfile_colnumbers', dest = 'len_db_columns', type = str, default  = '',
                          help = "Column numbers (starting from 1) with SID and sequence length (for time series: SID and start and finish time points) in the length database file (default: '1,2' or '1,2,3')")
input_format.add_argument('-sq', '--seqlenfile_quotes', dest = 'len_db_quotes', action = 'store_true', help = "Read in quotes from the length database file (default: ignore quotes)")
input_format.add_argument('-md', '--mapfile_delim', dest = 'group_map_delimiter', type = str, default = '\t', help = 'Delimiter in the group map file (default: <tab>)')
input_format.add_argument('-mh', '--mapfile_headers', dest = 'group_map_headers', type = int, default = 0, help = 'Number of header lines to discard in the group map file')
input_format.add_argument('-mc', '--mapfile_colnumbers', dest = 'group_map_columns', type = str, default  = '1,2', help = "Column numbers (starting from 1) with SID and GID in the group map file (default: '1,2')")
input_format.add_argument('-mq', '--mapfile_quotes', dest = 'group_map_quotes', action = 'store_true', help = "Read in quotes from the group map file (default: ignore quotes)")
input_format.add_argument('-a1d', '--anno1file_delim', dest = 'anno1_delimiter', type = str, default = '\t', help = 'Delimiter in the first annotation file (default: <tab>)')
input_format.add_argument('-a1h', '--anno1file_headers', dest = 'anno1_headers', type = int, default = 0, help = 'Number of header lines to discard in the first annotation file')
input_format.add_argument('-a1c', '--anno1file_colnumbers', dest = 'anno1_columns', type = str, default  = '',
                          help = "Column numbers (starting from 1) with start position, end position, SID, GID, site name (skip those not provided) in the first annotation file (default: '3,4,1,2,5' skip-adjusted)")
input_format.add_argument('-a1q', '--anno1file_quotes', dest = 'anno1_quotes', action = 'store_true', default = False, help = "Read in quotes from the first annotation file (default: ignore quotes)")
input_format.add_argument('-a2d', '--anno2file_delim', dest = 'anno2_delimiter', type = str, default = '\t', help = 'Delimiter in the second annotation file (default: <tab>)')
input_format.add_argument('-a2h', '--anno2file_headers', dest = 'anno2_headers', type = int, default = 0, help = 'Number of header lines to discard in the second annotation file')
input_format.add_argument('-a2c', '--anno2file_colnumbers', dest = 'anno2_columns', type = str, default  = '',
                          help = "Column numbers (starting from 1) with start position, end position, SID, GID, site name (skip those not provided) in the second annotation file (default: '3,4,1,2,5' skip-adjusted)")
input_format.add_argument('-a2q', '--anno2file_quotes', dest = 'anno2_quotes', action = 'store_true', help = "Read in quotes from the second annotation file (default: ignore quotes)")
input_alternatives.add_argument('-l', '--seqlen_value', dest = 'seq_len', type = int, default = 0, help = 'Length of all the sequences; sequence length database must not be provided')
input_alternatives.add_argument('-ss', '--single_sequence', dest = 'single_sequence', action = 'store_true', help = "Process single sequence; SIDs must not be provided")
input_alternatives.add_argument('-ts', '--timeseries_start', dest = 'series_start', type = str, default = '', help = 'Start of all the time series')
input_alternatives.add_argument('-tf', '--timeseries_finish', dest = 'series_finish', type = str, default = '', help = 'Finish of all the time series')
input_alternatives.add_argument('-a1as', '--anno1file_all_sequences', dest = 'anno1_all_sequences', action = 'store_true', help = "Treat all sites in the first annotation as belonging to all the SIDs")
input_alternatives.add_argument('-a1ag', '--anno1file_all_groups', dest = 'anno1_all_groups', action = 'store_true', help = "Treat all sites in the first annotation as belonging to all the GIDs")
input_alternatives.add_argument('-a2as', '--anno2file_all_sequences', dest = 'anno2_all_sequences', action = 'store_true', help = "Treat all sites in the second annotation as belonging to all the SIDs")
input_alternatives.add_argument('-a2ag', '--anno2file_all_groups', dest = 'anno2_all_groups', action = 'store_true', help = "Treat all sites in the second annotation as belonging to all the GIDs")
input_alternatives.add_argument('-sg', '--sequences_as_groups', dest = 'sequences_as_groups', action = 'store_true', help = "Treat all the SIDs also as GIDs (i.e., form one-sequnce groups)")
input_alternatives.add_argument('-nOg', '--non_overlapping_groups', dest = 'non_overlapping_groups', action = 'store_true',
                                help = "The group mapping contains only non-overlapping groups (GIDs in the annotation files must not be provided)")
input_controls.add_argument('-n', '--site_names', dest = 'site_names', action = 'store_true', help = "Read in SCE names in addition")
input_controls.add_argument('-t', '--time_unit', dest = 'time_unit', default = 'none', choices = ['none', 'sec', 'min', 'hour', 'day'], help = "Time unit if the sequences are time series (default: 'none')")
input_controls.add_argument('-a1r', '--anno1file_resolve', dest = 'anno1_resolve_overlaps', default = 'all', choices = ['all', 'first', 'last', 'merge'],
                            help = "Resolve overlaps within the first annotation: leave {all} sites, only the {first} one, only the {last} one, or {merge} overlapping sites (default: 'all')")
input_controls.add_argument('-a1bs', '--anno1file_begin_shift', dest = 'anno1_begin_shift', type = int, default = 0, help = "Constant shift of site begin positions in the first annotation (symbols)")
input_controls.add_argument('-a1es', '--anno1file_end_shift', dest = 'anno1_end_shift', type = int, default = 0, help = "Constant shift of site end positions in the first annotation (symbols)")
input_controls.add_argument('-a2r', '--anno2file_resolve', dest = 'anno2_resolve_overlaps', default = 'all', choices = ['all', 'first', 'last', 'merge'],
                            help = "Resolve overlaps within the second annotation: leave {all} sites, only the {first} one, only the {last} one, or {merge} overlapping sites (default: 'all')")
input_controls.add_argument('-a2bs', '--anno2file_begin_shift', dest = 'anno2_begin_shift', type = int, default = 0, help = "Constant shift of site begin positions in the second annotation (symbols)")
input_controls.add_argument('-a2es', '--anno2file_end_shift', dest = 'anno2_end_shift', type = int, default = 0, help = "Constant shift of site end positions in the second annotation (symbols)")
input_controls.add_argument('-e', '--end_overflow_policy', dest = 'end_overflow_policy', default = 'forbid', choices = ['forbid', 'trim', 'ignore', 'circular'],
                            help = "Policy on overflowing site ends: {forbid}, {ignore} site, {trim} it or treat sequences as {circular} (default: 'forbid')")
input_controls.add_argument('-z', '--zero_for_na', dest = 'na_zeros', action = 'store_true', help = 'Treat NA values as zeros while calculating averages')
input_controls.add_argument('-min', '--min_group_size', dest = 'min_group_size', type = int, default = 1, help = 'Minimal size (number of sequences) of a group')
input_controls.add_argument('-max', '--max_group_size', dest = 'max_group_size', type = int, default = 0, help = 'Maximal size (number of sequences) of a group (0=infinity)')
input_controls.add_argument('-d', '--detect', dest = 'detect', default = 'none', choices = ['none', 'strand', 'frame'], help = "Detect {strand} or reading {frame} (GenBank of BED input only) (default; 'none')")
output_files_extra.add_argument('-od', '--outfile_detailed', dest = 'output_file_detailed', type = str, default = '', help = 'Detailed output file')
output_files_extra.add_argument('-os', '--outfile_sites', dest = 'output_file_site', type = str, default = '', help = 'Output TSV file with site-wise statistics')
output_files_extra.add_argument('-ou', '--outfile_union', dest = 'output_file_union', type = str, default = '', help = 'Output TSV file with the union of two annotations')
output_files_extra.add_argument('-oi', '--outfile_intersection', dest = 'output_file_intersection', type = str, default = '', help = 'Output TSV file with the intersection of two annotations')
output_files_extra.add_argument('-oc1', '--outfile_complement_1', dest = 'output_file_complement1', type = str, default = '', help = 'Output TSV file with the complement of the first annotation')
output_files_extra.add_argument('-oc2', '--outfile_complement_2', dest = 'output_file_complement2', type = str, default = '', help = 'Output TSV file with the complement of the second annotation')
output_files_extra.add_argument('-ore1', '--outfile_rel_enrichment_1', dest = 'output_file_re1', type = str, default = '', help = 'Output TSV file with the sites of relative enrichement in the first annotation')
output_files_extra.add_argument('-ore2', '--outfile_rel_enrichment_2', dest = 'output_file_re2', type = str, default = '', help = 'Output TSV file with the sites of relative enrichement in the second annotation')
output_controls.add_argument('-osd', '--outfile_sites_diff', dest = 'site_difference', default = 'all', choices = ['all', 'matched', 'unmatched', 'discrepant'],
                             help = "Limit site-wise statics to {matched}, {unmatched} or {discrepant} sites (default: 'all')")
output_controls.add_argument('-c', '--clean', dest = 'clean', action = 'store_true', help = 'Produce cleaned output TSV (without comments and averaged values)')
output_controls.add_argument('-sort', '--sort_output', dest = 'sort_output', action = 'store_true', help = 'Sort the main output table by GID')
output_controls.add_argument('-sum', '--calculate_sums', dest = 'calculate_sums', action = 'store_true', help = 'Calculate sums in addition to averages for counts')
other_options.add_argument('-preparse', '--preparse_mapfile', dest = 'preparse_group_map', action = 'store_true', help = 'Preparse the group mapping before parsing the sequence length database')
other_options.add_argument('-w', '--warning_level', dest = 'warnings', type = int, default = 1, help = 'Warnings level: 0 - no warnings, 1- standard')
other_options.add_argument('-q', '--quiet', dest = 'quiet', action = 'store_true', help = 'Quiet run: do not print progress')
arg_processor = ArgumentProcessor(arg_parser)
opt = arg_processor.prepare_input_options()
global_state = GlobalState(opt)

#Parsing input files
file_parser = CSVParser(opt, global_state)
input_file_processor = InputFileProcessor(opt, file_parser)
input_data = input_file_processor.process_input_files()

#Processing data
data_processor = DataProcessor(opt, global_state, input_data)
data_processor.process()

if not opt.quiet:
    print('Finished!')