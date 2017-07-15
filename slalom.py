import sys, argparse
from slalom_structures import GlobalState
from slalom_auxiliar import ArgumentProcessor, CSVParser, InputFileProcessor, DataProcessor

if __name__ != '__main__':
    sys.exit()

#Parsing input arguments
usage = '%(prog)s [options] -s SEQ_LEN_DB_FILE [-m GROUP_MAP_FILE] -a1 ANNO_1_FILE -a2 ANNO_2_FILE -o OUTPUT_FILE'
version = '%(prog)s version 1.0'
arg_parser = argparse.ArgumentParser(usage = usage)
arg_parser.add_argument('--version', action = 'version', version = version)
arg_parser.add_argument('-s', metavar = 'SEQ_LEN_DB_FILE', dest = 'len_db', default = '', help = 'Input file with the database of sequence lengths')
arg_parser.add_argument('-l', dest = 'seq_len', type = int, default = 0, help = 'Length of all sequences')
arg_parser.add_argument('-m', metavar = 'GROUP_MAP_FILE', dest = 'group_map', default = '', help = 'Input file with the sequence group mapping')
arg_parser.add_argument('-names', dest = 'site_names', action = 'store_true', default = False, help = "Read in site names in addition")
arg_parser.add_argument('-a1', metavar = 'ANNO_1_FILE', dest = 'anno1', required = True, help = 'Input file with the first annotation')
arg_parser.add_argument('-a2', metavar = 'ANNO_2_FILE', dest = 'anno2', required = True, help = 'Input file with the second annotation')
arg_parser.add_argument('-o', dest = 'output_file', required = True, help = 'Output TSV file with calculated performance measures')
arg_parser.add_argument('-od', dest = 'output_file_detailed', default = '', help = 'Detailed output file')
arg_parser.add_argument('-os', dest = 'output_file_site', default = '', help = 'Output TSV file with site-wise statistics')
arg_parser.add_argument('-os_diff', dest = 'site_difference', default = 'all', choices = ['all', 'discrepant', 'unmatched'], help = "Limit site-wise statics to (unmatched) or (discrepant) sites (default: 'all')")
arg_parser.add_argument('-ou', dest = 'output_file_union', default = '', help = 'Output TSV file with the union of two annotations')
arg_parser.add_argument('-oi', dest = 'output_file_intersection', default = '', help = 'Output TSV file with the intersection of two annotations')
arg_parser.add_argument('-oc1', dest = 'output_file_complement1', default = '', help = 'Output TSV file with the complement of the first annotations')
arg_parser.add_argument('-oc2', dest = 'output_file_complement2', default = '', help = 'Output TSV file with the complement of the second annotations')
arg_parser.add_argument('-ore1', dest = 'output_file_re1', default = '', help = 'Output TSV file with the sites of relative enriched in the first annotations')
arg_parser.add_argument('-ore2', dest = 'output_file_re2', default = '', help = 'Output TSV file with the sites of relative enriched in the second annotations')
arg_parser.add_argument('-ss', dest = 'single_sequence', action = 'store_true', default = False, help = "Process single sequence; do not ask for SIDs")
arg_parser.add_argument('-t', dest = 'time_unit', default = 'none', choices = ['none', 'sec', 'min', 'hour', 'day'], help = "Time unit if the sequences are time series (default: 'none')")
arg_parser.add_argument('-ts', dest = 'series_start', default = '', help = 'Start of all time series')
arg_parser.add_argument('-tf', dest = 'series_finish', default = '', help = 'Finish of all time series')
arg_parser.add_argument('-sd', dest = 'len_db_delimiter', default = '\t', help = 'Delimiter in the length database file (default: <tab>)')
arg_parser.add_argument('-sh', dest = 'len_db_headers', type = int, default = 0, help = 'Number of header lines to discard in the length database file')
arg_parser.add_argument('-sc', dest = 'len_db_columns', default  = '1,2',
                        help = "Column numbers (starting from 1) with SID and sequence length (for time series: SID and start and finish time points) in the length database file (default: '1,2')")
arg_parser.add_argument('-sq', dest = 'len_db_quotes', action = 'store_true', default = False, help = "Read in quotes from the length database file (default: ignore quotes)")
arg_parser.add_argument('-md', dest = 'group_map_delimiter', default = '\t', help = 'Delimiter in the group map file (default: <tab>)')
arg_parser.add_argument('-mh', dest = 'group_map_headers', type = int, default = 0, help = 'Number of header lines to discard in the group map file')
arg_parser.add_argument('-mc', dest = 'group_map_columns', default  = '1,2', help = "Column numbers (starting from 1) with SID and GID in the group map file (default: '1,2')")
arg_parser.add_argument('-mq', dest = 'group_map_quotes', action = 'store_true', default = False, help = "Read in quotes from the group map file (default: ignore quotes)")
arg_parser.add_argument('-a1d', dest = 'anno1_delimiter', default = '\t', help = 'Delimiter in the first annotation file (default: <tab>)')
arg_parser.add_argument('-a1h', dest = 'anno1_headers', type = int, default = 0, help = 'Number of header lines to discard in the first annotation file')
arg_parser.add_argument('-a1c', dest = 'anno1_columns', default  = '2,3,1',
                        help = "Column numbers (starting from 1) with begin position, end position, SID[, GID][, site name] in the first annotation file (default: '2,3,1')")
arg_parser.add_argument('-a1r', dest = 'anno1_resolve_overlaps', default = 'all', choices = ['all', 'first', 'last', 'merge'],
                        help = "Resolve overlaps within the first annotation: leave (all) sites, only the (first) one, only the (last) one, or (merge) overlapping sites (default: 'all')")
arg_parser.add_argument('-a1q', dest = 'anno1_quotes', action = 'store_true', default = False, help = "Read in quotes from the first annotation file (default: ignore quotes)")
arg_parser.add_argument('-a1bs', dest = 'anno1_begin_shift', type = int, default = 0, help = "Constant shift of site begin positions in the first annotation (symbols)")
arg_parser.add_argument('-a1es', dest = 'anno1_end_shift', type = int, default = 0, help = "Constant shift of site end positions in the first annotation (symbols)")
arg_parser.add_argument('-a1as', dest = 'anno1_all_sequences', action = 'store_true', default = False, help = "Treat all sites in the first annotation as belonging to all the SIDs")
arg_parser.add_argument('-a1ag', dest = 'anno1_all_groups', action = 'store_true', default = False, help = "Treat all sites in the first annotation as belonging to all the GIDs")
arg_parser.add_argument('-a2d', dest = 'anno2_delimiter', default = '\t', help = 'Delimiter in the second annotation file (default: <tab>)')
arg_parser.add_argument('-a2h', dest = 'anno2_headers', type = int, default = 0, help = 'Number of header lines to discard in the second annotation file')
arg_parser.add_argument('-a2c', dest = 'anno2_columns', default  = '2,3,1',
                        help = "Column numbers (starting from 1) with begin position, end position, SID[, GID][, site name in the second annotation file (default: '2,3,1')")
arg_parser.add_argument('-a2r', dest = 'anno2_resolve_overlaps', default = 'all', choices = ['all', 'first', 'last', 'merge'],
                        help = "Resolve overlaps within the second annotation: leave (all) sites, only the (first) one, only the (last) one, or (merge) overlapping sites (default: 'all')")
arg_parser.add_argument('-a2q', dest = 'anno2_quotes', action = 'store_true', default = False, help = "Read in quotes from the second annotation file (default: ignore quotes)")
arg_parser.add_argument('-a2bs', dest = 'anno2_begin_shift', type = int, default = 0, help = "Constant shift of site begin positions in the second annotation (symbols)")
arg_parser.add_argument('-a2es', dest = 'anno2_end_shift', type = int, default = 0, help = "Constant shift of site end positions in the second annotation (symbols)")
arg_parser.add_argument('-a2as', dest = 'anno2_all_sequences', action = 'store_true', default = False, help = "Treat all sites in the second annotation as belonging to all the SIDs")
arg_parser.add_argument('-a2ag', dest = 'anno2_all_groups', action = 'store_true', default = False, help = "Treat all sites in the second annotation as belonging to all the GIDs")
arg_parser.add_argument('-e', dest = 'end_overflow_policy', default = 'error', choices = ['error', 'trim', 'ignore'], help = "Policy on overflowing site ends: raise (error), (ignore) site or (trim) it (default: 'error')")
arg_parser.add_argument('-b', dest = 'benchmark', action = 'store_true', default = False, help = 'Treat the first annotation as benchmark (default: the annotations are equal)')
arg_parser.add_argument('-nature', dest = 'predictor_nature', default = 'neutral', choices = ['lagging', 'neutral', 'leading'], help = "Predictor nature: (lagging), (neutral) or (leading) (default: 'neutral')")
arg_parser.add_argument('-ovs', dest = 'overlap_symbols', type = int, default = 1, help = 'Minimal overlaping symbols of a site to be considered as match, [1,inf.)')
arg_parser.add_argument('-ovp', dest = 'overlap_part', type = float, default = 0.0, help = 'Minimal overlaping part of a site to be considered as match, [0,1]')
arg_parser.add_argument('-ova', dest = 'overlap_apply', default = 'shortest', choices = ['shortest', 'longest', 'current', 'patched'],
                        help = "Apply the overlapping criteia to the (shortest), of two site, the (longest), the (current) or the curent in (patched) mode (default: 'shortest')")
arg_parser.add_argument('-n', dest = 'enrichment_count', type = int, default = 0, help = 'Minimal count to consider a symbol enriched')
arg_parser.add_argument('-gross', dest = 'gross', action = 'store_true', default = False, help = 'Count symbols gross in case of overlaps (applicable only in Boolean mode)')
arg_parser.add_argument('-g', dest = 'groupwise', action = 'store_true', default = False, help = 'Average basic measures group-wise (default: sequence-wise)')
arg_parser.add_argument('-z', dest = 'na_zeros', action = 'store_true', default = False, help = 'Treat NA values as zeros while calculating averages')
arg_parser.add_argument('-minsize', dest = 'min_group_size', type = int, default = 1, help = 'Minimal size (number of sequences) of a group')
arg_parser.add_argument('-maxsize', dest = 'max_group_size', type = int, default = 0, help = 'Maximal size (number of sequences) of a group (0=infinity)')
arg_parser.add_argument('-c', dest = 'clean', action = 'store_true', default = False, help = 'Produce cleaned output TSV (without comments and averaged values)')
arg_parser.add_argument('-preparse', dest = 'preparse_group_map', action = 'store_true', default = False, help = 'Preparse the group mapping before parsing the sequence length database')
arg_parser.add_argument('-w', dest = 'warnings', type = int, default = 1, help = 'Warnings level: 0 - no warnings, 1- standard')
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

print('Finished!')