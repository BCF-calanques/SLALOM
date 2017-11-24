"""SLALOM (StatisticaL Analysis of Locus Overlap Method)
Copyright (C) 2017  Roman Prytuliak

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see https://www.gnu.org/licenses/."""

import os, sys, re, math, datetime, time, copy, linecache, argparse
import numpy as np
from operator import itemgetter
from slalom_structures import DefaultOrderedDict, InputData, CurrentSequence, BasicBooleanMeasures, BasicEnrichmentMeasures, PerformanceMeasures, FileHandlers, EnrichmentCountType

def error(message):
    """Function for error reporting"""
    sys.stderr.write('Error: {}\n'.format(message))
    sys.stderr.flush()
    sys.exit(1)

class CustomHelpFormatter(argparse.HelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings:
            metavar, = self._metavar_formatter(action, action.dest)(1)
            return metavar
        else:
            parts = []
            if action.nargs == 0:
                parts.extend(action.option_strings)
            else:
                default = action.dest.upper()
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append(option_string)
                if args_string.startswith('{'):
                    ending = 'enum' + args_string
                elif action.type is None:
                    ending = ''
                elif action.type == EnrichmentCountType:
                    ending = "int or 'gross'"
                elif action.type == str:
                    ending = 'string'
                else:
                    ending = action.type.__name__
                parts[-1] += '    ' + ending
            return ', '.join(parts)

class GenBankMethods:
    """Class to keep methods for working with GenBank input"""
    def get_seq_len(fname):
        """Method to get the sequence length from a specified GenBank file"""
        try:
            return int(linecache.getline(fname, 1).split()[2])
        except:
            error('Cannor read the sequence length from the GenBank file "{}"'.format(fname))
    def gen_record(file_obj, detect_strand, detect_frame, seq_length = 0):
        """Generator for getting the next CDS record from a GenBank file preparsed as TSV"""
        expect_name = False
        line_idx = -1
        for line in file_obj:
            line_idx += 1
            if line.startswith('FEATURES'):
                break
        for line in file_obj:
            line_idx += 1
            try:
                if expect_name:
                    site_name = line.split('"')[1]
                    yield line_idx, '{}\t{}\t{}\t{}'.format(begin, end, frame, site_name)
                    expect_name = False
                    continue
                fields = line.split()
                if fields and (fields[0] == 'CDS'):
                    if fields[1].startswith('c'):
                        fields[1] = fields[1].split('(')[1][: -1]
                        strand_forward = False
                    else:
                        strand_forward = True
                    begin, end = fields[1].split('..')
                    begin = begin.strip('<')
                    end = end.strip('>')
                    if detect_frame:
                        reminder = (int(begin) if strand_forward else (seq_length - int(end) + 1)) % 3
                        if reminder == 0:
                            reminder = 3
                    frame = ('{:+d}'.format(reminder * (1 if strand_forward else -1)) if detect_frame else ('+' if strand_forward else '-')) if detect_strand else '*'
                    expect_name = True
            except (IndexError, ValueError):
                error('Error while parsing the line {} of the file "{}". The format is not GenBank'.format(line_idx + 1, file_obj.name))

class BEDMethods:
    """Class to keep methods for working with BED input"""
    frame_regex = re.compile('[-+][123]$')
    def files_have_names(anno1_fname, anno2_fname2):
        """Method to detect if provided input BED files contain site names"""
        for fname in (anno1_fname, anno2_fname2):
            if len(linecache.getline(fname, 1).split('\t')) < 4:
                return False
        return True
    def save_seq_len_record(seq_len_db, SID, seq_length, detect_strand, detect_frame):
        """Method to expand SIDs if strand or frame detection is required"""
        if not detect_strand:
            seq_len_db[SID] = seq_length
            return
        if SID.endswith('+') or SID.endswith('-'):
            error('If strand detection in BED files is enabled, chromosome names ending with plus or minus signs are not allowed')
        for sign in ('+', '-'):
            if detect_frame:
                if BEDMethods.frame_regex.search(SID):
                    error('If strand detection in BED files is enabled, chromosome names ending with "{}" are not allowed'.format(SID[-2: ]))
                for frame in range(1, 4):
                    seq_len_db[SID + sign + frame] = seq_length
            else:
                seq_len_db[SID + sign] = seq_length
    def gen_record(file_obj, detect_strand, detect_frame, site_names, seq_len_map = None):
        """Generator for getting the next record from a BED file converted according to the user request and file structure"""
        for line_idx, line in enumerate(file_obj):
            line = line.strip('\r\n')
            if not line:
                continue
            fields = line.split('\t')
            try:
                SID = fields[0]
                if detect_strand:
                    strand = fields[5]
                    if strand not in ('+', '-'):
                        error('Error while parsing the line {} of the file "{}". The strand must be eithe \'+\' or \'-\''.format(line_idx + 1, file_obj.name))
                    if detect_frame:
                        if strand == '+':
                            begin_idx = int(fields[1])
                        else:
                            end_res = int(fields[2])
                            seq_length = seq_len_map[SID + '+1']
                        reminder = ((int(begin_idx) if strand_forward else (seq_length - int(end_res))) + 1) % 3
                        if reminder == 0:
                            reminder = 3
                    SID += '{:+d}'.format(reminder * (1 if strand == '+' else -1)) if detect_frame else strand
                yield line_idx, '{}\t{}\t{}{}'.format(SID, fields[1], fields[2], ('\t' + fields[3]) if site_names else '')
            except IndexError:
                errot('Error while parsing the line {} of the file "{}". Not enough columns'.format(line_idx + 1, file_obj.name))

class ArgumentValidator:
    """Class that contain means to command line argument validation"""
    prefixes = {'s': 'len_db', 'm': 'group_map', 'a1': 'anno1', 'a2': 'anno2'}
    suffixes = {'d': 'delimiter', 'h': 'headers', 'c': 'columns', 'q': 'quotes', 'bs': 'begin_shift', 'es': 'end_shift'}
    misc_keys = {'-Os': 'overlap_symbols', '-Op': 'overlap_part', '-max': 'max_group_size', '-w': 'warnings', '-l': 'seq_len'}
    nonnegative_int_regex = re.compile('^\+?\d+$')
    def __init__(self, opt):
        self.opt = opt
        self.file_control_regex = '-({})({})'.format('|'.join(self.prefixes), '|'.join(self.suffixes))
    def _get_file_control_option_value(self, key):
        """Method to retrieve values of file control command line options by their keys"""
        regex_search = re.search(self.file_control_regex, key)
        if not regex_search:
            return None
        dest = '{}_{}'.format(self.prefixes[regex_search.group(1)], self.suffixes[regex_search.group(2)])
        return getattr(self.opt, dest)
    def _set_file_control_option_value(self, key, value):
        """Method to set values of file control command line options by their keys"""
        regex_search = re.search(self.file_control_regex, key)
        if not regex_search:
            error('The key "{}" does not exist'.format(key))
        dest = '{}_{}'.format(self.prefixes[regex_search.group(1)], self.suffixes[regex_search.group(2)])
        setattr(self.opt, dest, value)
    def preliminary_validate(self):
        """Method for ensuring the default values of parameters to be set internally"""
        if (self.opt.genbank or self.opt.bed) and (self.opt.anno1_columns + self.opt.anno2_columns != ''):
            error('Column numbers in the annotation files cannot be changed in simplified (GenBank or BED) modes')
        if self.opt.genbank and ((getattr(self.opt, self.misc_keys['-l']) > 0) or (self.opt.len_db)):
            error('Sequence lenght must not be provided in the simplified GenBank mode. This information will be read from the annotation files')
    def preliminary_set_the_internal_parameters(self):
        """Method for setting the iternal options relevant for the simplified modes on the basis of user input"""
        self.opt.detect_strand = True if self.opt.detect != 'none' else False
        self.opt.detect_frame = True if self.opt.detect == 'frame' else False
    def set_the_simplified_mode(self):
        """Method for settting the relevant default options for a simplified mode"""
        if self.opt.genbank:
            seqlen = GenBankMethods.get_seq_len(self.opt.anno1)
            if seqlen != GenBankMethods.get_seq_len(self.opt.anno1):
                error('Sequence lenght between the two GenBank files does not match')
            setattr(self.opt, self.misc_keys['-l'], seqlen)
            self.opt.group_map = ''
            self._set_file_control_option_value('-a1c', '1,2,3,4')
            self._set_file_control_option_value('-a2c', '1,2,3,4')
            self.opt.site_names = True
            if self.opt.detect_strand:
                self.opt.sort_output = True
        elif self.opt.bed:
            self.opt.site_names = BEDMethods.files_have_names(self.opt.anno1, self.opt.anno2)
            self._set_file_control_option_value('-a1bs', self._get_file_control_option_value('-a1bs') + 1)
            self._set_file_control_option_value('-a2bs', self._get_file_control_option_value('-a2bs') + 1)
    def set_the_internal_parameters(self):
        """Method for setting the iternal options on the basis of user input"""
        self.opt.gross = False
        self.opt.enrichment_count = str(self.opt.enrichment_count)
        if self.opt.enrichment_count == 'gross':
            self.opt.gross = True
            self.opt.enrichment_count = 0
        elif not ArgumentValidator.nonnegative_int_regex.search(self.opt.enrichment_count):
            error("invalid value for the enrichment count. It must be a non-negative integer of 'gross'")
        else:
            self.opt.enrichment_count = int(self.opt.enrichment_count)
        if self.opt.detect_strand:
            self.opt.sequences_as_groups = True
        self.opt.grouped = self.opt.group_map or self.opt.sequences_as_groups
        self.opt.circular = True if self.opt.end_overflow_policy == 'circular' else False
    def validate_file_paths(self):
        """Method to check for validity of given input and output file paths"""
        if (self.opt.seq_len == 0) and (not self.opt.series_start) and (not os.path.isfile(self.opt.len_db)):
            if not self.opt.len_db:
                error('Neither input file with the sequence length table nor the common sequence length is not provided')
            error('The sequence length table file does not exist')
        if self.opt.group_map:
            if not os.path.isfile(self.opt.group_map):
                error('Input file with the sequence group mapping does not exist')
        if not os.path.isfile(self.opt.anno1):
            error('Input file with the first annotation does not exist')
        if not os.path.isfile(self.opt.anno2):
            error('Input file with the second annotation does not exist')
    def validate_file_column_numbers(self):
        """Method to check for validity the options listing column numbers in the files"""
        for key in ('-sc', '-mc', '-a1c', '-a2c'):
            value = self._get_file_control_option_value(key)
            default_ = True if not value else False
            if key == '-mc':
                if not self.opt.group_map:
                    continue
                n = 2
            elif key == '-sc':
                if not self.opt.len_db:
                    continue
                if self.opt.time_unit == 'none':
                    n = 2
                    if default_:
                        value = '1,2'
                else:
                    n = 3
                    if default_:
                        value = '1,2,3'
            else:
                if self.opt.single_sequence or getattr(self.opt, 'anno' + key[2] + '_all_sequences'):
                    n = 2
                    if default_:
                        value = '1,2'
                elif (not self.opt.group_map) or self.opt.non_overlapping_groups or getattr(self.opt, 'anno' + key[2] + '_all_groups'):
                    n = 3
                    if default_:
                        value = '2,3,1'
                else:
                    n = 4
                    if default_:
                        value = '3,4,1,2'
                if self.opt.site_names:
                    n += 1
                    if default_:
                        value += ',{}'.format(n)
            regex = '^([1-9]\d*,){{{}}}[1-9]\d*$'.format(n - 1)
            if not re.search(regex, value):
                if n > 1:
                    error("Invalid format for the option '{}'. Expected a list of {} comma-delimited positive integers".format(key, n))
                else:
                    error("Invalid format for the option '{}'. Expected a positive integer".format(key))
            if default_:
                self._set_file_control_option_value(key, value)
    def validate_delimiters(self):
        """Method to check for validity the delimiters for the input files"""
        regex = re.compile('''^[ \t,;.:/]?$''')
        for key in ('-sd', '-md', '-a1d', '-a2d'):
            value = self._get_file_control_option_value(key)
            if not regex.match(value):
                error("Invalid value for the option '{}'. Expected a character from the set ' \t,;.:/' or empty string".format(key))
    def validate_numerical_options_boundaries(self):
        """Method to check if numerical option values lie in correct boundaries"""
        for key in ('-sh', '-mh', '-a1h', '-a2h'):
            if self._get_file_control_option_value(key) < 0:
                error("Invalid value for the option '{}'. Expected a non-negative integer".format(key))
        for key in ('-l', '-max'):
            if getattr(self.opt, self.misc_keys[key]) < 0:
                error("Invalid value for the option '{}'. Expected a non-negative integer".format(key))
        for key in ('-Os', ):
            if getattr(self.opt, self.misc_keys[key]) < 1:
                error("Invalid value for the option '{}'. Expected a positive integer".format(key))
        for key in ('-Op', ):
            value = getattr(self.opt, self.misc_keys[key])
            if (value < 0.0) or (value > 1.0):
                error("Invalid value for the option '{}'. Expected a value in range [0,1]".format(key))
        for key in ('-w', ):
            value = getattr(self.opt, self.misc_keys[key])
            if (value < 0) or (value > 1):
                error("Invalid value for the option '{}'. Expected an integer in range [0,1]".format(key))
    def validate_logic(self):
        """Method to validate the logic of the interplay of different parameters"""
        present = False
        for i in ('1', '2'):
            if getattr(self.opt, 'output_file_re' + i):
                present = True
        if present and (not getattr(self.opt, 'enrichment_count')):
            error('Relative enrichment output files can be written only in enrichment mode')
        if self.opt.len_db and (self.opt.seq_len > 0):
            error('Positive lengrh of all sequence and the file with the sequence length table cannot be provided simultaneously')
        if (self.opt.len_db == 0) and self.opt.single_sequence:
            error('If single sequence is analyzed, its length must be provided as the corresponding argument')
        if ((self.opt.anno1_resolve_overlaps != 'all') or (self.opt.anno2_resolve_overlaps != 'all')) and getattr(self.opt, 'enrichment_count'):
            error('Overlap resolving within a given annotation is not possible only in enrichment mode')
        if self.opt.preparse_group_map and (not self.opt.group_map):
            error('To preprse the group mapping it must be provided')
        if self.opt.non_overlapping_groups and (not self.opt.group_map):
            error('Non-overlapping groups may be specified only if the group mapping file is provided')
        if self.opt.sequences_as_groups and self.opt.group_map:
            error('If the group mapping is provided, treating sequences as groups is not possible')
        if self.opt.single_sequence and self.opt.sequences_as_groups:
            error('If SIDs are not provided, treating sequences as groups is not possible')
        if self.opt.seq_len and (self.opt.time_unit != 'none'):
            error('Positive length of all sequences is not compatible with time series')
        if (self.opt.series_start or self.opt.series_finish) and (self.opt.time_unit == 'none'):
            error('Start and finish of all sequences can be specified only for time series')
        if bool(self.opt.series_start) != bool(self.opt.series_finish):
            error('Start and finish of all sequences must be provided together')
        if (not self.opt.len_db) and self.opt.group_map and (not self.opt.preparse_group_map):
            if self.opt.warnings:
                print('The sequence length table file is not provided. The preparsing of the group mapping is activated automatically')
        if self.opt.single_sequence and self.opt.group_map:
            error('Group mapping is not compatible with the option of processing single sequence')
        if self.opt.sequences_as_groups and (self.opt.min_group_size > 1):
            error('The option of treating sequences as groups is not compatible with the minimal group size greater than 1')
        if self.opt.max_group_size and (self.opt.min_group_size > self.opt.max_group_size):
            error('Minimal group size cannot not be greater as the maximal size')
        if self.opt.anno1_all_sequences and self.opt.anno2_all_sequences:
            error('The option of treating all sites as belonging to all sequences cannot be activated for both annotations')
        if (self.opt.anno1_all_sequences and self.opt.anno1_all_groups) or (self.opt.anno2_all_sequences and self.opt.anno2_all_groups):
            error('The options of treating all sites as belonging to all sequences and all groups are not compatible for the same annotation')
        if (not self.opt.len_db) and (not self.opt.group_map) and (self.opt.anno1_all_sequences or self.opt.anno2_all_sequences):
            error('For the options of treating all sites as belonging to all sequences either the sequence length table or the group mapping is required')
        if (self.opt.anno1_all_groups or self.opt.anno2_all_groups) and self.opt.sequences_as_groups:
            error('The option of treating of sequences as groups is not compatible with an option of treating all sites as belonging to all groups')
        if (self.opt.anno1_all_groups or self.opt.anno2_all_groups) and (not self.opt.group_map):
            error('The options of treating all sites as belonging to all groups can be activated only if the group mapping is provided')
        if (not self.opt.benchmark) and (self.opt.predictor_nature != 'neutral'):
            error('The predictor nature can be only neutral if neither annotation is set as benchmark')
        if (self.opt.averaging == 'sequence') and self.opt.len_adjust:
            error('Sequence length adjustment is not applicable under sequence-wise averaging')
        if ((self.opt.anno1_resolve_overlaps == 'merge') or (self.opt.anno2_resolve_overlaps == 'merge')) and self.opt.site_names:
            error('Named sites are not compatible with merging while resolving overlaps')
        if (self.opt.overlap_apply == 'patched') and (self.opt.site_difference == 'discrepant'):
            error('Showing only discreoant in the site-wise output is not compatible with the patched overlap logic')
        if self.opt.sort_output and (not self.opt.grouped) and (self.opt.detect == 'none'):
            if self.opt.warnings:
                print('As no sequence groups are defined, output sorting is meaningless')
        if self.opt.calculate_sums and (not self.opt.grouped):
            error('Calculating sums is not applicable to non-grouped data')
        if (self.opt.detect != 'none') and (not self.opt.genbank) and (not self.opt.bed):
            error('Strand/frame detection is supported only in a simplified GenBank or BED mode')
        if self.opt.group_map and (self.opt.genbank or self.opt.bed):
            error('Group mapping files are not supported in simplified (GenBank or BED) modes')
        if (self.opt.genbank or self.opt.bed) and (self.opt.anno1_delimiter + self.opt.anno2_delimiter != '\t\t'):
            error('Delimiter in the annotation files cannot be changed in simplified (GenBank or BED) modes')
        if (self.opt.genbank or self.opt.bed) and (self.opt.anno1_headers + self.opt.anno2_headers != 0):
            error('Number of header rows in the annotation files cannot be changed in simplified (GenBank or BED) modes')
        if self.opt.circular and (self.opt.predictor_nature != 'neutral'):
            error('Lagging or leading predictor nature is not compatible with circular sequences')
        if self.opt.circular and (self.opt.time_unit != 'none'):
            error('Time series cannot be circular')
        if self.opt.circular and ((self.opt.anno1_resolve_overlaps in ('first', 'last')) or (self.opt.anno2_resolve_overlaps in ('first', 'last'))):
            error("For circular sequences, only the following choices for overlap resolving are supported: 'all', 'merge'")

class ArgumentProcessor:
    """Class to coordinate command line argument parsing"""
    def __init__(self, arg_parser):
        self.arg_parser = arg_parser
    def _check_if_empty(self):
        """Method to check if the argument list is empty"""
        if len(sys.argv) == 1:
            self.arg_parser.print_help()
            sys.exit()
    def prepare_input_options(self):
        """Method to genarae validated object with the command line options"""
        self._check_if_empty()
        opt = self.arg_parser.parse_args()
        validator = ArgumentValidator(opt)
        validator.preliminary_set_the_internal_parameters()
        validator.set_the_simplified_mode()
        validator.set_the_internal_parameters()
        validator.validate_file_paths()
        validator.validate_file_column_numbers()
        validator.validate_delimiters()
        validator.validate_numerical_options_boundaries()
        validator.validate_logic()
        return opt
        
class CSVParser:
    """Class to parse the input CSV files"""
    field_regex_quoted = '''((?:[^{0}"']|"[^"]*(?:"|$)|'[^']*(?:'|$))+|(?={0}{0})|(?={0}$)|(?=^{0}))'''
    field_regex_simple = '(?:[^{0}]+|(?={0}{0})|(?={0}$)|(?=^{0}))'
    quote_compiled = re.compile('''['"]''')
    int_regex = re.compile('^[+-]?\d+$')
    pos_int_regex = re.compile('^\+?[1-9]\d*$')
    time_formats_compiled = [re.compile(x) for x in ('^\d\d/\d\d/\d{4} \d\d:\d\d:\d\d$', '^\d\d/\d\d/\d{4} \d\d:\d\d$', '^\d\d\.\d\d\.\d{4} \d\d:\d\d:\d\d$', '^\d\d\.\d\d\.\d{4} \d\d:\d\d$')]
    time_formats = ['%m/%d/%Y %H:%M:%S', '%m/%d/%Y %H:%M', '%d.%m.%Y %H:%M:%S', '%d.%m.%Y %H:%M']
    def __init__(self, opt, global_state):
        self.opt = opt
        self.global_state = global_state
        self.input_data = InputData()
    def _parse_input_file(self, opt_prefix, preliminary = False):
        """Method to parse an input file"""
        column_indices = tuple(int(x) - 1 for x in getattr(self.opt, opt_prefix + '_columns').split(','))
        filename = getattr(self.opt, opt_prefix)
        delimiter = getattr(self.opt, opt_prefix + '_delimiter')
        collapse_spaces = True if not delimiter else False
        delimiter = delimiter if delimiter else ' '
        quotes_as_escaped = getattr(self.opt, opt_prefix + '_quotes')
        file_field = re.compile(getattr(CSVParser, 'field_regex_simple' if quotes_as_escaped else 'field_regex_quoted').format(delimiter))
        with open(filename, 'r') as ifile:
            for i in range(getattr(self.opt, opt_prefix + '_headers')):
                next(ifile)
            line_generator = enumerate(ifile)
            if opt_prefix.startswith('a'):
                if self.opt.genbank:
                    line_generator = GenBankMethods.gen_record(ifile, self.opt.detect_strand, self.opt.detect_frame, self.auto_seq_len)
                elif self.opt.bed:
                    line_generator = BEDMethods.gen_record(ifile, self.opt.detect_strand, self.opt.detect_frame, self.opt.site_names, self.input_data.seq_len)
            for line_idx, line in line_generator:
                if collapse_spaces:
                    line = re.sub(' +', ' ', line).strip()
                try:
                    values = itemgetter(*column_indices)(file_field.findall(line.strip('\n')))
                except IndexError:
                    error('Error while parsing the line {} of the file "{}". Not enough columns delimited by "{}" identified'.format(line_idx + 1, filename, delimiter))
                if not quotes_as_escaped:
                    values = [CSVParser.quote_compiled.sub('', el) for el in values]
                try:
                    self._save_record(opt_prefix, values, preliminary)
                except RuntimeError as e:
                    error('Error while parsing the line {} of the file "{}". {}'.format(line_idx + 1, filename, str(e)))
    def _duration_in_units(self, start_time_point, finish_time_point):
        """Method to calculate the distance in time, measured in speciied by the user units, between two points"""
        return math.floor((time.mktime(finish_time_point) - time.mktime(start_time_point)) / self.global_state.time_unit_seconds)
    def _convert_interval_to_time_structs(self, interval):
        """Mathod to convert time strings in an interval to time_struct objects"""
        for interval_idx, time_str in enumerate(interval):
            if ' ' not in time_str:
                time_str += ' 00:00'
            recognized = False
            for time_format_idx, time_format_compiled in enumerate(CSVParser.time_formats_compiled):
                if time_format_compiled.search(time_str):
                    recognized = True
                    try:
                        interval[interval_idx] = time.strptime(time_str, CSVParser.time_formats[time_format_idx])
                    except ValueError:
                        raise RuntimeError('Time format was not recognized. Supported formats: "mm/dd/yyyy HH:MM[:SS]" and "dd.mm.yyyy HH:MM[:SS]"') from None
                    break
            if not recognized:
                raise RuntimeError('Time format was not recognized. Supported formats: "mm/dd/yyyy HH:MM[:SS]" and "dd.mm.yyyy HH:MM[:SS]"')
    def _save_seq_len_db_record(self, values, not_first_to_check):
        """Method to save a sequence length table record"""
        def _save_record():
            """Closure to write the values to the corresponding dictionaries"""
            nonlocal seq_length
            if self.opt.bed:
                BEDMethods.save_seq_len_record(self.input_data.seq_len, SID, seq_length, self.opt.detect_strand, self.opt.detect_frame)
            else:
                self.input_data.seq_len[SID] = seq_length
            if self.global_state.time_unit_seconds:
                nonlocal interval
                self.input_data.time_series_starts[SID] = interval[0]
        if not self.global_state.time_unit_seconds:
            if self.opt.single_sequence:
                seq_length = values[0]
                SID = ''
            else:
                SID, seq_length = values
            if not self.pos_int_regex.search(seq_length):
                raise RuntimeError('Sequence length must be a positive integer')
        else:
            if self.opt.single_sequence:
                start, finish = values
                SID = ''
            else:
                SID, start, finish = values
            interval = [start, finish]
            self._convert_interval_to_time_structs(interval)
            seq_length = self._duration_in_units(interval[0], interval[1])
            if seq_length < 1:
                raise RuntimeError('The time interval must contain at least 1 time unit')
        seq_length = int(seq_length)
        if not_first_to_check:
            if SID not in self.input_data.seq_len.keys():
                return
            if self.input_data.seq_len[SID] is None:
                _save_record()
        else:
            _save_record()
        if (self.opt.single_sequence) and (not SID):
            raise RuntimeError('An SID cannot be empty')
        if '"' in SID:
            raise RuntimeError('An SID cannot contain double quotes')
        if self.input_data.seq_len[SID + (('+1' if self.opt.detect_frame else '+') if self.opt.detect_strand else '')] != seq_length:
            raise RuntimeError('Inconsistency in the sequence length table. Different length/duration values for a duplicating SID "{}".'.format(SID))
        if self.global_state.time_unit_seconds and (self.input_data.time_series_starts[SID] != interval[0]):
            raise RuntimeError('Inconsistency in the sequence length table. Different start values for a duplicating SID "{}".'.format(SID))
    def _save_group_map_record(self, values, preliminary = False):
        """Method to save a group mapping record"""
        SID, GID = values
        if not GID:
            raise RuntimeError('A GID cannot be empty')
        if '"' in GID:
            raise RuntimeError('A GID cannot contain double quotes')
        if preliminary:
            self.input_data.seq_len[SID] = self.auto_seq_len
            if self.global_state.time_unit_seconds:
                self.input_data.time_series_starts[SID] = self.auto_series_start
        elif SID not in self.input_data.seq_len.keys():
            if self.opt.warnings:
                print('Warning: SID "{}" is not the sequence length table. The group mapping record is ignored'.format(SID))
            return
        if SID not in self.input_data.group_map[GID]:
            self.input_data.group_map[GID].append(SID)
            if self.opt.non_overlapping_groups:
                if SID in self.reverse_group_map.keys():
                    raise RuntimeError('SID "{}" belongs to multiple groups in the group mapping, despite the non-overlapping groups were claimed')
                self.reverse_group_map[SID] = GID
    def _save_annotation_record(self, opt_prefix, values):
        """Method to save an annotation record"""
        if self.opt.site_names:
            site_name = values[-1]
            if not site_name:
                raise RuntimeError('A site name cannot be empty')
            if '"' in site_name:
                raise RuntimeError('A site name cannot contain double quotes')
            values = values[: -1]
        if self.opt.single_sequence or getattr(self.opt, opt_prefix + '_all_sequences'):
            begin, end = values
            GID = ''
            SID = ''
        elif self.opt.non_overlapping_groups:
            begin, end, SID = values
            GID = self.reverse_group_map[SID]
        elif self.opt.group_map and (not getattr(self.opt, opt_prefix + '_all_groups')):
            begin, end, SID, GID = values
        else:
            begin, end, SID = values
            GID = ''
        if self.opt.sequences_as_groups:
            GID = SID
        GID_list = [GID] if GID else list(self.input_data.group_map.keys())
        for GID_ in GID_list:
            SID_list = [SID] if SID else self.input_data.group_map[GID_]
            for SID_ in SID_list:
                if SID_ not in self.input_data.seq_len.keys():
                    if (not self.opt.len_db) and (not self.opt.group_map):
                        self.input_data.group_map[GID_].append(SID_)
                        self.input_data.seq_len[SID_] = self.auto_seq_len
                        if self.global_state.time_unit_seconds:
                            self.input_data.time_series_starts[SID_] = self.auto_series_start
                    elif self.opt.warnings:
                        print('Warning: SID "{}" is not in the sequence length table. The annotation record is ignored'.format(SID_))
                        return
                if self.opt.sequences_as_groups and self.opt.len_db:
                    self.input_data.group_map[GID_] = [SID_]
                elif GID_ not in self.input_data.group_map.keys():
                    if self.opt.warnings:
                        print('Warning: GID "{}" is not in the group mapping. The annotation record is ignored'.format(GID_))
                    return
                if SID_ not in self.input_data.group_map[GID_]:
                    if getattr(self.opt, opt_prefix + '_all_groups'):
                        continue
                    if self.opt.warnings:
                        print('Warning: SID "{}" does not belong to the group "{}" in the group mapping. The annotation record is ignored'.format(SID_, GID_))
                    return
                if not self.global_state.time_unit_seconds:
                    if (not self.int_regex.search(begin)) or (not self.int_regex.search(end)):
                        raise RuntimeError('Site begin and end position must be integers')
                    begin_ = int(begin)
                    end_ = int(end)
                else:
                    interval = [begin, end]
                    self._convert_interval_to_time_structs(interval)
                    begin_ = self._duration_in_units(self.input_data.time_series_starts[SID_], interval[0]) + 1
                    end_ = self._duration_in_units(self.input_data.time_series_starts[SID_], interval[1])
                begin_ += getattr(self.opt, opt_prefix + '_begin_shift')
                end_ += getattr(self.opt, opt_prefix + '_end_shift')
                if self.opt.circular:
                    if end_ - begin_ >= self.input_data.seq_len[SID_]:
                        raise RuntimeError('Site length cannot exceed the sequence length')
                    if begin_ > end_:
                        raise RuntimeError('Site begin position cannot exceed the end position. For circular sequences, end positions exceeding the sequence length should be used')
                    if (begin_ < 1) or (end_ > self.input_data.seq_len[SID_]):
                        begin_ = begin_ % self.input_data.seq_len[SID_]
                        end_ = end_ % self.input_data.seq_len[SID_] + self.input_data.seq_len[SID_]
                else:
                    if begin_ < 1:
                        if self.opt.end_overflow_policy == 'forbid':
                            raise RuntimeError('Site begin position must be positive')
                        elif self.opt.end_overflow_policy == 'trim':
                            if end_ < 1:
                                return
                            begin_ = 1
                        elif self.opt.end_overflow_policy == 'ignore':
                            return
                    if begin_ > end_:
                        raise RuntimeError('Site begin position cannot exceed the end position')
                    if end_ > self.input_data.seq_len[SID_]:
                        if self.opt.end_overflow_policy == 'forbid':
                            raise RuntimeError('Site end position cannot exceed the sequence length')
                        elif self.opt.end_overflow_policy == 'trim':
                            end_ = self.input_data.seq_len[SID_]
                            if begin_ > end_:
                                return
                        elif self.opt.end_overflow_policy == 'ignore':
                            return
                no = int(opt_prefix[-1])
                self.input_data.sites[no][GID_][SID_].append([begin_, end_])
                if self.opt.site_names:
                    self.input_data.sites[no][GID_][SID_][-1].append(site_name)
    def _save_record(self, opt_prefix, values, preliminary = False):
        """Method to save a record from an input file"""
        if opt_prefix == 'len_db':
            self._save_seq_len_db_record(values, self.opt.preparse_group_map)
        elif opt_prefix == 'group_map':
            self._save_group_map_record(values, preliminary)
        elif opt_prefix in ('anno1', 'anno2'):   
            self._save_annotation_record(opt_prefix, values)
    def _sort_annotations(self):
        """Method to sort the annotated sites for every sequence by begin symbol number"""
        for i in (1, 2):
            for group in self.input_data.sites[i].values():
                for sites in group.values():
                    sites.sort(key = lambda x: x[0])
    def _resolve_overlaps_within_annotations(self):
        """Method to resolve groups of overlapping sites within a given annotation according to the user-defined policy"""
        for i in (1, 2):
            policy = getattr(self.opt, 'anno{}_resolve_overlaps'.format(i))
            if policy == 'all':
                continue
            for group in self.input_data.sites[i].values():
                for sites in group.values():
                    sites_new = []
                    if policy == 'first':
                        last_end = 0
                        for site in sites:
                            if site[0] > last_end:
                                sites_new.append(site)
                            last_end = site[1]
                    elif policy == 'last':
                        next_begin = float('inf')
                        for site in reversed(sites):
                            if site[1] < next_begin:
                                sites_new.insert(0, site)
                            next_begin = site[0]
                    elif policy == 'merge':
                        last_end = 0
                        new_begin = 0
                        for site in sites:
                            if site[0] > last_end:
                                if new_begin > 0:
                                    sites_new.append([new_begin, last_end])
                                new_begin = site[0]
                            last_end = site[1]
                        if new_begin > 0:
                            sites_new.append([new_begin, last_end])
                    sites[: ] = sites_new       
    def calc_and_set_auto_seq_len(self):
        """Method to calculate the sequence length and, if applicable, the start of time series, on the basis of the input parameters if the sequence length table is not provided"""
        if self.opt.len_db:
            self.auto_seq_len = None
        else:
            if not self.global_state.time_unit_seconds:
                self.auto_seq_len = self.opt.seq_len
            else:
                interval = [self.opt.series_start, self.opt.series_finish]
                try:
                    self._convert_interval_to_time_structs(interval)
                except RuntimeError as e:
                    raise RuntimeError('Time series start and end: ' + e.args[0]) from None
                self.auto_seq_len = self._duration_in_units(interval[0], interval[1])
                self.auto_series_start = interval[0]
    def parse_sequence_length_db(self):
        """Method to parse the input sequence length table file"""
        self._parse_input_file('len_db')
        if self.opt.preparse_group_map:
            for SID in list(self.input_data.seq_len.keys()):
                if self.input_data.seq_len[SID] is None:
                    del self.input_data.seq_len[SID]
                    if self.opt.warnings:
                        print('Warning: SID "{}" is not the sequence length table. The group mapping record is ignored'.format(SID))
            self.input_data.group_map = DefaultOrderedDict(list)
        if not self.input_data.seq_len:
            error('The sequence length table does not contain any SIDs that can be retained')
        if not self.opt.quiet:
            print('The sequence length table has been read from "{}"'.format(getattr(self.opt, 'len_db')))
    def parse_group_map(self, preliminary = False):
        """Method to parse the input group map file"""
        if (not preliminary) and (not self.input_data.seq_len):
            error('The sequence length table must be parsed before the group mapping')
        if not self.opt.group_map:
            if not self.opt.sequences_as_groups:
                self.input_data.group_map[''].extend(self.input_data.seq_len.keys())
                if self.opt.single_sequence:
                    self.input_data.group_map[''] = ['']
                    self.input_data.seq_len[''] = self.auto_seq_len
            return
        if self.opt.non_overlapping_groups:
            self.reverse_group_map = {}
        self._parse_input_file('group_map', preliminary)
        if (self.opt.min_group_size > 1) or (self.opt.max_group_size > 0):
            for GID in list(self.input_data.group_map.keys()):
                SID_list = self.input_data.group_map[GID]
                if (len(SID_list) < self.opt.min_group_size) or (self.opt.max_group_size and (len(SID_list) > self.opt.max_group_size)):
                    del self.input_data.group_map[GID]
        if not self.input_data.group_map:
            error('The sequence length table does not contain any SIDs that can be retained')
        if not self.opt.quiet:
            print('The group mapping has been{} read from "{}"'.format(' preliminary' if self.auto_seq_len is None else '', getattr(self.opt, 'group_map')))
    def parse_annotations(self):
        """Method to parse the input annotation files"""
        if (self.opt.len_db and (not self.input_data.seq_len)) or ((not self.input_data.group_map) and (not self.opt.sequences_as_groups)):
            error('The annotation files must be parsed after the sequence length table and the group mapping')
        self._parse_input_file('anno1')
        if not self.opt.quiet:
            print('The first annotation has been read from "{}"'.format(getattr(self.opt, 'anno1')))
        self._parse_input_file('anno2')
        if not self.opt.quiet:
            print('The second annotation has been read from "{}"'.format(getattr(self.opt, 'anno2')))
        if (not self.input_data.sites[1]) or (not self.input_data.sites[2]):
            print(self.input_data.sites[2].keys())##/
            error('An annotation must not be empty')
        self._sort_annotations()
        self._resolve_overlaps_within_annotations()
    def get_data(self):
        return self.input_data

class InputFileProcessor:
    """Class for coordinating the input file processing"""
    def __init__(self, opt, file_parser):
        self.opt = opt
        self.file_parser = file_parser
    def process_input_files(self):
        """Method to coordinate processing of the input files"""
        self.file_parser.calc_and_set_auto_seq_len()
        if self.opt.preparse_group_map or (not self.opt.len_db):
            self.file_parser.parse_group_map(preliminary = True)
        if self.opt.len_db:
            self.file_parser.parse_sequence_length_db()
            self.file_parser.parse_group_map()
        self.file_parser.parse_annotations()
        return self.file_parser.get_data()

class BasicSequenceCalculator:
    """Abstract class for calculating basic measures and write into files required output annotations in a paricular sequence"""
    def __init__(self, global_state, opt, current_seq):
        self.global_state = global_state
        self.opt = opt
        self.current_seq = current_seq
        self.results = None
    def _in_union(self, idx):
        """Method to check if given symbol is in the annotation union"""
        raise NotImplementedError("Method '_in_union' is not implemented")
    def _in_intersection(self, idx):
        """Method to check if given symbol is in the annotation intersection"""
        raise NotImplementedError("Method '_in_intersection' is not implemented")
    def _in_complement1(self, idx):
        """Method to check if given symbol is in the annotation complement of the first"""
        raise NotImplementedError("Method '_in_complement1' is not implemented")
    def _in_complement2(self, idx):
        """Method to check if given symbol is in the annotation complement of the second"""
        raise NotImplementedError("Method '_in_complement2' is not implemented")
    def _in_re1(self, idx):
        """Method to check if given symbol is in the annotation of relative enrichment for the first annotation"""
        raise NotImplementedError("Method '_in_re1' is not implemented")
    def _in_re2(self, idx):
        """Method to check if given symbol is in the annotation of relative enrichment for the first annotation"""
        raise NotImplementedError("Method '_in_re2' is not implemented")
    def _write_site(self, file_handler, begin_idx, end_res):
        """Auxiliary method to write a site to an output annotation file"""
        group = (self.current_seq.GID + '\t' if self.current_seq.GID else '')
        file_handler.write('{}{}\t{}\t{}\n'.format(group, self.current_seq.SID, begin_idx + 1, end_res))
    def write_to_files(self, file_handlers):
        """Method to write the required output annotations"""
        for type_ in FileHandlers.output_file_types:
            file_handler = getattr(file_handlers, type_)
            if file_handler is None:
                continue
            if type_ in ('detailed', 'site'):
                continue
            else:
                in_site = False
                begin_idx = None
                last_end_res = None
                for idx in range(self.current_seq.length):
                    idx_ = (idx - self.current_seq.length) if (self.opt.circular and (idx >= self.current_seq.length)) else idx
                    if getattr(self, '_in_' + type_)(idx_):
                        if not in_site:
                            in_site = True
                            begin_idx = idx
                    elif in_site:
                        if self.opt.circular and (begin_idx == 0) and getattr(self, '_in_' + type_)(self.current_seq.length - 1):
                            in_site = False
                            last_end_res = idx + self.current_seq.length
                            continue
                        self._write_site(file_handler, begin_idx, idx)
                        in_site = False
                if in_site:
                    self._write_site(file_handler, begin_idx, last_end_res if last_end_res else self.current_seq.length)
    def calculate_residue_wise(self):
        """Method to calculate residue-wise measures for a given sequence"""
        raise NotImplementedError("Method 'calculate_residue_wise' is not implemented")
    def get_results(self):
        return self.results

class BasicBooleanSequenceCalculator(BasicSequenceCalculator):
    """Class for calculating basic Boolean measures and write into files required output annotations in a particular sequence"""
    def __init__(self, global_state, opt, current_seq):
        BasicSequenceCalculator.__init__(self, global_state, opt, current_seq)
        self.seq_length = current_seq.length
        self.seq = np.zeros(shape = current_seq.length, dtype = 'i1')
        self.results = BasicBooleanMeasures()
        self._classify_symbols()
    def _classify_symbols(self):
        """Method to classify symbols in the sequence by their occurrence in the annotations"""
        for site in self.current_seq.sites[1]:
            for idx in range(site[0] - 1, site[1]):
                if self.opt.circular and (idx >= self.seq_length):
                    idx -= self.seq_length
                self.seq[idx] = 1
        for site in self.current_seq.sites[2]:
            for idx in range(site[0] - 1, site[1]):
                if self.opt.circular and (idx >= self.seq_length):
                    idx -= self.seq_length
                if self.seq[idx] == 0:
                    self.seq[idx] = 2
                elif self.seq[idx] == 1:
                    self.seq[idx] = 3
    def _in_union(self, idx):
        """Method to check if given symbol is in the Bollean annotation union"""
        return True if self.seq[idx] >= 1 else False
    def _in_intersection(self, idx):
        """Method to check if given symbol is in the Bollean annotation intersection"""
        return True if self.seq[idx] == 3 else False
    def _in_complement1(self, idx):
        """Method to check if given symbol is in the Bollean annotation complement of the first"""
        return True if self.seq[idx] == 2 else False
    def _in_complement2(self, idx):
        """Method to check if given symbol is in the Bollean annotation complement of the second"""
        return True if self.seq[idx] == 1 else False
    def _get_overlapped_symbols_raw(self, begin, end, begin_, end_):
        """Method to calculate number of shared symbols between two sites in the same sequence, without further restrictions"""
        return max(min(end - begin_, end_ - begin, end - begin, end_ - begin_) + 1, 0)
    def _get_overlapped_symbols(self, site, site_, annotation_idx):
        """Method to calculate number of shared symbols between two sites in the same sequence"""
        if self.opt.predictor_nature != 'neutral':
            if (self.opt.predictor_nature == 'lagging') == (annotation_idx == 1):
                if site_[0] < site[0]:
                    return 0
            elif site_[0] > site[0]:
                return 0
        elif self.opt.circular:
            alt_0 = self._get_overlapped_symbols_raw(site[0], site[1], site_[0], site_[1]) 
            alt_1 = self._get_overlapped_symbols_raw(site[0], site[1], site_[0] + self.seq_length, site_[1] + self.seq_length) if site[1] > self.seq_length else 0
            alt_2 = self._get_overlapped_symbols_raw(site[0] + self.seq_length, site[1] + self.seq_length, site_[0], site_[1]) if site_[1] > self.seq_length else 0
            return max(alt_0, alt_1, alt_2)
        return self._get_overlapped_symbols_raw(site[0], site[1], site_[0], site_[1])
    def _get_site_length(self, site, site_):
        """Method to calculate the effective site length according to the input settings"""
        if self.opt.overlap_apply == 'shortest':
            return min(site[1] - site[0], site_[1] - site_[0]) + 1
        elif self.opt.overlap_apply == 'longest':
            return max(site[1] - site[0], site_[1] - site_[0]) + 1
        elif self.opt.overlap_apply in ('current', 'patched'):
            return site[1] - site[0] + 1
    def _check_overlap_sufficiency(self, overlapped_symbols, site_length):
        """Method to check if a goven overlap between sites satisfies the input overlap criteria"""
        return (overlapped_symbols >= self.opt.overlap_symbols) and (overlapped_symbols / site_length >= self.opt.overlap_part)
    def _write_measure_info_to_detailed_file(self, type_, description, detailed_file_h):
        """Writing the information on sybol counts in a specific category"""
        symbols_n = getattr(self.results, type_)
        category_name = getattr(self.global_state, type_ + '_name')
        ending = '{} is' if symbols_n == 1 else 's{} are'
        ending = ending.format(' gross' if (self.opt.gross and (type_ != 'aa')) else '')
        detailed_file_h.write('{}{} symbol{} {}{}'.format(self.global_state.indent_site, symbols_n, ending, description, category_name) + os.linesep)
    def calculate_residue_wise(self, detailed_file_h):
        """Method to calculate Boolean residue-wise measures for a given sequence and write the information to the detailed output file"""
        if self.opt.gross:
            self.results.pa = 0
            self.results.ap = 0
            for i in (1, 2):
                j = 3 - i
                for site in self.current_seq.sites[i]:
                    matched_symbols_n = np.sum(self.seq[site[0] - 1: site[1]] == 3)
                    site_length = site[1] - site[0] + 1
                    unmatched_symbols_n = site_length - matched_symbols_n
                    self.results.pp_[i] += matched_symbols_n
                    attr_name = 'pa' if i == 1 else 'ap'
                    setattr(self.results, attr_name, getattr(self.results, attr_name) + unmatched_symbols_n)
                if detailed_file_h:
                    ending = ('' if self.results.pp_[i] == 1 else 's') + ' gross'
                    message = '{}{} symbol{} present in the {} are also present in the {}'
                    detailed_file_h.write(message.format(self.global_state.indent_site, self.results.pp_[i], ending, self.global_state.anno_name[i], self.global_state.anno_name[j]) + os.linesep)
        else:
            self.results.pp = np.sum(self.seq == 3)
            self.results.pp_[1] = self.results.pp
            self.results.pp_[2] = self.results.pp
            if detailed_file_h:
                self._write_measure_info_to_detailed_file('pp', 'present in both annotations', detailed_file_h)
            self.results.pa = np.sum(self.seq == 1)
            self.results.ap = np.sum(self.seq == 2)
        self.results.aa = np.sum(self.seq == 0)
        if detailed_file_h:
            description = 'present exclusively in the ' + self.global_state.anno_name[1]
            self._write_measure_info_to_detailed_file('pa', description, detailed_file_h)
            description = 'present exclusively in the ' + self.global_state.anno_name[2]
            self._write_measure_info_to_detailed_file('ap', description, detailed_file_h)
            self._write_measure_info_to_detailed_file('aa', 'absent in both annotations', detailed_file_h)
    def calculate_site_wise(self, detailed_file_h, site_file_h):
        """Method to calculate site-wise measures and write the site-wise information to the detailed output file"""
        if self.opt.overlap_apply in ('shortest', 'longest', 'current'):
            for i in (1, 2):
                j = 3 - i
                for site in self.current_seq.sites[i]:
                    if self.opt.gross:
                        self.results.site_len[i] += site[1] - site[0] + 1
                    found_match = False
                    for site_ in self.current_seq.sites[j]:
                        if (not self.opt.circular) and (site_[0] > site[1]):
                            break
                        overlapped_symbols = self._get_overlapped_symbols(site, site_, i)
                        site_length_effective = self._get_site_length(site, site_)
                        if self._check_overlap_sufficiency(overlapped_symbols, site_length_effective):
                            self.results.site_m[i] += 1
                            found_match = True
                            break
                    if not found_match:
                        self.results.site_nm[i] += 1
                    if (detailed_file_h is not None) or (site_file_h is not None):
                        if found_match:
                            length_perc_1 = round(100 * overlapped_symbols / (site[1] - site[0] + 1))
                            length_perc_2 = round(100 * overlapped_symbols / (site_[1] - site_[0] + 1))
                            begin_ = site_[0]
                            end_ = site_[1]
                        else:
                            overlapped_symbols = length_perc_1 = length_perc_2 = 0
                            begin_ = end_ = '-'
                        if detailed_file_h is not None:
                            site_name_addition = ' ("{}")'.format(site[2]) if self.opt.site_names else ''
                            message = '{}Site {}-{}{} of the {}: '.format(self.global_state.indent_site, site[0], site[1], site_name_addition, self.global_state.anno_name[i])
                            if found_match:
                                ending = '' if overlapped_symbols == 1 else 's'
                                site_name_addition_ = ' ("{}")'.format(site_[2]) if self.opt.site_names else ''
                                message += 'overlaps with site {}-{}{} of the {} by {} symbol{} ({}% and {}% of the site lengths respectively)'
                                message = message.format(site_[0], site_[1], site_name_addition_, self.global_state.anno_name[j], overlapped_symbols, ending, length_perc_1, length_perc_2)
                            else:
                                message += 'no sufficient overlap found'
                            detailed_file_h.write(message + os.linesep)
                        if site_file_h is not None:
                            if self.opt.site_difference == 'matched':
                                if not found_match:
                                    continue
                            elif self.opt.site_difference == 'unmatched':
                                if found_match:
                                    continue
                            elif self.opt.site_difference == 'discrepant':
                                if (length_perc_1 == 100) and (length_perc_2 == 100):
                                    continue
                            list_ = [self.current_seq.GID] if self.opt.group_map else []
                            list_.extend([self.current_seq.SID, self.global_state.anno_short_name[i], site[0], site[1]])
                            if self.opt.site_names:
                                list_.append(site[2])
                            list_.extend([overlapped_symbols, length_perc_1, length_perc_2, begin_, end_])
                            if self.opt.site_names:
                                list_.append(site_[2] if found_match else '')
                            message = ('{}\t' * (len(list_) - 1) + '{}').format(*list_)
                            site_file_h.write(message + os.linesep)
        elif self.opt.overlap_apply == 'patched':
            for i in (1, 2):
                for site in self.current_seq.sites[i]:
                    if self.opt.circular and (site[1] > self.seq_length):
                        matched_symbols = np.sum(self.seq[site[0] - 1: ] == 3) + np.sum(self.seq[: site[1] - self.seq_length] == 3)
                    else:
                        matched_symbols = np.sum(self.seq[site[0] - 1: site[1]] == 3)
                    site_length = site[1] - site[0] + 1
                    found_match = self._check_overlap_sufficiency(matched_symbols, site_length)
                    if found_match:
                        self.results.site_m[i] += 1
                    else:
                        self.results.site_nm[i] += 1
                    if (detailed_file_h is not None) or (site_file_h is not None):
                        if found_match:
                            length_perc = round(100 * matched_symbols / site_length)
                        else:
                            length_perc = 0
                        if detailed_file_h is not None:
                            message = '{}Site {}-{} of the {}: '.format(self.global_state.indent_site, site[0], site[1], self.global_state.anno_name[i])
                            if found_match:
                                ending = '' if matched_symbols == 1 else 's'
                                message += 'overlaps by total {} symbol{} ({}% of the site length) with sites from the {}'.format(matched_symbols, ending, length_perc, self.global_state.anno_name[j])
                            else:
                                message += 'no sufficient overlap found'
                            detailed_file_h.write(message + os.linesep)
                    if site_file_h is not None:
                        if found_match and (self.opt.site_difference == 'unmatched'):
                            continue
                        list_ = [self.current_seq.GID] if self.opt.group_map else []
                        list_.extend([self.current_seq.SID, self.global_state.anno_short_name[i], site[0], site[1]])
                        if self.opt.site_names:
                            list_.append(site[2])
                        list_.extend([overlapped_symbols, length_perc])
                        message = ('{}\t' * (len(list_) - 1) + '{}').format(*list_)
                        site_file_h.write(message + os.linesep)
        else:
            error('Unknown overlap apply method')
        for i in (1, 2):
            if not self.opt.gross:
                self.results.site_len[i] = int(np.sum(self.seq == i) + np.sum(self.seq == 3))
            sites_n = self.results.site_m[i] + self.results.site_nm[i]
            if sites_n == 0:
                if detailed_file_h:
                    detailed_file_h.write('{}There are no sites in the {}'.format(self.global_state.indent_site, self.global_state.anno_name[i]) + os.linesep)
                continue
            else:
                ending = '' if sites_n == 1 else 's'
                verb = 'is' if sites_n == 1 else 'are'
                symbols = self.results.site_len[i]
                ending1 = '' if symbols == 1 else 's'
                insert0 = '' if self.opt.gross else ' unique'
                insert1 = ' gross' if self.opt.gross else ''
                message_base = '{}There {} {} site{} in the {} with total length {}{} symbol{}{}'
                if detailed_file_h:
                    detailed_file_h.write(message_base.format(self.global_state.indent_site, verb, sites_n, ending, self.global_state.anno_name[i], symbols, insert0, ending1, insert1) + os.linesep)
            if detailed_file_h:
                sites_n = self.results.site_m[i]
                ending = ' is' if sites_n == 1 else 's are'
                detailed_file_h.write('{}{} {} site{} matched in the {}'.format(self.global_state.indent_site, sites_n, self.global_state.anno_name[i], ending, self.global_state.anno_name[j]) + os.linesep)
                sites_n = self.results.site_nm[i]
                ending = ' has' if sites_n == 1 else 's have'
                detailed_file_h.write('{}{} {} site{} no match in the {}'.format(self.global_state.indent_site, sites_n, self.global_state.anno_name[i], ending, self.global_state.anno_name[j]) + os.linesep)
        
class BasicEnrichmentSequenceCalculator(BasicSequenceCalculator):
    """Class for calculating basic enrichment measures and write into files required output annotations in a particular sequence"""
    def __init__(self, global_state, opt, current_seq):
        BasicSequenceCalculator.__init__(self, global_state, opt, current_seq)
        self.n = self.opt.enrichment_count
        bytes_required = self._estimate_required_precison()
        self.seq_length = current_seq.length
        self.seq = [None] + [np.zeros(shape = current_seq.length, dtype = 'i' + str(bytes_required)) for x in range(2)]
        self.results = BasicEnrichmentMeasures()
        self._count_occurrences()
    def _count_occurrences(self):
        """Method to count the occurences in the annotations for each symbol in the sequence"""
        for i in (1, 2):
            for site in self.current_seq.sites[i]:
                for idx in range(site[0] - 1, site[1]):
                    if self.opt.circular and (idx >= self.seq_length):
                        idx -= self.seq_length
                    self.seq[i][idx] += 1
    def _in_union(self, idx):
        """Method to check if given symbol is in the enrichment annotation union"""
        for i in (1, 2):
            if self.seq[i][idx] >= self.n:
                return True
        return False
    def _in_intersection(self, idx):
        """Method to check if given symbol is in the enrichment annotation intersection"""
        check = 0
        for i in (1, 2):
            if self.seq[i][idx] >= self.n:
                check += 1
        return True if check == 2 else False
    def _in_complement1(self, idx):
        """Method to check if given symbol is in the enrichment annotation complement of the first"""
        return True if (self.seq[1][idx] < self.n) and (self.seq[2][idx] >= self.n) else False
    def _in_complement2(self, idx):
        """Method to check if given symbol is in the enrichment annotation complement of the second"""
        return True if (self.seq[1][idx] >= self.n) and (self.seq[2][idx] < self.n) else False
    def _in_re1(self, idx):
        """Method to check if given symbol is in the annotation of relative enrichment for the first annotation"""
        return True if self.seq[1][idx] - self.seq[2][idx] >= self.n else False
    def _in_re2(self, idx):
        """Method to check if given symbol is in the annotation of relative enrichment for the first annotation"""
        return True if self.seq[2][idx] - self.seq[1][idx] >= self.n else False
    def _estimate_required_precison(self):
        """Method to provide a higher estimation for the reeqired integer precision of the counts"""
        count_limit = max(len(self.current_seq.sites[1]), len(self.current_seq.sites[2])) + 1
        bytes_required = math.ceil((math.ceil(math.log(count_limit, 2)) + 1) / 8)
        if bytes_required > 8:
            error("Too  many sites annotated. Maximal number is {}".format(2 ** 63 - 1))
        bytes_required = [1, 1, 2, 4, 4, 8, 8, 8, 8][bytes_required]
        return bytes_required
    def _write_measure_info_to_detailed_file(self, symbols_n, description, file_handler):
        """Writing the information on symbol counts in a specific category"""
        ending = ' is' if symbols_n == 1 else 's are'
        file_handler.write('{}{} symbol{} enriched in {}'.format(self.global_state.indent_site, symbols_n, ending, description) + os.linesep)
    def calculate_residue_wise(self, detailed_file_h):
        """Method to calculate count residue-wise measures for a given sequence"""
        for i in (1, 2):
            j = 2 if i == 1 else 1
            self.results.e[i] = np.sum(self.seq[i] >= self.n)
            description = 'the ' + self.global_state.anno_name[i]
            if detailed_file_h:
                self._write_measure_info_to_detailed_file(self.results.e[i], description, detailed_file_h)
            self.results.re[i] = np.sum((self.seq[i] - self.seq[j]) >= self.n)
        self.results.ee = np.sum((self.seq[1] >= self.n) * (self.seq[2] >= self.n))
        if detailed_file_h:
            self._write_measure_info_to_detailed_file(self.results.ee, 'both annotations', detailed_file_h)
        self.results.ne = np.sum((self.seq[1] < self.n) * (self.seq[2] < self.n))
        if detailed_file_h:
            self._write_measure_info_to_detailed_file(self.results.ne, 'neither annotations', detailed_file_h)
        self.results.nre = len(self.seq) - self.results.re[1] - self.results.re[2]
        
class PerformanceCalculator:
    """Class to calculate all selected performance measures on the basis of the basic measures"""
    def __init__(self, basic_measures, performance_maeasures):
        self.basic_measures = basic_measures
        self.performance_measures = performance_maeasures
    def _calc_p1(self):
        """Method to calculate fraction of symbols present in the first annotation"""
        p1 = self.basic_measures.pp + self.basic_measures.pa
        self.performance_measures.set_value('p1', p1)
    def _calc_p2(self):
        """Method to calculate fraction of symbols present in the second annotation"""
        p2 = self.basic_measures.pp + self.basic_measures.ap
        self.performance_measures.set_value('p2', p2)
    def _calc_pp(self):
        """Method to copy fraction of symbols present in both annotations"""
        pp = self.basic_measures.pp
        self.performance_measures.set_value('pp', pp)
    def _calc_pp1(self):
        """Method to copy fraction of symbols gross present in the first annotation that are also present in the second"""
        pp1 = self.basic_measures.pp_[1]
        self.performance_measures.set_value('pp1', pp1)
    def _calc_pp2(self):
        """Method to copy fraction of symbols gross present in the second annotation that are also present in the first"""
        pp2 = self.basic_measures.pp_[2]
        self.performance_measures.set_value('pp2', pp2)
    def _calc_pa(self):
        """Method to copy fraction of symbols present exclusively in the first annotation"""
        pa = self.basic_measures.pa
        self.performance_measures.set_value('pa', pa)
    def _calc_ap(self):
        """Method to copy fraction of symbols present exclusively in the second annotation"""
        ap = self.basic_measures.ap
        self.performance_measures.set_value('ap', ap)
    def _calc_aa(self):
        """Method to copy fraction of symbols absent in both annotation"""
        aa = self.basic_measures.aa
        self.performance_measures.set_value('aa', aa)
    def _calc_rc2(self):
        """Method to calculate symbol-wise recall for the second annotation"""
        denominator = self.basic_measures.pp_[1] + self.basic_measures.pa
        rc2 = self.basic_measures.pp_[1] / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('rc2', rc2)
    def _calc_pr2(self):
        """Method to calculate symbol-wise precision for the second annotation"""
        denominator = self.basic_measures.pp_[2] + self.basic_measures.ap
        pr2 = self.basic_measures.pp_[2] / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('pr2', pr2)
    def _calc_sp2(self):
        """Method to calculate symbol-wise specificity for the second annotation"""
        denominator = self.basic_measures.aa + self.basic_measures.ap
        sp2 = self.basic_measures.aa / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('sp2', sp2)
    def _calc_npv2(self):
        """Method to calculate symbol-wise negative predictive value for the second annotation"""
        denominator = self.basic_measures.aa + self.basic_measures.pa
        npv2 = self.basic_measures.aa / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('npv2', npv2)
    def _calc_in2(self):
        """Method to calculate symbol-wise informedness for the second annotation"""
        in2 = self.performance_measures.get_value('rc2') + self.performance_measures.get_value('sp2') - 1
        self.performance_measures.set_value('in2', in2)
    def _calc_mk2(self):
        """Method to calculate symbol-wise markedness for the second annotation"""
        mk2 = self.performance_measures.get_value('pr2') + self.performance_measures.get_value('npv2') - 1
        self.performance_measures.set_value('mk2', mk2)
    def _calc_pc(self):
        """Method to calculate symbol-wise performance coefficient"""
        denominator = self.basic_measures.pp_[2] + self.basic_measures.ap + self.basic_measures.pa
        pc = self.basic_measures.pp_[2] / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('pc', pc)
    def _calc_acc(self):
        """Method to calculate symbol-wise accuracy ACC"""
        acc = self.basic_measures.pp + self.basic_measures.aa
        self.performance_measures.set_value('acc', acc)
    def _calc_mcc(self):
        """Method to calculate symbol-wise Matthews correlation coefficient"""
        numerator = self.basic_measures.pp * self.basic_measures.aa + self.basic_measures.ap * self.basic_measures.pa
        temp = (self.basic_measures.pp + self.basic_measures.pa) * (self.basic_measures.pp + self.basic_measures.ap) * (self.basic_measures.aa + self.basic_measures.ap) * (self.basic_measures.aa + self.basic_measures.pa)
        mcc = numerator / math.sqrt(temp) if temp > 0.0 else float('nan')
        self.performance_measures.set_value('mcc', mcc)
    def _calc_f1(self):
        """Method to calculate symbol-wise F1 score"""
        denominator = 2 * self.basic_measures.pp_[1] * self.basic_measures.pp_[2] + self.basic_measures.pp_[2] * self.basic_measures.pa + self.basic_measures.pp_[1] * self.basic_measures.ap
        f1 = 2 * self.basic_measures.pp_[1] * self.basic_measures.pp_[2] / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('f1', f1)
    def _calc_site_n1(self):
        """Method to calculate number of sites in the first annotation"""
        site_n1 = self.basic_measures.site_m[1] + self.basic_measures.site_nm[1]
        self.performance_measures.set_value('site_n1', site_n1)
    def _calc_site_n2(self):
        """Method to calculate number of sites in the second annotation"""
        site_n2 = self.basic_measures.site_m[2] + self.basic_measures.site_nm[2]
        self.performance_measures.set_value('site_n2', site_n2)
    def _calc_site_n1_m(self):
        """Method to copy number of matched sites in the first annotation"""
        site_n1_m = self.basic_measures.site_m[1]
        self.performance_measures.set_value('site_n1_m', site_n1_m)
    def _calc_site_n2_m(self):
        """Method to copy number of matched sites in the second annotation"""
        site_n2_m = self.basic_measures.site_m[2]
        self.performance_measures.set_value('site_n2_m', site_n2_m)
    def _calc_site_len1(self):
        """Method to copy total length of sites in the first annotation"""
        site_len1 = self.basic_measures.site_len[1]
        self.performance_measures.set_value('site_len1', site_len1)
    def _calc_site_len2(self):
        """Method to copy total length of sites in the second annotation"""
        site_len2 = self.basic_measures.site_len[2]
        self.performance_measures.set_value('site_len2', site_len2)
    def _calc_site_rc2(self):
        """Method to calculate site-wise recall for the second annotation"""
        denominator = self.basic_measures.site_m[1] + self.basic_measures.site_nm[1]
        site_rc2 = self.basic_measures.site_m[1] / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('site_rc2', site_rc2)
    def _calc_site_pr2(self):
        """Method to calculate site-wise recall for the second annotation"""
        denominator = self.basic_measures.site_m[2] + self.basic_measures.site_nm[2]
        site_pr2 = self.basic_measures.site_m[2] / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('site_pr2', site_pr2)
    def _calc_site_pc2(self):
        """Method to calculate site-wise performance coefficient for the second annotation"""
        denominator = self.basic_measures.site_m[2] + self.basic_measures.site_nm[2] + self.basic_measures.site_nm[1]
        site_pc2 = self.basic_measures.site_m[2] / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('site_pc2', site_pc2)
    def _calc_site_f1(self):
        """Method to calculate site-wise F1 score"""
        temp = 2 * self.basic_measures.site_m[1] * self.basic_measures.site_m[2]
        denominator = temp + self.basic_measures.site_nm[1] * self.basic_measures.site_m[2] + self.basic_measures.site_nm[2] * self.basic_measures.site_m[1]
        site_f1 = temp / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('site_f1', site_f1)
    def _calc_site_pcv(self):
        """Method to calculate site-wise positive correlation value"""
        temp = self.basic_measures.site_m[1] + self.basic_measures.site_m[2]
        denominator = temp + self.basic_measures.site_nm[1] + self.basic_measures.site_nm[2]
        site_pcv = temp / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('site_pcv', site_pcv)
    def _calc_e_p1(self):
        """Method to copy fraction of symbols enriched in the first annotation"""
        e_p1 = self.basic_measures.e[1]
        self.performance_measures.set_value('e_p1', e_p1)
    def _calc_e_p2(self):
        """Method to copy fraction of symbols enriched in the second annotation"""
        e_p2 = self.basic_measures.e[2]
        self.performance_measures.set_value('e_p2', e_p2)
    def _calc_e_pp(self):
        """Method to copy fraction of symbols enriched in both annotations"""
        e_pp = self.basic_measures.ee
        self.performance_measures.set_value('e_pp', e_pp)
    def _calc_e_pa(self):
        """Method to calculate fraction of symbols enriched exclusively in the first annotation"""
        e_pa = self.basic_measures.e[1] - self.basic_measures.ee
        self.performance_measures.set_value('e_pa', e_pa)
    def _calc_e_ap(self):
        """Method to calculate fraction of symbols enriched exclusively in the second annotation"""
        e_ap = self.basic_measures.e[2] - self.basic_measures.ee
        self.performance_measures.set_value('e_ap', e_ap)
    def _calc_e_aa(self):
        """Method to copy fraction of symbols enriched in neither annotation"""
        e_aa = self.basic_measures.ne
        self.performance_measures.set_value('e_aa', e_aa)
    def _calc_e_rc2(self):
        """Method to calculate symbol-wise enrichment recall for the second annotation"""
        e_rc2 = self.basic_measures.ee / self.basic_measures.e[1] if self.basic_measures.e[1] > 0.0 else float('nan')
        self.performance_measures.set_value('e_rc2', e_rc2)
    def _calc_e_pr2(self):
        """Method to calculate symbol-wise enrichment precision for the second annotation"""
        e_pr2 = self.basic_measures.ee / self.basic_measures.e[2] if self.basic_measures.e[2] > 0.0 else float('nan')
        self.performance_measures.set_value('e_pr2', e_pr2)
    def _calc_e_sp2(self):
        """Method to calculate symbol-wise enrichment specificity for the second annotation"""
        denominator = 1 - self.basic_measures.e[1]
        e_sp2 = self.basic_measures.ne / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('e_sp2', e_sp2)
    def _calc_e_npv2(self):
        """Method to calculate symbol-wise enrichment negative predictive value for the second annotation"""
        denominator = 1 - self.basic_measures.e[2]
        e_npv2 = self.basic_measures.ne / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('e_npv2', e_npv2)
    def _calc_e_in2(self):
        """Method to calculate symbol-wise enrichment informedness for the second annotation"""
        e_in2 = self.performance_measures.get_value('e_rc2') + self.performance_measures.get_value('e_sp2') - 1
        self.performance_measures.set_value('e_in2', e_in2)
    def _calc_e_mk2(self):
        """Method to calculate symbol-wise enrichment markedness for the second annotation"""
        e_mk2 = self.performance_measures.get_value('e_pr2') + self.performance_measures.get_value('e_npv2') - 1
        self.performance_measures.set_value('e_mk2', e_mk2)
    def _calc_e_pc(self):
        """Method to calculate symbol-wise enrichment performance coefficient"""
        denominator = 1 - self.basic_measures.ne
        e_pc = self.basic_measures.ee / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('e_pc', e_pc)
    def _calc_e_acc(self):
        """Method to calculate symbol-wise enrichment accuracy ACC"""
        e_acc = self.basic_measures.ee + self.basic_measures.ne
        self.performance_measures.set_value('e_acc', e_acc)
    def _calc_e_mcc(self):
        """Method to calculate symbol-wise enrichment Matthews correlation coefficient"""
        only_in_1 = self.basic_measures.e[1] - self.basic_measures.ee
        only_in_2 = self.basic_measures.e[2] - self.basic_measures.ee
        numerator = self.basic_measures.ee * self.basic_measures.ne + only_in_1 * only_in_2
        denominator = math.sqrt((self.basic_measures.ee + only_in_1) * (self.basic_measures.ee + only_in_2) * (self.basic_measures.ne + only_in_1) * (self.basic_measures.ne + only_in_2))
        e_mcc = numerator / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('e_mcc', e_mcc)
    def _calc_e_f1(self):
        """Method to calculate symbol-wise enrichment F1 score"""
        denominator = self.basic_measures.e[1] + self.basic_measures.e[2]
        e_f1 = 2 * self.basic_measures.ee / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('e_f1', e_f1)
    def _calc_e_eac(self):
        """Method to calculate enrichment asymmetry coefficient"""
        denominator = self.basic_measures.e[1] + self.basic_measures.e[2] - self.basic_measures.ee
        e_eac = (self.basic_measures.re[1] + self.basic_measures.re[2]) / denominator if denominator > 0.0 else float('nan')
        self.performance_measures.set_value('e_eac', e_eac)
    def _calc_seq_n(self):
        """Method to copy the number of sequences in the group"""
        seq_n = self.basic_measures.seq_n
        self.performance_measures.set_value('seq_n', seq_n)
    def _calc_e_seq_n(self):
        """Method to copy the number of sequences in the group for the enrichment mode"""
        e_seq_n = self.basic_measures.seq_n
        self.performance_measures.set_value('e_seq_n', e_seq_n)
    def calculate_performance_measures(self):
        """Method to run all relevant measure calculating methods"""
        for measure in self.performance_measures.name_map:
            getattr(self, '_calc_' + measure.var_name)()

class CalculationCoordinator():
    """Class to coordinate the process of performance measures calculation in accordance with the given averaging approach"""
    def __init__(self, global_state, opt, input_data, file_handlers):
        self.global_state = global_state
        self.opt = opt
        self.input_data = input_data
        self.file_handlers = file_handlers
    def _process_sequence(self, current_seq):
        """Method to calculate basic measures for annotatopns of sites in a particular sequence in a particular group"""
        args = (self.global_state, self.opt, current_seq)
        if self.file_handlers.detailed is not None:
            ending = '' if current_seq.length == 1 else 's'
            seq_description = 'sequence "{}"'.format(current_seq.SID) if current_seq.SID else 'unnamed sequence'
            self.file_handlers.detailed.write('{}Information on the {} (length {} symbol{}):'.format(self.global_state.indent_seq, seq_description, current_seq.length, ending) + os.linesep)
        basic_sequence_calculator = BasicBooleanSequenceCalculator(*args) if self.opt.enrichment_count == 0 else BasicEnrichmentSequenceCalculator(*args)
        basic_sequence_calculator.calculate_residue_wise(self.file_handlers.detailed)
        if self.opt.enrichment_count == 0:
            basic_sequence_calculator.calculate_site_wise(self.file_handlers.detailed, self.file_handlers.site)
        basic_sequence_calculator.write_to_files(self.file_handlers)
        return basic_sequence_calculator.get_results()
    def process_group(self, GID):
        """Method to calculate all relevant performance measures for a giben sequence group"""
        if self.opt.averaging == 'sequence':
            group_performance_measures = None
        else:
            group_counts = None
        group_counts_ = None
        seq_length_sum = 0
        if self.opt.grouped and (self.file_handlers.detailed is not None):
            group_len = len(self.input_data.group_map[GID])
            self.file_handlers.detailed.write('Information on the group "{}" (contains {} sequence{}):'.format(GID, group_len, ('s' if group_len > 1 else '')) + os.linesep)
        for SID in self.input_data.group_map[GID]:
            seq_length = self.input_data.seq_len[SID]
            sites = [None] + [self.input_data.sites[i][GID][SID] for i in (1, 2)]
            current_seq = CurrentSequence(GID, SID, seq_length, sites)
            seq_counts = self._process_sequence(current_seq)
            if self.opt.averaging == 'sequence':
                seq_performance_measures = PerformanceMeasures(self.opt.enrichment_count, self.opt.benchmark, self.opt.gross)
                PerformanceCalculator(seq_counts, seq_performance_measures).calculate_performance_measures()
                seq_performance_measures.set_value('seq_n', 1)
                if group_performance_measures is None:
                    group_performance_measures = seq_performance_measures
                else:
                    group_performance_measures += seq_performance_measures
            else:
                if self.opt.len_adjust:
                    seq_counts /= seq_length
                if group_counts is None:
                    group_counts = seq_counts
                else:
                    group_counts += seq_counts
            seq_length_sum += seq_length
        group_seq_n = len(self.input_data.group_map[GID])
        if self.opt.averaging == 'sequence':
            for measure in group_performance_measures.name_map:
                if measure.basic:
                    group_performance_measures.set_count(measure.var_name, seq_length_sum)
            group_performance_measures.average(group_seq_n if self.opt.na_zeros else 0)
        else:
            group_counts.seq_n = group_seq_n
            if self.opt.averaging == 'dataset':
                group_counts_ = copy.deepcopy(group_counts)
            group_counts /= (group_seq_n if self.opt.len_adjust else seq_length_sum)
            group_performance_measures = PerformanceMeasures(self.opt.enrichment_count, self.opt.benchmark, self.opt.gross)
            PerformanceCalculator(group_counts, group_performance_measures).calculate_performance_measures()
        return group_performance_measures, group_counts_, seq_length_sum

class DataProcessor:
    """Class to calculate and save into corresponding files performance measures as well as output annotations for each group and the whole dataset"""
    def __init__(self, opt, global_state, input_data):
        self.opt = opt
        self.input_data = input_data
        self.file_handlers = FileHandlers()
        self.global_state = global_state
        self.calculator = CalculationCoordinator(global_state, opt, input_data, self.file_handlers)
        self.dataset_performance_measures = PerformanceMeasures(self.opt.enrichment_count, self.opt.benchmark, self.opt.gross)
    def _float_to_fixed_width_str(value, width):
        """Method to make the best attempt to represent a float as fixed-width string"""
        for i in range(width - 2, -1, -1):
            str0 = '{:.{}f}'.format(value, i)
            if len(str0) <= width:
                break
        return str0
    def _positive_int_to_fixed_width_str(value, width):
        """Method to make the best attempt to represent a posititve integer as fixed-width string"""
        str0 = str(value)
        if len(str0) <= width:
            return str0
        for i in range(width - 4, -1, -1):
            str0 = '{:.{}E}'.format(value, i)
            if len(str0) <= width + 1:
                break
            if (len(str0) == width + 2) and '+0' in str0:
                break
        return str0.replace('+', '').replace('E0', 'E')
    def _open_output_files(self):
        """Method to open required output annotation files"""
        for type_ in FileHandlers.output_file_types:
            filepath = getattr(self.opt, 'output_file_' + type_)
            if not filepath:
                continue
            file_handler = open(filepath, 'w')
            if type_ == 'detailed':
                pass
            elif type_ == 'site':
                header = ('Group\t' if self.opt.group_map else '') + 'Sequence\tAnnotation\tSite begin\tSite end\t'
                if self.opt.site_names:
                    header += 'Site name\t'
                if self.opt.overlap_apply == 'patched':
                    header += 'Matched symbols\tMatched perc.\n'
                else:
                    header += 'Overlapped symbols\tOverlapped perc.\tPartner overlapped perc.\tPartner begin\tPartner end'
                if self.opt.site_names:
                    header += '\tPartner name'
                file_handler.write(header + os.linesep)
            else:
                header = ('Group\t' if self.opt.grouped else '') + 'Sequence\tbegin\tend\n'
                file_handler.write(header)
            setattr(self.file_handlers, type_, file_handler)
    def _close_output_files(self):
        """Method to close ouptput annotation files"""
        for type_ in FileHandlers.output_file_types:
            handler = getattr(self.file_handlers, type_)
            if handler is None:
                continue
            handler.close()
            filepath = getattr(self.opt, 'output_file_' + type_)
            description = type_
            description.replace('complement', 'complement of')
            description = re.sub('$re', 'relative enrichment for', description)
            description.replace('1', ' the first')
            description.replace('1', ' the second')
            if type_ == 'detailed':
                output_type = 'output'
            elif type_ == 'site':
                output_type = 'statistics'
            else:
                output_type = 'annotatiion'
            description += ' ' + output_type
            if not self.opt.quiet:
                print("The {} file '{}' has been written".format(description, filepath))
    def _generate_header(self, grouped):
        """Method to form the header with basic launch information as well as relevant column names"""
        header = ''
        if not self.opt.clean:
            header += "# This file was generated at {} with SLALOM".format(str(datetime.datetime.now())[: -7]) + os.linesep
            header += '# Command line options (unquoted and unescaped): ' + ' '.join(sys.argv[1: ]) + os.linesep
            header += '# The following statistics have been calculated:' + os.linesep
        column_names = ((('Frame' if self.opt.genbank else 'Seq.') if self.opt.sequences_as_groups else 'Group') + '\t') if grouped else ''
        for measure in self.dataset_performance_measures.name_map:
            if not self.opt.clean:
                header += '#    {}{}: {}'.format(measure.displayed_name, '' if (measure.force_avg or self.opt.averaging != 'dataset') else '*', measure.description) + os.linesep
            column_names += measure.displayed_name + '\t'
        if not self.opt.clean:
            if self.opt.averaging == 'dataset':
                message = '# The dataset-wide averages are micro-averages of the group-wide metrics for the measures marked with asterisk and are macro-averages for the rest.'
            else:
                message = '# All the provided dataset-wide averages are macro-averages of the group-wide metrics'
            header += message + os.linesep
        header += column_names + os.linesep
        return header
    def _produce_bottom_lines_string(self, results, attr_names, groups_n):
        """Method to calculate averages and sums of relevant performance measures and save them to the string"""
        avg_string = 'Average\t' if groups_n else ''
        if self.opt.calculate_sums:
            sum_string = 'Sum\t'
        for attr_name in attr_names:
            value = results.get_value(attr_name)
            if self.opt.calculate_sums:
                sum_string += ((DataProcessor._positive_int_to_fixed_width_str(value, 7) if type(value) == int else '') + '\t')
            if not self.opt.na_zeros:
                if not math.isnan(value):
                    if (self.opt.averaging == 'dataset') and (type(value) == int):
                        value /= groups_n
                    elif groups_n:
                        value /= results.get_count(attr_name)
            else:
                if math.isnan(value):
                    value = 0.0
                else:
                    value = value / groups_n if (groups_n and self.opt.averaging != 'dataset') else value
            avg_string += (DataProcessor._positive_int_to_fixed_width_str(value, 7) if type(value)==int else DataProcessor._float_to_fixed_width_str(value, 6)) + '\t'
        return avg_string[: -1] + (os.linesep + sum_string[: -1] if self.opt.calculate_sums else '')
    def process(self):
        """Method to coordinate the input data processing and outputting"""
        self._open_output_files()
        with open(self.opt.output_file, 'w') as ofile:
            header = self._generate_header(self.opt.grouped)
            ofile.write(header)
            attr_names = [x.var_name for x in self.dataset_performance_measures.name_map]
            if self.opt.grouped:
                if self.opt.averaging == 'dataset':
                    dataset_counts = None
                    seq_length_sum_dataset = 0
                for GID in (sorted(self.input_data.group_map.keys()) if self.opt.sort_output else self.input_data.group_map.keys()):
                    group_performance_measures, group_counts, seq_length_sum_group = self.calculator.process_group(GID)
                    row = GID
                    for attr_name in attr_names:
                        value = group_performance_measures.get_value(attr_name)
                        value = '{:.4f}'.format(value) if type(value) != int else str(value)
                        row += '\t' + value
                    ofile.write(row + os.linesep)
                    if self.opt.averaging == 'dataset':
                        if dataset_counts is None:
                            dataset_counts = group_counts
                        else:
                            dataset_counts += group_counts
                        seq_length_sum_dataset += seq_length_sum_group
                    else:
                        self.dataset_performance_measures += group_performance_measures
                groups_n = len(self.input_data.group_map)
                if self.opt.averaging == 'dataset':
                    dataset_counts /= seq_length_sum_dataset
                    PerformanceCalculator(dataset_counts, self.dataset_performance_measures).calculate_performance_measures()
                if not self.opt.clean:
                    ofile.write('#' + '-' * (8 * (len(attr_names) + 1) - 1) + os.linesep)
            else:
                self.dataset_performance_measures = self.calculator.process_group('')[0]
                groups_n = 0
            ofile.write(self._produce_bottom_lines_string(self.dataset_performance_measures, attr_names, groups_n))
        if not self.opt.quiet:
            print("The output file '{}' with performance measures has been written".format(self.opt.output_file))
        self._close_output_files()
            