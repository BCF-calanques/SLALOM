import math, copy
from collections import defaultdict, OrderedDict, Callable

class DefaultOrderedDict(OrderedDict):
    """Default ordered dictionary"""
    # Source: http://stackoverflow.com/a/6190500/562769
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None) and (not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory
    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)
    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value
    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()
    def copy(self):
        return self.__copy__()
    def __copy__(self):
        return type(self)(self.default_factory, self)
    def __deepcopy__(self, memo):
        return type(self)(self.default_factory, copy.deepcopy(self.items()))
    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory, OrderedDict.__repr__(self))

class EnrichmentCountType:
    """Class to mark the type of the enrichment count command line parameter"""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return self.value

class InputData:
    """Class to hold the mapping and annotation data"""
    def __init__(self):
        self.seq_len = {}
        self.time_series_starts = {}
        self.group_map = DefaultOrderedDict(list)
        self.sites = (None, defaultdict(lambda: defaultdict(list)), defaultdict(lambda: defaultdict(list)))

class GlobalState:
    """Class to hold the global state of the program"""
    def __init__(self, opt):
        self.anno_name = [None, None, None]
        self.anno_name[1] = 'benchmark' if opt.benchmark else 'first annotation'
        self.anno_name[2] = 'prediction' if opt.benchmark else 'second annotation'
        self.anno_short_name = [None, None, None]
        self.anno_short_name[1] = 'b' if opt.benchmark else '1'
        self.anno_short_name[2] = 'p' if opt.benchmark else '2'
        self.pp_name = ' (true positives)' if opt.benchmark else ''
        self.pa_name = ' (false negatives)' if opt.benchmark else ''
        self.ap_name = ' (false positives)' if opt.benchmark else ''
        self.aa_name = ' (true negatives)' if opt.benchmark else ''
        self.indent_site = '        ' if opt.group_map else '    '
        self.indent_seq = '    ' if opt.group_map else ''
        self.anno_read_SIDs = [None, None, None]
        self.anno_read_SIDs[1] = True if (opt.single_sequence or opt.anno1_all_sequences) else False
        self.anno_read_SIDs[2] = True if (opt.single_sequence or opt.anno2_all_sequences) else False
        if opt.time_unit == 'sec':
            self.time_unit_seconds = 1
        elif opt.time_unit == 'min':
            self.time_unit_seconds = 60
        elif opt.time_unit == 'hour':
            self.time_unit_seconds = 3600
        elif opt.time_unit == 'day':
            self.time_unit_seconds = 3600 * 24
        else:
            self.time_unit_seconds = 0

class CurrentSequence:
    """Class to hold the information about the current sequence being processed"""
    def __init__(self, GID, SID, length, sites):
        self.GID = GID
        self.SID = SID
        self.length = length
        self.sites = sites

class BasicMeasures:
    """Class to hold required basic measures for a sequence of a group"""
    def __init__(self):
        self._attr_list = list(filter(lambda x: not x.startswith('_'), dir(self)))
    def __iadd__(self, other):
        for attr in self._attr_list:
            content = getattr(self, attr)
            if content is None:
                continue
            if not isinstance(content, list):
                setattr(self, attr, content + getattr(other, attr))
            else:
                for i in (1, 2):
                    content[i] += getattr(other, attr)[i]
        return self
    def __itruediv__(self, seq_length):
        for attr in self._attr_list:
            if attr.startswith('s'):
                continue
            content = getattr(self, attr)
            if content is None:
                continue
            if not isinstance(content, list):
                setattr(self, attr, content / seq_length)
            else:
                for i in (1, 2):
                    content[i] /= seq_length
        return self
    def __truediv__(self, seq_length):
        other = copy.deepcopy(self)
        other /= seq_length
        return other

class BasicBooleanMeasures(BasicMeasures):
    """Class to hold basic Boolean measures for a sequence of a group"""
    def __init__(self):
        self.pp = None
        self.pp_ = [None, 0, 0]
        self.pa = None
        self.ap = None
        self.aa = None
        self.site_m = [None, 0, 0]
        self.site_nm = [None, 0, 0]
        self.site_len = [None, 0, 0]
        self.seq_n = None
        BasicMeasures.__init__(self)

class BasicEnrichmentMeasures(BasicMeasures):
    """Class to hold basic count measures for a sequence of a group"""
    def __init__(self):
        self.e = [None, 0, 0]
        self.ee = None
        self.ne = None
        self.re = [None, 0, 0]
        self.nre = None
        BasicMeasures.__init__(self)

class MeasureType:
    """Class to store the information abour a specific measure displayed in the output"""
    def __init__(self, displayed_name, var_name, description, type_, force_avg, basic):
        self.displayed_name = displayed_name
        self.var_name = var_name
        self.description = description
        self.type_ = type_
        self.force_avg = force_avg
        self.basic = basic

class MeasureFullType(MeasureType):
    """Class to store the information abour a potentially displayed measure"""
    def __init__(self, displayed_name, var_name, description, type_, mode_Bs = False, mode_Bg = False, mode_En = False, mode_Eq = False, mode_Bn = False, force_avg = False, basic = False):
        MeasureType.__init__(self, displayed_name, var_name, description, type_, force_avg, basic)
        self.mode_Bs = mode_Bs
        self.mode_Bg = mode_Bg
        self.mode_En = mode_En
        self.mode_Eq = mode_Eq
        self.mode_Bn = mode_Bn 

class PerformanceMeasures:
    """Class to hold performance measures for a grpup or database-wide"""
    name_map_full = (
        MeasureFullType('Nseq', 'seq_n', 'Number of sequences', 'int', mode_Bs = True, mode_Bg = True, mode_En = True, mode_Eq = True, mode_Bn = True, force_avg = True),
        MeasureFullType('Rel.', 'p1', 'Share of relevant symbols (=TP+FN)', 'float', mode_Bs = True, mode_En = True, mode_Bn = True, basic = True),
        MeasureFullType('Sel.', 'p2', 'Share of selected/predicted symbols (=TP+FP)', 'float', mode_Bs = True, mode_En = True, mode_Bn = True, basic = True),
        MeasureFullType('P1', 'p1', 'Share of symbols present in the first annotation', 'float', mode_Bs = True, mode_Eq = True, basic = True),
        MeasureFullType('P2', 'p2', 'Share of symbols present in the second annotation', 'float', mode_Bs = True, mode_Eq = True, basic = True),
        MeasureFullType('E1', 'p1', 'Share of symbols enriched in the first annotation', 'float', mode_En = True, mode_Eq = True, basic = True),
        MeasureFullType('E2', 'p2', 'Share of symbols enriched in the second annotation', 'float', mode_En = True, mode_Eq = True, basic = True),
        MeasureFullType('TP', 'pp', 'Share of true positive symbols', 'float', mode_Bs = True, mode_En = True, mode_Bn = True, basic = True),
        MeasureFullType('TP_b', 'pp1', 'Ratio of true positive symbols gross in the benchmark to the overall symbol count', 'float', mode_Bg = True, mode_Bn = True, basic = True),
        MeasureFullType('TP_p', 'pp2', 'Ratio of true positive symbols gross in the prediction to the overall symbol count', 'float', mode_Bg = True, mode_Bn = True, basic = True),
        MeasureFullType('PP', 'pp', 'Share of symbols present in both annotations', 'float', mode_Bs = True, mode_Eq = True, basic = True),
        MeasureFullType('P1P2', 'pp1', 'Ratio of symbols gross present in the first annotation that are also present in the second to the overall symbol count', 'float', mode_Bg = True, mode_Eq = True, basic = True),
        MeasureFullType('P2P1', 'pp2', 'Ratio of symbols gross present in the second annotation that are also present in the first to the overall symbol count', 'float', mode_Bg = True, mode_Eq = True, basic = True),
        MeasureFullType('FP', 'ap', 'Share of false positive symbols', 'float', mode_Bs = True, mode_En = True, mode_Bn = True, basic = True),
        MeasureFullType('FP', 'ap', 'Ratio of false positive symbols gross to the overall symbol count', 'float', mode_Bg = True, mode_Bn = True, basic = True),
        MeasureFullType('AP', 'ap', 'Share of symbols absent in the first annotation but present in the second', 'float', mode_Bs = True, mode_Eq = True, basic = True),
        MeasureFullType('AP', 'ap', 'Ratio of symbols gross absent in the first annotation but present in the second to the overall symbol count', 'float', mode_Bg = True, mode_Eq = True, basic = True),
        MeasureFullType('FN', 'pa', 'Share of false negative symbols', 'float', mode_Bs = True, mode_En = True, mode_Bn = True, basic = True),
        MeasureFullType('FN', 'pa', 'Ratio of false negative symbols gross to the overall symbol count', 'float', mode_Bg = True, mode_Bn = True, basic = True),
        MeasureFullType('PA', 'pa', 'Share of symbols present in the first annotation but absent in the second', 'float', mode_Bs = True, mode_Eq = True, basic = True),
        MeasureFullType('PA', 'pa', 'Ratio of symbols gross present in the first annotation but absent in the second to the overall symbol count', 'float', mode_Bg = True, mode_Eq = True, basic = True),
        MeasureFullType('TN', 'aa', 'Share of true negative symbols', 'float', mode_Bs = True, mode_Bg = True, mode_En = True, mode_Bn = True, basic = True),
        MeasureFullType('AA', 'aa', 'Share of symbols absent in both annotations', 'float', mode_Bs = True, mode_Bg = True, mode_Eq = True, basic = True),
        MeasureFullType('EE', 'pp', 'Share of symbols enriched in both annotations', 'float', mode_En = True, mode_Eq = True, basic = True),
        MeasureFullType('XE1', 'pa', 'Share of symbols enriched exclusively in the first annotation', 'float', mode_En = True, mode_Eq = True, basic = True),
        MeasureFullType('XE2', 'ap', 'Share of symbols enriched exclusively in the second annotation', 'float', mode_En = True, mode_Eq = True, basic = True),
        MeasureFullType('NE', 'aa', 'Share of symbols enriched in neither annotation', 'float', mode_En = True, mode_Eq = True),
        MeasureFullType('TPR', 'rc2', 'Symbol-wise true positive rate (a.k.a. sensitivity or recall)', 'float', mode_Bs = True, mode_Bg = True, mode_En = True, mode_Bn = True),
        MeasureFullType('PPV', 'pr2', 'Symbol-wise positive predictive value (a.k.a. precision)', 'float', mode_Bs = True, mode_Bg = True, mode_En = True, mode_Bn = True),
        MeasureFullType('SPC', 'sp2', 'Symbol-wise specificity', 'float', mode_Bs = True, mode_Bg = True, mode_En = True, mode_Bn = True),
        MeasureFullType('NPV', 'npv2', 'Symbol-wise negative predictive value', 'float', mode_Bs = True, mode_Bg = True, mode_En = True, mode_Bn = True),
        MeasureFullType('I-ness.', 'in2', 'Symbol-wise informedness', 'float', mode_Bs = True, mode_Bg = True, mode_En = True, mode_Bn = True),
        MeasureFullType('M-ness.', 'mk2', 'Symbol-wise markedness', 'float', mode_Bs = True, mode_Bg = True, mode_En = True, mode_Bn = True),
        MeasureFullType('PC', 'pc', 'Symbol-wise performance coefficient', 'float', mode_Bs = True, mode_Bg = True, mode_En = True, mode_Bn = True),
        MeasureFullType('ACC', 'acc', 'Symbol-wise accuracy', 'float', mode_Bs = True, mode_En = True, mode_Eq = True, mode_Bn = True),
        MeasureFullType('MCC', 'mcc', 'Symbol-wise Matthews correlation coefficient', 'float', mode_Bs = True, mode_En = True, mode_Eq = True, mode_Bn = True),
        MeasureFullType('F1', 'f1', 'Symbol-wise F1 score', 'float', mode_Bs = True, mode_Bg = True, mode_En = True, mode_Eq = True, mode_Bn = True),
        MeasureFullType('EAC', 'eac', 'Enrichment asymmetry coefficient', 'float', mode_En = True, mode_Eq = True, mode_Bn = True),
        MeasureFullType('SiteN1', 'site_n1', 'Number of sites in the first annotation', 'int', mode_Bs = True, mode_Bg = True, mode_Eq = True, force_avg = True),
        MeasureFullType('SiteN2', 'site_n2', 'Number of sites in the second annotation', 'int', mode_Bs = True, mode_Bg = True, mode_Eq = True, force_avg = True),
        MeasureFullType('SiteNB', 'site_n1', 'Number of sites in the benchmark annotation', 'int', mode_Bs = True, mode_Bg = True, mode_Bn = True, force_avg = True),
        MeasureFullType('SiteNP', 'site_n2', 'Number of sites in the prediction annotation', 'int', mode_Bs = True, mode_Bg = True, mode_Bn = True, force_avg = True),
        MeasureFullType('SiteN1m', 'site_n1_m', 'Number of matched sites in the first annotation', 'int', mode_Bs = True, mode_Bg = True, mode_Eq = True, force_avg = True),
        MeasureFullType('SiteN2m', 'site_n2_m', 'Number of matched sites in the second annotation', 'int', mode_Bs = True, mode_Bg = True, mode_Eq = True, force_avg = True),
        MeasureFullType('SiteNBm', 'site_n1_m', 'Number of matched sites in the benchmark annotation', 'int', mode_Bs = True, mode_Bg = True, mode_Bn = True, force_avg = True),
        MeasureFullType('SiteNPm', 'site_n2_m', 'Number of matched sites in the prediction annotation', 'int', mode_Bs = True, mode_Bg = True, mode_Bn = True, force_avg = True),
        MeasureFullType('SiteL1', 'site_len1', 'Total length of sites in the first annotation', 'int', mode_Bs = True, mode_Bg = True, mode_Eq = True, force_avg = True),
        MeasureFullType('SiteL2', 'site_len2', 'Total length of sites in the second annotation', 'int', mode_Bs = True, mode_Bg = True, mode_Eq = True, force_avg = True),
        MeasureFullType('SiteLB', 'site_len1', 'Total length of sites in the benchmark annotation', 'int', mode_Bs = True, mode_Bg = True, mode_Bn = True, force_avg = True),
        MeasureFullType('SiteLP', 'site_len2', 'Total length of sites in the prediction annotation', 'int', mode_Bs = True, mode_Bg = True, mode_Bn = True, force_avg = True),
        MeasureFullType('SiteTPR', 'site_rc2', 'Site-wise true positive rate (a.k.a. sensitivity or recall)', 'float', mode_Bs = True, mode_Bg = True, mode_Bn = True),
        MeasureFullType('SitePPV', 'site_pr2', 'Site-wise positive predictive value (a.k.a. precision)', 'float', mode_Bs = True, mode_Bg = True, mode_Bn = True),
        MeasureFullType('SitePC', 'site_pc2', 'Site-wise performance coefficient', 'float', mode_Bs = True, mode_Bg = True, mode_Bn = True),
        MeasureFullType('SiteF1', 'site_f1', 'Site-wise F1 score', 'float', mode_Bs = True, mode_Bg = True, mode_Eq = True, mode_Bn = True),
        MeasureFullType('SitePCV', 'site_pcv', 'Site-wise positive correlation value', 'float', mode_Bs = True, mode_Bg = True, mode_Eq = True, mode_Bn = True)
    )
    def __init__(self, enrichment_count, benchmark, gross):
        self.name_map = []
        for measure in PerformanceMeasures.name_map_full:
            if enrichment_count and (not measure.mode_En):
                continue
            if (not enrichment_count) and (not (measure.mode_Bs or measure.mode_Bg)):
                continue
            if benchmark and (not measure.mode_Bn):
                continue
            if (not benchmark) and (not measure.mode_Eq):
                continue
            if gross and (not measure.mode_Bg):
                continue
            if (not gross) and (not (measure.mode_Bs or measure.mode_En)):
                continue
            self.name_map.append(MeasureType(measure.displayed_name, measure.var_name, measure.description, measure.type_, measure.force_avg, measure.basic))
        self._attr_list = []
        for measure in self.name_map:
            if enrichment_count:
                measure.var_name = 'e_' + measure.var_name
            setattr(self, measure.var_name, [float('nan') if measure.type_ == 'float' else 0, 0])
            self._attr_list.append(measure.var_name)
    def __iadd__(self, other):
        for attr in self._attr_list:
            content = getattr(self, attr)
            content_ = getattr(other, attr)
            if not math.isnan(content_[0]):
                if math.isnan(content[0]):
                    content[0] = content_[0]
                else:
                    content[0] += content_[0]
                content[1] += content_[1]
        return self
    def set_value(self, attr, value):
        """Method to set value of a measure"""
        getattr(self, attr)[0] = value
        getattr(self, attr)[1] = 1
    def set_count(self, attr, count):
        """Method to set count for a measure"""
        getattr(self, attr)[1] = count
    def get_value(self, attr):
        """Method to get value of a measure"""
        return getattr(self, attr)[0]
    def get_count(self, attr):
        """Method to get count for a measure"""
        return getattr(self, attr)[1]
    def average(self, divisor):
        """Method for averaging all the non-integer values on the basis of respective counts"""
        for measure in self.name_map:
            if measure.type_ == 'int':
                continue
            value = self.get_value(measure.var_name)
            if measure.basic:
                self.set_value(measure.var_name, value / self.get_count(measure.var_name))
            elif not math.isnan(value):
                self.set_value(measure.var_name, value / (divisor if divisor else self.get_count(measure.var_name)))
            elif divisor:
                self.set_value(measure.var_name, 0.0)

class FileHandlers:
    """Class to hold the handlers of the output annotation files"""
    output_file_types = ('union', 'intersection', 'complement1', 'complement2', 're1', 're2', 'detailed', 'site')
    def __init__(self):
        for attr_name in FileHandlers.output_file_types:
            setattr(self, attr_name, None)