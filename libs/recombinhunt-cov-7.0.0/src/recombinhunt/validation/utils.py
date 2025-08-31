import json
from collections import Counter
from itertools import zip_longest
from logging import warning
from typing import List
import inflect
import pandas as pd
from recombinhunt.core.environment import PangoLineageHierarchy


infl = inflect.engine()


class AssessedContributingLin:
    def __init__(self, alias_key_file_path):
        with open(alias_key_file_path, "r") as alias_key_file:
            self.alias_key = json.load(alias_key_file)
        # remove mappings to ""
        self.alias_key = {x:y for x,y in self.alias_key.items() if isinstance(y, list)}
        self.fixes_to_original_file()

    def contributing_to(self, lineage: str) -> list:
        return self.alias_key.get(lineage, [lineage])

    def fixes_to_original_file(self):
        # corrections to the alias_key.josn file suggested by Pango issues
        replacements = dict()
        replacements['XA'] = ["B.1.177", "B.1.1.7"]
        replacements['XD'] = ["B.1.617.2*", "BA.1*", "B.1.617.2*"]
        replacements['XAW'] = ["AY.122", "BA.2*", "AY.122"]
        replacements['XBU'] = ["BA.2.75*", "BQ.1*", "BA.2.75*"]
        replacements['XBD'] = ["BA.2.75.2","BA.5.2.1"]
        replacements['XBE'] = ['BA.5.2*', 'BE.4.1']
        replacements['XBF'] = ['BA.5.2', 'CJ.1']
        replacements['XBJ'] = ["BA.2.3.20","BA.5.2*"]
        if len(replacements):
            warning(f"Corrections applied to alias_key.json suggested by Pango issues:")
        for k, v in replacements.items():
            warning(f"{k} {self.alias_key.get(k, '?')} -> {replacements[k]}")
        self.alias_key.update(replacements)


class BreakpointsLocation:
    """
    Breakpoint locations from ground truth (PANGO issues). The safest choice is to consider the coordinates from source
    as edge-included, i.e., giving the interpretation that covers the widest possible region.
    For reasons of consistency with the rest of the code, we convert to the edge-exclude notation
    by adding -1 to the left and + 1 to the right. It's just a matter of encoding, so the breakpoint region doesn't
    change if the convention is respected,
    All the breakpoint regions obtained from this class methods are therefore edge-excluded.
    """

    breakpoints_location_map = {
        # CONVERT ORIGINAL COORDINATES TO EDGE-EXCLUDED
        k: tuple((br_pos_start - 1, br_pos_end + 1) for br_pos_start, br_pos_end in v)  # convert to edge-excluded
        for k, v in {
            # BELOW: COORDINATES AS FROM THE SOURCE (CONSIDERED EDGE-INCLUDED)
            "XA": ((21255, 21764),),
            "XB": ((21000, 22500),),
            "XC": ((26000, 26001),),
            "XD": ((22076, 22235), (25000, 25480)),
            "XE": ((10448, 11287),),
            "XF": ((5386, 6511),),
            "XG": ((5927, 6511),),
            "XH": ((10448, 11287),),
            "XJ": ((13200, 17400),),
            "XK": ((18163, 18164),),
            "XL": ((6518, 8392),),
            "XM": ((17410, 19995),),
            "XN": ((2834, 4183),),
            "XP": ((27385, 29509),),
            "XQ": ((4322, 5385),),
            "XR": ((4322, 4891),),
            "XS": ((9055, 10447),),
            "XT": ((26062, 26528),),
            "XU": ((6518, 9343),),
            "XV": ((13196, 15713),),
            "XW": ((2834, 4183),),
            "XY": ((11540, 12879),),
            "XZ": ((26062, 26250),),
            "XAA": ((8939, 9343),),
            "XAB": ((6516, 8392),),
            "XAC": ((25813, 26059), (27385, 29509)),
            "XAD": ((26063, 26528),),
            "XAE": ((24506, 26048),),
            "XAF": ((10448, 11287),),
            "XAG": ((6516, 8392),),
            "XAH": ((26859, 27381),),
            "XAJ": ((1, 2),),   # should be 2BP but it is unclear where they are
            "XAK": ((13195, 15240), (21618, 21762)),
            "XAL": ((17413, 19953),),
            "XAM": ((8088, 9191),),     # approximated 4,2cm - 5,7cm- over 6,8cm (==5K-10K) in total
            "XAN": ((17823, 21765),),   # accept most distant values
            "XAP": ((26063, 26528),),
            "XAQ": ((18166, 19953),),
            "XAR": ((2834, 4183),),
            "XAS": ((23040, 27787),),   # accept most distant values
            "XAT": ((26062, 26528),),
            "XAU": ((2835, 4183),),
            "XAV": ((15960, 17278),),
            "XAW": ((22036, 22587), (28271, 28311)),
            "XAY": ((1,2),(2,3), (3,4), (4,5)),   # 3/4BP unclear breakpoints. Combination of three different lineages.
            "XAZ": ((3359, 9865), (27385, 27386)),
            # XBA   # XBA == XAY, 3/4 BP
            "XBB": ((22892, 22934),),   # approximated by translating AA positions
            "XBC": ((2790, 4184), (22578, 22674), (25000, 25584)), # 4BP missing one
            "XBD": ((23013, 24615),),   # approximated by translating AA positions
            "XBE": ((22592, 23608),),   # approximated by translating AA positions (and parsing text)
            "XBF": ((9865, 9866),),     # approximated (just one point)
            "XBG": ((22601, 22915),),
            "XBH": ((15451, 22000),),
            "XBJ": ((23015, 25809),),
            "XBK": ((1, 2),),   # 1BP unclear breakpoint
            "XBL": ((406, 3795), (5184, 12443)),
            "XBM": ((22601, 22915),),
            "XBN": ((1, 2), ), # 1BP unclear breakpoint
            "XBP": ((22193, 22330),),
            "XBQ":   ((1,2),), # 1BP unclear breakpoint and unclear parent lineages
            "XBR": ((22034, 22189),),
            "XBS":   ((1,2),), # 1BP unknown breakpoint
            "XBT": ((5183, 9766), (22577, 22898)),
            "XBU": ((22577, 22893), (25416, 26275)),
            "XBV":   ((1,2),), # 1BP unknonw breakpoint
            "XBW": ((25417, 26274),),
            "XCA": ((1,2),) # 1BP unknown breakpoint
        }.items()
    }

    @staticmethod
    def known_lineages(breakpoints_num):
        return [k for k, v in BreakpointsLocation.breakpoints_location_map.items() if len(v) == breakpoints_num]

    @staticmethod
    def get(lineage, ith=None):
        try:
            br = BreakpointsLocation.breakpoints_location_map[lineage]
            if ith is None:
                return br
            else:
                try:
                    return br[ith]
                except IndexError:
                    return None
        except KeyError:
            return None

    @staticmethod
    def breakpoints_num(lineage):
        try:
            return len(BreakpointsLocation.breakpoints_location_map[lineage])
        except KeyError:
            return 0

    @staticmethod
    def all_breakpoints(lineage):
        return BreakpointsLocation.breakpoints_location_map.get(lineage)

    @staticmethod
    def ith_breakpoint(lineage, idx: int):
        try:
            return BreakpointsLocation.breakpoints_location_map[lineage][idx]
        except IndexError:
            return None
        except KeyError:
            return None

    @staticmethod
    def to_target_pos(breakpoint: tuple, target: list):
        """
        Converts a pair of breakpoint coordinates in 0-based target position.
        Breakpoint positions are intended as edge-excluded, i.e. with the breakpoint region
        standing in between the extremes and excluding both edges.
        Examples of conversion:

        - breakpoint within mutations 3 and 5 (3<br<5) of t -> returns indexes 3,5

        - breakpoint after last mutation of t (mutation x) -> returns indexes x,x+1

        - breakpoint before first mutation of t -> returns 0,0

        :param breakpoint:
        :param target:
        :return:
        """
        seq_change_pos = sorted([int(change.split('_')[0]) for change in target])
        t_pos_before_br = len([c for c in seq_change_pos if c <= breakpoint[0]]) - 1
        t_pos_before_br = max(0, t_pos_before_br)

        # if breakpoint is after last change of t (nth change), it returns
        t_pos_after_br = len(seq_change_pos) - len([c for c in seq_change_pos if c >= breakpoint[1]])

        return t_pos_before_br, t_pos_after_br

    @staticmethod
    def is_unknown(gt_genomic_breakpoints = None, lineage = None):
        if (gt_genomic_breakpoints is None and lineage is None) or \
                (gt_genomic_breakpoints is not None and lineage is not None):
            raise ValueError("Specify exactly one between gt_genomic_breakpoints or lineage arguments")
        elif lineage is not None:
            return (0,3) in BreakpointsLocation.all_breakpoints(lineage)
        else:
            return (0,3) in gt_genomic_breakpoints



class RhBreakpointsLocation:
    """
    Breakpoint locations from Recombinhunt experiements on consensus sequence. Breakpoints are given
    as edge-excluded 0-based positions on the target mutations of the consensus sequence.
    """

    breakpoints_location_map = {
            # BELOW: BR-POSITIONS FROM RECOMBINHUNT ON TARGET CONSENSUS
            "XA": ((12, 13),),
            # "XB": ((, ),),
            "XC": ((27, 28),),
            "XD": ((22, 23), (51, 52)),
            "XE": ((8, 9),),
            "XF": ((5, 6),),
            "XG": ((5, 6),),
            "XH": ((9, 10),),
            "XJ": ((12, 13),),
            "XK": ((14, 15),),
            "XL": ((7, 8),),
            "XM": ((13, 14),),
            # "XN": ((, ),),
            # "XP": ((, ),),
            "XQ": ((3, 4),),
            "XR": ((3, 4),),
            "XS": ((9, 10),),
            "XT": ((52, 53),),
            "XU": ((4, 5),),
            "XV": ((11, 12),),
            "XW": ((2, 3),),
            "XY": ((13, 14),),
            "XZ": ((55, 56),),
            "XAA": ((6, 7),),
            "XAB": ((5, 6),),
            "XAC": ((55, 56), (61, 62)),
            "XAD": ((54, 55),),
            "XAE": ((54, 55),),
            "XAF": ((8, 9),),
            "XAG": ((5, 6),),
            "XAH": ((55, 56),),
            # "XAK": ((, ), (, )),
            "XAL": ((14, 15),),
            "XAM": ((5, 6),),
            "XAN": ((10, 11),),
            "XAP": ((54, 55),),
            "XAQ": ((15, 16),),
            # "XAR": ((, ),),
            "XAS": ((60, 61),),
            # "XAT": ((, ),),
            "XAU": ((2, 3),),
            # "XAV": ((, ),),
            "XAW": ((50, 51), (96, 97)),
            "XAZ": ((2, 3),),
            "XBB": ((50, 51),),
            # "XBC": ((, ), (, ), (, )),
            "XBD": ((64, 65),),
            "XBE": ((62, 63),),
            # "XBF": ((, ),),
            "XBG": ((39, 40),),
            "XBH": ((3, 4), (18, 19)),
            "XBJ": ((71, 72),),
            "XBL": ((4, 5), (11, 12)),
            "XBM": ((5, 6), (32, 33)),
            # "XBT": ((, ), (, )),
            # "XBU": ((, ), (, ))
    }

    @staticmethod
    def breakpoints_num(lineage):
        try:
            return len(RhBreakpointsLocation.breakpoints_location_map[lineage])
        except KeyError:
            return 0

    @staticmethod
    def all_breakpoints(lineage):
        return RhBreakpointsLocation.breakpoints_location_map.get(lineage)

    @staticmethod
    def ith_breakpoint(lineage, idx: int):
        try:
            return RhBreakpointsLocation.breakpoints_location_map[lineage][idx]
        except IndexError:
            return None
        except KeyError:
            return None
        
def compute_X_perc_characterization(lists=None, strings=None, threshold=0.5) -> list:
    if lists is not None:
        assert isinstance(lists, list), "Malformed input"
        assert isinstance(lists[0], list), "Malformed input"
        sequences = lists
    elif strings is not None:
        assert isinstance(strings, list), "Malformed input"
        assert isinstance(strings[0], str), "Malformed input"
        sequences = [n.split(',') for n in strings]
    else:
        raise ValueError("One of arguments lists or strings must be initialized")

    change_frequency = Counter()
    for seq in sequences:
        change_frequency.update(seq)

    threshold_num = threshold * len(sequences)

    above_threshold = [change for change, counter in change_frequency.items() if counter >= threshold_num]
    return sorted(above_threshold, key=lambda x: int(x.split('_')[0]))

def compute_75_perc_characterization(lists=None, strings=None) -> list:
    if lists is not None:
        assert isinstance(lists, list), "Malformed input"
        assert isinstance(lists[0], list), "Malformed input"
        sequences = lists
    elif strings is not None:
        assert isinstance(strings, list), "Malformed input"
        assert isinstance(strings[0], str), "Malformed input"
        sequences = [n.split(',') for n in strings]
    else:
        raise ValueError("One of arguments lists or strings must be initialized")

    change_frequency = Counter()
    for seq in sequences:
        change_frequency.update(seq)

    threshold_num = 0.75 * len(sequences)

    above_threshold = [change for change, counter in change_frequency.items() if counter >= threshold_num]
    return sorted(above_threshold, key=lambda x: int(x.split('_')[0]))

def all_candidates_matching(regions: list, contributing_lineages: list, lh: PangoLineageHierarchy) -> bool:
    if len(regions) != len(contributing_lineages):
        return False
    else:
        OK_KO_values = []
        for reg, gt_lin in zip_longest(regions, contributing_lineages):
            BC_evaluation = lh.hierarchy_distance(reg.designated, gt_lin)
            if BC_evaluation == 0:  # if candidate == gt skip checks on alternative candidates
                OK_KO_values.append(0)
            else:
                alt_c_evaluation = [lh.hierarchy_distance(c, gt_lin) for c in
                                    reg.get_master_region().get_good_alternative_candidates()]
                OK_KO_values.append(min([BC_evaluation, *alt_c_evaluation]))
        if all([ev <= 1 for ev in OK_KO_values]):  # sub-lineages / equal
            return True
        elif any([ev >= 2 for ev in OK_KO_values]):  # super-lineages / different
            return False

def candidates_matching(regions: list, contributing_lineages: list, lh: PangoLineageHierarchy) -> List[bool]:
    OK_KO_values = []
    for r_idx, (gt_lin, reg) in enumerate(zip_longest(contributing_lineages, regions)):
        if reg is None or gt_lin is None:
            OK_KO_values.append(3)
            continue
        # get master region for double-breakpoints
        if len(regions) == 3 and r_idx % 2 == 0:
            reg = reg.master_region()
        BC_evaluation = lh.hierarchy_distance(reg.designated, gt_lin)
        if BC_evaluation == 0:  # if candidate == gt skip checks on alternative candidates
            OK_KO_values.append(0)
        else:
            alt_c_evaluation = [lh.hierarchy_distance(c, gt_lin)
                                for c in reg.get_good_alternative_candidates()]
            OK_KO_values.append(min([BC_evaluation, *alt_c_evaluation]))
    return [True if ev <= 1 else False for ev in OK_KO_values]

def rank_gt_in_regions(regions: list, contributing_lineages: list, lh: PangoLineageHierarchy) -> list:
    ranks = []
    for reg, gt_lin in zip_longest(regions, contributing_lineages):
        if gt_lin is None or reg is None:
            ranks.append("-")
        else:
            try:
                equivalent_lin = [True if lh.is_same_lineage(c, gt_lin) else False for c in reg.candidates]
                ranks.append(str(equivalent_lin.index(True) + 1))
            except ValueError:
                ranks.append("11")
    return ranks

def show_batch_sequence_private_mutations(environment, consensus_seq: list, test_sequence_dataset: list) -> pd.DataFrame:
        def filter_private_mutations(all_mut_list):
            return [f"{m} {environment.change_probability(m):.2e}" for m in all_mut_list if m not in consensus_seq]

        def sorting_key(mutation):
            return int(mutation.split("_")[0])

        def sorted_mutations(all_mut_list):
            return sorted(all_mut_list, key=sorting_key)


        private_mut_series = pd.Series(index=[x[0] for x in test_sequence_dataset],
                                       data=[x[2].split(',') for x in test_sequence_dataset], name="private mutations")
        private_mut_series = private_mut_series.apply(sorted_mutations)
        private_mut_series = private_mut_series.apply(filter_private_mutations)
        return pd.DataFrame(private_mut_series)

def hot_reload():
    print('Starting...')