import logging
import math
from typing import Optional, Tuple, List
import numpy as np
from pprint import pformat
from copy import deepcopy
from enum import Enum, auto


class NoCandidatesFound(Exception):
    pass
    # can be raised by search_L, search_L_fixed_direction
    # handled in model_1BP, model_2BP, but not in search_L for searching L1


class SingleCandidateGenome(Exception):
    pass


class BadOppositeRegion(Exception):
    pass


class Region:
    def __init__(self, merged_df):
        self.seq_change = merged_df.seq_change
        self.merg_pos = merged_df.merg_pos
        self.num_genome_positions = merged_df.shape[0]

        # region boundaries
        self.pos_start = None   # left end included
        self.pos_end = None     # right end excluded
        self.t_pos_start = None     # left end included
        self.t_pos_end = None       # right end excluded
        self.genomic_start = None   # left end included
        self.genomic_end = None     # right end excluded

        # candidates likelihood
        self.logp_values = None
        self.max_likelihood_value = None  # list of maximum likelihood values of all candidates
        self.pos_max = None  # list of positions of maximum likelihood of all candidates
        self.pos_max_t = None  # list of positions of maximum likelihood of all candidates in target coord.

        # candidates
        self.designated = None
        self.candidates = None  # list of all top-likelihood candidates (including designated)
        self.good_alternative_candidates = None

        # 4 filtering candidates on p_value
        self.candidates_similarity_stats = None  # dict{candidate -> [likelihood in pos_max designated,aic,p-value/None]

        # 4 filtering candidates on phylogenetic analysis
        self._alternative_candidates_in_same_branch = None  # list of alt. cand. on sane phylogen. branch as designated
        self._alternative_candidates_in_same_tree = None    # list of alt. cand. on sane phylogen. tree as designated

        # misc
        self.search_dir = None              # 0 is >> , 1 is <<
        self.master_region = None           # pointer to region connected to this one in case of 2BP

    def __str__(self):
        output = pformat(self.describe(), sort_dicts=False)
        chars_to_remove = ["{", "}", "'", '"', ","]
        for char in chars_to_remove:
            output = output.replace(char, " ")
        return output

    def describe(self):
        output = {
            "pos_start_in_t": self.t_pos_start + 1, # +1 converts to 1-based edge included
            "pos_end_in_t": self.t_pos_end,         # +0 converts to 1-based edge included
            "designated": self.designated,
            "good alternative candidates": self.good_alternative_candidates
            #"search_direction": self.search_dir
            #"all_candidates": self.candidates[0] if self.candidates is not None else None
        }
        return output

    def __repr__(self):
        return pformat({
            "designated": self.designated,
            "candidates": self.candidates,
            "candidates_similarity_stats": self.candidates_similarity_stats,
            "pos_start": self.pos_start,
            "t_pos_start": self.t_pos_start,
            "t_pos_end": self.t_pos_end,
            "pos_end": self.pos_end,
            "search_dir": self.search_dir,
            "logp_values": self.logp_values
        }, sort_dicts=False)

    # CONSTRUCTORS
    @staticmethod
    def from_candidate_list(candidates, direction, merged_df, l_edge_pos=None, r_edge_pos=None):
        obj = Region(merged_df)
        names, max_likelihood_v, max_pos, logp_values = candidates
        if direction == 0:
            obj.set_pos_start(l_edge_pos or 0)
            obj.set_pos_end(max_pos[0] + 1)
        else:
            obj.set_pos_start(max_pos[0])
            obj.set_pos_end(r_edge_pos or obj.num_genome_positions)
        obj.designated = names[0]
        obj.candidates = names
        obj.pos_max = max_pos
        obj.pos_max_t = [Region.convert_pos_to_t(p, obj.seq_change) for p in max_pos]
        obj.max_likelihood_value = max_likelihood_v
        obj.search_dir = direction
        obj.logp_values = logp_values

        return obj

    @staticmethod
    def from_single_candidate(name, max_l, max_pos, logp_values, direction, merged_df, l_edge_pos=None, r_edge_pos=None):
        obj = Region(merged_df)
        if direction == 0:
            obj.set_pos_start(l_edge_pos or 0)
            obj.set_pos_end(max_pos + 1)
        else:
            obj.set_pos_start(max_pos)
            obj.set_pos_end(r_edge_pos or obj.num_genome_positions)
        obj.designated = name
        obj.candidates = [name]
        obj.pos_max = [max_pos]
        obj.pos_max_t = [Region.convert_pos_to_t(max_pos, obj.seq_change)]
        obj.max_likelihood_value = [max_l]
        obj.search_dir = direction
        obj.logp_values = logp_values.reshape(len(logp_values),1)
        return obj

    def set_master_region(self, reg):
        """
        Make one region dependent on another. This region behaviour is not directly affected by this method. However,
        users can recall the master region to apply special logics. If not explicityl set, the master region is self.
        :param reg:
        :return:
        """
        self.master_region = reg

    def get_master_region(self):
        """
        Normally, every region is independent of other regions and self is returned upon call. However, when this
        region depends on another, the other region is returned. Region bounds are set explicilty through
        set_master_region().
        :return: self or the depending Region.
        """
        return self.master_region or self

    #   REGION BOUNDARIES
    def mask_where_occupied(self):
        region_mask = np.zeros(self.num_genome_positions, dtype=bool)
        region_mask[self.pos_start:self.pos_end] = True
        return region_mask

    def extend(self):
        self.set_pos_start(0)
        self.set_pos_end(self.num_genome_positions)

    def length_in_t(self):
        return np.sum(self.seq_change[self.pos_start:self.pos_end])

    def set_pos_start(self, new_pos):
        self.pos_start = new_pos

        # update pos in t
        t_change_pos_idx_start = 0  # smallest index of a change of t located >= pos_start (if not found is 0) as position of merged_df
        for i in range(self.pos_start, self.num_genome_positions):
            if self.seq_change.iloc[i]:
                t_change_pos_idx_start = i
                break
        self.t_pos_start = np.sum(
            self.seq_change[:t_change_pos_idx_start])  # smallest index of change of t >= pos_start as position of t
        self.t_pos_start = Region.convert_pos_to_t(self.pos_start, self.seq_change)

        self.genomic_start = self.merg_pos.iloc[new_pos]

    def set_pos_end(self, new_pos):
        self.pos_end = new_pos

        # update pos in t
        t_change_pos_idx_stop = self.num_genome_positions  # smallest index of a change of t located >= pos_end (if not found is num_genome_positions) as position of merged_df
        for j in range(self.pos_end, self.num_genome_positions):
            if self.seq_change.iloc[j]:
                t_change_pos_idx_stop = j
                break
        self.t_pos_end = np.sum(
            self.seq_change[:t_change_pos_idx_stop])  # smallest index of change of t >= pos_end as position of t

        if new_pos == self.num_genome_positions:
            self.genomic_end = self.merg_pos.iloc[self.num_genome_positions-1] + 1
        else:
            self.genomic_end = self.merg_pos.iloc[new_pos]

    # CANDIDATES LIKELIHOOD
    def likelihood_values(self, alternative_candidate=None):
        idx = self.__index_of_candidate(alternative_candidate)
        v = self.logp_values[:, idx]
        return v.reshape(self.num_genome_positions)

    def candidate_max_likelihood(self, alternative_candidate=None):
        candidate_idx = self.__index_of_candidate(alternative_candidate)
        return self.max_likelihood_value[candidate_idx]

    def candidate_max_likelihood_pos_t(self, candidate=None, include_master_region=False):
        try:
            index = self.__index_of_candidate(candidate)
            return self.pos_max_t[index]
        except IndexError as e:
            if include_master_region and self.get_master_region() is not None:
                return self.get_master_region().candidate_max_likelihood_pos_t(candidate)
            else:
                raise e

    def candidate_likelihood_at_dir_end(self, alternative_candidate=None):
        self.__check_is_valid_candidate(alternative_candidate)
        candidate = alternative_candidate or self.designated
        return self.candidates_similarity_stats[candidate][0]

    # 4 FILTERING CANDIDATES ON MAX LIKELIHOOD POS
    def candidate_max_likelihood_pos_t_distance2designated(self, candidate, include_master_region=False):
        try:
            index = self.__index_of_candidate(candidate)
            return self.pos_max_t[0] - self.pos_max_t[index]
        except IndexError as e:
            if include_master_region and self.get_master_region() is not None:
                return self.get_master_region().candidate_max_likelihood_pos_t_distance2designated(candidate)
            else:
                raise e

    def alternative_candidates_with_max_likelihood_pos_t_distance_below(self, threshold_value):
        reference_region = self  # self.master_region or self
        return [c for c in reference_region.alternative_candidates()
                if abs(reference_region.candidate_max_likelihood_pos_t_distance2designated(c)) <= threshold_value]

    # 4 FILTERING CANDIDATES ON P-VALUE
    def candidate_aic(self, alternative_candidate):
        self.__check_is_valid_candidate(alternative_candidate)
        candidate = alternative_candidate or self.designated
        return self.candidates_similarity_stats[candidate][1]

    def candidate_p_value(self, alternative_candidate):
        self.__check_is_valid_candidate(alternative_candidate)
        candidate = alternative_candidate or self.designated
        return self.candidates_similarity_stats[candidate][2]

    def alternative_candidates_with_p_value_above(self, threshold_value):
        reference_region = self#self.master_region or self
        return [c for c in reference_region.alternative_candidates() if reference_region.candidate_p_value(c) >= threshold_value]

    # 4 FILTERING CANDIDATES ON PHYLOGENETIC ANALYSIS
    def save_alternative_candidates_in_same_branch(self, lineage_hierarchy):
        #reference_region = self  # self.master_region or self
        self._alternative_candidates_in_same_branch = lineage_hierarchy.filter_same_hierarchy_as_first(self.candidates)[1:]

    def alternative_candidates_in_same_branch(self, lineage_hierarchy=None):
        if self._alternative_candidates_in_same_branch is None and not lineage_hierarchy:
            raise ValueError("Non initialized variable. First, either call save_alternative_candidates_in_same_branch()"
                             " or pass lineage_hierarchy as argument")
        elif lineage_hierarchy and self._alternative_candidates_in_same_branch is None:
            self.save_alternative_candidates_in_same_branch(lineage_hierarchy)
        return self._alternative_candidates_in_same_branch

    def save_alternative_candidates_in_same_tree(self, lineage_hierarchy):
        #reference_region = self  # self.master_region or self
        ancestor = lineage_hierarchy.farthest_ancestor_in_list(self.designated, self.alternative_candidates())
        self._alternative_candidates_in_same_tree = lineage_hierarchy.filter_same_hierarchy_as(self.alternative_candidates(), ancestor)

    def alternative_candidates_in_same_tree(self, lineage_hierarchy=None):
        if self._alternative_candidates_in_same_tree is None and not lineage_hierarchy:
            raise ValueError("Non initialized variable. First, either call save_alternative_candidates_in_same_branch()"
                             " or pass lineage_hierarchy as argument")
        elif lineage_hierarchy and self._alternative_candidates_in_same_tree is None:
            self.save_alternative_candidates_in_same_tree(lineage_hierarchy)
        return self._alternative_candidates_in_same_tree

    # MISC
    def alternative_candidates(self):
        return self.candidates[1:]

    def set_good_alternative_candidates(self, good_alt_c_list=None, multiple_lists=None):
        if good_alt_c_list is None and len(multiple_lists) == 0:
            raise ValueError("Function called without arguments")
        if good_alt_c_list is not None:
            if any([g not in self.alternative_candidates() for g in good_alt_c_list]):
                raise ValueError(f"Proposed 'good' alternative candidate(s) "
                                 f"{', '.join([g not in self.alternative_candidates() for g in good_alt_c_list])} "
                                 f"are not alternative candidates of this region.")
            self.good_alternative_candidates = good_alt_c_list
        else:
            self.good_alternative_candidates = [c for c in self.alternative_candidates()
                                                if all([c in fl for fl in multiple_lists])]

    def get_good_alternative_candidates(self):
        return self.good_alternative_candidates

    def __index_of_candidate(self, alternative_candidate=None):
        if not self.candidates:
            raise IndexError("Region doesn't have any assigned candidate.")
        else:
            if not alternative_candidate:
                return 0
            else:
                try:
                    return self.candidates.index(alternative_candidate)
                except IndexError:
                    raise IndexError(f"The argument ({alternative_candidate}) is not a candidate for this region. Valid"
                                     f" candidates are: {self.candidates}.")

    def __check_is_valid_candidate(self, candidate):
        if not self.candidates:
            raise ValueError("Region doesn't have any assigned candidate.")
        elif candidate and not candidate in self.candidates:
            raise ValueError(f"The argument ({candidate}) is not a candidate for this region. Valid "
                             f"candidates are: {self.candidates}.")

    @staticmethod
    def convert_pos_to_t(pos_start, seq_change: np.array):
        """
        0-based position of t preceding pos_start
        """
        # update pos in t
        t_change_pos_idx_start = 0  # smallest index of a change of t located >= pos_start (if not found is 0) as position of merged_df
        for i in range(pos_start, len(seq_change)):
            if seq_change.iloc[i]:
                t_change_pos_idx_start = i
                break
        t_pos_start = np.sum(
            seq_change[:t_change_pos_idx_start])  # smallest index of change of t >= pos_start as position of t
        return t_pos_start

# class dotdict(dict):
#     """dot.notation access to dictionary attributes"""
#     __getattr__ = dict.get
#     __setattr__ = dict.__setitem__
#     __delattr__ = dict.__delitem__


class GenomeView:
    def __init__(self, merged_df):
        self.merged_df = merged_df
        self.num_genome_positions: int = merged_df.shape[0]
        self.len_of_t = np.sum(merged_df.seq_change)
        self.regions: List[Region] = []

        # debug attributes
        self.gap_history_pos_t = []
        self.region_lengths_pre_gap_resolution = []

    def add_region(self, region):
        self.regions = sorted(self.regions + [region], key=lambda i: i.pos_start)
        self.region_lengths_pre_gap_resolution = sorted(self.region_lengths_pre_gap_resolution + [[region.t_pos_start, region.t_pos_end]], key=lambda i: i[0])
        if len(self.regions) > 1:
            self.gap_history_pos_t.append(self.breakpoints_in_t())

    def len_free_region(self) -> int:
        mask_free = np.ones(self.num_genome_positions, dtype=bool)
        for r in self.regions:
            mask_free[r.pos_start:r.pos_end] = False
        return mask_free.sum()

    def len_free_region_of_t(self) -> int:
        uncovered_genome_mask = self.merged_df.seq_change.copy()
        for reg in self.regions:
            uncovered_genome_mask[reg.pos_start:reg.pos_end] = False
        return uncovered_genome_mask.sum()

    def free_region(self) -> Optional[Tuple[int, int]]:
        mask_free = np.ones(self.num_genome_positions, dtype=bool)
        for r in self.regions:
            mask_free[r.pos_start:r.pos_end] = False
        if mask_free.sum() == 0:
            return None
        else:
            free_pos_idx = np.where(mask_free)[0]
            return free_pos_idx[0], free_pos_idx[-1] + 1

    def free_changes_of_t(self) -> list:
        """
        Returns a (possibly empty) list of free changes of target
        :return:
        """
        mask_free = np.ones(self.num_genome_positions, dtype=bool)
        for r in self.regions:
            mask_free[r.pos_start:r.pos_end] = False
        return self.merged_df[(mask_free) & (self.merged_df.seq_change)].index

    def extendable_regions(self):
        if not self.regions:
            raise ValueError("This genome is empty.")
        elif self.len_free_region() > 0:
            if len(self.regions) == 1:
                return self.regions[0]
            else:
                free_reg = self.free_region()
                if free_reg is not None:
                    l_edge, r_edge = free_reg
                    # return regions surrounding the free region
                    l_region = [r for r in self.regions if r.pos_end == l_edge][0]
                    r_region = [r for r in self.regions if r.pos_start == r_edge][0]
                    return l_region, r_region
                else:
                    raise ValueError("Current regions already cover the entire genome.")
        else:
            raise ValueError("Current regions already cover the entire genome.")

    def number_of_regions(self) -> int:
        return len(self.regions)

    def contributing_lineages(self) -> list:
        return [r.designated for r in self.regions]

    def list_good_alternative_candidates(self) -> list:
        """
        :return: a list of good alternative candidates for each unique region (unique designated candidate)
        """
        unique_regions = sorted([r.get_master_region() for r in self.regions[:2]], key=lambda r: r.pos_start)
        return [r.get_good_alternative_candidates() for r in unique_regions]

    def breakpoints(self):
        """
        Returns a list of pairs of genomic coordinates indicating the breakpoint as an interval in genomic coordinates.
        Breakpoint intervals are expressed as edges excluded; e.g. (22450,22789) indicates the breakpoint covers the
        region between 22450 and 22789, but does not include the coordinates 22450 and 22789.
        If the genome has only one region, an empty list is returned.
        :return: list
        """
        if self.number_of_breakpoints() > 0:
            return [(r1.genomic_end - 1, r2.genomic_start) for r1, r2 in zip(self.regions[:-1], self.regions[1:])]
        else:
            return []

    def breakpoints_in_t(self):
        """
        Returns a list of pairs of 0-based positions indicating the breakpoint as an interval in the positions of the
        target. Breakpoint intervals are expressed as edges excluded; e.g. (3,4) indicates the breakpoint sits in
        between the target positions 3 and 4 but does not include the positions 3 and 4.
        If the genome has only one region, an empty list is returned.
        :return: list
        """
        if self.number_of_breakpoints() > 0:
            return [(r1.t_pos_end - 1, r2.t_pos_start) for r1,r2 in zip(self.regions[:-1], self.regions[1:])]
        else:
            return []

    def number_of_breakpoints(self):
        return len(self.regions) - 1

    def positions_of_t(self):
        return np.where(self.merged_df.seq_change)[0]

    def update_true_genomic_end(self, genome_length):
        """
        If any region currently covering the rightmost part of the genome exists, the genome_length value is set as
        genomic end for that region.

        Use this method when regions in the genome has no free regions, but the maximum genomic region coordinate
        differs from the expected maximum coordinatye for this genome.
        :param genome_length:
        :return: None
        """
        for r in self.regions:
            if r.pos_end == self.num_genome_positions:
                r.genomic_end = genome_length

    def __str__(self):
        output = self.describe()
        output = pformat(output, sort_dicts=False)
        chars_to_remove = ["{", "}", "'", '"', "(", ")", "]", "[", ","]
        for char in chars_to_remove:
            output = output.replace(char, " ")
        return output

    def describe(self):
        return {
            "target length": self.len_of_t,
            "designated candidates": " + ".join([r.designated for r in self.regions]),
            "region details": [(idx + 1, r.describe()) for idx, r in enumerate(self.regions)]
        }

    def __repr__(self):
        return pformat({
            "num_genome_positions": self.num_genome_positions,
            "target length": self.len_of_t,
            "num_regions": len(self.regions)
        }, sort_dicts=False)


class DefaultParams:
    # defaults
    MIN_SEARCHABLE_REGION_LENGTH_TO = 3
    MIN_CANDIDATE_REGION_LENGTH = 3
    MIN_L2_ENCLOSED_REGION_LENGTH = 2
    ALT_CANDIDATE_P_VALUE_DIFFERENCE = 1e-05
    ALT_CANDIDATE_MAX_POS_DISTANCE_T = 1


class Experiment:

    def __init__(self, environment, candidates_hierarchy=None,
                 min_searchable_region_length=DefaultParams.MIN_SEARCHABLE_REGION_LENGTH_TO,
                 min_candidate_region_length=DefaultParams.MIN_CANDIDATE_REGION_LENGTH,
                 min_l2_enclosed_region_length=DefaultParams.MIN_L2_ENCLOSED_REGION_LENGTH,
                 alt_candidate_p_value_difference=DefaultParams.ALT_CANDIDATE_P_VALUE_DIFFERENCE,
                 alt_candidate_max_pos_distance_t=DefaultParams.ALT_CANDIDATE_MAX_POS_DISTANCE_T):
        self.env = environment
        if candidates_hierarchy is not None:
            self.lh = candidates_hierarchy
            logging.warn("Experiment argument candidates_hierarchy is kept for backward compatibility. You should provide an Environment with built in CandidatesHierarchy instead.")
        else:
            self.lh = environment.ch

        self.MIN_SEARCHABLE_REGION_LENGTH_TO = min_searchable_region_length
        self.MIN_CANDIDATE_REGION_LENGTH = min_candidate_region_length
        self.MIN_L2_ENCLOSED_REGION_LENGTH = min_l2_enclosed_region_length
        self.ALT_CANDIDATE_P_VALUE_DIFFERENCE = alt_candidate_p_value_difference
        self.ALT_CANDIDATE_MAX_POS_DISTANCE_T = alt_candidate_max_pos_distance_t

        self._reset_internal_variables()
    
    def _reset_internal_variables(self):
        self.p_merged_df = None
        self.merged_df = None
        self.change_probabilities = None
        self.genome_view = None

        self.aik = None
        self.p_values = None

        self.L1_dir = None

        self.p_val_partial_models = dict()
        self.best_model_key = None
        self.best_model_label = None
        self.discarded_model = None

        self.flags = []

        self.cache = dict()

    def __str__(self):
        experiment_output = {
            "AIK": self.aik,
            "p_values": {k1+" vs "+k2: Experiment.__format_pvalue(v) for (k1, k2), v in self.p_values.items()},
            "flags": ", ".join(self.get_flags() or "-")
        }
        experiment_output = pformat(experiment_output, sort_dicts=False)
        experiment_output = experiment_output.replace("recombinant_model", " + ".join([r.designated for r in self.genome_view.regions]))
        chars_to_remove = ["{", "}", "'", '"', ","]
        for char in chars_to_remove:
            experiment_output = experiment_output.replace(char, " ")

        output = str(self.genome_view) + "\n" + experiment_output
        return output

    @staticmethod
    def __format_pvalue(val):
        return f"{val:.2e}"

    def __repr__(self):
        return pformat({
            "aik": self.aik,
            "p_values": self.p_values
        }, sort_dicts=False)

    def set_target(self, nuc_changes: list):
        self._reset_internal_variables()
        seq_df = self.env.sequence_nuc_mutations2df(nuc_changes)
        self.merged_df = self.env.make_merged_df(seq_df)
        self.p_merged_df, self.change_probabilities = self.env.probabilities(seq_df)

    def model_1BP(self, L1_reg):
        """
        Can raise SingleCanididateGenome if no candidates for L2 are found.
        :param L1_reg:
        :return:
        """
        gen_1BP = GenomeView(self.merged_df)
        gen_1BP.add_region(L1_reg)

        # find L2 (hypothesis model is L1 >> + << L2)
        l_edge, r_edge = gen_1BP.free_region()
        L2_search_dir = abs(1 - L1_reg.search_dir)
        try:
            L2_reg = search_L_fixed_direction(self.merged_df, self.p_merged_df, self.change_probabilities,
                                              left_to_right=L2_search_dir == 0,
                                              l_edge_pos_idx=l_edge, r_edge_pos_idx=r_edge,
                                              min_target_len=self.MIN_CANDIDATE_REGION_LENGTH, cache=self.cache)
            gen_1BP.add_region(L2_reg)
        except NoCandidatesFound:
            self.flags.append(Experiment.Flags.Model_1BP_NoL2)
            raise SingleCandidateGenome

        if not self.lh.is_totally_different(L1_reg.designated, L2_reg.designated):
            self.flags.append(Experiment.Flags.Model_1BP_L1eqL2)

        # remove gap between L1 and L2
        try:
            extendable_regions = gen_1BP.extendable_regions()
        except ValueError:
            pass
        else:
            fill_gap(self.merged_df, self.p_merged_df, self.change_probabilities, *extendable_regions)

        # range of added region
        added_range = (L2_reg.pos_start, L2_reg.pos_end)
        return gen_1BP, added_range

    def model_2BP(self, L1_reg):
        """
        Can raise BadOppositeRegion if the length of the oppotiste region is < min(MIN_CANDIDATE_REGION_LENGTH
         , free_region_of_t).
         Can raise SingleCandidateGenome if there is a problem with the enclosed L2 region: no candidates are found
         for the enclosed L2 region or if it has length < MIN_CANDIDATE_REGION_LENGTH.
        :param L1_reg:
        :return:
        """
        gen_2BP = GenomeView(self.merged_df)
        gen_2BP.add_region(L1_reg)

        # find opposite region
        l_edge, r_edge = gen_2BP.free_region()
        opposite_region = find_dual_region(self.merged_df, self.p_merged_df, self.change_probabilities, L1_reg,
                                           l_edge, r_edge, cache=self.cache)
        if opposite_region.length_in_t() < min(self.MIN_CANDIDATE_REGION_LENGTH, gen_2BP.len_free_region_of_t()):
            self.flags.append(Experiment.Flags.Model_2BP_Bad_L1_opp)
            raise BadOppositeRegion
        else:
            gen_2BP.add_region(opposite_region)

        # find L2 enclosed region (hypothesis model is L1 >> + L2 + << L1)
        l_edge, r_edge = gen_2BP.free_region()
        try:
            L2_enclosed_reg = search_L(self.merged_df, self.p_merged_df, self.change_probabilities,
                              l_edge_pos_idx=l_edge, r_edge_pos_idx=r_edge, cache=self.cache)
        except NoCandidatesFound:
            self.flags.append(Experiment.Flags.Model_2BP_NoL2)
            raise SingleCandidateGenome
        else:
            if L2_enclosed_reg.length_in_t() < self.MIN_L2_ENCLOSED_REGION_LENGTH:
                if gen_2BP.len_free_region_of_t() < self.MIN_L2_ENCLOSED_REGION_LENGTH:
                    self.flags.append(Experiment.Flags.Model_2BP_NotEnoughSpace_ForL2)
                self.flags.append(Experiment.Flags.Model_2BP_Bad_L2)
                raise SingleCandidateGenome
            else:
                gen_2BP.add_region(L2_enclosed_reg)

        if not self.lh.is_totally_different(L1_reg.designated, L2_enclosed_reg.designated):
            self.flags.append(Experiment.Flags.Model_2BP_L1eqL2)

        # remove gap between L2_enclosed and L1_5'/L1_3'
        try:
            extendable_regions = gen_2BP.extendable_regions()
        except ValueError:
            pass
        else:
            fill_gap(self.merged_df, self.p_merged_df, self.change_probabilities, *extendable_regions)

        # range of added regions
        added_range = (min([r.pos_start for r in (opposite_region, L2_enclosed_reg)]), max([r.pos_end for r in (opposite_region, L2_enclosed_reg)]))
        return gen_2BP, added_range

    def run(self, merge_when_L1eqL2=False, save_discarded_model=False):
        # find L1 (can raise NoCandidatesFound if min_target_len is > 0 && cases L1+L2+L1 with breakpoints at edges)
        L1_reg: Region = search_L(self.merged_df, self.p_merged_df, self.change_probabilities, min_target_len=self.MIN_CANDIDATE_REGION_LENGTH, cache=self.cache)
        self.L1_dir = L1_reg.search_dir

        single_candidate_genome = GenomeView(self.merged_df)
        single_candidate_genome.add_region(L1_reg)
        try:
            if single_candidate_genome.len_free_region_of_t() < self.MIN_SEARCHABLE_REGION_LENGTH_TO:
                self.flags.append(Experiment.Flags.NotEnoughSpaceAfterL1)
                raise SingleCandidateGenome
            elif all(self.merged_df.iloc[:,:-2].loc[single_candidate_genome.free_changes_of_t()].sum(axis=0) < self.MIN_CANDIDATE_REGION_LENGTH):
                self.flags.append(Experiment.Flags.NoCandidateWith3TargetMutationsAfterL1)
                raise SingleCandidateGenome
            else:
                genome_1BP = genome_2BP = None

                try:
                    genome_1BP, new_reg_1BP_range = self.model_1BP(deepcopy(L1_reg))
                except SingleCandidateGenome:
                    raise # let if pass (it cannot find neither L1 designated, neither L2)

                try:
                    genome_2BP, new_reg_2BP_range = self.model_2BP(deepcopy(L1_reg))
                except BadOppositeRegion:   # raised if L1 opposite is not a good candidate -> thus, L1_5' + ... + L1_3' cannot be found
                    genome_2BP = new_reg_2BP_range = None
                except SingleCandidateGenome:
                    raise   # let if pass (means it found L1 opposite but couldn't find enclosed L2)

                if not genome_1BP and not genome_2BP:   # never verifies (ready for future)
                    raise SingleCandidateGenome
                elif not genome_1BP and genome_2BP:     # never verifies (ready for future)
                    self.genome_view = genome_2BP
                elif genome_1BP and not genome_2BP:
                    self.genome_view = genome_1BP
                else:
                    # compare new regions (all regions except L1) on p-value
                    widest_range = list(zip(new_reg_1BP_range, new_reg_2BP_range))
                    new_reg_l_edge = min(widest_range[0])
                    new_reg_r_edge = max(widest_range[1])
                    aic_partial_models = {
                        '1BP': aic_on_range(self.merged_df, self.p_merged_df, self.change_probabilities,
                                            genome_1BP.regions, new_reg_l_edge, new_reg_r_edge, self.cache),
                        '2BP': aic_on_range(self.merged_df, self.p_merged_df, self.change_probabilities,
                                            genome_2BP.regions, new_reg_l_edge, new_reg_r_edge, self.cache)
                    }

                    self.p_val_partial_models = {
                        "1BP vs 2BP": compare_aik(aic_partial_models['1BP'], aic_partial_models['2BP']),
                        "2BP vs 1BP": compare_aik(aic_partial_models['2BP'], aic_partial_models['1BP'])
                    }
                    if self.p_val_partial_models["2BP vs 1BP"] < self.p_val_partial_models["1BP vs 2BP"]:   # L1 + L2 + L1 is better
                        best_model = genome_2BP
                        self.flags.append(Experiment.Flags.Model_2BP_Best)
                        self.best_model_key = "2BP vs 1BP"
                        if save_discarded_model:
                            self.discarded_model = genome_1BP
                    else:    # L1 + L2 is better
                        self.flags.append(Experiment.Flags.Model_1BP_Best)
                        best_model = genome_1BP
                        self.best_model_key = "1BP vs 2BP"
                        if save_discarded_model:
                            self.discarded_model = genome_2BP
                    self.genome_view = best_model

            if merge_when_L1eqL2 and designated_candidates_in_same_hierarchy_branch(self.genome_view, self.lh):
                raise SingleCandidateGenome

        except SingleCandidateGenome:
            self.flags.append(Experiment.Flags.SingleCandidateGenome)
            L1_reg.extend()  # target not recombinant
            self.genome_view = single_candidate_genome

        # evaluate alternative candidates
        self.assertion_error_alt_c = update_region_p_values(self.merged_df, self.p_merged_df, self.change_probabilities, self.genome_view.regions, self.cache)
        [r.save_alternative_candidates_in_same_branch(self.lh) for r in self.genome_view.regions]
        # candidate max pos is evaluated internally within each region
        for r in self.genome_view.regions:
            r.set_good_alternative_candidates(multiple_lists=(
                r.alternative_candidates_with_p_value_above(self.ALT_CANDIDATE_P_VALUE_DIFFERENCE),
                r.alternative_candidates_with_max_likelihood_pos_t_distance_below(self.ALT_CANDIDATE_MAX_POS_DISTANCE_T),
                r.alternative_candidates_in_same_branch()
            ))

        # p_value of the model for designated candidates
        self.aik, self.p_values = aik_p_values(self.merged_df, self.p_merged_df, self.change_probabilities, self.genome_view.regions)

        # TODO self.genome_view.update_true_genomic_end(genome_length) # TODO read gneome_length from environment

        return self

    def get_flags(self):
        return [f.name for f in self.flags]

    class Flags(Enum):
        NotEnoughSpaceAfterL1 = auto()
        NoCandidateWith3TargetMutationsAfterL1 = auto()
        SingleCandidateGenome = auto()

        Model_1BP_NoL2 = auto()
        Model_1BP_Best = auto()
        Model_1BP_L1eqL2 = auto()

        Model_2BP_Bad_L1_opp = auto()
        Model_2BP_NoL2 = auto()
        Model_2BP_Bad_L2 = auto()
        Model_2BP_Best = auto()
        Model_2BP_NotEnoughSpace_ForL2 = auto()
        Model_2BP_L1eqL2 = auto()


def designated_candidates_in_same_hierarchy_branch(gen: GenomeView, lh):
    return any([not lh.is_totally_different(r1.designated, r2.designated)
                for r1,r2 in zip(gen.regions[1:], gen.regions[:-1])])


def logp4all(p_merged_df, change_probabilities, cache=dict()):
    """
    Returns the likelihood ratio scores for all the lineages in p_merged_df. 

    The cache may be used to avoid expensive re-computations. 
    The cache may be used also to speedup the method 'logp4lin'.

    :param p_merged_df:
    :param change_probabilities:
    :param cache:
    :return:
    """
    key = ("logp4all", p_merged_df.shape)
    try:
        _logp = cache[key]
    except (KeyError,TypeError):
        # ratio p_merged_df / change_probabilities
        c2lp_2darray = p_merged_df.iloc[:,:-2].values
        cp_2darray = change_probabilities.reshape((-1,1))
        out = np.ones(c2lp_2darray.shape)
        with np.errstate(divide='ignore', invalid='ignore'):
            p_ratio = np.divide(c2lp_2darray, cp_2darray, out=out, where=cp_2darray != 0) # >= 0 / >= 0
            # when change_probability is 0  -> return 1 instead of nan, but calculus is done anyway so ignore warning "divide"
            # when c2l probability is 0     -> return 0 but np raises a warning so ignore warning "invalid"
            # all special cases and coping strategies:
            # n/0 = inf and raises warnign divide by 0 -> (never happens)
            # 0/0 = nan and raises warning divide by 0 -> replace with 1, ignore warning "divide"
            # 0/n = 0   and raises warning invalid value -> keep 0, ignore warning "invalid"

        # logarithm
        # compute log(x) only for x > TOO_SMALL4LOG_THRESHOLD, else use MINUS_INF_LIKE
        TOO_SMALL4LOG_THRESHOLD = 0.001
        MINUS_INF_LIKE = -10
        _logp = np.where(p_ratio > TOO_SMALL4LOG_THRESHOLD, p_ratio, MINUS_INF_LIKE)
        _logp = np.log(_logp, out=_logp, where=_logp > 0)
        cache[key] = _logp
    return _logp


def logp4lin(l, p_merged_df, change_probabilities, cache=dict()):
    """
    Returns the likelihood ratio scores for the given lineage. 

    The cache may be used to read the output of logp4all (likelihood ratios for all the candidates) and skip the computation. 
    It is advisable to not save the computation of likelihood ratios for a single candidate since logp4all already saves the likelihood ratios for all.

    :param l:
    :param p_merged_df:
    :param change_probabilities:
    :param cache:
    :return:
    """
    # ln(x) results is MINUS_INF_LIKE for x < TOO_SMALL4LOG_THRESHOLD
    key = ("logp4all", p_merged_df.shape)
    try:
        logp4l = cache[key][:,p_merged_df.columns.get_loc(l)]
    except KeyError:
        TOO_SMALL4LOG_THRESHOLD = 0.001
        MINUS_INF_LIKE = -10

        c2l1_probabilities = p_merged_df[l].values
        with np.errstate(divide='ignore', invalid='ignore'):
            # when change_probability is 0 -> return 1 instead of nan, but calculus is done anyway so ignore divide warning
            # when c2l1_probability is 0 -> return 0 but np raises a warning so ignore invalid warning
            p_ratio = np.where(change_probabilities == 0, 1.0, c2l1_probabilities / change_probabilities) # >= 0 / >= 0

        logp4l = np.where(p_ratio > TOO_SMALL4LOG_THRESHOLD, p_ratio, MINUS_INF_LIKE)
        logp4l = np.log(logp4l, out=logp4l, where=logp4l > 0)
    return logp4l


def clogp4all(merged_df, p_merged_df, change_probabilities, left_to_right=True, cache=dict()):
    """
    Returns the cumulative sum of likelihood contributions for all the lineages in merged_df. 
    For a lineage L, likelihood values contribute positively to the cumulative sum if the mutation is held by the target sequence, 
    negatively if only by the lineage L, do not contribute at all otherwise.

    Computing the cumulative sum of values is a generally fast operation. 
    It is advised to avoid caching them to avoid large memory consumption for a small performance advantage.
    The cache may be used to fetch likelihood ratio scores instead.


    :param merged_df:
    :param p_merged_df:
    :param change_probabilities:
    :param left_to_right:
    :return:
    """
    key = ("clogp4all", merged_df.shape, left_to_right)
    try:
        result = cache[key]
    except KeyError: # then compute
        _logp = logp4all(p_merged_df, change_probabilities, cache)
        len_changes = merged_df.shape[0]
        changes_of_target = merged_df['seq_change'].values

        def lin_cumulative_sum(l, l_idx):
            logp4l = _logp[:,l_idx]
            # positive_contribs
            contribs = np.where(changes_of_target, logp4l, np.zeros(len_changes))       # +logp if change of target, +0 otherwise (neutral contrib for changes not in target and lin)
            # negative contribs
            changes_of_l = merged_df[l].values
            contribs = np.where((~changes_of_target & changes_of_l), -logp4l, contribs) # -logp if change of lin but not of target
            return np.cumsum(contribs[::1 if left_to_right else -1])

        result = np.array([lin_cumulative_sum(lin, lin_idx) for lin_idx, lin in enumerate(merged_df.columns[:-2])]).T
        cache[key] = result
    return result


def clogp4lin(l, merged_df, logp4l, left_to_right=True):
    """
    Returns the cumulative sum of likelihood contributions for the given lineage. 
    Likelihood values contribute positively to the cumulative sum if the mutation is held by the target sequence, 
    negatively if only by the given lineage, and do not contribute at all otherwise.
    :param l:
    :param merged_df:
    :param logp4l: the likelihood ratio values for lienage l
    :param left_to_right:
    :return:
    """
    changes_of_target = merged_df['seq_change'].values
    len_changes = merged_df.shape[0]

    # positive_contribs
    contribs = np.where(changes_of_target, logp4l, np.zeros(len_changes))       # +logp if change of target, +0 otherwise (neutral contrib for changes not in target and lin)
    # negative contribs
    changes_of_l = merged_df[l].values
    contribs = np.where((~changes_of_target & changes_of_l), -logp4l, contribs) # -logp if change of lin but not of target
    return np.cumsum(contribs[::1 if left_to_right else -1])


def top_candidates(merged_df, p_merged_df, change_probabilities, l_edge_pos_idx=None, r_edge_pos_idx=None,
                   left_to_right=True, force_include_lineage=None, min_target_len=0,
                   cache = dict()
                   ) -> Tuple[List[str], np.array, List[int], np.array]:
    """
    Computes and compares the likelihood of all the possible candidates in the region defined between l_edge_pos_idx and
    r_edge_pos_idx. Returns the top 10 candidates with the highest likelihoods in any of the points between the given
    interval. If more than 10 candidates have the same top likelihood value, all of them are returned.
    The candidates are described through a set of sorted lists/np.array (best candidate first). The lists/np.array tell
    respectively: the candidate names, the maximum value, the position of the maximum value and
    all the likelihood values (covering the entire genome, not only the given interval).
    :param merged_df:
    :param p_merged_df:
    :param change_probabilities:
    :param l_edge_pos_idx:
    :param r_edge_pos_idx:
    :param left_to_right:
    :param force_include_lineage:
    :param min_target_len: minimum number of target changes that acceptable candidates must cover
    :param cache: (optional) a dictionary for memoization
    :return:
    """
    ## cumulative likelihood for all lineages
    # search in cache
    key = ("clogp4all", merged_df.shape, left_to_right)
    try:
        all_likelihoods = cache[key]
    except KeyError:
        # or compute
        all_likelihoods = clogp4all(merged_df, p_merged_df, change_probabilities, left_to_right, cache=cache)
        if not left_to_right:
            all_likelihoods = np.flip(all_likelihoods, axis=0)
        cache[key] = all_likelihoods

    # limit search
    all_likelihoods_window = all_likelihoods[l_edge_pos_idx:r_edge_pos_idx]

    # shift likelihood to reset accumulation before window's position
    if (left_to_right and l_edge_pos_idx) or (not left_to_right and r_edge_pos_idx is not None):
        shift_pos_idx = (l_edge_pos_idx - 1) if left_to_right else min(r_edge_pos_idx, merged_df.shape[0]-1)
        all_likelihoods_window_shift = all_likelihoods[shift_pos_idx,:]
        all_likelihoods_window = all_likelihoods_window - all_likelihoods_window_shift

    # maximum value of each lineage
    try:
        max4lin = np.max(all_likelihoods_window, axis=0)
    except ValueError:  # raised if `all_likelihoods_window` is empty.
        raise NoCandidatesFound(f"Empty region between l_edge_pos_idx and r_edge_pos_idx. Likelihood of candidates cannot be computed.")

    # find position of max value relative to all positions
    pos_max4lin_window = np.array([np.where(all_likelihoods_window[:, idx] == v)[0][0 if left_to_right else -1] for idx, v in enumerate(max4lin)])
    pos_max4lin = pos_max4lin_window if not l_edge_pos_idx else pos_max4lin_window + l_edge_pos_idx

    # keep lineages whose max value location is >= min_target_len changes of target from edge
    if left_to_right:
        t_changes_count_to_max_pos = np.array([np.sum(merged_df.seq_change[l_edge_pos_idx:i]) for i in (pos_max4lin + 1)])
    else:
        t_changes_count_to_max_pos = np.array([np.sum(merged_df.seq_change[i:r_edge_pos_idx]) for i in pos_max4lin])
    lin_max_pos_above_thresh_idx = np.where(t_changes_count_to_max_pos >= min_target_len)[0]
    if lin_max_pos_above_thresh_idx.size == 0:
        raise NoCandidatesFound("No candidate covering min_target_len changes has been found.")

    # pick first 10
    max4lin_above_thresh = max4lin[lin_max_pos_above_thresh_idx]    # sub-array of max4lin holding only max values of lin that passed previous filter
    mask_lin_at_max_value = max4lin_above_thresh == np.max(max4lin_above_thresh)   # tells which elements have maximum value
    n_lineages_to_pickup = np.sum(mask_lin_at_max_value)                            # tells how many have maximum value
    if n_lineages_to_pickup > 10:   # lots of even candidates in the chart... pointless sort -> return all top
        idx_top_lin = lin_max_pos_above_thresh_idx[mask_lin_at_max_value]
    else:
        n_lineages_to_pickup = min(max4lin_above_thresh.size, 10)
        partition_idx = max4lin_above_thresh.size - n_lineages_to_pickup
        idx_top_lin = lin_max_pos_above_thresh_idx[np.argsort(max4lin_above_thresh)[partition_idx:][::-1]]


    # force include one lineage in the top10
    if force_include_lineage:
        added_lin_idx = merged_df.columns.get_loc(force_include_lineage.upper())
        if added_lin_idx not in idx_top_lin:
            idx_top_lin = np.append(idx_top_lin, [added_lin_idx])
    names_top_lin = list(merged_df.columns[idx_top_lin])

    # likelihood values of the top candidates
    max_top_lin = max4lin[idx_top_lin]

    # position of the maximum value of the top candidates
    pos_max_top_lin = pos_max4lin[idx_top_lin]

    if len(names_top_lin) == 0 or len(max_top_lin) == 0 or len(pos_max_top_lin) == 0:
        raise AssertionError("names_top_lin, max_top_lin, pos_max_top_lin should never be empty. Exception exist for that condition.")

    return names_top_lin, max_top_lin, pos_max_top_lin, all_likelihoods[:,idx_top_lin]


def search_L(merged_df, p_merged_df, change_probabilities, l_edge_pos_idx=None, r_edge_pos_idx=None, min_target_len=0, cache=dict()):
    """
    Search any likely candidate for the region defined in the range <l_edge_pos_idx, r_edge_pos_idx> having maximum
    likelihood such that it coverts at least min_target_len changes of positions of the target genome.
    Search is performed by looking at increasing likelihoods in both directions, but only the list of the candidate
    that has the highest value is returned.
    May raise NoCandidatesFound.
    :param merged_df:
    :param p_merged_df:
    :param change_probabilities:
    :param l_edge_pos_idx:
    :param r_edge_pos_idx:
    :param min_target_len:
    :param cache: (optional) a dictionary for memoization
    :return:
    """
    try:
        candidates2dx = top_candidates(merged_df, p_merged_df, change_probabilities, l_edge_pos_idx=l_edge_pos_idx,
                                       r_edge_pos_idx=r_edge_pos_idx, left_to_right=True, min_target_len=min_target_len, cache=cache)
    except NoCandidatesFound:
        candidates2dx = None
    try:
        candidates2sx = top_candidates(merged_df, p_merged_df, change_probabilities, l_edge_pos_idx=l_edge_pos_idx,
                                       r_edge_pos_idx=r_edge_pos_idx, left_to_right=False, min_target_len=min_target_len, cache=cache)
    except NoCandidatesFound:
        candidates2sx = None
    if candidates2dx and candidates2sx:
        direction = int(candidates2dx[1][0] < candidates2sx[1][0])   # 0 is >> , 1 is <<
    elif candidates2dx:
        direction = 0
    elif candidates2sx:
        direction = 1
    else:
        raise NoCandidatesFound("No candidates found searching in both directions.")
    candidates = [candidates2dx, candidates2sx][direction]
    interval = Region.from_candidate_list(candidates, direction, merged_df, l_edge_pos_idx, r_edge_pos_idx)
    return interval


def search_L_fixed_direction(merged_df, p_merged_df, change_probabilities, left_to_right,
                             l_edge_pos_idx=None, r_edge_pos_idx=None, min_target_len=0, cache=dict()):
    """
    Like search_L it returns any good candidate lineage satisfying the requirement min_target_len in the given
    positional range, but the search direction is fixed.
    May raise NoCandidatesFound.
    :param merged_df:
    :param p_merged_df:
    :param change_probabilities:
    :param left_to_right:
    :param l_edge_pos_idx:
    :param r_edge_pos_idx:
    :param min_target_len:
    :param cache: (optional) a dictionary for memoization
    :return:
    """

    candidates = top_candidates(merged_df, p_merged_df, change_probabilities,
                                l_edge_pos_idx=l_edge_pos_idx,r_edge_pos_idx=r_edge_pos_idx,
                                left_to_right=left_to_right, min_target_len=min_target_len, cache=cache)
    # don't catch exception
    return Region.from_candidate_list(candidates, left_to_right == 0, merged_df, l_edge_pos_idx, r_edge_pos_idx)


def find_dual_region(merged_df, p_merged_df, change_probabilities, L1_reg, l_edge, r_edge, cache=dict()):
    """
    This function is a wrapper around compute_region_for_candidate that decides the direction in which to search the
    candidate of L1_reg and once the new region is found, L1_reg is set as master_region for it.

    :param merged_df: lineage characterizations and definition of the target seqeunce.
    :param p_merged_df: change probability for each lineage.
    :param change_probabilities: global change probability.
    :param L1_reg: the region of the first candidate identified in the sequence.
    :param l_edge: the left boundary of the region where to search the candidate lineage.
    :param r_edge: the right boundary of the region where to search the candidate lineage.
    :return: Region.
    """
    new_direction = abs(1 - L1_reg.search_dir)

    new_region = compute_region_for_candidate(L1_reg.designated, merged_df, p_merged_df, change_probabilities,
                                              l_edge, r_edge, new_direction == 0, cache=cache)
    new_region.set_master_region(L1_reg)
    return new_region


def compute_region_for_candidate(l, merged_df, p_merged_df, change_probabilities,
                                 l_edge_pos_idx=None, r_edge_pos_idx=None,
                                 left_to_right=True, cache=dict()):
    """
    This function computes a Region for the given candidate within the given positional range. The returned region may
    also 0 length if the candidate is not suited for the region/target genome.
    :param l:
    :param merged_df:
    :param p_merged_df:
    :param change_probabilities:
    :param l_edge_pos_idx:
    :param r_edge_pos_idx:
    :param left_to_right:
    :return:
    """
    cl = clogp4lin(l, merged_df, logp4lin(l, p_merged_df, change_probabilities, cache), left_to_right)
    cl = cl if left_to_right else cl[::-1]

    cl_window = cl[l_edge_pos_idx:r_edge_pos_idx]
    max_l = cl_window.max()
    max_l_pos_idx = np.where(cl_window == max_l)[0][0 if left_to_right else -1]

    if l_edge_pos_idx:
        max_l_pos_idx = max_l_pos_idx + l_edge_pos_idx

    direction = 0 if left_to_right else 1
    return Region.from_single_candidate(l, max_l, max_l_pos_idx, cl, direction, merged_df, l_edge_pos_idx, r_edge_pos_idx)


# def fill_gap(merged_df, p_merged_df, change_probabilities, L1_reg: Region, L2_reg: Region):
#     """
#     Sum the cumulative likelihoods of the regions within the gap. Finds the position that maximises the sum.
#     Divides the gap in intervals between positions of target (left edge is a target position, right edge not).
#     Finds the interval where the position of maximum is in. Assigns the target mutation at the left edge of the
#     interval to the left-adjacent region, the rest to the right-adjacent region.
#     """
#     l_gap_edge = L1_reg.pos_end
#     r_gap_edge = L2_reg.pos_start
#
#     logp4L_SX = L1_reg.likelihood_values(L1_reg.designated)
#     logp4L_DX = L2_reg.likelihood_values(L2_reg.designated)
#
#     sum_CL_gap_l1l2 = np.sum([logp4L_SX[l_gap_edge:r_gap_edge], logp4L_DX[l_gap_edge:r_gap_edge]], axis=0)
#
#     # shift_sx = logp4L_SX[max(1,l_gap_edge) - 1]
#     # shift_dx = logp4L_DX[min(r_gap_edge, merged_df.shape[0]-1)]
#     # sum_CL_gap_l1l2 = sum_CL_gap_l1l2 - (shift_sx + shift_dx)
#
#     # find pos of max likelihood
#     max_sum = sum_CL_gap_l1l2.max()
#     max_pos_idx = np.where(sum_CL_gap_l1l2 == max_sum)[0][0] + l_gap_edge
#
#     positions_of_t = np.where(merged_df.seq_change)[0]
#     interval_between_positions_of_t = list(zip(positions_of_t[:-1], positions_of_t[1:]))
#     max_interval = [(l, r) for l, r in interval_between_positions_of_t if l <= max_pos_idx < r][0]
#
#     L1_reg_pos_end = max_interval[0] + 1  # +1 because we want the start of the interval to be included
#     L2_reg_pos_start = L1_reg_pos_end
#
#     # Inspection figure (debug purpose)
#     # print(l_gap_edge, r_gap_edge)
#     # fig = go.Figure()
#     # fig.add_trace(
#     #     go.Scatter(x=list(range(l_gap_edge, r_gap_edge)), y=logp4L_SX[l_gap_edge:r_gap_edge], name="region SX"))
#     # fig.add_trace(
#     #     go.Scatter(x=list(range(l_gap_edge, r_gap_edge)), y=logp4L_DX[l_gap_edge:r_gap_edge], name="region DX"))
#     # fig.add_trace(go.Scatter(x=list(range(l_gap_edge, r_gap_edge)), y=sum_CL_gap_l1l2, name="sum"))
#     # fig.add_trace(go.Scatter(x=np.where(merged_df.seq_change[l_gap_edge:r_gap_edge])[0] + l_gap_edge,
#     #                          y=[1 for p in positions_of_t], name="positions of t"))
#     # fig.add_trace(go.Scatter(x=[max_pos_idx], y=[1], name="max of sum"))
#     # fig.add_trace(go.Scatter(x=[L2_reg_pos_start], y=[1], name="new start of region DX"))
#     # fig.show("png")
#
#     L1_reg.set_pos_end(L1_reg_pos_end)
#     L2_reg.set_pos_start(L2_reg_pos_start)
#
#     return logp4L_SX, logp4L_DX, sum_CL_gap_l1l2


def fill_gap(merged_df, p_merged_df, change_probabilities, left_region: Region, right_region: Region):
    """
    Sum the cumulative likelihoods of the regions within the gap. Finds the position that maximises the sum within the
    gap (max_pos) but ignoring the positions of target changes. max_pos can be at the left or right of a change of the
    candidate regions. The sum of the likelihoods is lower at the left of max_pos, and is lower or equal to
    the right of max_pos. Everything to the left of max_pos is assigned to the region at left, max_pos and everything to
    the right is assigned to the region at right.
    :param merged_df:
    :param p_merged_df:
    :param change_probabilities:
    :param left_region:
    :param right_region:
    :return:
    """
    l_gap_edge = left_region.pos_end
    r_gap_edge = right_region.pos_start

    logp4L_SX = left_region.likelihood_values(left_region.designated)
    logp4L_DX = right_region.likelihood_values(right_region.designated)

    sum_CL_gap_l1l2 = np.sum([logp4L_SX[l_gap_edge:r_gap_edge], logp4L_DX[l_gap_edge:r_gap_edge]], axis=0)

    positions_of_t_gap_idx = np.where(merged_df.seq_change[l_gap_edge:r_gap_edge])[0]
    sum_CL_gap_l1l2[positions_of_t_gap_idx] = np.nan

    # find pos of max likelihood
    if all(np.isnan(sum_CL_gap_l1l2)):
        logging.warning("Gap made of only target mutations (How can it be? It's an almost impossible condition). Target mutation that minimizes the loss has been assigned to the right region).")
        sum_CL_gap_l1l2 = np.sum([logp4L_SX[l_gap_edge:r_gap_edge], logp4L_DX[l_gap_edge:r_gap_edge]], axis=0)
        max_pos_idx = np.argmax(sum_CL_gap_l1l2) + l_gap_edge
    else:
        max_pos_idx = np.nanargmax(sum_CL_gap_l1l2) + l_gap_edge  # returns max value ignoring any nan

    # Inspection figure (debug purpose)
    # fig = go.Figure()
    # print(l_gap_edge, r_gap_edge)
    # fig.add_trace(go.Scatter(x=list(range(l_gap_edge,r_gap_edge)), y=logp4L_SX[l_gap_edge:r_gap_edge], name="region SX"))
    # fig.add_trace(go.Scatter(x=list(range(l_gap_edge,r_gap_edge)), y=logp4L_DX[l_gap_edge:r_gap_edge], name="region DX"))
    # fig.add_trace(go.Scatter(x=list(range(l_gap_edge,r_gap_edge)), y=sum_CL_gap_l1l2, name="sum"))
    # fig.add_trace(go.Scatter(x=positions_of_t_gap_idx + l_gap_edge, y=[1 for p in positions_of_t_gap_idx], name="positions of t"))
    # fig.add_trace(go.Scatter(x=[max_pos_idx], y=[1], name="max of sum"))
    # fig.show("png")

    left_region.set_pos_end(max_pos_idx)
    right_region.set_pos_start(max_pos_idx)

    return logp4L_SX, logp4L_DX, sum_CL_gap_l1l2


# def fill_gap(merged_df, p_merged_df, change_probabilities, L1_reg: Region, L2_reg: Region):
#     """
#     Sum the cumulative likelihoods of the regions within the gap. Finds the position that maximises the sum within the
#     gap (max_pos). Assigns max_pos to the region that is increasing its likelihood around max_pos. Progressively
#     increasing windows centered around max_pos are used to understand the trend of the likelihood of the two regions.
#     Assigns everything to the right of max_pos to the region at right, and everyting to the left to the region at left.
#     :param merged_df:
#     :param p_merged_df:
#     :param change_probabilities:
#     :param L1_reg:
#     :param L2_reg:
#     :return:
#     """
#     l_gap_edge = L1_reg.pos_end
#     r_gap_edge = L2_reg.pos_start
#
#     # TODO replace logp4lin with log values saved in the two regions
#     logp4L_SX = logp4lin(L1_reg.designated, merged_df, p_merged_df, change_probabilities, left_to_right=True)
#     logp4L_DX = logp4lin(L2_reg.designated, merged_df, p_merged_df, change_probabilities, left_to_right=False)[::-1]
#
#     sum_CL_gap_l1l2 = np.sum([logp4L_SX[l_gap_edge:r_gap_edge], logp4L_DX[l_gap_edge:r_gap_edge]], axis=0)
#
#     # shift_sx = logp4L_SX[max(1,l_gap_edge) - 1]
#     # shift_dx = logp4L_DX[min(r_gap_edge, merged_df.shape[0]-1)]
#     # sum_CL_gap_l1l2 = sum_CL_gap_l1l2 - (shift_sx + shift_dx)
#
#     # find pos of max likelihood
#     max_sum = sum_CL_gap_l1l2.max()
#     max_pos_idx = np.where(sum_CL_gap_l1l2 == max_sum)[0][0]
#
#     # assign the position to the lineage that was increasing his likelihood before the maximum position (almost always L2)
#     L1_non_decreasing = L2_non_decreasing = False
#     for window_size in range(2, max_pos_idx + 1):
#         L1_non_decreasing = non_decreasing(logp4L_SX[l_gap_edge + max_pos_idx - window_size:l_gap_edge + max_pos_idx])
#         L2_non_decreasing = non_decreasing(logp4L_DX[l_gap_edge + max_pos_idx - window_size:l_gap_edge + max_pos_idx])
#         if L1_non_decreasing and not L2_non_decreasing:
#             l_max_pos_idx = r_max_pos_idx = max_pos_idx + 1  # assigns max pos to L1
#             break
#         elif L2_non_decreasing and not L1_non_decreasing:
#             l_max_pos_idx = r_max_pos_idx = max_pos_idx  # assigns max pos to L2
#             break
#     if L1_non_decreasing == L2_non_decreasing:  # fallback case
#         l_max_pos_idx = r_max_pos_idx = max_pos_idx  # assigns max pos to the most likely situation (L2)
#
#     # try:
#     #     print("assignment of max of gap to L1/L2 concluded between", max_pos_idx - window_size, max_pos_idx, "with",
#     #           L1_non_decreasing, L2_non_decreasing)
#     # except UnboundLocalError:
#     #     pass
#
#     # change reference of max_pos_idx from gap to merged_df
#     L1_reg.set_pos_end(l_max_pos_idx + l_gap_edge)
#     L2_reg.set_pos_start(r_max_pos_idx + l_gap_edge)
#
#     return logp4L_SX, logp4L_DX, sum_CL_gap_l1l2
#
#
# def non_decreasing(L):
#     return all(x <= y for x, y in zip(L, L[1:]))


def compare_aik(v1, against_v2):
    try:
        ans = math.exp((round(v1, 2) - round(against_v2, 2)) / 2)
    except OverflowError:
        if (v1 - against_v2) < 0:
            ans = np.NINF
        else:
            ans = np.inf
    return ans


def log_values_pc2l(l, merged_df, p_merged_df, left_to_right=True, change_mask=None):
    TOO_SMALL4LOG_THRESHOLD = 1e-09
    MINUS_INF_LIKE = -20

    if change_mask is not None:
        changes_of_target = merged_df['seq_change'][change_mask].values
        changes_of_l1 = merged_df[l][change_mask].values
        len_changes = np.sum(change_mask)

        c2l1_probabilities = p_merged_df[[l]][change_mask].values.reshape(len_changes)  # >= 0

    else:
        changes_of_target = merged_df['seq_change'].values
        changes_of_l1 = merged_df[l].values
        len_changes = merged_df.shape[0]

        c2l1_probabilities = p_merged_df[[l]].values.reshape(len_changes)  # >= 0

    logp = np.where(c2l1_probabilities > TOO_SMALL4LOG_THRESHOLD, c2l1_probabilities, MINUS_INF_LIKE)
    logp = np.log(logp, out=logp, where=logp > 0)

    l1_running_logp = []
    l1_logp = 0
    for i in range(len_changes) if left_to_right else range(len_changes - 1, -1, -1):
        if changes_of_target[i]:
            l1_logp += logp[i]
        else:
            l1_logp -= logp[i]
        l1_running_logp.append(l1_logp)
    return l1_running_logp


def log_values_pc(merged_df, change_probabilities, left_to_right=True, change_mask=None):
    TOO_SMALL4LOG_THRESHOLD = 1e-09
    MINUS_INF_LIKE = -20

    if change_mask is not None:
        changes_of_target = merged_df['seq_change'][change_mask].values
        len_changes = np.sum(change_mask)
        c_probabilities = change_probabilities[change_mask]  # > 0

    else:
        changes_of_target = merged_df['seq_change'].values
        len_changes = merged_df.shape[0]
        c_probabilities = change_probabilities  # > 0

    logp = np.where(c_probabilities > TOO_SMALL4LOG_THRESHOLD, c_probabilities, MINUS_INF_LIKE)
    logp = np.log(logp, out=logp, where=logp > 0)

    l1_running_logp = []
    l1_logp = 0
    for i in range(len_changes) if left_to_right else range(len_changes - 1, -1, -1):
        if changes_of_target[i]:
            l1_logp += logp[i]
        else:
            l1_logp -= logp[i]
        l1_running_logp.append(l1_logp)
    return l1_running_logp


def aik_p_values(merged_df, p_merged_df, change_probabilities, regions: list) -> Tuple[dict, dict]:
    if len(regions) is None:
        raise ValueError
    num_genome_positions = merged_df.shape[0]

    candidates = [r.designated for r in regions]
    candidates_regions = [(r.pos_start, r.pos_end) for r in regions]

    candidates_change_mask = np.zeros(num_genome_positions, dtype=bool)  # i.e. parameters
    for c in candidates:
        candidates_change_mask = candidates_change_mask | merged_df[c]
    candidates_change_mask = candidates_change_mask | merged_df.seq_change
    n_parameters = candidates_change_mask.sum()

    # AIK candidates
    aik_candidates = []
    for c in dict.fromkeys(candidates):
        c_log_pc2l = log_values_pc2l(c, merged_df, p_merged_df, change_mask=candidates_change_mask)[-1]
        aik = 2 * n_parameters - 2 * c_log_pc2l
        aik_candidates.append((c, aik))
    aik = dict(aik_candidates)

    # AIK SC2
    c_log_pc = log_values_pc(merged_df, change_probabilities, change_mask=candidates_change_mask)[-1]
    aik_sc2 = 2 * n_parameters - 2 * c_log_pc
    aik["general_sc2_model"] = aik_sc2

    # AIK recombination
    if len(candidates) > 1:
        c_log_candidate_regions = []
        for c, (s, e) in zip(candidates, candidates_regions):
            c_region_change_mask = candidates_change_mask.copy()
            c_region_change_mask[:s] = False
            c_region_change_mask[e:] = False
            c_region_log_pc2l = log_values_pc2l(c, merged_df, p_merged_df, left_to_right=True,
                                                change_mask=c_region_change_mask)[-1]
            c_log_candidate_regions.append(c_region_log_pc2l)
        aik_recombination = 2 * (n_parameters + 1) - 2 * (sum(c_log_candidate_regions))
        aik["recombinant_model"] = aik_recombination

    # P-VALUES
    p_values = []
    if len(candidates) > 1:    # compare recombinant model wrt single assignments
        for c in dict.fromkeys(candidates):
            p = compare_aik(aik["recombinant_model"], aik[c])
            p_values.append((("recombinant_model", c), p))
    else:                      # or compare non-recombinant model against general sc2 model
        c = candidates[0]
        p = compare_aik(aik[c], aik["general_sc2_model"])
        p_values.append(((c, "general_sc2_model"), p))
    p_values = dict(p_values)

    return aik, p_values

#
# def p_values_comparing_candidates(merged_df, p_merged_df, region: Region) -> dict:
#     if not region.all_other_candidates:
#         raise ValueError("No alternative candidates to compare with.")
#     num_genome_positions = merged_df.shape[0]
#
#     candidates = [region.designated] + region.all_other_candidates
#
#     # mask of the mutations of all the candidates + target
#     candidates_change_mask = np.zeros(num_genome_positions, dtype=bool)
#     for c in candidates:
#         candidates_change_mask = candidates_change_mask | merged_df[c]
#     candidates_change_mask = candidates_change_mask | merged_df.seq_change
#     # remove mutations of candidates outside the region of interest
#     candidates_change_mask[:region.pos_start] = False
#     candidates_change_mask[region.pos_end:] = False
#     # count parameters
#     n_parameters = candidates_change_mask.sum() # i.e. parameters
#
#     # AIK candidates
#     aik_candidates = []
#     for c in dict.fromkeys(candidates):
#         c_log_pc2l = log_values_pc2l(c, merged_df, p_merged_df, change_mask=candidates_change_mask)[-1]
#         aik = 2 * n_parameters - 2 * c_log_pc2l
#         aik_candidates.append((c, aik))
#     aik = dict(aik_candidates)
#
#     # P-VALUES
#     p_values = []
#     for alt_c in dict.fromkeys(region.all_other_candidates):
#         p = compare_aik(aik[region.designated], aik[alt_c])
#         p_values.append(((region.designated, alt_c), p))
#     p_values = dict(p_values)
#
#     return p_values


def update_region_p_values(merged_df, p_merged_df, change_probabilities, regions: list, cache=dict()):
    """
    For every region, computes the likelihood, aic and p-value of every candidate and saves them inside the attribute
    Region.candidates_similarity_stats
    :param merged_df:
    :param p_merged_df:
    :param change_probabilities:
    :param regions:
    :return:
    """
    is_2bp = len(regions) == 3
    if is_2bp:
        L1_5p, L2, L1_3p = regions

        # ALTERNATIVE CANDIDATES L1
        L1_changes = (
                mask_changes_of_candidates_or_target(merged_df, L1_5p) |
                mask_changes_of_candidates_or_target(merged_df, L1_3p)
        )
        all_other_candidates = L1_5p.alternative_candidates() + L1_3p.alternative_candidates()
        all_candidates = [L1_5p.designated] + all_other_candidates

        ## AIC
        n_parameters = L1_changes.sum()
        likelihood_candidates = dict()
        aic_candidates = []
        for c in dict.fromkeys(all_candidates):
            logp4l = logp4lin(c, p_merged_df, change_probabilities, cache)[L1_changes]
            c_log = clogp4lin(c, merged_df[L1_changes], logp4l)[-1]
            likelihood_candidates[c] = c_log
            aic = 2 * n_parameters - 2 * c_log
            aic_candidates.append((c, aic))
        aic_candidates = dict(aic_candidates)

        ## P-VALUES
        p_values = []
        for c in dict.fromkeys(all_other_candidates):
            p = compare_aik(aic_candidates[L1_5p.designated], aic_candidates[c])
            p_values.append((c, p))
        p_values = dict(p_values)

        output = {
            L1_5p.designated: (
                likelihood_candidates[L1_5p.designated],
                aic_candidates[L1_5p.designated],
                1
            )
        }
        output.update({
            c: (
                likelihood_candidates[c],
                aic_candidates[c],
                p_values[c]
            ) for c in all_other_candidates
        })
        L1_5p.candidates_similarity_stats = L1_3p.candidates_similarity_stats = output

        # ALTERNATIVE CANDIDATES L2
        L2_changes = mask_changes_of_candidates_or_target(merged_df, L2)

        ## AIC
        n_parameters = L2_changes.sum()
        likelihood_candidates = dict()
        aic_candidates = []
        for c in dict.fromkeys(L2.candidates):
            c_log = clogp4lin(c, merged_df[L2_changes], logp4lin(c, p_merged_df[L2_changes], change_probabilities[L2_changes], cache))[-1]
            likelihood_candidates[c] = c_log
            aic = 2 * n_parameters - 2 * c_log
            aic_candidates.append((c, aic))
        aic_candidates = dict(aic_candidates)

        ## P-VALUES
        p_values = []
        for c in dict.fromkeys(L2.alternative_candidates()):
            p = compare_aik(aic_candidates[L2.designated], aic_candidates[c])
            p_values.append((c, p))
        p_values = dict(p_values)

        output = {
            L2.designated: (
                likelihood_candidates[L2.designated],
                aic_candidates[L2.designated],
                1
            )
        }
        output.update({
            c: (
                likelihood_candidates[c],
                aic_candidates[c],
                p_values[c]
            ) for c in L2.alternative_candidates()
        })
        L2.candidates_similarity_stats = output
    else:
        for r in regions:
            r_changes = mask_changes_of_candidates_or_target(merged_df, r)

            ## AIC
            n_parameters = r_changes.sum()
            likelihood_candidates = dict()
            aic_candidates = []
            for c in dict.fromkeys(r.candidates):
                c_log = clogp4lin(c, merged_df[r_changes], logp4lin(c, p_merged_df[r_changes], change_probabilities[r_changes], cache))[-1]
                likelihood_candidates[c] = c_log
                aic = 2 * n_parameters - 2 * c_log
                aic_candidates.append((c, aic))
            aic_candidates = dict(aic_candidates)

            ## P-VALUES
            p_values = []
            for c in dict.fromkeys(r.alternative_candidates()):
                p = compare_aik(aic_candidates[r.designated], aic_candidates[c])
                p_values.append((c, p))
            p_values = dict(p_values)

            output = {
                r.designated: (
                    likelihood_candidates[r.designated],
                    aic_candidates[r.designated],
                    1
                )
            }
            output.update({
                c: (
                    likelihood_candidates[c],
                    aic_candidates[c],
                    p_values[c]
                ) for c in r.alternative_candidates()
            })
            r.candidates_similarity_stats = output


def mask_changes_of_candidates_or_target(merged_df, region: Region) -> np.array:
    candidates_change_mask = np.zeros(merged_df.shape[0], dtype=bool)
    for c in region.candidates:
        candidates_change_mask = candidates_change_mask | merged_df[c]
    candidates_change_mask = candidates_change_mask | merged_df.seq_change

    # remove mutations of candidates outside the region of interest
    candidates_change_mask[:region.pos_start] = False
    candidates_change_mask[region.pos_end:] = False
    return candidates_change_mask


def aic_on_range(merged_df, p_merged_df, change_probabilities, regions: list, l_edge_pos_idx, r_edge_pos_idx, cache=dict()):
    """
    Returns the AIC for the positional range <l_edge_pos_idx, r_edge_pos_idx>. Only the regions that are partially or
    completely overlappin with the given range are used and only for the portion included in the range.
    :param merged_df:
    :param p_merged_df:
    :param change_probabilities:
    :param regions:
    :param l_edge_pos_idx:
    :param r_edge_pos_idx:
    :return:
    """
    if len(regions) == 0:
        raise ValueError("regions list is empty.")
    num_genome_positions = merged_df.shape[0]

    interested_regions = [r for r in regions if not (r.pos_end <= l_edge_pos_idx or r.pos_start >= r_edge_pos_idx)]
    interested_candidates = [r.designated for r in interested_regions]

    # change mask of the changes of interested candidates + target
    candidates_change_mask = np.zeros(num_genome_positions, dtype=bool)  # i.e. parameters
    for c in interested_candidates:
        candidates_change_mask = candidates_change_mask | merged_df[c]
    candidates_change_mask = candidates_change_mask | merged_df.seq_change
    candidates_change_mask[:l_edge_pos_idx] = False
    candidates_change_mask[r_edge_pos_idx:] = False
    n_parameters = candidates_change_mask.sum()

    # AIK
    pos_interested_regions = [(r.pos_start, r.pos_end) for r in interested_regions]
    c_log_candidate_regions = []
    for c, (s, e) in zip(interested_candidates, pos_interested_regions):
        ss = max(s, l_edge_pos_idx)
        ee = min(e, r_edge_pos_idx)
        # compute logp4lin for candidate c only on changes of t and c, withing the region boundaries
        logp4l = logp4lin(c, p_merged_df, change_probabilities, cache)[ss:ee]
        c_region_clog = clogp4lin(c, merged_df[ss:ee], logp4l)[-1]
        # left_to_right / right_to_left we don't care
        c_log_candidate_regions.append(c_region_clog)
    # pp({
    #     "n_parameters": n_parameters,
    #     "candidates": interested_candidates,
    #     "candidates_clog": c_log_candidate_regions,
    #     "sum_clog": sum(c_log_candidate_regions)
    # })
    return 2 * (n_parameters + len(interested_candidates) - 1) - 2 * (sum(c_log_candidate_regions))



# def cluster_alternative_candidates_old(merged_df, p_merged_df, change_probabilities, regions: list, is_2bp, max_p_value_diff):
#     if is_2bp:
#         L1_5p, L2, L1_3p = regions
#
#         L1_changes = (
#                 mask_changes_of_candidates_or_target(merged_df, L1_5p) |
#                 mask_changes_of_candidates_or_target(merged_df, L1_3p)
#         )
#         alt_cand_p_values = alternative_candidate_p_values(merged_df, p_merged_df, change_probabilities,
#                                                            L1_5p.designated, L1_5p.all_other_candidates, L1_changes)
#         good_alternative_candidates = dict(
#             [(cand, p_val) for cand, p_val in alt_cand_p_values.items() if
#              p_val <= max_p_value_diff])
#         # assign good alternative candidates to both regions L1 5' and L1 3'
#         L1_5p.alternative_designations = L1_3p.alternative_designations = good_alternative_candidates
#
#         # find good alternative candidates for the region L2
#         L2_changes = mask_changes_of_candidates_or_target(merged_df, L2)
#         alt_cand_p_values = alternative_candidate_p_values(merged_df, p_merged_df, change_probabilities,
#                                                            L2.designated, L2.all_other_candidates, L2_changes)
#         L2.alternative_designations = dict(
#             [(cand, p_val) for cand, p_val in alt_cand_p_values.items() if
#              p_val <= max_p_value_diff])
#     else:
#         for r in regions:
#             r_changes = mask_changes_of_candidates_or_target(merged_df, r)
#             alt_cand_p_values = alternative_candidate_p_values(merged_df, p_merged_df, change_probabilities,
#                                                                r.designated, r.all_other_candidates, r_changes)
#             r.alternative_designations = dict(
#                 [(cand, p_val) for cand, p_val in alt_cand_p_values.items() if
#                  p_val <= max_p_value_diff])


# def alternative_candidate_p_values(merged_df, p_merged_df, change_probabilities, designated: str, alternative_candidates: list,
#                                    candidates_change_mask):
#     # AIK candidates
#     n_parameters = candidates_change_mask.sum()
#     aik_candidates = []
#     for c in dict.fromkeys([designated] + alternative_candidates):
#         c_log_pc2l = log_values_pc2l(c, merged_df, p_merged_df, change_mask=candidates_change_mask)[-1]
#         aik = 2 * n_parameters - 2 * c_log_pc2l
#         aik_candidates.append((c, aik))
#     aik = dict(aik_candidates)
#
#     # AIK SC2
#     c_log_pc = log_values_pc(merged_df, change_probabilities, change_mask=candidates_change_mask)[-1]
#     aik_sc2 = 2 * n_parameters - 2 * c_log_pc
#     aik["general_sc2_model"] = aik_sc2
#
#     # P-VALUES
#     p_values = []
#     for c in dict.fromkeys([designated] + alternative_candidates):
#         p = compare_aik(aik[c], aik["general_sc2_model"])
#         p_values.append(((c, "general_sc2_model"), p))
#     p_values = dict(p_values)
#
#     return p_values

# TODO delete old and commented functions
