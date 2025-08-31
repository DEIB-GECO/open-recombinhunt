import ast
import logging
import math
import random
from pprint import pformat
from typing import Any

import numpy as np
from plotly.subplots import make_subplots
from plotly import graph_objects as go
from tqdm import tqdm
from recombinhunt.core.environment import Environment, PangoLineageHierarchy
from recombinhunt.core.method import Experiment
from recombinhunt.validation.utils import *


class GroundTruth:
    P_VAL_ACCEPT_THRESHOLD = 1e-05

    def __init__(self, lineage_hierarchy: PangoLineageHierarchy,
                 environment, test_sequence_dataset=None, max_samples=None, shuffle_test_dataset=False,
                 experiment_output_csv=None):
        self.lh = lineage_hierarchy

        assert any([test_sequence_dataset, experiment_output_csv]) \
               and not all([test_sequence_dataset, experiment_output_csv]), \
            "One in (test_sequence_dataset,experiment_output_csv) should not be None"

        # get experiment's output
        if experiment_output_csv is not None:
            experiment_output = self.read_batch_experiment_output_from_file(experiment_output_csv)
        else:
            experiment_results = self.run_batch_experiment(environment, test_sequence_dataset, max_samples,
                                                           shuffle_test_dataset)
            experiment_output = self.format_experiment_output_as_df(experiment_results)
            experiment_output.to_csv(self.experiment_output_file_name())
            print(f"Batch run of experiment saved to: {self.experiment_output_file_name()}")
        self.batch_experiment_output = experiment_output

    # COMPUTE EXPERIMENTS
    def run_batch_experiment(self, env: Environment, test_sequence_dataset, max_samples=None, shuffle=False) -> list:
        if shuffle:
            if max_samples:
                test_sequence_dataset = random.sample(test_sequence_dataset, min(max_samples, len(test_sequence_dataset)))
            else:
                random.shuffle(test_sequence_dataset)
        elif max_samples:
            test_sequence_dataset = test_sequence_dataset[:max_samples]
        # else do nothing more


        experiment_results = []
        for exp_input_num, experiment_input in tqdm(enumerate(test_sequence_dataset)):
            try:
                name, true_lin, nuc_changes = experiment_input
            except:
                logging.exception(f"Failed attempt to parse input at line #{exp_input_num+1}: {experiment_input}\n"
                                  f"Due to the following reason:")
            else:
                try:
                    nuc_changes = nuc_changes.split(',')  # convert string to list
                    exp = Experiment(env, self.lh)
                    exp.set_target(nuc_changes)
                    exp.run()
                    experiment_results.append(self.format_experiment_output(name, true_lin, nuc_changes, exp))
                except:
                    logging.exception(f"Sequence {name} skipped due to the following reason:")
        return experiment_results

    def gt_contributing_lineages(self, *args) -> list:
        raise NotImplementedError

    # EXPERIMENT OUTPUT FORMATTING
    def format_experiment_output(self, name: str, true_lin: str, nuc_changes: list, exp: Experiment) -> list:
        is_recombinant = exp.genome_view.number_of_regions() > 1
        matching_candidates = all_candidates_matching(
            exp.genome_view.regions, self.gt_contributing_lineages(true_lin), self.lh)
        designated_candidates = exp.genome_view.contributing_lineages()
        good_alt_candidates = [r.get_good_alternative_candidates() for r in exp.genome_view.regions]
        breakpoints = exp.genome_view.breakpoints()
        breakpoints_in_t = exp.genome_view.breakpoints_in_t()
        number_of_breakpoints = exp.genome_view.number_of_breakpoints()
        p_values = [i for i in exp.p_values.values()]
        seq_len = exp.genome_view.len_of_t
        mutations = exp.merged_df[exp.merged_df.seq_change].seq_change.index.tolist()

        return [name, is_recombinant, matching_candidates, designated_candidates, good_alt_candidates,
                number_of_breakpoints, breakpoints, breakpoints_in_t, p_values, seq_len, mutations]

    # WRITE/READ EXPERIMENT OUTPUT
    @staticmethod
    def format_experiment_output_as_df(experiment_output: list) -> pd.DataFrame:
        return (
            pd.DataFrame(experiment_output,
                         columns=['name', 'is_recombinant', 'matching_candidates', 'designated_candidates',
                                  'good_alt_candidates', 'number_of_breakpoints', 'breakpoints', 'breakpoints_in_t',
                                  'p_values', 'seq_len', 'mutations'])
            .set_index('name')
        )

    @staticmethod
    def read_batch_experiment_output_from_file(file_name) -> pd.DataFrame:
        try:
            experiment_output_1 = pd.read_csv(file_name, index_col="name")
            experiment_output_1.designated_candidates = experiment_output_1.designated_candidates.map(ast.literal_eval).values.tolist()
            experiment_output_1.good_alt_candidates = experiment_output_1.good_alt_candidates.map(ast.literal_eval).values.tolist()
            experiment_output_1.breakpoints = experiment_output_1.breakpoints.map(ast.literal_eval).values.tolist()
            experiment_output_1.breakpoints_in_t = experiment_output_1.breakpoints_in_t.map(ast.literal_eval).values.tolist()
            experiment_output_1.p_values = experiment_output_1.p_values.map(ast.literal_eval).values.tolist()
            experiment_output_1.mutations = experiment_output_1.mutations.map(ast.literal_eval).values.tolist()
        except:
            logging.exception("Error parsing the given experiment output due to the following reason")
            raise
        return experiment_output_1

    def experiment_output_file_name(self) -> str:
        raise NotImplementedError

    def get_batch_experiment_output(self) -> pd.DataFrame:
        return self.batch_experiment_output

    # COMPARISON OUTPUT
    def validation_output_file_name(self) -> str:
        raise NotImplementedError

    def get_full_figure_output(self) -> go.Figure:
        raise NotImplementedError

    def get_table_row(self) -> list:
        raise NotImplementedError

    def get_breakpoint_histogram(self) -> go.Figure:
        raise NotImplementedError


class RecombinantGroundTruth(GroundTruth):
    def __init__(self, lineage, gt_contributing_lineages, gt_breakpoint_ranges, lh: PangoLineageHierarchy,
                 environment, test_sequence_dataset=None, max_samples=None, shuffle_test_dataset=False,
                 experiment_output_csv=None):
        self.lineage = lineage
        self._gt_contributing_lineages = gt_contributing_lineages
        self.gt_breakpoint_ranges = gt_breakpoint_ranges or []
        self.gt_number_of_breakpoints = len(self.gt_breakpoint_ranges)
        super().__init__(lh, environment, test_sequence_dataset, max_samples, shuffle_test_dataset,
                         experiment_output_csv)
        # following code requires self.batch_experiment_output
        try:
            self.predefined_consensus_seq = environment.x_characterizing_nuc_mutations(self.lineage)
        except KeyError:
            logging.warning(f"target sequence of {lineage} unavailable. Using consensus on dataset for breakpoint location.")
            self.predefined_consensus_seq = self.samples_consensus_sequence()
        # following code requires self.predefined consensus_seq
        self.gt_breakpoint_ranges_t = [BreakpointsLocation.to_target_pos(x, self.predefined_consensus_seq) for x in
                                       self.gt_breakpoint_ranges]

    def __str__(self):
        return pformat(self.__dict__, sort_dicts=False)

    def experiment_output_file_name(self):
        return f'experiment_output_{self.lineage}.csv'

    def validation_output_file_name(self):
        return f'validation_output_{self.lineage}'

    def gt_contributing_lineages(self, *args) -> list:
        return self._gt_contributing_lineages

    def get_full_figure_output(self):
        true_lineage = self.lineage
        number_of_sequences = self.batch_experiment_output.shape[0]

        # find the number of plot to draw
        number_of_breakpoints = max(self.batch_experiment_output.groupby("number_of_breakpoints", sort=False).size().index.tolist())
        box_plots_of_breakpoints_titles = [f'{infl.ordinal(br_idx+1)} breakpoint (target pos.)' for br_idx in range(number_of_breakpoints)]
        # box_plots_of_breakpoints_titles = [x if y%2==0 else x[:-13] + "(target length)" for x in box_plots_of_breakpoints_titles for y in range(2)]
        # box_plots_of_breakpoints_titles_coord = [f'{infl.ordinal(br_idx+1)} breakpoint (genomic pos)' for br_idx in range(number_of_breakpoints)]
        # box_plots_of_breakpoints_titles_coord = [x if y%2==0 else x[:-13] + "(genomic length)" for x in box_plots_of_breakpoints_titles_coord for y in range(2)]
        # box_plots_of_breakpoints_titles += box_plots_of_breakpoints_titles_coord
        number_of_rows = int(math.ceil((4+len(box_plots_of_breakpoints_titles))/3))
        subplot_titles = ["Recombinants detected", "Correctness of candidates", "N° Breakpoints"]+box_plots_of_breakpoints_titles+["P-values < 1e-05"]
        number_of_plots = len(subplot_titles)
        fig = make_subplots(rows=number_of_rows, cols=3, subplot_titles=subplot_titles)
        plot_ordinal_number = 0

        # MODEL CLASS evaluation
        plot_ordinal_number += 1
        is_recombinant_number = np.sum(self.batch_experiment_output.is_recombinant)
        is_not_recombinant_number = self.batch_experiment_output.shape[0] - is_recombinant_number
        fig.add_trace(go.Bar(x=["True", "False"], y=[is_recombinant_number, is_not_recombinant_number], name='recombination detected'), row=1, col=1)
        fig.update_layout(dict1={
            "yaxis": dict(tickformat=',d', title_text="number of sequences"),
            "xaxis": dict(title_text="is recombinant?")
        })


        # CANDIDATES evaluation
        plot_ordinal_number += 1
        number_exact = np.sum(self.batch_experiment_output.matching_candidates)
        number_diff = len(self.batch_experiment_output.matching_candidates) - number_exact
        fig.add_trace(go.Bar(x=["as ground truth", "different"], y=[number_exact, number_diff],
                             name='candidates found vs true candidates'), row=1, col=2)
        fig.update_layout(dict1={
            "yaxis2": dict(tickformat=',d', title_text="number of sequences"),
            "xaxis2": dict(title_text="candidate lineages")
        })

        # N° BREAKPOINTS
        plot_ordinal_number += 1
        breakpoint_number_value_counts = self.batch_experiment_output.groupby("number_of_breakpoints", sort=False).size()
        # Series counting the num. of sequences having 0,1,2 breakpoints (actually may include less than all three options)
        br_ = dict.fromkeys([0,1,2], 0)     # dict{0:num.seq.0-br, 1:num.seq.1-br, 2:num.seq.2-br}
        for k,v in breakpoint_number_value_counts.items():
            br_[k] = v
        fig.add_trace(go.Bar(x=list(br_.keys()), y=list(br_.values()), name='breakpoints quantity'), row=1, col=3)
        fig.update_layout(dict1={
            "yaxis3": dict(tickformat=',d', title_text="number of sequences"),
            "xaxis3": dict(title_text="breakpoints nr.")
        })

        # BREAKPOINTS STABILITY IN T
        row_n = 2
        col_n = 1
        max_seq_len = np.max(self.batch_experiment_output.seq_len)
        highest_breakpoint_num = max(breakpoint_number_value_counts.index.tolist())
        for br_idx in range(highest_breakpoint_num):  # sequences with 1-br -> br_idx=0, 2-br -> br_idx=1, 0 -> self.batch_experiment_output.breakpoints_in_t.values is empty list
            plot_ordinal_number += 1
            # find 1st breakpoint pair of positions of all sequences having at least one breakpoint
            # (otherwise self.batch_experiment_output.breakpoints_in_t is an empty list and is skipped)
            ith_breakpoint_values = []
            for breakpoints_of_sequence in self.batch_experiment_output.breakpoints_in_t.values:
                try:
                    ith_breakpoint_values.append(breakpoints_of_sequence[br_idx])
                except IndexError:
                    pass
            # unzip pairs
            first_pos_t = [i[0]+1 for i in ith_breakpoint_values]   # +1 for 1-based conversion
            last_pos_t = [i[1]+1 for i in ith_breakpoint_values]    # +1 for 1-based conversion
            breakpoint_length_t = [last_pos_t[i] - first_pos_t[i] for i in range(len(last_pos_t))]
            # make x/y values for histogram
            all_pos_of_t = dict.fromkeys(list(range(1,max_seq_len+1)), 0)
            for p in first_pos_t:
                all_pos_of_t[p] += 1

            fig.add_trace(go.Bar(x=list(all_pos_of_t.keys()), y=list(all_pos_of_t.values()), name=f"{infl.ordinal(br_idx+1)} breakpoint start pos"), row=row_n, col=col_n)
            fig.update_layout(dict1={
                f"yaxis{plot_ordinal_number}": dict(tickformat=',d', title_text="number of sequences"),
                f"xaxis{plot_ordinal_number}": dict(title_text="breakpoint coordinate on target")
            })
            # fig.add_trace(go.Violin(y=first_pos_t, name=f"{infl.ordinal(br_idx+1)} breakpoint length", box_visible=True, line_color='black',
            #                        meanline_visible=True, fillcolor='lightseagreen'), row=row_n, col=col_n)

            col_n += 1
            if col_n > 3:
                row_n += 1
                col_n = 1
            # fig.add_trace(go.Box(y=breakpoint_length_t, name=f"{infl.ordinal(br_idx+1)} breakpoint length", boxpoints='all'), row=row_n, col=col_n)
            # # fig.add_trace(go.Violin(y=breakpoint_length_t, name=f"{infl.ordinal(br_idx+1)} breakpoint length", box_visible=True, line_color='black',
            # #                        meanline_visible=True, fillcolor='lightseagreen'), row=row_n, col=col_n)
            # col_n += 1
            # if col_n > 3:
            #     row_n += 1
            #     col_n = 1

        # # BREAKPOINTS STABILITY IN GENOMIC COORDINATES
        # # row_n = 2
        # # col_n = 1
        # for br_idx in range(max(breakpoint_number_value_counts.index.tolist())):
        #     # plot 2 distinct box plot (for start/end coordinates) for each breakpoint
        #     ith_breakpoint_values = []
        #     for breakpoints_of_sequence in self.batch_experiment_output.breakpoints.values:
        #         try:
        #             ith_breakpoint_values.append(breakpoints_of_sequence[br_idx])
        #         except IndexError:
        #             pass
        #     # unzip pairs
        #     first_pos = [i[0] for i in ith_breakpoint_values]
        #     last_pos = [i[1] for i in ith_breakpoint_values]
        #     breakpoint_length = [last_pos[i] - first_pos[i] for i in range(len(last_pos))]
        #
        #     fig.add_trace(go.Histogram(x=first_pos, name=f"{infl.ordinal(br_idx+1)} breakpoint start pos", xbins_size=100), row=row_n, col=col_n)
        #     # fig.add_trace(go.Violin(y=first_pos_t, name=f"{infl.ordinal(br_idx+1)} breakpoint length", box_visible=True, line_color='black',
        #     #                        meanline_visible=True, fillcolor='lightseagreen'), row=row_n, col=col_n)
        #
        #     col_n += 1
        #     if col_n > 3:
        #         row_n += 1
        #         col_n = 1
        #     fig.add_trace(go.Box(y=breakpoint_length, name=f"{infl.ordinal(br_idx+1)} breakpoint length", boxpoints='all'), row=row_n, col=col_n)
        #     # fig.add_trace(go.Violin(y=breakpoint_length, name=f"{infl.ordinal(br_idx+1)} breakpoint length", box_visible=True, line_color='black',
        #     #                        meanline_visible=True, fillcolor='lightseagreen'), row=row_n, col=col_n)
        #     col_n += 1
        #     if col_n > 3:
        #         row_n += 1
        #         col_n = 1


        # P-VALUE test
        plot_ordinal_number += 1
        below_threshold = np.sum(self.batch_experiment_output.p_values.apply(lambda pvaluelist: all([x <= GroundTruth.P_VAL_ACCEPT_THRESHOLD for x in pvaluelist])))
        above_threshold = self.batch_experiment_output.shape[0] - below_threshold
        fig.add_trace(go.Bar(x=[True, False], y=[below_threshold, above_threshold], name="p-value test"), row=row_n, col=col_n)
        fig.update_layout(dict1={
            f"yaxis{plot_ordinal_number}": dict(tickformat=',d', title_text="number of sequences"),
            f"xaxis{plot_ordinal_number}": dict(title_text="p-value of model <= threshold")
        })


        fig.update_layout(height=900, width=1300,
                          title=dict(text=f"<i><b>Case {true_lineage}</b></i><br>"
                                     f"<sup>Contributing lineages {' + '.join(self.gt_contributing_lineages())}<br>"
                                     f"N° breakpoints {self.gt_number_of_breakpoints} "
                                     f"{' '.join(['-'.join([str(u) for u in r]) for r in self.gt_breakpoint_ranges])}"
                                     f"<br>"
                                     f"N° sequences {number_of_sequences}</sup>"),
                          margin=dict(l=50, r=50, t=200, b=50))
        return fig

    def get_table_row(self) -> list:
        # Lineage, n. seq [0-100], % correct recombinant, % correct lineage candidates, % n. bp, %p-value<1e-5
        number_of_sequences = self.batch_experiment_output.shape[0]

        number_is_recombinant = np.sum(self.batch_experiment_output.is_recombinant)
        perc_is_recombinant = number_is_recombinant / number_of_sequences

        number_exact_candidates = np.sum(self.batch_experiment_output.matching_candidates)
        perc_exact_candidates = number_exact_candidates / number_of_sequences

        breakpoint_number_value_counts = self.batch_experiment_output.groupby("number_of_breakpoints", sort=False).size()
        # Series counting the num. of sequences having 0,1,2 breakpoints (actually may include less than all three options)
        br_ = dict.fromkeys([0, 1, 2], 0)  # dict{0:num.seq.0-br, 1:num.seq.1-br, 2:num.seq.2-br}
        for k, v in breakpoint_number_value_counts.items():
            br_[k] = v
        perc_number_0_breakpoints = br_[0] / number_of_sequences
        perc_number_1_breakpoints = br_[1] / number_of_sequences
        perc_number_2_breakpoints = br_[2] / number_of_sequences

        pval_below_threshold = sum(self.batch_experiment_output.p_values.apply(
            lambda pvaluelist: all([x <= GroundTruth.P_VAL_ACCEPT_THRESHOLD for x in pvaluelist])))
        #pval_above_threshold = self.batch_experiment_output.shape[0] - pval_below_threshold
        pval_below_threshold /= number_of_sequences

        return [self.lineage, number_of_sequences, round(perc_is_recombinant, 2), round(perc_exact_candidates, 2),
                round(perc_number_0_breakpoints, 2), round(perc_number_1_breakpoints, 2), round(perc_number_2_breakpoints, 2),
                round(pval_below_threshold, 2)]

    def get_breakpoint_histogram(self, include_case_details=True, Y_AXIS_RESCALED_MAX_VALUE=79,
                                 GT_ADDITIONAL_HEIGHT=4) -> go.Figure:
        colors = {
            "1 BP": ['rgb(0, 99, 197)', 'rgb(199,236,255)'],
            "2 BP": ['rgb(255, 127, 14)', 'rgb(244,220,171)']
        }
        fig = make_subplots(shared_xaxes=True, specs=[[{"secondary_y": True}]])

        consensus_seq_len = len(self.consensus_sequence())

        single_br = list()
        double_br = list()
        for br_tuples in self.batch_experiment_output.breakpoints_in_t.values:
            if len(br_tuples) == 1:
                single_br.append(br_tuples[0][0] + 1)
            elif len(br_tuples) > 1:
                double_br.extend([br_tuple[0] + 1 for br_tuple in br_tuples])

        data = dict()
        c = Counter()
        if single_br:
            data["1 BP"] = single_br
            c.update(single_br)
        if double_br:
            data["2 BP"] = double_br
            c.update(double_br)
        highest_bin_value = np.max(list(c.values()))

        multiplier = int(round(Y_AXIS_RESCALED_MAX_VALUE / highest_bin_value, 0))
        multiplied_highest_bin_value = highest_bin_value * multiplier

        # rh breakpoints on consensus sequence
        gt_br_color = (colors["1 BP"] if RhBreakpointsLocation.breakpoints_num(self.lineage) == 1 else colors["2 BP"])[
            1]
        for br_idx, (s, e) in enumerate(RhBreakpointsLocation.all_breakpoints(self.lineage) or []):
            br_point = s + 1
            fig.add_trace(go.Histogram(
                x=(multiplied_highest_bin_value + GT_ADDITIONAL_HEIGHT) * [br_point], marker_color=gt_br_color,
                xbins=dict(start=1 - 0.5, end=consensus_seq_len + 0.5, size=1)
            ),
                secondary_y=False)
        if not RhBreakpointsLocation.all_breakpoints(self.lineage):
            # plot somethig out of visible box
            fig.add_trace(go.Histogram(
                x=[-10], marker_color=gt_br_color,
                xbins=dict(start=1 - 0.5, end=consensus_seq_len + 0.5, size=1)
            ),
                secondary_y=False)

        for d in data:  # histogram
            fig.add_trace(
                go.Histogram(x=multiplier * data[d], histnorm=None, name=d, marker_color=colors[d][0], nbinsx=10,
                             xbins=dict(start=1 - 0.5, end=consensus_seq_len + 0.5, size=1)
                             ), secondary_y=True)
        # histogram options
        fig.update_layout(bargap=0.5)
        fig.update_layout(barmode='stack')
        fig.update_layout(scattermode="overlay")

        # axes options
        x_axis_plot_range = [0, consensus_seq_len + 1]
        fig.update_layout(xaxis=dict(range=x_axis_plot_range, autorange=False))
        fig.update_xaxes(tickvals=[consensus_seq_len])
        fig.update_xaxes(ticks="outside", tickwidth=2, tickcolor='black', ticklen=8)
        fig.update_yaxes(tickvals=[0, multiplied_highest_bin_value])
        fig.update_yaxes(labelalias={f"{multiplied_highest_bin_value}": f"{highest_bin_value}"})
        fig.update_yaxes(ticks="outside", ticklabelposition="outside bottom", tickwidth=4, tickcolor='black', ticklen=8)
        fig.update_xaxes(tickfont=dict(size=36, color="black"))
        fig.update_yaxes(tickfont=dict(size=36, color="black"))
        fig.update_layout(yaxis=dict(range=[0, multiplied_highest_bin_value + GT_ADDITIONAL_HEIGHT], autorange=False))
        fig.update_layout(yaxis2=dict(range=[0, multiplied_highest_bin_value + GT_ADDITIONAL_HEIGHT], autorange=False))
        fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
        fig.update_yaxes(showline=True, linewidth=2, linecolor='black')
        fig.update_layout(xaxis=dict(showgrid=False), yaxis=dict(showgrid=False))
        fig.update_layout(yaxis2=dict(visible=False))

        # title, background color
        fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))
        fig.update_layout(title_x=0.04, title_y=0.98)
        fig.update_layout(title_font_color="black", title_font_size=44)
        fig.update_layout(title_font_family="Verdana")
        fig.update_layout(showlegend=False)
        fig.update_layout(plot_bgcolor='rgb(255,255,255)')
        fig.update_layout(paper_bgcolor='rgb(255,255,255)')

        if include_case_details:
            plot_title = (f"<i><b>{self.lineage}</b></i><br>"
                          f"<sup>Contributing lineages {' + '.join(self.gt_contributing_lineages())}<br>"
                          f"N° breakpoints {self.gt_number_of_breakpoints} "
                          f"{' '.join(['-'.join([str(u) for u in r]) for r in self.gt_breakpoint_ranges])}"
                          f"<br>"
                          f"N° sequences {self.batch_experiment_output.shape[0]}</sup>")
            fig.update_layout(title=dict(text=plot_title))
            fig.update_layout(margin=dict(l=50, r=50, t=200, b=50))
            fig.update_layout(dict1={
                f"yaxis": dict(tickformat=',d', title_text="number of sequences"),
                f"xaxis": dict(title_text="breakpoint coordinate on target")
            })
            fig.update_layout(showlegend=True)

        return fig

    # ADDITIONAL CUSTOM FEATURES
    def samples_consensus_sequence(self):
        return compute_75_perc_characterization(lists=self.batch_experiment_output.mutations.tolist())

    def consensus_sequence(self):
        return self.predefined_consensus_seq

    def show_batch_sequence_private_mutations(self, environment) -> pd.DataFrame:
        try:
            consensus_seq = set(self.samples_consensus_sequence())
        except KeyError:
            private_mut_series = pd.Series(index=self.batch_experiment_output.index,
                                           data="no consensus sequence", name="private mutations")
            return pd.DataFrame(private_mut_series)
        else:
            def filter_private_mutations(all_mut_list):
                return [f"{m} {environment.change_probability(m):.2e}" for m in all_mut_list if m not in consensus_seq]

            private_mut_series = (
                self.batch_experiment_output.mutations
                .apply(filter_private_mutations))
            return pd.DataFrame(private_mut_series).rename(columns={"mutations": "private mutations"})



class NotRecombinantMixedSequence(GroundTruth):
    def __init__(self, lineage_hierarchy: PangoLineageHierarchy,
                 environment, test_sequence_dataset=None, max_samples=None, shuffle_test_dataset=False,
                 experiment_output_csv=None):
        super().__init__(lineage_hierarchy, environment, test_sequence_dataset, max_samples, shuffle_test_dataset,
                         experiment_output_csv)

    def experiment_output_file_name(self):
        return f'experiment_output_NOT_RECOMBINANT.csv'

    def validation_output_file_name(self):
        return f"validation_output_NOT_RECOMBINANT"

    def gt_contributing_lineages(self, *args) -> list:
        assert not (args[0].startswith('X') and '.' not in args[0]), \
            f"Input test dataset contains recombinant lineage {args[0]}"
        return args[0]

    def get_full_figure_output(self):
        # TODO include check on correctness of contrib lineages
        number_of_sequences = self.batch_experiment_output.shape[0]
        fig = make_subplots(rows=1, cols=2, subplot_titles=["Recombinants detected", "P-values < 1e-05"])

        # MODEL CLASS evaluation
        is_recombinant_number = np.sum(self.batch_experiment_output.is_recombinant)
        is_not_recombinant_number = self.batch_experiment_output.shape[0] - is_recombinant_number

        fig.add_trace(go.Bar(x=["True", "False"], y=[is_recombinant_number, is_not_recombinant_number], name='recombination detected'), row=1, col=1)


        # P-VALUE test
        below_threshold = np.sum(self.batch_experiment_output.p_values.apply(
            lambda pvaluelist: all([x <= GroundTruth.P_VAL_ACCEPT_THRESHOLD for x in pvaluelist])))
        above_threshold = self.batch_experiment_output.shape[0] - below_threshold
        fig.add_trace(go.Bar(x=[True, False], y=[below_threshold, above_threshold], name="p-value test"), row=1, col=2)


        fig.update_layout(height=500, width=900,
                          title=dict(text=f"<i><b>Case NOT RECOMBINANT</b></i><br>"
                                     f"<sup>N° sequences {number_of_sequences}</sup>"),
                          yaxis=dict(tickformat=',d', title_text="number of sequences"),
                          yaxis2=dict(tickformat=',d', title_text="number of sequences"),
                          xaxis=dict(title_text="is recombinant?"),
                          xaxis2=dict(title_text="p-value of model <= threshold"),
                          margin=dict(l=50, r=50, t=200, b=50))
        return fig

    def get_table_row(self) -> list:
        # TODO (1 row for every lineage)
        raise NotImplementedError

    def get_breakpoint_histogram(self) -> go.Figure:
        # even if gt is not recombinant, experiment can be !
        fig = go.Figure()
        breakpoint_number_value_counts = self.batch_experiment_output.groupby("number_of_breakpoints", sort=False).size()
        max_seq_len = np.max(self.batch_experiment_output.seq_len)
        highest_breakpoint_num = max(breakpoint_number_value_counts.index.tolist())
        for br_idx in range(
                highest_breakpoint_num):  # sequences with 1-br -> br_idx=0, 2-br -> br_idx=1, 0 -> self.batch_experiment_output.breakpoints_in_t.values is empty list
            # find 1st breakpoint pair of positions of all sequences having at least one breakpoint
            # (otherwise self.batch_experiment_output.breakpoints_in_t is an empty list and is skipped)
            ith_breakpoint_values = []
            for breakpoints_of_sequence in self.batch_experiment_output.breakpoints_in_t.values:
                try:
                    ith_breakpoint_values.append(breakpoints_of_sequence[br_idx])
                except IndexError:
                    pass
            # unzip pairs
            first_pos_t = [i[0] + 1 for i in ith_breakpoint_values]  # +1 for 1-based conversion
            last_pos_t = [i[1] + 1 for i in ith_breakpoint_values]  # +1 for 1-based conversion
            breakpoint_length_t = [last_pos_t[i] - first_pos_t[i] for i in range(len(last_pos_t))]
            # make x/y values for histogram
            all_pos_of_t = dict.fromkeys(list(range(1, max_seq_len + 1)), 0)
            for p in first_pos_t:
                all_pos_of_t[p] += 1

            fig.add_trace(go.Bar(x=list(all_pos_of_t.keys()), y=list(all_pos_of_t.values()),
                                 name=f"{infl.ordinal(br_idx + 1)} breakpoint start pos"))
        fig.update_layout(dict1={
            f"yaxis": dict(tickformat=',d', title_text="number of sequences"),
            f"xaxis": dict(title_text="breakpoint coordinate on target")
        })
        return fig
