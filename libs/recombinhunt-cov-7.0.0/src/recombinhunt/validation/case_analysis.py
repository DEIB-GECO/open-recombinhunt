import sys
from itertools import cycle, chain
from typing import List, Optional
from IPython.core.display import display_html
from tabulate import tabulate
import json
from pathlib import Path

from recombinhunt.core.environment import PangoLineageHierarchy
from recombinhunt.core.method import Experiment
from recombinhunt.validation.utils import *
from recombinhunt.core.graphics import *
from recombinhunt.core.method import Region


class CaseAnalysis:
    def __init__(self, experiment: Experiment,
             case_name: str, number_of_sequences: int, case_number: int, case_group_name: str,
             lineage_hierarchy: PangoLineageHierarchy, acl: Optional[AssessedContributingLin] = None,
             gt_candidates: Optional[List[str]] = None, gt_genomic_breakpoints: Optional[tuple] = None):
        self.exp = experiment
        self.case_name = case_name
        self.n_seq = number_of_sequences
        self.lh = lineage_hierarchy
        self.cl = acl

        # Description
        self.group_name = case_group_name
        self.case_number = case_number

        self.gt_genomic_breakpoints = gt_genomic_breakpoints  # breakpoints of GT in genomic coordinates
        self.gt_n_breakpoints = None
        self.gt_genomic_breakpoints_str = None

        self.gt_target_breakpoints = None  # 1-based breakpoints of GT in T
        self.gt_target_breakpoints_str = None  # string version

        self.gt_candidates = gt_candidates
        self.gt_candidates_str = None

        self.target_nuc_changes = None
        self.changes_number = None

        # Experiment
        self.best_candidates_str = None
        self.alternative_candidates_str = None
        self.dir_l1 = None
        self.rank_L1_L2 = None
        self.bc_target_breakpoints_str = None  # 1-based string version of best candidate breakpoints in t
        self.bc_genomic_breakpoints_str = None  # 1-based string version of best candidate genomic breakpoints
        self.OK_KO = None
        self.gap = None
        self.initial_R_span = None
        self.PV_candidato1 = None
        self.PV_candidato2 = None

        self._analyse()

    def _analyse(self):
        try:
            # --- GT Independent Analysis (from experiment results) ---
            self.target_nuc_changes = target_series = self.exp.merged_df.seq_change
            self.target_nuc_changes = target_series[target_series].index.tolist()
            self.changes_number = len(self.target_nuc_changes)

            self.best_candidates_str = " + ".join(self.exp.genome_view.contributing_lineages())

            alternative_candidates = self.exp.genome_view.list_good_alternative_candidates()
            alternative_candidates = [f'[{", ".join(alt_list)}]' for alt_list in alternative_candidates]
            self.alternative_candidates_str = ", ".join(alternative_candidates)

            self.dir_l1 = '<<' if self.exp.L1_dir == 1 else '>>'

            Br_BC_tup = self.exp.genome_view.breakpoints_in_t()
            Br_BC = [self.format_as_range_1_based(x) for x in Br_BC_tup]
            self.bc_target_breakpoints_str = ", ".join(Br_BC)
            self.bc_genomic_breakpoints_str = ', '.join([' - '.join([str(x) for x in tup]) for tup in self.exp.genome_view.breakpoints()])

            Gap_tup = []
            for gaps in self.exp.genome_view.gap_history_pos_t:
                Gap_tup.append(", ".join([self.format_as_range_1_based(x) for x in gaps if not x[1] == x[0] + 1]))
            self.gap = ' -> '.join(Gap_tup)

            initial_R_span = [[x[0], x[1] - 1] for x in
                              self.exp.genome_view.region_lengths_pre_gap_resolution]  # make edge-included
            initial_R_span = [self.format_as_range_1_based(x) for x in initial_R_span]
            self.initial_R_span = ','.join(initial_R_span)

            try:
                self.PV_candidato1 = self.format_pvalue(
                    self.exp.p_values.get(("recombinant_model", self.exp.genome_view.contributing_lineages()[0])))
            except IndexError:
                self.PV_candidato1 = "-"
            try:
                self.PV_candidato2 = self.format_pvalue(
                    self.exp.p_values.get(("recombinant_model", self.exp.genome_view.contributing_lineages()[1])))
            except IndexError:
                self.PV_candidato2 = "-"

            

            # --- GT Dependent Analysis (from GT information) CONDITIONAL---
            if self.gt_candidates is not None and self.gt_genomic_breakpoints is not None:
                # This block runs only if Ground Truth is provided
                true_parents = self.cl.contributing_to(self.case_name) # Uses case_name
                self.gt_candidates = ' + '.join(true_parents)
                self.gt_n_breakpoints = len(self.gt_genomic_breakpoints)
                self.rank_L1_L2 = ' '.join(rank_gt_in_regions(self.exp.genome_view.regions, true_parents, self.lh))
                self.OK_KO = "OK" if all_candidates_matching(self.exp.genome_view.regions, true_parents, self.lh) else "KO"

                if BreakpointsLocation.is_unknown(self.gt_genomic_breakpoints):
                    self.gt_target_breakpoints_str = "-"
                    self.gt_genomic_breakpoints_str = "-"
                else:
                    self.gt_target_breakpoints = [BreakpointsLocation.to_target_pos(br, self.target_nuc_changes) for br in self.gt_genomic_breakpoints]
                    self.gt_target_breakpoints = [self.format_as_range_1_based(x) for x in self.gt_target_breakpoints]
                    self.gt_target_breakpoints_str = ", ".join(self.gt_target_breakpoints)
                    self.gt_genomic_breakpoints_str = ', '.join(
                        [' - '.join([str(x) for x in tup]) for tup in self.gt_genomic_breakpoints])
            else:
                # This block runs if Ground Truth is NOT provided
                self.gt_candidates = "N/A"
                self.gt_n_breakpoints = None # or 0
                self.rank_L1_L2 = "N/A"
                self.OK_KO = "N/A"
                self.gt_target_breakpoints_str = "N/A"
                self.gt_genomic_breakpoints_str = "N/A"

        except Exception as e:
            print(f"Error while analysing case number {self.case_number} {self.gt_lineage}")
            raise e
        


    def get_issue(self):
        
        if self.exp.genome_view.number_of_breakpoints() == 0:
            return "issue_0BP"
        elif self.gt_n_breakpoints > self.exp.genome_view.number_of_breakpoints():
            return "issue_1BP"
        elif self.gt_n_breakpoints < self.exp.genome_view.number_of_breakpoints():
            return "issue_2BP"
        elif self.OK_KO == "OK":
            return "non_issue_OK"
        else:
            return "issue_KO"

    def analysis_table_row(self):
        return [self.case_number, self.case_name, self.n_seq, self.changes_number, self.group_name,
                self.gt_candidates, self.OK_KO, self.initial_R_span, self.gap,
                self.best_candidates_str, self.dir_l1, self.rank_L1_L2,
                self.bc_target_breakpoints_str, self.gt_target_breakpoints_str,
                self.PV_candidato1, self.PV_candidato2,
                self.alternative_candidates_str]

    def region_table(self, region: Region):
        num_seq = [self.exp.env.number_of_sequences_of_lineage(c) for c in region.candidates]

        max_likelihood_pos = [region.candidate_max_likelihood_pos_t(c) + 1 for c in region.candidates]
        max_likelihood = region.max_likelihood_value
        likelihood_in_bc_max_pos = [region.candidate_likelihood_at_dir_end(c) for c in region.candidates]
        aic_candidate = [region.candidate_aic(c) for c in region.candidates]
        p_val_candidate = [region.candidate_p_value(c) for c in region.candidates]
        pval_in_threshold = ["*" if p >= self.exp.ALT_CANDIDATE_P_VALUE_DIFFERENCE else "" for p in p_val_candidate]

        p_val_candidate[0] = None
        likelihood_in_bc_max_pos[0] = None
        aic_candidate[0] = None

        close_alternative_candidates = region.alternative_candidates_with_max_likelihood_pos_t_distance_below(self.exp.ALT_CANDIDATE_MAX_POS_DISTANCE_T)
        close_max_likelihood_pos = ["*" if c in close_alternative_candidates or c == region.designated else "" for c in
                                    region.candidates]

        phylogenetic_ok_candidates = region.alternative_candidates_in_same_branch()
        same_phylogenetic_tree_of_BC = ["*" if c in phylogenetic_ok_candidates or c == region.designated else "" for c
                                        in region.candidates]
        return pd.DataFrame(
            {
                "num_seq": num_seq,
                "t_ch_MAX": max_likelihood_pos,
                "max_CL": max_likelihood,
                "CL@BC_t_ch_MAX": likelihood_in_bc_max_pos,
                "aic": aic_candidate,
                "PV": p_val_candidate,
                "PV_OK": pval_in_threshold,
                "t_ch_MAX_OK": close_max_likelihood_pos,
                "phyl_OK": same_phylogenetic_tree_of_BC

            },
            index=[c for c in region.candidates]
        )

    def print_case_details(self, out):
        # Intestazione
        out.write('<pre style="font-size:1.5em;">')
        
        # Conditionally display the 'test' status
        test_status_str = f"\t\t\t\ttest: {self.OK_KO}" if self.OK_KO != "N/A" else ""
        out.write(f"<h3> Case {self.case_number} ({self.group_name}): {self.case_name}{test_status_str}</h3><br>")
        out.write("</pre>")

        # render information body as table
        out.write('<pre style="font-size:1.5em;">')
        exp_1BPvs2BPmodel_comparison = '-' if not self.exp.p_val_partial_models else f"{self.exp.best_model_key}: {self.format_pvalue(self.exp.p_val_partial_models[self.exp.best_model_key])}"
        
        ### previous version of table
        # table = [
        #     [f"<u>Target</u>: {self.n_seq} samples", f"<u>Number of changes</u>: {self.changes_number}"],
        #     [f"<u>GT</u>: {self.gt_candidates}", f"<u>GT BR</u>: {self.gt_target_breakpoints_str}", f"<u>GT BR coord</u>: {self.gt_genomic_breakpoints_str}", f"<u>Rank_L1_L2</u>: {self.rank_L1_L2}"],
        #     [f"<u>BC</u>: {self.best_candidates_str}", f"<u>BC BR</u>: {self.bc_target_breakpoints_str}", f"<u>BC BR coord</u>: {self.bc_genomic_breakpoints_str}"],
        #     [f"<u>Direction L1:</u> {self.dir_l1}", f"<u>Initial region span</u>: {self.initial_R_span}", f"<u>Gap history (edge excluded)</u>: {self.gap}"],
        #     [f"<u>Alt. candidates</u>: {self.alternative_candidates_str}"],
        #     [f"<u>Model 1BP/2BP comparison</u>: ", exp_1BPvs2BPmodel_comparison],
        #     [f"<u>Rec. model vs L1</u>: {self.PV_candidato1}", f"<u>Rec. model vs L2</u>: {self.PV_candidato2}"],
        #     [f"<u>Flags</u>: {', '.join(self.exp.get_flags())}"]
        # ]

        table = []
        # Always add the Target row
        table.append([f"<u>Target</u>: {self.n_seq} samples", f"<u>Number of changes</u>: {self.changes_number}"])
        
        # Conditionally add the Ground Truth row
        if self.gt_candidates != "N/A":
            table.append([f"<u>GT</u>: {self.gt_candidates}", f"<u>GT BR</u>: {self.gt_target_breakpoints_str}", f"<u>GT BR coord</u>: {self.gt_genomic_breakpoints_str}", f"<u>Rank_L1_L2</u>: {self.rank_L1_L2}"])

        # Always add the Best Candidate (BC) row and subsequent rows
        table.append([f"<u>BC</u>: {self.best_candidates_str}", f"<u>BC BR</u>: {self.bc_target_breakpoints_str}", f"<u>BC BR coord</u>: {self.bc_genomic_breakpoints_str}"])
        table.append([f"<u>Direction L1:</u> {self.dir_l1}", f"<u>Initial region span</u>: {self.initial_R_span}", f"<u>Gap history (edge excluded)</u>: {self.gap}"])
        table.append([f"<u>Alt. candidates</u>: {self.alternative_candidates_str}"])
        table.append([f"<u>Model 1BP/2BP comparison</u>: ", exp_1BPvs2BPmodel_comparison])
        table.append([f"<u>Rec. model vs L1</u>: {self.PV_candidato1}", f"<u>Rec. model vs L2</u>: {self.PV_candidato2}"])
        table.append([f"<u>Flags</u>: {', '.join(self.exp.get_flags())}"])

        out.write(tabulate(table, tablefmt='unsafehtml'))
        out.write("</pre>")

        # out.write('<pre style="font-size:1.5em;">')
        # out.write(f"{self.n_seq} samples\t\t\t\tNumber of changes: {self.changes_number}<br>")
        # out.write(f"GT: {self.gt_candidates}\t\t\t\t{self.gt_target_breakpoints_str}\t\t\t\t{self.gt_genomic_breakpoints_str}\t\t\t\tRank_L1_L2 {self.rank_L1_L2}<br>")
        # out.write(f"BC: {self.best_candidates_str}\t\t\t\t{self.bc_target_breakpoints_str}\t\t\t\t{self.bc_genomic_breakpoints_str}\t\t\t\tL1_dir {self.dir_l1}\t\t\t\t<br>")
        # out.write(f"Initial_R_span {self.initial_R_span}\t\t\t\tGap history (edge excluded) {self.gap}<br>")
        # out.write(f"Alt. candidates: {self.alternative_candidates_str}<br>")
        # exp_1BPvs2BPmodel_comparison = '\t' + '\t\t\t\t'.join(
        #     [f'{"* " if mod_key == self.exp.best_model_key else ""}{mod_key}: {mod_pval}' for mod_key, mod_pval in
        #      self.exp.p_val_partial_models.items()])
        # out.write(f"1BP vs 2BP test: {'-' if not self.exp.best_model_key else exp_1BPvs2BPmodel_comparison}<br>")
        # out.write(f"P_L1 {self.PV_candidato1}\t\t\t\tP_L2 {self.PV_candidato2}<br>")
        # out.write(f"Flags: {', '.join(self.exp.get_flags())}<br>")
        # out.write("</pre>")

        html_tables_side_by_side(
            *[self.region_table(r) for r in self.exp.genome_view.regions],
            titles=[f"{r.designated} {'>>' if r.search_dir == 0 else '<<'}" for r in self.exp.genome_view.regions],
            out=out
        )

        fig = plot_likelihood(self.exp.genome_view, xaxis="changes", changes="target")
        # reduce plot size
        fig.update_layout(
            title="Cumulative Likelihood per-region",
            autosize=False,
            width=1600,
            height=600,
            xaxis=dict(showgrid=False),
            yaxis=dict(showgrid=False)
            # margin=dict(
            #     l=50,
            #     r=50,
            #     b=100,
            #     t=100,
            #     pad=4
            # ),
            # paper_bgcolor="LightSteelBlue",
        )
        fig.write_html(out, include_plotlyjs="cdn", include_mathjax="cdn", auto_open=False)
        # fig.write_image(out, format='svg')
        out.write("<br>")
        fig = plot_likelihood_whole_genome(self.exp.genome_view, xaxis="changes", changes="target")
        # reduce plot size
        fig.update_layout(
            title="Cumulative Likelihood whole genome",
            autosize=False,
            width=1600,
            height=600,
            # margin=dict(
            #     l=50,
            #     r=50,
            #     b=100,
            #     t=100,
            #     pad=4
            # ),
            # paper_bgcolor="LightSteelBlue",
            xaxis=dict(showgrid=False),
            yaxis=dict(showgrid=False)
        )
        fig.write_html(out, include_plotlyjs="cdn", include_mathjax="cdn", auto_open=False)
        # fig.write_image(out, format='svg')
        out.write("<br>")

        # Target sequence
        out.write("<h5>Target sequence</h5><br>")
        target_nuc_changes_str = ""
        last_line_length = 0
        for _change in self.target_nuc_changes:
            temp_str = target_nuc_changes_str + ", " + _change
            if len(temp_str) - last_line_length > 190:
                target_nuc_changes_str += ", <br>"
                last_line_length = len(target_nuc_changes_str)
                target_nuc_changes_str += _change
            else:
                target_nuc_changes_str = temp_str
        out.write(target_nuc_changes_str + "<br>")

    def save_structured_report(self, output_dir_path: Path):
        """
        Saves the key results of the case analysis into separate, machine-readable files.
        """
        try:
            # 1. Create the output directory
            output_dir_path.mkdir(parents=True, exist_ok=True)

            # 2. Save the main summary table as a JSON file
            exp_1BPvs2BPmodel_comparison = '-' if not self.exp.p_val_partial_models else f"{self.exp.best_model_key}: {self.format_pvalue(self.exp.p_val_partial_models[self.exp.best_model_key])}"
            summary_data = {
                "case_number": self.case_number,
                "case_name": self.case_name,
                "group_name": self.group_name,
                "number_of_sequences": self.n_seq,
                "number_of_changes": self.changes_number,
                "best_candidates": self.best_candidates_str,
                "best_candidates_breakpoints_target": self.bc_target_breakpoints_str,
                "best_candidates_breakpoints_genomic": self.bc_genomic_breakpoints_str,
                "direction_L1": self.dir_l1,
                "initial_region_span": self.initial_R_span,
                "gap_history": self.gap,
                "alternative_candidates": self.alternative_candidates_str,
                "model_comparison": exp_1BPvs2BPmodel_comparison,
                "p_value_vs_L1": self.PV_candidato1,
                "p_value_vs_L2": self.PV_candidato2,
                "flags": ', '.join(self.exp.get_flags())
            }
            with open(output_dir_path / "summary.json", 'w') as f:
                json.dump(summary_data, f, indent=4)

            # 3. Save each region table as a separate CSV file
            for i, region in enumerate(self.exp.genome_view.regions):
                region_df = self.region_table(region)
                region_df.to_csv(output_dir_path / f"region_{i+1}_table.csv")

            # 4. Save the Plotly graphs as JSON data
            fig_per_region = plot_likelihood(self.exp.genome_view, xaxis="changes", changes="target")
            fig_whole_genome = plot_likelihood_whole_genome(self.exp.genome_view, xaxis="changes", changes="target")
            
            with open(output_dir_path / "plot_per_region.json", 'w') as f:
                f.write(fig_per_region.to_json())
            
            with open(output_dir_path / "plot_whole_genome.json", 'w') as f:
                f.write(fig_whole_genome.to_json())

            # 5. Save the target mutations to a text file
            with open(output_dir_path / "target_mutations.txt", 'w') as f:
                f.write(",".join(self.target_nuc_changes))

        except Exception as e:
            print(f"Error saving structured report for case {self.case_number} ({self.case_name}): {e}")

    @staticmethod
    def format_pvalue(val):
        if val is not None:
            return f"{val:.2e}"
        else:
            return "-"

    @staticmethod
    def format_as_range_1_based(values, parentheses_on_close_values=False):
        one_based_values = [x + 1 for x in values]
        string_values = '-'.join([str(x) for x in one_based_values])
        if parentheses_on_close_values and values[1] == (values[0] + 1):
            return f"({string_values})"
        else:
            return string_values


def html_tables_side_by_side(*args, titles=cycle(['']), out=sys.stdout):
    items = []
    for df, title in zip(args, chain(titles, cycle(['</br>']))):
        html_str = '<div style="display:inline-block">'
        html_str += f'<h2 style="text-align: center;">{title}</h2>'
        html_str += '<th style="text-align:center"><td style="vertical-align:top">'
        html_str += df.to_html()  # .replace('table','table style="display:inline-block"')
        html_str += '</td></th></div>'
        items.append(html_str)
    html_str = '<div style="display:inline-block"><span style="opacity:0;">|||||</span></div>'.join(items)  # +spacing between items
    html_str += '</br>'
    if out == sys.stdout or out is None:
        display_html(html_str, raw=True)
    else:
        print(html_str, file=out)


# if __name__ ==  '__main__':
#     cl = AssessedContributingLin("../validation_data/alias_key.json")
#     lh = PangoLineageHierarchy("../validation_data/alias_key.json")
#     env = Environment("../environments/env_2023_04_11")
#
#     exp = Experiment(env, lh)
#     exp.set_target(env.x_characterizing_nuc_mutations("XBB"))
#     exp.run()
#
#     ca = CaseAnalysis(exp, "XBB", lh, cl)




    """
            ###### PREVIOUS VERSION ######
            # Description
            self.target_nuc_changes = target_series = self.exp.merged_df.seq_change
            self.target_nuc_changes = target_series[target_series].index.tolist()
            self.changes_number = len(self.target_nuc_changes)
            #self.gt_n_breakpoints = len(self.gt_genomic_breakpoints)
            #self.gt_candidates = ' + '.join(self.cl.contributing_to(self.case_name))

            # Experiment
            self.best_candidates_str = " + ".join(self.exp.genome_view.contributing_lineages())

            alternative_candidates = self.exp.genome_view.list_good_alternative_candidates()
            alternative_candidates = [f'[{", ".join(alt_list)}]' for alt_list in alternative_candidates]
            self.alternative_candidates_str = ", ".join(alternative_candidates)

            self.dir_l1 = '<<' if self.exp.L1_dir == 1 else '>>'
            #self.rank_L1_L2 = rank_gt_in_regions(self.exp.genome_view.regions, self.cl.contributing_to(self.case_name), self.lh)
            #self.rank_L1_L2 = ' '.join(self.rank_L1_L2)

            Br_BC_tup = self.exp.genome_view.breakpoints_in_t()
            Br_BC = [self.format_as_range_1_based(x) for x in Br_BC_tup]
            self.bc_target_breakpoints_str = ", ".join(Br_BC)
            self.bc_genomic_breakpoints_str = ', '.join([' - '.join([str(x) for x in tup]) for tup in self.exp.genome_view.breakpoints()])

            if BreakpointsLocation.is_unknown(self.gt_genomic_breakpoints):
                self.gt_target_breakpoints_str = "-"
                self.gt_genomic_breakpoints_str = "-"
            else:
                self.gt_target_breakpoints = [BreakpointsLocation.to_target_pos(br, self.target_nuc_changes) for br in self.gt_genomic_breakpoints]
                self.gt_target_breakpoints = [self.format_as_range_1_based(x) for x in self.gt_target_breakpoints]
                self.gt_target_breakpoints_str = ", ".join(self.gt_target_breakpoints)
                self.gt_genomic_breakpoints_str = ', '.join(
                    [' - '.join([str(x) for x in tup]) for tup in self.gt_genomic_breakpoints])

            self.OK_KO = "OK" if all_candidates_matching(self.exp.genome_view.regions,
                                                         self.cl.contributing_to(self.case_name),
                                                         self.lh) else "KO"

            # Gap
            Gap_tup = []
            for gaps in self.exp.genome_view.gap_history_pos_t:
                Gap_tup.append(", ".join([self.format_as_range_1_based(x) for x in gaps if not x[1] == x[0] + 1]))
            self.gap = ' -> '.join(Gap_tup)

            # Region lengths
            initial_R_span = [[x[0], x[1] - 1] for x in
                              self.exp.genome_view.region_lengths_pre_gap_resolution]  # make edge-included
            initial_R_span = [self.format_as_range_1_based(x) for x in initial_R_span]
            self.initial_R_span = ','.join(initial_R_span)

            try:
                self.PV_candidato1 = self.format_pvalue(
                    self.exp.p_values.get(("recombinant_model", self.exp.genome_view.contributing_lineages()[0])))
            except IndexError:
                self.PV_candidato1 = "-"
            try:
                self.PV_candidato2 = self.format_pvalue(
                    self.exp.p_values.get(("recombinant_model", self.exp.genome_view.contributing_lineages()[1])))
            except IndexError:
                self.PV_candidato2 = "-"
        except Exception as e:
            print(f"Error while analysing case number {self.case_number} {self.gt_lineage}")
            raise e
    """