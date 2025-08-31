import logging
from pprint import pformat

from recombinhunt.core.environment import PangoLineageHierarchy
from recombinhunt.core.method import Experiment
from recombinhunt.validation.utils import *


class NotRecombinantMixedSequence:
    def __init__(self, lineage_hierarchy: PangoLineageHierarchy):
        self.lh = lineage_hierarchy
        self.batch_experiment_output = None
        self.experiment = None

    def __str__(self):
        return pformat(self.__dict__, sort_dicts=False)

    def consensus_output_file_name(self):
        return f'consensus_output_NOT_REC'

    #  COMPARISON WITH GT
    def validate_consensus_sequence(self, environment, test_sequence_dataset=None, computed_average_75_percent=None) -> str:
        """
        Computes or reads an averaged sequence using all the nucleotide changes with frequency >= 75%. A single experiment is run and the
        outcome is compared against the given ground truth. The cvomparison is printed alongside the likelihood plot.
        :param environment:
        :param test_sequence_dataset: input dataset
        :param computed_average_75_percent:  pre-computed averaged sequence
        :return: Nothing. The function prints the comparison output.
        """
        assert any([test_sequence_dataset, computed_average_75_percent]) and not all([test_sequence_dataset, computed_average_75_percent]), "One in (test_sequence_dataset,computed_average_75_percent) should not be None"
        if computed_average_75_percent:
            self.consensus_seq = computed_average_75_percent
        else:
            self.consensus_seq = compute_75_perc_characterization(strings=[nuc_changes for name, _, nuc_changes in test_sequence_dataset])

        try:
            exp = Experiment(environment, self.lh)
            exp.set_target(self.consensus_seq)
            exp.run()
        except:
            logging.exception(f"Error while running experiment due to the following reason:")
            raise
        self.experiment = exp

        return self.consensus_seq_comparison(exp)

    def gt_contributing_lineages(self, *args) -> list:
        raise NotImplementedError

    def consensus_seq_comparison(self, exp: Experiment):
        output = [f"Recombination", f"\tdetected: {len(exp.genome_view.regions) > 1}"]
        # model class

        return '\n'.join(output)


