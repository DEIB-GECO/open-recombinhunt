import pandas as pd
import numpy as np
import sys
from os.path import sep
import json
import fastjsonschema
import warnings
import copy


class CandidatesHierarchy:
    def __init__(self, environment):
        pass

    valid_schema = fastjsonschema.compile({
        "type": "object",
        "properties": {
            "class_name": {
                "type": "string",
            },
            "arguments": {
                "type": "array",
            },
        },
        "required": ["class_name", "arguments"],
        "additionalProperties": False,
    })

    def is_totally_different(self, cand1: str, cand2):
        return True
    
    def filter_same_hierarchy_as_first(self, candidate_list):
        return []
        

class Environment:

    ## ENVIRONMENT CONSTRUCTION

    def __init__(self, data_dir_path, ignore_lineage=tuple()):
        self.data_dir_path = data_dir_path
        self.lc_df = self._read_parquet_or_pickle_file("lc_df")
        self.c2lp_df = self._read_parquet_or_pickle_file("change2lineage_probability")
        self.cp_df = self._read_parquet_or_pickle_file("change_probability")
        self.x_characterization = None      # is initialized on request
        self.lc_quality_df = self._read_parquet_or_pickle_file("lc_quality_df")

        self.lc2lp_df = self._read_or_compute_lc2lp_df()
        self.lc2p_df = self._read_or_compute_lc2p_df()
        self.lc_mutations = set(self.lc_df.index)
        self._c2lp_mutations = set(self.c2lp_df.index)
        self._cp_mutations = set(self.cp_df.index)
        
        for l in ignore_lineage:
            try:
                # remove lineages
                self.lc_df.drop(columns=l.upper(), inplace=True)
                self.c2lp_df.drop(columns=l.upper(), inplace=True)
                self.lc2lp_df.drop(columns=l.upper(), inplace=True)
                # TODO remove orphan characteristic changes
            except KeyError:
                pass

        self.ch = self._init_candidates_hierarchy(self.data_dir_path + sep + "candidates_hierarchy.json")     

    def copy_with_exclusions(self, ignore_lineage: tuple):
        ignore_lineage = [l.upper() for l in ignore_lineage if l in self.included_lineages()]
        env = copy.copy(self)
        env.lc_df = self.lc_df.drop(columns=ignore_lineage)
        env.c2lp_df = self.c2lp_df.drop(columns=ignore_lineage)
        env.lc2lp_df = self.lc2lp_df.drop(columns=ignore_lineage)
        env.lc_quality_df = self.lc_quality_df.drop(columns=ignore_lineage)
        return env
    
    def embed_candidates_hierarchy(self, hierarchy: CandidatesHierarchy):
        """
        Method making Environment of RecombinHunt 4.+ backward compatible. 
        If the candidate hierarchy file is not provided in the environment, 
        call this method passing a CandidatesHierarchy object.
        """
        self.ch = hierarchy
    
    ## ENVIRONMENT USAGE
        
    @staticmethod
    def sequence_nuc_mutations2df(sequence_changes):
        temp_df = pd.DataFrame(index=sequence_changes)
        temp_df['seq_pos'] = temp_df.index.str.split('_').str[0]
        temp_df = temp_df.astype({'seq_pos': int})
        return temp_df

    def make_merged_df(self, seq_df):
        """
        Merge lineage characterization DF and sequence characterization DF in a new DF containing changes from both
        as index and sorted by genomic coordinates. Includes, in order, one boolean column for every lineage and two 
        more columns without relative order guaranteed: one boolean column for the sequence (seq_change), a column 
        of genomic coordinates (merg_pos).
        :param seq_df: sequence DF as returned from Data.sequence_nuc_mutations2df
        :return: a DF
        """
        target_mutations = set(seq_df.index)
        target_and_not_lc_mutations = list(target_mutations - self.lc_mutations)

        target_and_not_lc_df = pd.concat([
            pd.DataFrame(
                np.zeros((len(target_and_not_lc_mutations), len(self.included_lineages())), dtype=bool),
                index=target_and_not_lc_mutations,
                columns=self.included_lineages()), 
            seq_df.loc[target_and_not_lc_mutations].rename(columns={'seq_pos':'merg_pos'})],
            axis=1)

        merged_df = pd.concat([self.lc_df.rename(columns={'lc_pos':'merg_pos'}), target_and_not_lc_df])
        merged_df['seq_change'] = merged_df.index.isin(list(target_mutations))
        merged_df = merged_df.sort_index().sort_values(by='merg_pos', kind='stable')             # sort by position but preserve relative index order
        return merged_df

    def probabilities(self, seq_df):
        target_mutations = set(seq_df.index)
        
        # P_MERGED_DF
        # df of seq_df mutations that are not in lc_df (a consequence of the following condition) and not in c2lp_df
        target_not_lc__c2lp_unkown__mutations = list(target_mutations - self._c2lp_mutations)
        target_not_lc__c2lp_unkown__df = pd.concat([
            pd.DataFrame(
                np.zeros((len(target_not_lc__c2lp_unkown__mutations), len(self.included_lineages())), dtype=np.float64),
                index=target_not_lc__c2lp_unkown__mutations,
                columns=self.included_lineages()), 
            seq_df.loc[target_not_lc__c2lp_unkown__mutations].rename(columns={'seq_pos':'merg_pos'})],
            axis=1)
        
        # mutations of seq_df that are not in lc_df but in c2lp_df
        target_not_lc__c2lp_known__mutations = list((target_mutations - self.lc_mutations) & self._c2lp_mutations)
        
        p_merged_df = pd.concat([
            self.lc2lp_df,                                                                                   # lc_changes having a probability in c2lp
            self.c2lp_df.loc[target_not_lc__c2lp_known__mutations].rename(columns={'pos':'merg_pos'}),  # changes in (target but not characteristic) having a probability in c2lp
            target_not_lc__c2lp_unkown__df                                                              # changes not having a probability in c2lp
            ])
        p_merged_df['seq_change'] = p_merged_df.index.isin(list(target_mutations))      # set seq_change column
        p_merged_df = p_merged_df.sort_index()                                          # sort (merg_pos,index)
        p_merged_df = p_merged_df.sort_values(by='merg_pos', kind='stable')             # 'stable' preserve relative index order (faster than reindex(merged_df.index))

        # CHANGE_PROBABILITIES
        # df of seq_df mutations that are not in lc_df (a consequence of the following condition) and not in cp_df
        target_not_lc__cp_unknown__mutations = list(target_mutations - self._cp_mutations)
        target_not_lc__cp_unknown__df = seq_df.loc[target_not_lc__cp_unknown__mutations]
        target_not_lc__cp_unknown__df = target_not_lc__cp_unknown__df.rename(columns={'seq_pos':'pos'})
        target_not_lc__cp_unknown__df['probability'] = 0.0

        # mutations of seq_df that are not in lc_df but in cp_df
        target_not_lc__cp_known__mutations = list((target_mutations - self.lc_mutations) & self._cp_mutations)

        change_probabilities = pd.concat([
            self.lc2p_df,                                              # lc_changes having a probability in c2lp
            self.cp_df.loc[target_not_lc__cp_known__mutations],   # changes in (target but not characteristic) having a probability in cp_df
            target_not_lc__cp_unknown__df,                        # changes not having a probability in cp_df
        ])
        change_probabilities = change_probabilities.sort_index()                            # sort (merg_pos,index)
        change_probabilities = change_probabilities.sort_values(by='pos', kind='stable')    # 'stable' preserve relative index order (faster than reindex(merged_df.index))
        change_probabilities = change_probabilities['probability'].values 
        
        return p_merged_df, change_probabilities

    def change_probability(self, change):
        try:
            return self.cp_df.at[change, 'probability']
        except KeyError:
            return 0.0
        
    def included_lineages(self):
        return self.lc_df.columns[:-1].tolist()
    
    def number_of_sequences_of_lineage(self, lineage_name):
        try:
            if "NmC" in lineage_name: return self.lc_quality_df.at['num', lineage_name]
            else: return self.lc_quality_df.at['num', lineage_name.upper()]
        except KeyError:
            return 0

    ## ENVIRONEMNT CREATION HELPER METHODS   
          
    def _read_parquet_or_pickle_file(self, file_name_without_extension):
        mypath = self.data_dir_path + sep + file_name_without_extension
        try:
            file = pd.read_parquet(mypath + ".parquet")
        except FileNotFoundError:
            file = pd.read_pickle(mypath + ".pickle")
        return file

    def _read_or_compute_lc2lp_df(self):
        mypath = self.data_dir_path + sep + "lc2lp_df"
        try:
            file = pd.read_parquet(mypath + ".parquet")
        except FileNotFoundError:
            try:
                file = pd.read_pickle(mypath + ".pickle")
            except FileNotFoundError:
                file = self.c2lp_df.loc[self.lc_df.index].rename(columns={'pos':'merg_pos'})
        return file
    
    def _read_or_compute_lc2p_df(self):
        mypath = self.data_dir_path + sep + "lc2p_df"
        try:
            file = pd.read_parquet(mypath + ".parquet")
        except FileNotFoundError:
            try:
                file = pd.read_pickle(mypath + ".pickle")
            except FileNotFoundError:
                file = self.cp_df.loc[self.lc_df.index]
        return file

    def __getattribute__(self, item):
        if item == "x_characterization":
            self.__setattr__("x_characterization", self.__get_x_characterization())
        return super(Environment, self).__getattribute__(item)

    def __get_x_characterization(self):
        mypath = self.data_dir_path + sep + "lc_df_with_X"
        try:
            try:
                xlc_df = pd.read_parquet(mypath + ".parquet")
            except FileNotFoundError:
                xlc_df = pd.read_pickle(mypath + ".pickle")
        except FileNotFoundError:
            raise FileNotFoundError(f"Characterization for recombinant lineages not available in the current "
                                    f"environment {self.data_dir_path}")
        xlc_df = xlc_df[[l for l in xlc_df.columns if l.startswith('X')]]
        return xlc_df

    def _init_candidates_hierarchy(self, path):
        try:
            with open(path, "r") as ch_instr_file:
                ch_instructions = json.loads(ch_instr_file.read())
            CandidatesHierarchy.valid_schema(ch_instructions)   # assert the file is properly formatted
            class_obj = getattr(sys.modules[__name__], ch_instructions["class_name"])       # find the referenced class
        except OSError:
            # if file not exists (OSError)
            warnings.warn("RecombinHunt 5.+ expects a candidates_hierarchy.json in the environment dir but the file is missing. Backward compatibility is possible by"
                    " calling Environment.embed_candidates_hierarchy() with a CandidatesHierarchy object. Resorting to default CandidatesHierarchy.")
            return CandidatesHierarchy(None)
        except fastjsonschema.JsonSchemaException:
            # if candidates_hierarchy.json does not respect the required schema
            warnings.warn("Invalid format of candidates_hierarchy.json. Check the correct definition in CandidatesHierarchy.valid_schema. Resorting to default CandidatesHierarchy.")
            return CandidatesHierarchy(None)
        except AttributeError:  
            # if unknown class name
            warnings.warn(f'Class {ch_instructions["class_name"]} not implemented. Resorting to default CandidatesHierarchy.')
            return CandidatesHierarchy(None)
        else:
            # instantiate the class
            init_args = ch_instructions["arguments"]
            return class_obj(*init_args, environment=self) if init_args else class_obj(environment=self)

    ## MISCELLANEA

    def approximate_memory_usage_bytes(self, pretty_print=True):
        """
        Returns the summed memory usage of the inner pandas and numpy objects in bytes. Objects generated after initialization and CandidatesHierarchy are not considered.

        :param pretty_print: prints the memory usage using the closest measurement unit (KB/MB/GB).
        """
        lc_df_mu: pd.Series = self.lc_df.memory_usage(index=True, deep=True).sum()
        c2lp_df_mu: pd.Series = self.c2lp_df.memory_usage(index=True, deep=True)
        cp_df_mu: pd.Series = self.cp_df.memory_usage(index=True, deep=True)
        lc_quality_df_mu: pd.Series = self.lc_quality_df.memory_usage(index=True, deep=True)
        total = sum([
            # lc_df          c2lp_df          cp_df           _c2lp_mutations     _cp_mutations     lc_quality_df
            lc_df_mu.sum(), c2lp_df_mu.sum(), cp_df_mu.sum(), c2lp_df_mu.iloc[0], cp_df_mu.iloc[0], lc_quality_df_mu.sum()
        ])
        if pretty_print:
            unit_measure = {
                0: "Bytes", 1: "KB", 2: "MB", 3: "GB", 4: "TB" 
            }
            mult = 0
            close_unit = 1000**mult
            while close_unit <= total:
                mult += 1
                close_unit = 1000**mult
            mult -= 1
            close_unit = 1000**mult
            clean_total = total / close_unit
            print(f"{clean_total:.3f}", unit_measure[mult])
        return total
    
    def x_characterizing_nuc_mutations(self, x_lineage_name):
        if x_lineage_name.upper() not in self.x_characterization.columns:
            raise KeyError(f"Lineage characterization for {x_lineage_name} is not available in the environment "
                           f"{self.data_dir_path}")
        return self.x_characterization.index[self.x_characterization[x_lineage_name.upper()]].values.tolist()

    ## VALIDITY
    def assert_valid(self):
        print('Starting formal checks of the environment...')
        """
        Verifies that tables 'lc_df', 'change2lineage_probability' and 'change2probability' contain the right number and type of columns. 
        Verifies that all values are defined.
        Verifies that mutations in table 'lc_df' are contained in those of 'change2lineage_probability'.
        Verifies that mutations in table 'change2lineage_probability' are contained in 'change2probability'.
        Warns if exist any mutations of table 'lc_df' that is not characteristic to any lineage. 
        """
        self.formal_check_lc_df()
        self.formal_check_change2lineage_probability()
        self.formal_check_change2probability()
        self.check_consistency_of_variants()
        self.check_consistency_of_probabilities()

    def formal_check_lc_df(self):
        """
        Check that table 'lc_df' contains a column 'lc_pos' of type 'int64' and that all preceeding columns are of type 'bool'. No other columns are contained in the table. 
        Check that all vaules are defined (not null).
        Issue a warning if not all mutations in table 'lc_df' are characteristic of at least one variant. 
        """
        # columns lc_pos existing and with correct type
        assert 'lc_pos' in self.lc_df.columns, "Column 'lc_pos' is missing from table 'lc_df'"
        assert pd.api.types.is_integer_dtype(self.lc_df.lc_pos.dtype) or np.issubdtype(self.lc_df.lc_pos.dtype, np.integer), "Column 'lc_pos' in table 'lc_df' is not of integer type"
        # No Null values
        assert not self.lc_df.isna().stack().any(), "Table 'lc_df' must not contain null values."
        # contains columns other than lc_pos
        assert self.lc_df.shape[1] > 1, "Table 'lc_df' should contain one column named 'lc_pos' and other additonal columns for the variants"
        # no columns after lc_pos
        assert self.lc_df.columns[-1] == 'lc_pos', "Table 'lc_df' should not contain other columns after the column 'lc_pos'"
        # columns preceeding lc_pos are of type bool
        lc_pos_col_index = self.lc_df.columns.get_indexer_for(['lc_pos'])[0]
        preceeding_col_types = self.lc_df.dtypes[:lc_pos_col_index]
        assert (preceeding_col_types == 'bool').all(), "Columns before 'lc_pos' in table 'lc_df' must be of type 'bool'."
        assert len(preceeding_col_types) > 0, "Columns named as the variants should precede the column 'lc_pos' in table 'lc_df'."
        # every mutation is characteristic to at least one lineage
        mutations_not_used = self.lc_df.index[self.lc_df[self.lc_df.columns[:lc_pos_col_index]].sum(axis=1) == 0]
        if len(mutations_not_used) > 0:
            warnings.warn(f"{len(mutations_not_used)} / {self.lc_df.shape[0]} mutations listed in table 'lc_df' are not characteristic of any variant. To improve efficiency, it is recommended to remove unused mutations from the table index.")

    def formal_check_change2lineage_probability(self):
        """
        Check that table 'change2lienage_probability' contains a column 'pos' of type 'int64' and that all preceeding columns are of type 'float64'. No other columns are contained in the table. 
        Check that all vaules are defined (not null).
        Check that probabilities are defined between 0 and 1. 
        """
        # columns pos existing and with correct type
        assert 'pos' in self.c2lp_df.columns, "Column 'pos' is missing from table 'change2lineage_probability'"
        assert pd.api.types.is_integer_dtype(self.c2lp_df.pos.dtype) or np.issubdtype(self.c2lp_df.pos.dtype, np.integer), "Column 'pos' in table 'change2lineage_probability' is not of intger type"
        # No Null values
        assert not self.c2lp_df.isna().stack().any(), "Table 'change2lineage_probability' must not contain null values."
        # contains columns other than pos
        assert self.c2lp_df.shape[1] > 1, "Table 'change2lineage_probability' should contain one column named 'pos' and other additonal columns for the variants"
        # no columns after pos
        assert self.c2lp_df.columns[-1] == 'pos', "Table 'change2lineage_probability' should not contain other columns after the column 'pos'."
        # columns preceeding pos are of type float64
        pos_col_index = self.c2lp_df.columns.get_indexer_for(['pos'])[0]
        preceeding_col_types = self.c2lp_df.dtypes[:pos_col_index]
        assert all([pd.api.types.is_float_dtype(col_dtype) or np.issubdtype(col_dtype, np.floating) for col_dtype in preceeding_col_types]), "Columns before 'pos' in table 'change2lineage_probability' must be of type float."
        assert len(preceeding_col_types) > 0, "Columns named as the variants should precede the column 'pos' in table 'change2lineage_probability'."
        # probability values are between 0 and 1
        assert self.c2lp_df[self.c2lp_df.columns[:pos_col_index]].stack().between(0, 1, inclusive='both').all(), "Table 'change2lineage_probability' contains invalid probability values (not enclosed between 0 and 1)."
         
    def formal_check_change2probability(self):
        """
        Check that table 'change2probability' contains exactly two columns: 'probability' and 'pos' of type float64 and int64 respectively.
        Check that all vaules are defined (not null).
        Check that probabilities are defined between 0 and 1. 
        """
        # columns pos existing and with correct type
        assert 'pos' in self.cp_df.columns, "Column 'pos' is missing from table 'change_probability'"
        assert pd.api.types.is_integer_dtype(self.cp_df.pos.dtype) or np.issubdtype(self.cp_df.pos.dtype, np.integer), "Column 'pos' in table 'change2_probability' is not of intger type"
        # column probability exisits with correct dtype
        assert 'probability' in self.cp_df.columns, "Column 'probability' is missing from table 'change_probability'"
        assert pd.api.types.is_float_dtype(self.cp_df.probability.dtype) or np.issubdtype(self.cp_df.probability.dtype, np.floating), "Column 'pos' in table 'change_probability' is not of type float"
        # no other columns
        assert self.cp_df.shape[1] == 2, "Table 'change_probability' should not contain columns other than 'probability' and 'pos'"
        # No Null values
        assert not self.cp_df.isna().stack().any(), "Table 'change_probability' must not contain null values."
        # probability values are between 0 and 1
        assert self.cp_df.probability.between(0, 1, inclusive='both').all(), "Table 'change_probability' contains invalid probability values (not enclosed between 0 and 1)."

    def check_consistency_of_variants(self):
        """
        Check that variants in tables 'lc_df' and 'change2lienage_probability' match exactly.
        """
        lc_variants = self.lc_df.columns[:-1]
        c2lp_variants = self.c2lp_df.columns[:-1]
        try:
            pd.testing.assert_index_equal(lc_variants, c2lp_variants, check_order=True, check_names=False)
        except AssertionError:
            print("Names and order of variants must match in tables 'lc_df', 'change2lineage_probability'.")
            raise

    def check_consistency_of_probabilities(self):
        """
        Check that mutations in table 'lc_df' are defined also in tables 'change2lienage_probability' and 'change2probability'.
        Check that mutations in table 'change2lienage_probability' are defined also in table 'change2probability'.
        """
        lc_df_mutations = set(self.lc_df.index)
        c2lp_mutations = set(self.c2lp_df.index)
        cp_mutations = set(self.cp_df.index)
        assert lc_df_mutations.issubset(c2lp_mutations) and lc_df_mutations.issubset(cp_mutations), "Mutations in table 'lc_df' must be defined also in tables 'change2lineage_probability' and 'change_probability'."
        assert c2lp_mutations.issubset(cp_mutations), "Mutations in table 'change2lineage_probability' must be defined also in table 'change_probability'."


class PangoLineageHierarchy(CandidatesHierarchy):
    def __init__(self, alias_key_file_path, environment=None):
        super().__init__(environment)
        if environment is not None:         # then it's created from Envrionment
            alias_key_file_path = environment.data_dir_path + sep + alias_key_file_path
        with open(alias_key_file_path, "r") as alias_key_file:
            self.alias_key = json.load(alias_key_file)
        # remove mappings to ""
        self.alias_key = {x:y for x,y in self.alias_key.items() if y != ""}
        #self.invert_alias_key = self.__prepare_invert_mapping()

    # def __prepare_invert_mapping(self) -> dict:
    #     # remove mappings to lists (recombinant lineage mappings)
    #     alias_key = {x:y for x,y in self.alias_key.items() if not isinstance(y, list)}
    #     # invert mapping
    #     invert_mapping = dict()
    #     for x, y in alias_key.items():
    #         mapped_aliases = invert_mapping.get(y, list())
    #         mapped_aliases.append(x)
    #         invert_mapping[y] = mapped_aliases
    #     return invert_mapping
    #
    # def aliases_of(self, lineage:str) -> list:
    #     ris = self.invert_alias_key.get(lineage.upper(), [])
    #     if not res1:
    #         res2 = self.alias_key.get(lineage.upper(), [])
    #         return res2

    def unfold_lineage(self, lineage: str):
        lineage = lineage.upper()
        if lineage.startswith("X"):
            return lineage
        prefix, sep, postfix = lineage.partition(".")
        # exception for lineages expressed as AY* (star not preceded by dot)
        prefix_no_star, star, _ = prefix.partition("*")
        return self.alias_key.get(prefix_no_star, prefix_no_star) + sep + postfix + star

    def is_sublineage(self, lin1: str, of_lin2: str) -> bool:
        """
        Returns True if the first is a child or nephew of the second lineage; False otherwise.
        Example:
        X vs X.1 -> False
        X.1 vs X -> True
        X.1 vs X.2 -> False
        X.1 vs X*/X.* -> True
        X.1 vs X.1*/X.1.* -> False (X.1 matches exactly X.1*/X.1.*  -- see is_same_lineage())
        :param lin1:
        :param of_lin2:
        :return: bool
        """
        child_name = self.unfold_lineage(lin1.rstrip(".*"))
        putative_parent_name = self.unfold_lineage(of_lin2.rstrip(".*"))
        return putative_parent_name in child_name and (child_name != putative_parent_name)

    def is_same_lineage(self, lin1: str, as_lin2: str) -> bool:
        """
        Returns True if the two lineages match exactly; False otherwise.
        Example:
        X vs X -> True
        X vs X.1 -> False
        X vs X*/X.* -> True
        X vs X.1*/X.1.* -> False
        :param lin1:
        :param as_lin2:
        :return:
        """
        l1 = self.unfold_lineage(lin1.rstrip(".*"))
        l2 = self.unfold_lineage(as_lin2.rstrip(".*"))
        return l1 == l2

    def is_superlineage(self, lin1: str, of_lin2: str) -> bool:
        """
        Returns True if the first lineage is an ancestor of the second lineage; False otherwise.
        Example:
        X vs X -> False
        X vs X.1 -> True
        X vs X*/X.* -> False
        X vs X.1*/X.1.* -> True
        :param lin1:
        :param of_lin2:
        :return:
        """
        child_name = self.unfold_lineage(of_lin2.rstrip(".*"))
        putative_parent_name = self.unfold_lineage(lin1.rstrip(".*"))
        return putative_parent_name in child_name and (child_name != putative_parent_name)

    def is_totally_different(self, lin1: str, lin2):
        lin1 = self.unfold_lineage(lin1.rstrip(".*"))
        lin2 = self.unfold_lineage(lin2.rstrip(".*"))
        return lin1 != lin2 and not (lin1 in lin2 or lin2 in lin1)

    def hierarchy(self, lin1: str, with_lin2: str):
        lin1 = self.unfold_lineage(lin1.rstrip(".*"))
        with_lin2 = self.unfold_lineage(with_lin2.rstrip(".*"))
        if lin1 == with_lin2:
            return "equal"
        elif lin1 in with_lin2:
            return "super-lineage"
        elif with_lin2 in lin1:
            return "sub-lineage"
        else:
            return "different"

    def hierarchy_distance(self, lin1: str, with_lin2: str):
        lin1 = self.unfold_lineage(lin1.rstrip(".*"))
        with_lin2 = self.unfold_lineage(with_lin2.rstrip(".*"))
        if lin1 == with_lin2:   # lin 1 == lin 2
            return 0
        elif lin1 in with_lin2: # lin 1 is super-lineage of lin 1
            return 2
        elif with_lin2 in lin1: # lin 1 is sub-lineage of lin 2
            return 1
        else:                   # otherwise
            return 3

    def is_matching_candidate(self, candidate_lineage: str, true_lineage: str) -> bool:
        # Notice: this method is equal to is_same_lineage or (is_sublineage and true_lineage[-1] = '*')
        # It has been rewritten for reasons of efficiency
        candidate_lineage = self.unfold_lineage(candidate_lineage.rstrip(".*"))
        true_lineage = self.unfold_lineage(true_lineage.rstrip(".*"))
        return (
            candidate_lineage == true_lineage                                                   # == is_same_lineage
            or (
                    (true_lineage in candidate_lineage and candidate_lineage != true_lineage)   # == is_sublineage
                    and true_lineage[-1] == "*"
            )
        )

    def filter_same_hierarchy_as_first(self, lineage_list):
        return [l for l in lineage_list if self.hierarchy_distance(lineage_list[0], l) < 3]

    def filter_same_hierarchy_as(self, lineage_list, as_lineage):
        return [l for l in lineage_list if self.hierarchy_distance(as_lineage, l) < 3]

    def farthest_ancestor_in_list(self, from_sub_lineage, lineage_list):
        ancestor = from_sub_lineage
        for lin in lineage_list:
            if self.is_superlineage(lin, ancestor):
                ancestor = lin
        return ancestor


# if __name__ == "__main__":
#     lh = PangoLineageHierarchy("../../validation_data/alias_key.json")

#     ##### TEST CASES
#     print("\t", "\t".join(["sub", "exact", "matching_sol", "super"]))
#     for l1 in ["x", "x.1", "x.2"]:
#         print(l1, "vs")
#         for l2 in ["x", "x.1", "x.2", "x*", "x.*", "x.1*", "x.1.*"]:
#             print("\t", l2, lh.is_sublineage(l1,l2), lh.is_same_lineage(l1, l2), lh.is_matching_candidate(l1, l2), lh.is_superlineage(l1,l2))

#     # for l1,l2 in [("x", "x"), ("x", "x.1"), ("x", "x.2"), ("x", "x*"), ("x", "x.*"), ("x", "x.1*"), ("x", "x.1.*"),
#     #               ("x.1", "x"), ("x.1", "x.1"), ("x.1", "x.2"), ("x.1", "x*"), ("x.1", "x.*"), ("x.1", "x.1*"), ("x.1", "x.1.*"),
#     #               ("x.2", "x"), ("x.2", "x.1"), ("x.2", "x.2"), ("x.2", "x*"), ("x.2", "x.*"), ("x.2", "x.1*"), ("x.2", "x.1.*")]:
#     #     print(l1,"vs", l2, "->", lh.is_sublineage(l1,l2), lh.is_same_lineage(l1, l2), is_matching_candidate(l1, l2))


# if __name__ == "__main__":
#     with open("/Users/tom/Developer/recombinhunt-dev/environments/env_2023_04_11/candidates_hierarchy.json") as ch_instr_file:
#         ch_instructions = json.loads(ch_instr_file.read())
#     print(ch_instructions)

#     CandidatesHierarchy.DECLARATION_SCHEMA(ch_instructions)

#     # fastjsonschema.validate({'type': 'array'}, ["a", "b"])
#     # fastjsonschema.validate({'type': 'array'}, 2)
#     # fastjsonschema.validate({'type': 'array'}, [])
