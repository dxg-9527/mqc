import cobra
import pandas as pd
import re

from mqc.defaults import *
from ..utils import add_rxn
from mqc.config import Config
from mqc.control.formula_charge_pairing import missing_input_charges, complete_reference_pair


class Preprocess:
    """
    Model preprocessing related.

    Attributes
    ----------
    diff_results : python.Dictionary
        The dictionary structure of memote.MemoteResult objects.
    configuration : memote.MemoteConfiguration
        A memote configuration structure.

    """

    def __init__(self, file_path: str, cfg: Config):
        """
        load initial model.

        Parameters
        ----------
        file_path : str
            model file

        """
        self.cfg = cfg
        self.pattern = r"\([^()]*\)n"
        try:
            if file_path.endswith(".json"):
                self.model = cobra.io.load_json_model(file_path)
                self.check_model = cobra.io.load_json_model(file_path)
            elif file_path.endswith(".yaml"):
                self.model = cobra.io.load_yaml_model(file_path)
                self.check_model = cobra.io.load_yaml_model(file_path)
            elif file_path.endswith(".mat"):
                self.model = cobra.io.load_matlab_model(file_path)
                self.check_model = cobra.io.load_matlab_model(file_path)
            else:
                self.model = cobra.io.read_sbml_model(file_path)
                self.check_model = cobra.io.read_sbml_model(file_path)
        except:
            self.model = ""
            return None
        self._input_charge_missing = missing_input_charges(file_path, self.model)
        missing_input_charges(file_path, self.check_model)

    def add_meta_formula(self, meta_met):
        """
        Add formula to meta.

        Parameters
        ----------
        model : cobra.Model

        Return
        ------
        cobra.Model
            Model after adding formula

        """
        df = pd.read_excel(self.cfg.meta_met)
        meta_formula = dict(zip(df["meta_id"], df["formula"]))
        for metId in meta_met:
            met = self.model.metabolites.get_by_id(metId)
            if metId.split("@")[0] in list(df["meta_id"]):
                if met.formula and met.formula != "X":
                    met.formula = re.sub(self.pattern, "", met.formula)
                else:
                    met.formula = str(meta_formula[metId.split("@")[0]])
            if not met.compartment:
                met.compartment = metId.split("@")[-1]

    def add_exchange_rxn(self, met_exchange):
        """
        Add formula to meta.

        Parameters
        ----------
        model : cobra.Model

        Return
        ------
        cobra.Model
            Model after adding formula

        """
        h_c, h2o_c, pi_c, nh4_c, o2_c, co2_c, so4_c = "", "", "", "", "", "", ""
        for met in self.model.metabolites:
            if met.id in H:
                h_c = met.id
            if met.id in H2O:
                h2o_c = met.id
            if met.id in PI:
                pi_c = met.id
            if met.id in NH4:
                nh4_c = met.id
            if met.id in O2:
                o2_c = met.id
            if met.id in CO2:
                co2_c = met.id
            if met.id in SO4:
                so4_c = met.id
        print(
            h_c,
            h2o_c,
            pi_c,
            nh4_c,
            o2_c,
            co2_c,
            so4_c,
            "add_exchange_rxn.................",
        )
        rxn_name = [
            "add_h2o_exchange",
            "add_nh4_exchange",
            "add_pi_exchange",
            "add_h_exchange",
            "add_o2_exchange",
            "add_co2_exchange",
            "add_so4_exchange",
        ]
        rxn_exp = [
            h2o_c + " <=>",
            nh4_c + " <=>",
            pi_c + " <=>",
            h_c + " <=>",
            o2_c + " <=>",
            co2_c + " <=>",
            so4_c + " <=>",
        ]
        for i in range(7):
            if rxn_name[i] == met_exchange:
                add_rxn(self.model, rxn_name[i], rxn_exp[i])
                add_rxn(self.check_model, rxn_name[i], rxn_exp[i])

    def add_seed_formula(self, seed_met):
        """
        Add formula to modelseed.

        Parameters
        ----------
        model : cobra.Model

        Returns
        -------
        cobra.Model
            Model after adding formula

        """
        df = pd.read_excel(self.cfg.modelseed_met)
        modelseed_formula = dict(zip(df["id"], df["formula"]))
        modelseed_name = dict(zip(df["id"], df["name"]))
        modelseed_charge = dict(zip(df["id"], df["charge"]))
        for metId in seed_met:
            met = self.model.metabolites.get_by_id(metId)
            if metId.split("_")[0] in list(df["id"]):
                met.name = modelseed_name[met.id.split("_")[0]]
                if met.formula and met.formula != "X":
                    met.formula = re.sub(self.pattern, "", met.formula)
                else:
                    met.formula = str(modelseed_formula[met.id.split("_")[0]])
                if met.charge:
                    pass
                else:
                    if not pd.isna(modelseed_charge[met.id.split("_")[0]]):
                        met.charge = int(modelseed_charge[met.id.split("_")[0]])
            if not met.compartment:
                met.compartment = metId.split("_")[-1]

    def add_virtual_formula(self, virtual_met):
        """
        Add formula to modelseed.

        Parameters
        ----------
        model : cobra.Model

        Returns
        -------
        cobra.Model
            Model after adding formula

        """
        df = pd.read_excel(self.cfg.virtual_met)
        virtual_formula = dict(zip(df["ID"], df["Formula"]))
        virtual_name = dict(zip(df["ID"], df["Name"]))
        virtual_charge = dict(zip(df["ID"], df["Charge"]))
        for metId in virtual_met:
            met = self.model.metabolites.get_by_id(metId)
            for ids in df["ID"]:
                if metId.split("[")[0] == ids.split("[")[0]:

                    if met.formula and met.formula != "X":
                        met.formula = re.sub(self.pattern, "", met.formula)
                    else:
                        met.formula = str(virtual_formula[ids])
                    if met.charge:
                        pass
                    else:
                        if not pd.isna(virtual_charge[ids]):
                            met.charge = int(virtual_charge[ids])
                    if met.name:
                        pass
                    else:
                        met.name = virtual_name[ids]
            if not met.compartment:
                met.compartment = metId.split("[")[-1].split("]")[0]

    def change_virtual_to_bigg(self, virtual_met):
        """Convert bracketed compartment suffixes to underscore-separated identifiers."""
        for target in (self.model, self.check_model):
            for metId in virtual_met:
                if metId in target.metabolites:
                    new_id = re.sub(r"\[(.*?)\]", r"_\1", metId)
                    if new_id in target.metabolites:
                        raise ValueError('Conflicting bracket/underscore identities: ' + metId)
        for target in (self.model, self.check_model):
            for metId in virtual_met:
                if metId in target.metabolites:
                    target.metabolites.get_by_id(metId).id = re.sub(r"\[(.*?)\]", r"_\1", metId)

    def add_kegg_formula(self, kegg_met):
        """Add formula to kegg"""
        table = pd.read_excel(self.cfg.kegg_met)
        kegg_formula = dict(zip(table["id"], table["formula"]))
        kegg_name = dict(zip(table["id"], table["name"]))
        kegg_charge = dict(zip(table["id"], table["charge"]))
        for metId in kegg_met:
            reference_id = metId.split('_')[0]
            if reference_id in kegg_formula:
                met = self.model.metabolites.get_by_id(metId)
                met.name = kegg_name[reference_id]
                if met.formula and met.formula != "X":
                    met.formula = re.sub(self.pattern, "", met.formula)
                else:
                    met.formula = str(kegg_formula[reference_id])
                if met.charge:
                    pass
                else:
                    if not pd.isna(kegg_charge[reference_id]):
                        met.charge = int(kegg_charge[reference_id])

    def add_other_formula(self, other_met):
        """Add formula to metacyc"""
        bigg_metacyc_met = pd.read_excel(self.cfg.bigg_metacyc)
        bigg_metacyc_formula = dict(zip(bigg_metacyc_met['id'], bigg_metacyc_met['formula']))
        bigg_metacyc_charge = dict(zip(bigg_metacyc_met['id'], bigg_metacyc_met['charge']))
        bigg_metacyc_name = dict(zip(bigg_metacyc_met['id'], bigg_metacyc_met['name']))
        for metId in other_met:
            met = self.model.metabolites.get_by_id(metId)
            charge_missing = metId in self._input_charge_missing
            if metId in list(bigg_metacyc_met['id']):      # 代谢物ID在总表里面，就可以对formula进行添加
                met.name = bigg_metacyc_name[metId]
                if met.formula and met.formula != 'X': # formula存在时检测是否有问题
                    met.formula = re.sub(self.pattern, "", met.formula)
                else:
                    met.formula = str(bigg_metacyc_formula[metId])
                if charge_missing:
                    charge_missing = not complete_reference_pair(
                        met, bigg_metacyc_formula[metId], bigg_metacyc_charge[metId])
            else:  # 代谢物ID不在总表里面，需要对formula进行检测，将'C10H17O8PR2(C5H8O5PR)n'这种不合法formula变为'C10H17O8PR2'
                if met.formula and met.formula != 'X':
                    met.formula = re.sub(self.pattern, "", met.formula)     
            if not met.formula:
                for ids in list(bigg_metacyc_met['id']):
                    if metId.replace('_', '') == ids.replace('_', '') and ids.count('_') == metId.count('_') + 1:
                        if '__' in ids:
                            met.formula = str(bigg_metacyc_formula[ids])
            if charge_missing:
                for ids in list(bigg_metacyc_met['id']):
                    if metId.replace('_', '') == ids.replace('_', '') and ids.count('_') == metId.count('_') + 1:
                        if '__' in ids:
                            if complete_reference_pair(met, bigg_metacyc_formula[ids], bigg_metacyc_charge[ids]):
                                break
            if not met.compartment:
                met.compartment = metId.split('_')[-1]
                print(met.id,'add_other_formula.................',met.compartment)
