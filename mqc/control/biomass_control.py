import cobra
import time
import json
import numpy as np
import cplex
from collections import OrderedDict
from copy import deepcopy
from mqc.config import Config
from mqc.utils import *
from mqc.control.rules_control import Rules
from mqc.control.preprocessing_control import Preprocess
from mqc.control.model_control import ModelPreprocess
from mqc.control.gapfill_chemistry import assess_gapfill_chemistry, record_gapfill_assessment
from mqc.control.gap_compartments import target_inventory, foreign_compartments
from mqc.control.proton_duplicate import native_proton_index, proton_duplicate_conflict
from mqc.control.boundary_diagnostics import unresolved_boundary_source
from mqc.control.auxiliary_supply import allocate_auxiliary_supply
from mqc.control.gap_rule_policy import (
    enabled as rule_policy_enabled, apply_native_rules, configure_candidates,
    candidate_copy, rule_pfba,
)


class Biomasses:

    def search_pathway(self, model, target, model_info):
        """Keep the native primary optimum; use rule costs for the pFBA tie-break."""
        if model.slim_optimize() > 0:
            solution = rule_pfba(model, model_info.get('gapfill_candidate_costs', {}))
            return list(solution.fluxes[abs(solution.fluxes) > 1e-6].keys())
        return []

    @staticmethod
    def is_uncoupled_redox_shortcut(reaction):
        """Recognize bare reduced/oxidized carrier conversions, not enzymes.

        Such two-metabolite conversions have no electron donor/acceptor. They
        must not enter a biological model as gap-filling reactions. Real coupled
        redox reactions and transport of the same carrier are not excluded.
        """
        if len(reaction.metabolites) != 2 or not reaction.reactants or not reaction.products:
            return False
        aliases = [
            (NADH_NAME + ['nadh'], NAD_NAME + ['nad', 'nad+']),
            (NADPH + ['nadph'], NADP + ['nadp', 'nadp+']),
            (FADH2 + ['fadh2'], FAD + ['fad']),
            (FMNH2 + ['fmnh2'], FMN + ['fmn']),
        ]
        metabolites = list(reaction.metabolites)
        labels = [{(met.name or '').lower(), met.id.lower(), met.id.split('_')[0].lower()} for met in metabolites]
        for reduced, oxidized in aliases:
            red, ox = {s.lower() for s in reduced}, {s.lower() for s in oxidized}
            if (labels[0] & red and labels[1] & ox) or (labels[1] & red and labels[0] & ox):
                return True
        return False

    def __init__(self, cfg: Config):
        """
        Get model information and check_model object.
        """
        self.cfg = cfg

    def recover_bio_objective(self, model_info, model):
        """"""
        model.objective = model_info["initial_rxn"]["biomass_rxn_id"]

    def test(self, model, check_model):
        """"""

    def get_biomass_info(self, model_info, model_control_info):
        """"""
        model_info["biomass_original"] = model_info["initial_rxn"]["biomass_rxn_flux"]
        temp = [
            f"{model_info['initial_rxn']['biomass_rxn_id']} : {model_info['initial_rxn']['biomass_rxn_exp']}\n"
        ]
        temp.extend([f"Initial growth rate : {model_info['biomass_original']}\n"])
        temp.extend([model_info["initial_rxn"]["check_biomass_composition"]])
        write_biomass(model_control_info, temp, "Biomass_information")

    def add_rely_biomass_rxn(
        self, model_info, model, check_model, model_control_info, initial_rxn_id
    ):
        """"""
        temp, rely_rxn = [f"biomass-dependent responses:"], []
        with model:
            for rxn in model_info["reactions"]:
                if (
                    "nadhs_modify" in rxn.keys()
                    or "atps_modify" in rxn.keys()
                    or "nets_modify" in rxn.keys()
                    or "yields_modify" in rxn.keys()
                ):
                    model.reactions.get_by_id(rxn["id"]).bounds = (
                        check_model.reactions.get_by_id(rxn["id"]).bounds
                    )
            for rxn in model_info["reactions"]:
                if (
                    "nadhs_modify" in rxn.keys()
                    or "atps_modify" in rxn.keys()
                    or "nets_modify" in rxn.keys()
                    or "yields_modify" in rxn.keys()
                ):
                    model.objective = initial_rxn_id
                    model.reactions.get_by_id(rxn["id"]).bounds = (0, 0)
                    print(
                        rxn["id"],
                        model.reactions.get_by_id(rxn["id"]).bounds,
                        model.slim_optimize(),
                        "//////////////////////////",
                    )
                    if model.slim_optimize() <= 1e-6:
                        rely_rxn.append(rxn["id"])
                        rxn["rely_biomass"] = "true"
                    model.reactions.get_by_id(rxn["id"]).bounds = (
                        check_model.reactions.get_by_id(rxn["id"]).bounds
                    )

        print("rely_biomass_rxn:............................")
        rely_rxn = list(set(rely_rxn))
        for rxnId in rely_rxn:
            model.reactions.get_by_id(rxnId).bounds = check_model.reactions.get_by_id(
                rxnId
            ).bounds
            print(rxnId, model.reactions.get_by_id(rxnId).bounds)

        if rely_rxn:
            temp.extend(rely_rxn)
            model_control_info["check_biomass_production"]["rely_biomass_rxn"] = (
                " ".join(temp)
            )

    def get_general_library(self, model):
        """"""
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':读取总库')
        meta_met = []
        # controler = Preprocess(self.cfg.general_library_path,self.cfg)
        # with controler.model:
        #     model_pre_process = ModelPreprocess(self.cfg)
        #     model_pre_process.preprocess_met(controler)
        #     for rxn in controler.model.reactions:
        #         if 'C' in rxn.check_mass_balance():
        #             ks.append(rxn.id)
        # print(ks,'...................')
        for met in model.metabolites:
            if "@" in met.id:  # meta
                meta_met.append(met.id)
        if len(meta_met)/len(model.metabolites) > 0.7:
            library_path = self.cfg.general_library_path_meta
            general_model = cobra.io.read_sbml_model(library_path)
            print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':读取meta库')
        else:
            library_path = self.cfg.general_library_path
            general_model = cobra.io.read_sbml_model(library_path)
            print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':读取其他库')
        from mqc.control.library_charge import repair_library_charges
        self.library_charge_pairs = repair_library_charges(
            general_model, library_path, getattr(self.cfg, 'bigg_metacyc', None))
        from mqc.control.metabolite_chemistry import repair_known_charges
        self.library_charge_errata = repair_known_charges(general_model)
        from mqc.control.flavodoxin_identity import repair_model as repair_flavodoxin
        self.library_flavodoxin_repairs = repair_flavodoxin(general_model)
        if getattr(self.cfg, 'repair_single_electron_templates', False):
            from mqc.control.electron_stoichiometry import repair_model
            # The library is an in-memory copy. Resolve shared identities as
            # COBRA will on insertion; never change the reference file.
            resolved = {m.id: m for m in general_model.metabolites}
            resolved.update({m.id: m for m in model.metabolites})
            self.electron_library_repairs = repair_model(general_model, resolved=resolved)
        else:
            self.electron_library_repairs = []
        # Exclude diagnostic-like redox shortcuts before any gap-filling LP.
        # Filtering only carbon imbalance would admit e.g. NADH <=> NAD+.
        self.rejected_gap_redox_shortcuts = [r.id for r in general_model.reactions
                                            if self.is_uncoupled_redox_shortcut(r)]
        general_model.remove_reactions(self.rejected_gap_redox_shortcuts)
        # for k in ks:
        #     model.reactions.get_by_id(k).remove_from_model()
        # cobra.io.write_sbml_model(model, "/home/dengxiao/mqc/mqc/summary/test.xml")
        # return cobra.io.load_json_model(self.cfg.general_library_path)
        return general_model

    def add_sink_rxn(self, big_model, model_info, check_model):
        """"""
        sink_rxn = []
        initial_rxn_id = model_info["initial_rxn"]["biomass_rxn_id"]
        initial_rxn = big_model.reactions.get_by_id(initial_rxn_id)
        for met in initial_rxn.products:
            if (
                met.name not in (ADP_NAME + H_NAME + PI_NAME + H2O_NAME + PPI_NAME)
                and "biomass" not in met.name.lower()
                and f"SK_{met.id}" not in big_model.reactions
            ):
                big_model.add_boundary(metabolite=met, type="sink")
                check_model.add_boundary(metabolite=met, type="sink")
                sink_rxn.append(f"SK_{met.id}")
                big_model.reactions.get_by_id(f"SK_{met.id}").bounds = (0, 1000)
                check_model.reactions.get_by_id(f"SK_{met.id}").bounds = (0, 1000)
                print(
                    "Metabolites receiving sink reactions:",
                    met,
                    "............................",
                )
        model_info["sink_rxn"] = sink_rxn

    def get_big_general_library(self, model_info, model, general_library):
        """"""
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':获取大模型')
        initial_rxn_id = model_info["initial_rxn"]["biomass_rxn_id"]
        model_info['electron_stoichiometry_library'] = getattr(self, 'electron_library_repairs', [])
        model_info['metabolite_charge_errata_library'] = getattr(self, 'library_charge_errata', [])
        model_info['reference_library_charge_pairs'] = getattr(self, 'library_charge_pairs', {})
        model_info['flavodoxin_identity_library'] = getattr(self, 'library_flavodoxin_repairs', [])
        model_file = model_info["model_file"]
        # The stored path is an output, not a valid cache key for model state.
        # Later repair calls may have changed bounds, chemistry or reactions.
        # Rebuild from that current state; retain the original S*1000 search.
        gap_model = model.copy()
        with gap_model as model:
            for rxns in model.reactions:  # 给系数扩大1000倍,添加的demand反应不改
                temp_dict = {}
                for met in rxns.metabolites:
                    temp_dict[met] = rxns.get_coefficient(met) * (1000)
                rxns.add_metabolites(temp_dict, combine=False) # 相同时会覆盖掉原代谢物。乘以1000是为了pfba计算
            # for rxns in model.reactions:
            #     if rxns.id not in general_library.reactions and len(rxns.products) != 0:
            #         general_library.add_reaction(rxns) 
            add_rxn, no_add_rxn = [], ['NO2t2r','NARK','NO3t2','NO2t3','NO3t','NAR_syn_2','NOR_syn_2','NO3t7','NO3t3','NO2t2rpp','NTRIRy_1','FNOR_1','HYDA','KGD2','GLUSfx','NDH_1_1_um','NDH_1_1_um','NDH_1_4_um','AMMQT8']
            chemistry = {"policy": "known_non_hydrogen_atoms_and_formula_coverage",
                         "decision_stage": "target_model_after_metabolite_resolution",
                         "construction": "target_resolved_prefilter_then_batch_add_and_recheck",
                         "experimental_rates_used": False, "assessments": 0,
                         "counts": {}, "records": [], "accepted_candidate_count": 0,
                         "limitations": "H/charge-only differences require separate review; no automatic proton, water or formula corrections."}
            model_info["gapfill_chemistry"] = chemistry
            # Match COBRA's ID resolution, including the first reference
            # definition for metabolites absent from the original target.
            # This read-only map also covers rejected candidates, preserving
            # the annotation semantics of add-all-then-filter construction.
            resolved = {met.id: met for met in model.metabolites}
            # Freeze BEFORE adding library metabolites: rejected candidates
            # must not authorize novel compartments for later candidates.
            allowed_compartments = set(target_inventory(model))
            compartment_audit = dict(allowed=sorted(allowed_compartments), excluded=[])
            model_info['gapfill_compartment_policy'] = compartment_audit
            proton_index = native_proton_index(model)
            model_info['gapfill_proton_duplicate_exclusions'] = []
            for rxns in general_library.reactions:
                if rxns.id not in model.reactions and rxns.id != initial_rxn_id and len(rxns.metabolites) != 1 and rxns.id not in no_add_rxn:
                    conflict = foreign_compartments(rxns, allowed_compartments, resolved)
                    if conflict:
                        compartment_audit['excluded'].append(conflict)
                        continue
                    assessment = assess_gapfill_chemistry(rxns)
                    record_gapfill_assessment(chemistry, assessment, "reference")
                    for met in rxns.metabolites:
                        resolved.setdefault(met.id, met)
                    assessment = assess_gapfill_chemistry(rxns, resolved_metabolites=resolved)
                    record_gapfill_assessment(chemistry, assessment, "target_model")
                    conflict = proton_duplicate_conflict(rxns, assessment, resolved, proton_index)
                    if conflict:
                        model_info['gapfill_proton_duplicate_exclusions'].append(conflict)
                        continue
                    if assessment["allowed"]:
                        add_rxn.append(rxns.copy())
            # Repeated one-reaction solver updates and undo records are costly
            # for a universal library. Batch only candidates already assessed
            # using target identities; retain a post-merge safety assertion.
            model.add_reactions(add_rxn)
            # self.convert_nan_to_null(model)
            delect_rxn = []
            for candidate in add_rxn:
                rxn = model.reactions.get_by_id(candidate.id)
                assessment = assess_gapfill_chemistry(rxn)
                record_gapfill_assessment(chemistry, assessment, "post_merge")
                if not assessment["allowed"]:
                    delect_rxn.append(rxn.id)
            model.remove_reactions(delect_rxn)
            if rule_policy_enabled(self.cfg):
                configure_candidates(model, [r.id for r in add_rxn if r.id in model.reactions], model_info)
            chemistry["accepted_candidate_count"] = len(add_rxn) - len(delect_rxn)
            # for rxnId in delect_rxn:
            #     model.reactions.get_by_id(rxnId).remove_from_model()
            convert_nan_to_null(model)
            if model_file.endswith(".json"):
                cobra.io.save_json_model(model, self.cfg.big_general_library_json)
                self.cfg.big_general_library = self.cfg.big_general_library_json
            elif model_file.endswith(".yaml"):
                cobra.io.save_yaml_model(model,self.cfg.big_general_library_yaml)
                self.cfg.big_general_library = self.cfg.big_general_library_yaml
            elif model_file.endswith(".mat"):
                cobra.io.save_matlab_model(model,self.cfg.big_general_library_mat)
                self.cfg.big_general_library = self.cfg.big_general_library_mat
            else:
                cobra.io.write_sbml_model(model, self.cfg.big_general_library_xml)
                self.cfg.big_general_library = self.cfg.big_general_library_xml
            # try:
            #     cobra.io.write_sbml_model(model, self.cfg.big_general_library_xml)
            #     # file_path = self.cfg.big_general_library_xml
            #     self.cfg.big_general_library = self.cfg.big_general_library_xml
            # except:
            #     cobra.io.save_json_model(model, self.cfg.big_general_library_json)
            #     # file_path = self.cfg.big_general_library_json
            #     self.cfg.big_general_library = self.cfg.big_general_library_json
            # print('大模型路径:',self.cfg.big_general_library)
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':获取大模型完毕')
        return model

    def get_big_model(self, model_info, model, general_library):
        """"""

        model = self.get_big_general_library(model_info, model, general_library)
        model_file = model_info["model_file"]
        print(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            ": Loading the merged model",
        )
        if model_file.endswith(".json"):
            big_model = cobra.io.load_json_model(self.cfg.big_general_library)
        elif model_file.endswith(".yaml"):
            big_model = cobra.io.load_yaml_model(self.cfg.big_general_library)
        elif model_file.endswith(".mat"):
            big_model = cobra.io.load_matlab_model(self.cfg.big_general_library)
        else:
            big_model = cobra.io.read_sbml_model(self.cfg.big_general_library)

        print(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            ": Merged model loaded",
        )
        return big_model

    def bio_precursor_gapfilling(self, model_info, big_model, no_synthesized_mets, model_control_info):
        """"""
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':看大模型中有多少个物质不能合成')
        to_gap, no_gap, model_info['no_gap_metid'] = [], [], []
        # test = ["cpd15665_c0","cpd00166_c0","cpd15668_c0","cpd15669_c0","cpd15667_c0"]
        temp = [f"After quality control, the initial rate of biomass is 0, the metabolite {','.join(no_synthesized_mets)} of biomass components cannot be synthesized"]
        for metId in no_synthesized_mets:
            with big_model:
                objectiveId = add_demand(model_info, big_model, metId, independent=True)
                big_model.objective = objectiveId
                ini_flux = big_model.slim_optimize()
                if pd.isna(ini_flux) or ini_flux <= 1e-6:
                    no_gap.append(metId)
                    continue
                # if metId in test:
                #     write_flux_file(model_info, big_model, objectiveId, metId)
                to_gap.append(metId)
        # if no_gap:
        #     temp.extend([f"In the original model, the precursor substance {','.join(no_synthesized_mets)} of biomass cannot be synthesized, but in the larger model,{len(no_gap)} metabolites cannot be synthesized by gapfilling: {','.join(no_gap)}"])
        # else:
        #     temp.extend([f"In the original model, the precursor substance {','.join(no_synthesized_mets)} of biomass cannot be synthesized, but in the larger model, all substances can be synthesized by gapfilling"])
        # temp.extend([f"Gap is applied to the substances that can be synthesized: {to_gap}"])
        write_biomass(model_control_info, temp, "Gapfilling_of_biomass_components")
        model_info['no_gap_metid'] = no_gap
        return to_gap

    def get_exchange_rxn(
        self, model_info, big_model, no_synthesized_mets, model_control_info
    ):
        """"""

        to_gap_metid = self.bio_precursor_gapfilling(
            model_info, big_model, no_synthesized_mets, model_control_info
        )
        print(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            ": Identifying exchanges for synthesizable metabolites",
        )
        if to_gap_metid:
            reaction = cobra.Reaction("EX_gap")
            for met in to_gap_metid:
                reaction.add_metabolites({big_model.metabolites.get_by_id(met): -1})
            reaction.bounds = (0, 1)
            big_model.add_reactions([reaction])
            return reaction.id
        else:
            return ""

    def get_new_way_rxn(
        self, model_info, big_model, no_synthesized_mets, model_control_info
    ):
        """"""
        exchange_rxnid = self.get_exchange_rxn(
            model_info, big_model, no_synthesized_mets, model_control_info
        )
        print(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            ": Finding a target pathway in the merged model for initially infeasible biomass",
        )
        new_way_rxn = []
        if exchange_rxnid:
            big_model.objective = exchange_rxnid
            big_model.reactions.get_by_id(exchange_rxnid).bounds = (
                0,
                0.1,
            )  # Set an upper bound on the demand reaction.
            if big_model.slim_optimize() > 0:
                pfba_solution = rule_pfba(big_model, model_info.get('gapfill_candidate_costs', {}))
                need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-6]

                new_way_rxn = list(need_fluxes.keys())
            return new_way_rxn
        else:
            return []

    def gapfilling(self, model, model_info, demand_id, general_library):
        """"""
        big_model = self.get_big_model(model_info, model, general_library)
        print(
            f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}: Pathway for {demand_id} in the merged model"
        )
        new_way_rxn = []
        with big_model:
            rxns = model.reactions.get_by_id(demand_id).copy()
            # A fresh merged network already includes this temporary query
            # at native S*1000 scale. COBRA ignores a duplicate ID on add;
            # reset its equation so the intended unscaled demand is used.
            # Keep the solver variables: remove/re-add of an objective ID
            # inside a COBRA context can break restoration of that objective.
            if demand_id in big_model.reactions:
                query = big_model.reactions.get_by_id(demand_id)
                query.subtract_metabolites(dict(query.metabolites))
                query.add_metabolites(rxns.metabolites)
            else:
                big_model.add_reactions([rxns])
            big_model.objective = demand_id
            big_model.reactions.get_by_id(demand_id).bounds = (
                0,
                0.1,
            )  # Set an upper bound on the demand reaction.
            if big_model.slim_optimize() > 0:
                pfba_solution = rule_pfba(big_model, model_info.get('gapfill_candidate_costs', {}))  # Optimize the target in the merged network and add the active synthesis route to the model.
                need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-6]

                new_way_rxn = list(need_fluxes.keys())
        print(
            f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}: Pathway for {demand_id} in the merged model was computed: {new_way_rxn}"
        )
        return new_way_rxn

    def biomass_gapfilling(self, model, model_info, initial_rxn_id, general_library):
        """"""
        big_model = self.get_big_model(model_info, model, general_library)
        print(
            f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}: Pathway for {initial_rxn_id} in the merged model"
        )
        new_way_rxn = []
        with big_model:
            print(initial_rxn_id, "....................")
            big_model.objective = initial_rxn_id
            big_model.reactions.get_by_id(initial_rxn_id).bounds = (
                0,
                0.1,
            )  # Set an upper bound on the demand reaction.

            print(big_model.slim_optimize(), "..............")
            if big_model.slim_optimize() > 0:
                pfba_solution = rule_pfba(big_model, model_info.get('gapfill_candidate_costs', {}))  # Optimize the target in the merged network and add the active synthesis route to the model.
                need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-6]
                print("........................")

                new_way_rxn = list(need_fluxes.keys())
        print(
            f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}: Pathway for {initial_rxn_id} in the merged model was computed"
        )
        return new_way_rxn

    def add_met_name(self, model, need_add_rxn, general_library, obj):
        """"""
        print(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            ": Annotating metabolites in gap-filling reactions",
        )
        seed_met, meta_met, virtual_met, kegg_met, other_met = [], [], [], [], []
        pattern = r"^C\d{5}$"
        for rxnId in need_add_rxn:

            rxns = general_library.reactions.get_by_id(rxnId)
            for met in rxns.metabolites:
                if met.id.startswith("cpd"):
                    seed_met.append(met.id)
                elif "@" in met.id:
                    meta_met.append(met.id)
                elif "[" in met.id and "]" in met.id:
                    virtual_met.append(met.id)
                elif re.match(pattern, met.id):
                    kegg_met.append(met.id)
                else:
                    other_met.append(met.id)
        if len(seed_met) != 0:
            obj.add_seed_formula(seed_met)
        if len(meta_met) != 0:
            obj.add_meta_formula(meta_met)

        if len(kegg_met) != 0:
            obj.add_kegg_formula(kegg_met)
        if len(other_met) != 0:
            obj.add_other_formula(other_met)

    def add_gapfilling_rxn(self, model, model_info, check_model, general_library, new_way_rxn, obj):
        """"""
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':开始添加gapfilling反应:')
        need_add_rxn = []
        if not new_way_rxn:
            gap_info = ["Can not able to gapfill !!!"]
        else:    
            t1 = time.time()
            gap_info = ["The gapfilling process added the reaction:"]
            for rxnId in new_way_rxn:
                if rxnId not in model.reactions and rxnId in general_library.reactions:
                    rxns = candidate_copy(general_library, rxnId, model_info)
                    model.add_reactions([rxns])
                    need_add_rxn.append(rxnId)
                    model_info["gap_num"] += 1
                    # rxn_exp = [f"{rxns}"]
                    # if model_info['model_identifier'] == 'modelseed' or model_info['model_identifier'] == 'metaNetx':
                    #     rxn_exp = [f"{rxns.id}: {rxns.build_reaction_string(use_metabolite_names=True)}"]
                    # gap_info.extend(rxn_exp)
            print('need_add_rxn: ',need_add_rxn)
            t2 = time.time()
            print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':正在gapfilling添加反应', 'time:',t2-t1)
            self.add_met_name(model, need_add_rxn, general_library, obj)
            for rxnId in need_add_rxn:
                # if rxnId not in model.reactions and rxnId in general_library.reactions:
                rxns = model.reactions.get_by_id(rxnId)
                rxn_exp = [f"{rxns.id}: {rxns.build_reaction_string(use_metabolite_names=True)}"]
                gap_info.extend(rxn_exp)
            if len(need_add_rxn) == 0:
                gap_info.extend(["no gap response,since it already exists!"])
        model_info["gaped_rxn"].extend(need_add_rxn)
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':添加gapfilling反应完毕')
        return need_add_rxn, gap_info

    def add_gap_rxn(self, model, need_add_rxn, general_library, check_model):
        """"""
        for rxnId in need_add_rxn:
            rxns = general_library.reactions.get_by_id(rxnId).copy()
            if rxnId not in model.reactions:
                if rule_policy_enabled(self.cfg):
                    from mqc.control.gap_rule_policy import hydrolysis_bounds
                    rxns.bounds = hydrolysis_bounds(rxns)
                    from mqc.control.biotin_direction import biotin_bounds
                    rxns.bounds = biotin_bounds(rxns)
                    from mqc.control.deacetylase_direction import deacetylase_bounds
                    rxns.bounds = deacetylase_bounds(rxns)
                model.add_reactions([rxns])
                check_model.add_reactions([rxns.copy()])

    def check_bio_component_is_zero(self, model_info, model, model_control_info, initial_rxn_id):
        """"""
        t1 = time.time()
        no_synthesized_mets, synthesized_mets, num = [], [], 0
        ATP_NADH = [ATP_NAME+CTP_NAME+GTP_NAME+UTP_NAME+ITP_NAME+NADH_NAME+NADPH+FADH2+Q8H2_NAME+FMNH2+MQL8+DMMQL8]
        for met in model.reactions.get_by_id(initial_rxn_id).reactants:
            num += 1
            if met.name in ATP_NADH:
                pass 
            with model:
                if met.name in ACCOA and met.compartment in C_COMPARTMENT:
                    objectiveId = add_accoa(model, str(met.id))
                else:
                    objectiveId = add_demand(model_info, model, str(met.id), independent=True)
                model.objective = objectiveId
                ini_flux = model.slim_optimize()
                synthesized_mets.extend([f"{str(met.id)} : {round(ini_flux,3)}"])
                # ini_flux = model.optimize().fluxes[objectiveId]
                if pd.isna(ini_flux) or ini_flux <= 1e-6:
                    no_synthesized_mets.append(str(met.id))
                print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':检查生物量组成物质哪些为0', objectiveId,ini_flux,model.reactions.get_by_id(objectiveId).bounds)
        temp = [f"The biomass equation consists of {num} total components"]
        if no_synthesized_mets:
            temp.extend([f"among which the synthesis rate of component {','.join(no_synthesized_mets)} was 0 "])
        else:
            temp.extend([f"none of these precursor species have a rate of 0, namely: {','.join(synthesized_mets)}"])
        # model_control_info["check_biomass_production"]["Check_synthesizability_of_biomass_components"] = ', '.join(temp)
        write_biomass(model_control_info, temp, "Check_synthesizability_of_biomass_components")
        t2 = time.time()
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':检查生物量组成物质哪些为0', t2-t1,'s')
        return no_synthesized_mets

    def gap_initial_zero_bio(
        self,
        model_info,
        model,
        model_control_info,
        general_library,
        no_synthesized_mets,
        big_model,
        obj,
        check_model,
    ):
        """"""
        temp = []
        new_way_rxn = self.get_new_way_rxn(
            model_info, big_model, no_synthesized_mets, model_control_info
        )
        need_add_rxn, gap_info = self.add_gapfilling_rxn(
            model, model_info, check_model, general_library, new_way_rxn, obj
        )
        model_info["bio_modified"].extend(need_add_rxn)
        temp.extend(gap_info)

        if getattr(self.cfg, 'native_reaction_coverage', False):
            from pathlib import Path
            from mqc.control.native_coverage import select_native_coverage
            reports = model_info.setdefault('native_reaction_coverage', [])
            directory = Path(self.cfg.output_dir) / ('native_coverage_%02d' % len(reports))
            # Actual physical network after ordinary precursor insertion;
            # original scaled search is still available in memory here.
            cobra.io.save_json_model(model, str(Path(self.cfg.output_dir) / 'coverage_entry_model.json'))
            original_bounds = {r.id: r.bounds for r in model.reactions}
            extra, report = select_native_coverage(model, big_model, model_info, directory)
            reports.append(report)
            if extra:
                added, messages = self.add_gapfilling_rxn(
                    model, model_info, check_model, general_library, extra, obj)
                if added != extra:
                    raise RuntimeError('Native coverage insertion differed from selection')
                for rid in added:
                    actual = model.reactions.get_by_id(rid)
                    expected = big_model.reactions.get_by_id(rid)
                    if actual.bounds != expected.bounds or actual.reaction != expected.reaction:
                        raise RuntimeError('Native coverage insertion changed candidate: ' + rid)
                model_info['bio_modified'].extend(added)
                report['added'] = added
                temp.extend(['Native reaction coverage within precursor repair:'] + messages)
            if any(model.reactions.get_by_id(r).bounds != bounds for r, bounds in original_bounds.items()):
                raise RuntimeError('Native coverage changed original bounds')
            cobra.io.save_json_model(model, str(Path(self.cfg.output_dir) / 'coverage_exit_model.json'))
            (directory / 'report.json').write_text(json.dumps(report, indent=2))

        write_biomass(model_control_info, temp, "Gapfilling_of_biomass_components")
        return no_synthesized_mets

    def get_no_synthesized_mets_flux(self, model, model_info, no_synthesized_mets):
        """"""
        no_synthesized_mets_dic = {}
        for metId in no_synthesized_mets:
            with model:
                objectiveId = add_demand(model_info, model, metId, independent=True)
                model.objective = objectiveId
                now_flux = model.slim_optimize()   
                no_synthesized_mets_dic[metId] = round(now_flux,3)
        return no_synthesized_mets_dic

    def reducing_carbon_source_input(self, model_info, model, check_model):
        """"""
        model_info["no_gap_met_input"] = []
        for rxn in check_model.reactions:
            if len(rxn.metabolites) == 1:
                r_met = list(rxn.metabolites.keys())[0]
                if r_met.id in model_info["no_gap_metid"]:

                    print(
                        "reducing_carbon_source_input111",
                        model.reactions.get_by_id(rxn.id).bounds,
                        rxn.bounds,
                        rxn.id,
                        r_met.id,
                        rxn,
                    )

                    model.reactions.get_by_id(rxn.id).bounds = (-1000, 1000)
                    model_info["no_gap_met_input"].append(r_met.id)
                    model_info["no_gap_rxn_input"].append(rxn.id)
                for metId in model_info["no_gap_metid"]:
                    mets1 = model.metabolites.get_by_id(r_met.id)
                    mets2 = model.metabolites.get_by_id(metId)

                    if mets1.name == mets2.name:
                        print(
                            "reducing_carbon_source_input222",
                            model.reactions.get_by_id(rxn.id).bounds,
                            rxn.bounds,
                            rxn.id,
                            r_met.id,
                            rxn,
                        )
                        model.reactions.get_by_id(rxn.id).bounds = (-1000, 1000)
                        model_info["no_gap_met_input"].append(r_met.id)
                        model_info["no_gap_rxn_input"].append(rxn.id)
        for rxn in model_info["reactions"]:
            if rxn["id"] in model_info["no_gap_rxn_input"]:
                rxn["bounds"] = (-1000, 1000)
        print("biomass", model.slim_optimize())

    def update_no_gap_metid(self, model_info, model, no_synthesized_mets, flag):
        """"""
        no_gap = []
        for metId in no_synthesized_mets:
            with model:
                if flag == 0:
                    objectiveId = add_demand(model_info, model, metId, independent=True)
                else:
                    mets = model.metabolites.get_by_id(metId)
                    try:
                        demand_rxn = model.add_boundary(mets, type = 'demand') # 添加demand反应
                        objectiveId = demand_rxn.id
                    except:
                        objectiveId = f"DM_{metId}"
                model.objective = objectiveId
                ini_flux = model.slim_optimize()     
                # print(metId,objectiveId,ini_flux,'ccccccccccccc') 
                if pd.isna(ini_flux) or ini_flux <= 1e-6:
                    no_gap.append(metId)
        model_info['no_gap_metid'] = no_gap
        return no_gap

    def find_min_carbon_source_set(
        self, model_info, model, check_model, model_control_info
    ):
        """"""
        carbon_source_rxn, flag = [], 0
        for rxn in model_info["reactions"]:
            if "carbon_source" in rxn.keys():
                carbon_source_rxn.append(rxn["id"])
        for rxnId in model_info["sink_rxn"]:
            rxns = model.reactions.get_by_id(rxnId)
            mets = list(rxns.metabolites.keys())[0]
            if len(rxns.metabolites) == 1 and "C" in mets.formula:
                carbon_source_rxn.append(rxnId)
                print(rxnId, "find_min_carbon_source_set")

        min_set_dict = get_min_carbon_source_set(
            model_info, model, check_model, carbon_source_rxn, flag
        )
        print("min_set_dict", min_set_dict, "sssssssssssssssss")
        if min_set_dict:
            for k, v in min_set_dict.items():
                min_set_key = k
                min_set_all = v
            glucose_rxnId = model_info["glucose_rxnId"]
            if min_set_all:
                if (
                    "min_set" == min_set_key
                ):  # Prefer glucose or sucrose among single sources that support growth.
                    if glucose_rxnId in min_set_all:
                        min_set = glucose_rxnId
                    else:
                        carbon_rxnId = find_other_carbon_source(model_info, model)
                        min_set = (
                            carbon_rxnId
                            if carbon_rxnId in min_set_all
                            else min_set_all[0]
                        )
                elif (
                    "min_set1" == min_set_key
                ):  # Otherwise select a pair of sources that supports growth.
                    min_set = min_set_all[0]
                else:  # Retain all members of the source set found by recursive removal.
                    min_set = min_set_all
            if type(min_set) == str:
                model_info["min_set"] = [min_set]
            else:
                model_info["min_set"] = min_set
        else:
            # An unsuccessful search is not a request to erase the current medium.
            return
        print(model_info["min_set"], 'model_info["min_set"]')
        model_info["automatic_build"] = 0
        get_carbon_proportion(model_info, model, check_model, model_control_info)
        set_c_source_supply(model_info, model, "bio", check_model)
        model_control_info["boundary_information"]["carbon_source_boundary"] = ",".join(
            model_info["carbon_source_boundary"]
        )

    def add_no_gap_sink_rxn(self, model_info, model, check_model):
        """"""
        sink_rxn = []
        initial_rxn_id = model_info["initial_rxn"]["biomass_rxn_id"]
        initial_rxn = model.reactions.get_by_id(initial_rxn_id)
        if len(model_info["no_gap_metid"]) <= (len(initial_rxn.reactants) * 0.2):
            print("...........................")
            print(
                model_info["no_gap_metid"],
                len(initial_rxn.reactants),
                len(initial_rxn.reactants) * 0.2,
            )
            print("...........................")
            for metId in model_info["no_gap_metid"]:
                mets = model.metabolites.get_by_id(metId)
                model.add_boundary(mets, type="sink", lb=-1000, ub=0)
                check_model.add_boundary(mets, type="sink", lb=-1000, ub=0)
                sink_rxn.append(f"SK_{metId}")
        if "sink_rxn" in model_info.keys():
            model_info["sink_rxn"].extend(sink_rxn)
        else:
            model_info["sink_rxn"] = sink_rxn

    def retry_joint_biomass(self, model_info, model, model_control_info,
                           general_library, obj, check_model):
        """Try the existing biomass search once before changing nutrient supply.

        Equal precursor demands do not establish feasibility of the complete
        biomass equation: simultaneous ratios, cofactors and byproducts can
        still block growth. Keep the first precursor search, native S*1000,
        current medium and rule costs. Validate any extra route in an unscaled
        copy before inserting it. Final MQC cycle/yield guards still apply.
        """
        biomass = model_info['initial_rxn']['biomass_rxn_id']
        model.objective = biomass
        before = model.optimize()
        if before.status != 'optimal' or before.objective_value > 1e-6:
            return False
        boundaries = {r.id: r.bounds for r in model.reactions if len(r.metabolites) == 1}
        report = dict(status='not_restored', growth_before=before.objective_value,
                      target=biomass, experimental_values_used=False,
                      supply_bounds_unchanged=True, attempted_candidates=[])
        model_info.setdefault('joint_biomass_retries', []).append(report)
        # The existing builder updates candidate metadata. Keep a failed
        # attempt from replacing the first search's accepted direction ledger.
        retry_info = deepcopy(model_info)
        selected = self.biomass_gapfilling(model, retry_info, biomass, general_library)
        extra = list(dict.fromkeys(rid for rid in selected
            if rid not in model.reactions and rid in general_library.reactions
            and len(general_library.reactions.get_by_id(rid).metabolites) > 1))
        report['attempted_candidates'] = extra
        if not extra:
            report['reason'] = 'no_additional_route_found'
            return False
        trial = model.copy()
        trial.add_reactions([candidate_copy(general_library, rid, retry_info) for rid in extra])
        result = trial.optimize()
        report.update(trial_solver_status=result.status,
                      trial_growth=result.objective_value if result.status == 'optimal' else None)
        if result.status != 'optimal' or not np.isfinite(result.objective_value) or result.objective_value <= 1e-6:
            report['reason'] = 'extra_route_failed_native_biomass_check'
            return False
        if boundaries != {r.id:r.bounds for r in trial.reactions if len(r.metabolites)==1}:
            raise RuntimeError('Joint biomass retry changed nutrient supply')
        bounds = model_info.setdefault('gapfill_candidate_bounds', {})
        costs = model_info.setdefault('gapfill_candidate_costs', {})
        for rid in extra:
            bounds[rid] = list(trial.reactions.get_by_id(rid).bounds)
            cost = retry_info.get('gapfill_candidate_costs', {}).get(rid, 1.)
            if cost == 1:
                costs.pop(rid, None)
            else:
                costs[rid] = cost
        added, messages = self.add_gapfilling_rxn(
            model, model_info, check_model, general_library, extra, obj)
        model_info['bio_modified'].extend(added)
        report.update(status='restored_in_recorded_medium', added=added,
                      reason='existing_whole_biomass_search_restored_joint_feasibility')
        write_biomass(model_control_info,
            ['Joint biomass feasibility was restored before changing nutrient supply.'] + messages,
            'Gapfilling_of_biomass_components')
        return True

    def preprocess_initial_zero_bio(
        self,
        model_info,
        model,
        model_control_info,
        general_library,
        no_synthesized_mets,
        obj,
        check_model,
    ):
        """"""
        print("bio", model.slim_optimize())
        if model.slim_optimize() <= 1e-6:
            model_info["bio_nogrow"] = "0"
            big_model = self.get_big_model(model_info, model, general_library)
            model_control_info["check_biomass_production"]["score"] = 0
            no_synthesized_mets = self.gap_initial_zero_bio(
                model_info,
                model,
                model_control_info,
                general_library,
                no_synthesized_mets,
                big_model,
                obj,
                check_model,
            )
            model.objective = model_info["initial_rxn"]["biomass_rxn_id"]
            print("bio22222", model.slim_optimize())
            if model.slim_optimize() <= 1e-6:
                self.retry_joint_biomass(model_info, model, model_control_info,
                                         general_library, obj, check_model)
            if model.slim_optimize() <= 1e-6:
                no_gap = self.update_no_gap_metid(
                    model_info, model, no_synthesized_mets, 0
                )
                self.reducing_carbon_source_input(
                    model_info, model, check_model
                )  # Open exchanges for the supplied metabolites.
                if model.slim_optimize() <= 1e-6:
                    new_no_gap = self.update_no_gap_metid(
                        model_info, model, no_synthesized_mets, 1
                    )
                    self.find_min_carbon_source_set(
                        model_info, model, check_model, model_control_info
                    )  # Open the selected carbon sources.
                    print(model.slim_optimize(), "sssssssssssssss")
                    if model.slim_optimize() <= 1e-6:
                        self.add_no_gap_sink_rxn(
                            model_info, model, check_model
                        )  # Add sink reactions for precursors that cannot be synthesized.
                        print(
                            model.slim_optimize(),
                            "aaaaaaaaaaaaaaa",
                            model_info["no_gap_metid"],
                            new_no_gap,
                            "dddddddd",
                            no_gap,
                        )
                        if model.slim_optimize() <= 1e-6:
                            temp = [
                                f"After gap, the model still cannot simulate growth, the biomass components {','.join(model_info['no_gap_metid'])} cannot be synthesized, please check the model boundary conditions and monomer synthesis pathway"
                            ]
                            write_biomass(
                                model_control_info,
                                temp,
                                "Gapfilling_of_biomass_components",
                            )
                            return -1
                        else:
                            temp = [
                                f"After gap, the model still cannot simulate growth, the biomass components {','.join(model_info['no_gap_metid'])} cannot be synthesized"
                            ]
                            temp.extend(
                                [
                                    f"But in the initial model, ({','.join(model_info['no_gap_met_input'])}) can be directly used as a carbon source input, so supply ({','.join(model_info['no_gap_met_input'])}) as carbon source,\n the model still cannot simulate growth, then open the minimum carbon source set ({','.join(model_info['min_set'])}), the model still cannot simulate growth;\n so add the sink reaction of the non synthesized components ({','.join(model_info['no_gap_metid'])}) yields a biomass value of {round(model.slim_optimize(),3)}"
                                ]
                            )
                            write_biomass(
                                model_control_info,
                                temp,
                                "Gapfilling_of_biomass_components",
                            )
                            model_info["initial_rxn"][
                                "biomass_rxn_flux"
                            ] = model.slim_optimize()
                    else:
                        temp = [
                            f"After gap, the model still cannot simulate growth, the biomass components {','.join(model_info['no_gap_metid'])} cannot be synthesized"
                        ]
                        temp.extend(
                            [
                                f"But in the initial model, ({','.join(model_info['no_gap_met_input'])}) can be directly used as a carbon source input, so supply ({','.join(model_info['no_gap_met_input'])}) as carbon source, the model still cannot simulate growth,\n then open the minimum carbon source set ({','.join(model_info['min_set'])}) yields a biomass value of {round(model.slim_optimize(),3)}"
                            ]
                        )
                        write_biomass(
                            model_control_info, temp, "Gapfilling_of_biomass_components"
                        )
                        model_info["initial_rxn"][
                            "biomass_rxn_flux"
                        ] = model.slim_optimize()
                else:
                    temp = [
                        f"After gap, the model still cannot simulate growth, the biomass components {','.join(no_gap)} cannot be synthesized"
                    ]
                    temp.extend(
                        [
                            f"But in the initial model, ({','.join(model_info['no_gap_met_input'])}) can be directly used as a carbon source input, so supply ({','.join(model_info['no_gap_met_input'])}) as carbon source yields a biomass value of {round(model.slim_optimize(),3)}"
                        ]
                    )
                    write_biomass(
                        model_control_info, temp, "Gapfilling_of_biomass_components"
                    )
                    model_info["initial_rxn"][
                        "biomass_rxn_flux"
                    ] = model.slim_optimize()
            else:
                no_synthesized_mets_dic = self.get_no_synthesized_mets_flux(
                    model, model_info, no_synthesized_mets
                )
                temp = [
                    f"After gapfilling, {','.join(no_synthesized_mets)} can be synthesized as follows: {no_synthesized_mets_dic}"
                ]
                temp.extend(
                    [
                        f"The growth rate of this model is {round(model.slim_optimize(),3)} after adding the above reaction"
                    ]
                )
                write_biomass(
                    model_control_info, temp, "Gapfilling_of_biomass_components"
                )
                model_info["initial_rxn"]["biomass_rxn_flux"] = model.slim_optimize()
        else:
            temp = ["All of the biomass components can be synthesized."]
            write_biomass(model_control_info, temp, "Gapfilling_of_biomass_components")
        return 0

    def get_adp_h_pi(self, model):
        """"""
        adp_identifier, h_identifier, pi_identifier = "", "", ""
        for met in model.metabolites:
            if met.name in ADP_NAME and met.compartment in C_COMPARTMENT:
                adp_identifier = met.id
            if met.name in H_NAME and met.compartment in C_COMPARTMENT:
                h_identifier = met.id
            if met.name in PI_NAME and met.compartment in C_COMPARTMENT:
                pi_identifier = met.id
        return adp_identifier, h_identifier, pi_identifier

    def is_missing_sustaining_energy(
        self, model_info, model, model_control_info, initial_rxn_id
    ):
        """"""
        initial_rxns = model.reactions.get_by_id(initial_rxn_id)

        _, bio_text = keep_decimal_places(initial_rxns)
        missing_identifier, sustain, is_missing = [], [], 0
        products_mets = [
            met.id for met in model.reactions.get_by_id(initial_rxn_id).products
        ]
        adp_identifier, h_identifier, pi_identifier = self.get_adp_h_pi(model)
        if len(set(ADP + H + PI) & set(products_mets)) == 0:
            sustain = [
                "Growth-related maintenance (h, adp, pi) in biomass is missing, please make up"
            ]
            is_missing += 1
        elif (
            len(set(ADP) & set(products_mets)) == 0
            or len(set(H) & set(products_mets)) == 0
            or len(set(PI) & set(products_mets)) == 0
        ):
            if len(set(ADP) & set(products_mets)) == 0:
                if len(set(H) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(h_identifier)
                if len(set(PI) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(pi_identifier)
                initial_rxns.reaction = bio_text + f" + {coefficient} {adp_identifier}"
                missing_identifier.append(adp_identifier)
            elif len(set(H) & set(products_mets)) == 0:
                if len(set(ADP) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(adp_identifier)
                if len(set(PI) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(pi_identifier)
                missing_identifier.append(h_identifier)
                initial_rxns.reaction = bio_text + f" + {coefficient} {h_identifier}"
            elif len(set(PI) & set(products_mets)) == 0:
                if len(set(H) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(h_identifier)
                if len(set(ADP) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(adp_identifier)
                missing_identifier.append(pi_identifier)
                initial_rxns.reaction = bio_text + f" + {coefficient} {pi_identifier}"
            sustain = [
                f"The growth-related maintenance energy in biomass is incorrect and missing, and {','.join(missing_identifier)} is added to the product side of biomass equation."
            ]
            is_missing += 1
        if len(sustain) != 0:
            write_biomass(
                model_control_info,
                sustain,
                "Check_synthesizability_of_biomass_components",
            )
        initial_bio_value = model.slim_optimize()

        _, initial_biomass_rxn = keep_decimal_places(initial_rxns)
        if is_missing != 0:
            model_control_info["check_biomass_production"]["score"] = 0
            ini_bio_value = model_info["initial_rxn"]["biomass_rxn_flux"]
            sustain_result = [
                f"Initial biomass growth rate (before improving Growth-Related Maintenance Energy): {round(ini_bio_value,5)}"
            ]
            sustain_result.extend(
                [
                    f"Initial biomass equation (before improving Growth-Related Maintenance Energy): {bio_text}"
                ]
            )
            sustain_result.extend(
                [
                    f"Initial biomass growth rate (after improving Growth-Related Maintenance Energy): {round(initial_bio_value,5)}"
                ]
            )
            sustain_result.extend(
                [
                    f"Initial biomass equation (after improving Growth-Related Maintenance Energy): {initial_biomass_rxn}"
                ]
            )
        else:
            sustain_result = [
                f"Initial biomass growth rate (growth related maintenance can be intact): {round(initial_bio_value,5)}"
            ]
            sustain_result.extend(
                [
                    f"Initial biomass equation (growth related maintenance energy intact): {initial_biomass_rxn}"
                ]
            )
        return sustain_result

    def get_macromolecule(self, model, initial_rxn_id, model_control_info):
        """"""
        macromolecule, special_macro = [], [
            "dna",
            "rna",
            "protein",
            "lipid",
            "peptido",
            "cellwall",
        ]
        for metId in model.reactions.get_by_id(initial_rxn_id).reactants:
            met = model.metabolites.get_by_id(str(metId))
            if (
                pd.isna(met.formula)
                or "R" in met.formula
                or "X" in met.formula
                or met.formula == "nan"
                or not met.formula
                or met.formula == "null"
            ):
                metName = met.name
                if any(
                    metName.lower() in macro for macro in MACROMOLECULE
                ):  # Identify macromolecules by their reference names.
                    macromolecule.append(str(metId))
                if any(
                    macro in metName.lower() or macro in str(metId).lower()
                    for macro in special_macro
                ):
                    macromolecule.append(str(metId))
        macromolecule = list(set(macromolecule))
        temp = [f"Macromolecules include the following components: {macromolecule}"]
        write_biomass(model_control_info, temp, "Check_macromolecules")
        return macromolecule

    def get_smallmolecule(
        self, model_info, model, initial_rxn_id, macromolecule, model_control_info
    ):
        """"""
        smallmolecule = []
        initial_rxn = model.reactions.get_by_id(initial_rxn_id)

        for met in initial_rxn.metabolites:
            if (
                str(met.id) not in macromolecule
                and met not in initial_rxn.products
                and met.name not in SMALLMOLECULE
            ):
                smallmolecule.append(str(met.id))
        temp = [f"Small molecule metabolites include : {smallmolecule}"]
        write_biomass(model_control_info, temp, "Check_small_molecules")
        return smallmolecule

    def set_atpm_bounds(self, model_info, model):
        """Temporarily relax maintenance without losing pre-relaxation bounds."""
        restore_bounds = model_info.setdefault('atpm_restore_bounds', {})
        reaction_info = {rxn['id']: rxn for rxn in model_info['reactions']}
        for atpm_id in ATPM:
            if atpm_id in model_info['all_rxn_id']:
                reaction = model.reactions.get_by_id(atpm_id)
                restore_bounds.setdefault(atpm_id, reaction.bounds)
                reaction.bounds = (0,1000)
                # Do not clear the flag of a previously visited ATPM alias.
                reaction_info[atpm_id]['atpms_modify'] = "true"

    def get_formula(self, model, ids, modelInfo, check_model, obj):
        """"""
        readict, prodict, resultdict, flag, delect_rxn, need_add_rxn = (
            {},
            {},
            {},
            0,
            "",
            [],
        )
        mets = model.metabolites.get_by_id(ids)
        if (
            pd.isna(mets.formula)
            or "R" in mets.formula
            or "X" in mets.formula
            or mets.formula == "nan"
            or not mets.formula
            or mets.formula == "null"
        ):
            with model:

                if f"SK_{ids}" in model.reactions:
                    model.reactions.get_by_id(f"SK_{ids}").bounds = (0, 0)
                demand_rxn = add_demand(modelInfo, model, ids)
                model.objective = demand_rxn
                print(mets.id, mets.formula, model.slim_optimize(), ".................")
                if model.slim_optimize() > 0:
                    pfba_solution = pfba(
                        model
                    )  # Retain reactions with an absolute flux above 1e-6.
                    need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-6]
                else:
                    need_fluxes = {}

                for r, v in need_fluxes.items():
                    if r == demand_rxn:
                        continue
                    rxns = model.reactions.get_by_id(r)
                    if len(rxns.metabolites) == 1:
                        all_met = [m for m in rxns.metabolites][0]
                        products_mets = [m.id for m in rxns.products]
                        reactants_mets = [m.id for m in rxns.reactants]
                        if (
                            pd.isna(all_met.formula)
                            or "R" in all_met.formula
                            or "X" in all_met.formula
                            or all_met.formula == "nan"
                            or not all_met.formula
                            or all_met.formula == "null"
                        ):
                            flag += 1
                        if flag > 1:
                            return 0

                        if (reactants_mets and v < 0) or (products_mets and v > 0):
                            for key, value in all_met.elements.items():
                                if (
                                    key in readict
                                ):  # Merge the coefficient dictionaries and sum shared entries.
                                    readict[key] += value
                                else:
                                    readict[key] = value
                        if (reactants_mets and v > 0) or (products_mets and v < 0):
                            print(
                                reactants_mets, products_mets, v, all_met.formula, "xxx"
                            )
                            for key, value in all_met.elements.items():
                                if key in prodict:
                                    prodict[key] += value
                                else:
                                    prodict[key] = value
            print(readict, prodict)

            for key in OrderedDict.fromkeys(
                list(readict.keys()) + list(prodict.keys())
            ):  # Preserve key order when combining coefficient dictionaries.
                resultdict[key] = readict.get(key, 0) - prodict.get(key, 0)
            print(resultdict)
            resultdict = {
                key: value for key, value in resultdict.items() if value != 0
            }  # Remove zero coefficients.
            print(resultdict)
            result_string = "".join(
                f"{key}" if value == 1 else f"{key}{value}"
                for key, value in resultdict.items()
            )
            print(result_string, "xxx")
            print(mets.id, mets.formula, model.slim_optimize(), ".................")
            if not result_string:
                return 0
            print(result_string)
            mets.formula = result_string
            print(mets.id, mets.formula, model.slim_optimize(), ".................")
            print(mets.id, mets.formula, result_string)

    def calculateResult(
        self, model, initial_rxn_id, model_info, sink_rxn, check_model, obj
    ):
        """"""
        result, final = 0, 0

        if model.slim_optimize() > 0:
            pfba_solution = pfba(model)
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-6]
        else:
            need_fluxes = {}

        for ids, v in need_fluxes.items():
            rxns = model.reactions.get_by_id(ids)
            all_mets = [met.id for met in rxns.metabolites]
            all_mets_name = [met.name for met in rxns.metabolites]
            if (
                len(all_mets) == 1
                and ids != initial_rxn_id
                and "biomass" not in all_mets_name[0].lower()
                and ids not in sink_rxn
                and ids not in model_info["no_gap_rxn_input"]
            ):
                mets = model.metabolites.get_by_id(all_mets[0])
                if (
                    pd.isna(mets.formula)
                    or "R" in mets.formula
                    or "X" in mets.formula
                    or mets.formula == "nan"
                    or not mets.formula
                    or mets.formula == "null"
                ):
                    met_formula = self.get_formula(
                        model, all_mets[0], model_info, check_model, obj
                    )
                    print(
                        model_info["no_gap_rxn_input"],
                        ids,
                        met_formula,
                        "ddddddddddddddddddddd",
                    )
                    if met_formula == 0:
                        return all_mets[0]
                result -= v * rxns.get_coefficient(all_mets[0]) * relative_molecular_mass(model, all_mets[0])

        try:
            final = abs(result / model.slim_optimize())
        except ZeroDivisionError:
            final = 0
        return final

    def get_initial_biomass(
        self, model_info, model, initial_rxn_id, model_control_info, check_model, obj
    ):
        """"""
        self.set_atpm_bounds(model_info, model)
        initial_biomass = self.calculateResult(
            model, initial_rxn_id, model_info, [], check_model, obj
        )
        model_info["biomass"]["initial_biomass"] = initial_biomass
        if type(initial_biomass) != str and (990 < initial_biomass < 1010):
            temp = [
                f"The biomass relative molecular mass is {round(initial_biomass,3)}mg, corrected"
            ]
            write_biomass(model_control_info, temp, "Gapfilling_of_biomass_components")

        else:
            model_info["bio_1g"] = "0"

    def smallmolecule_iscoupling(
        self,
        model,
        model_info,
        small_demand_id,
        initial_rxn_id,
        metId,
        general_library,
        model_control_info,
        obj,
        check_model,
    ):
        """"""
        iscoupling, close_coupling_rxn, need_add_rxn, temp_list = 0, [], [], []
        sink_rxn = model_info["sink_rxn"]
        initial_flux = model.slim_optimize()
        print(initial_flux, type(initial_flux))
        if pd.isna(initial_flux) or initial_flux <= 1e-6:
            temp_list = [f"{metId} rate is 0, so gap:"]
            new_way_rxn = self.gapfilling(
                model, model_info, small_demand_id, general_library
            )
            need_add_rxn, gap_info = self.add_gapfilling_rxn(
                model, model_info, check_model, general_library, new_way_rxn, obj
            )
            model_info["bio_modified"].extend(need_add_rxn)
            temp_list.extend(gap_info)
            write_biomass(model_control_info, temp_list, "Check_small_molecules")
            if len(need_add_rxn) == 0:
                temp_list = [
                    f"In biomass generation, {metId} must be coupled with macromolecules, the gapfilling cannot be performed, and there is an error! ! ! Please correct and try again"
                ]
                write_biomass(model_control_info, temp_list, "Check_small_molecules")
                return -1
            elif pd.isna(model.slim_optimize()):
                temp_list = [
                    f"The compositing rate is nan after {metId} gapfilling, there is an error! ! ! Please correct and try again"
                ]
                write_biomass(model_control_info, temp_list, "Check_small_molecules")
                return -1
            else:
                final_flux = model.slim_optimize()
                temp_list = [
                    f"The composite rate after {metId} gapfilling is: {round(final_flux,5)}"
                ]
                write_biomass(model_control_info, temp_list, "Check_small_molecules")
                return need_add_rxn, close_coupling_rxn
        while True:
            status = 0  # Exit when no biomass or macromolecule coupling is detected in this iteration.
            if model.slim_optimize() > 0:
                pfba_solution = pfba(model)
                need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-6]

            else:
                need_fluxes = {}
            for r in list(
                need_fluxes.keys()
            ):  # Resolve macromolecule synthesis that depends on the biomass reaction.
                rxns = model.reactions.get_by_id(r)
                all_mets = [m.id for m in rxns.metabolites]
                met_formula = list(rxns.metabolites.keys())[0].formula
                # Unknown material entering the model can supply a precursor;
                # an outward biomass drain cannot. Keep the native biomass
                # coupling check, but never mistake secretion for supply.
                if r == initial_rxn_id or unresolved_boundary_source(rxns, need_fluxes[r]):
                    if (
                        r != small_demand_id and r not in sink_rxn
                    ):  # Check for coupling through biomass or macromolecular exchange reactions.
                        print(r, sink_rxn, "xxxx")
                        met_formula = self.get_formula(
                            model, all_mets[0], model_info, check_model, obj
                        )
                        if met_formula != 0:
                            continue
                        status += 1
                        iscoupling += 1
                        temp_list = [
                            f"{metId} synthesis rate is: {round(model.slim_optimize(),5)}"
                        ]
                        model.reactions.get_by_id(r).bounds = (0, 0)
                        close_coupling_rxn.append(r)
                        temp_list.extend(
                            [
                                f"{small_demand_id} synthesis is coupled with {r}, close {r}"
                            ]
                        )
                        write_biomass(
                            model_control_info, temp_list, "Check_small_molecules"
                        )
            if status == 0:
                break
        if iscoupling == 0 and initial_flux > 1e-6:
            return 0
        else:
            final_flux2 = model.slim_optimize()
            if final_flux2 <= 1e-6 or pd.isna(final_flux2):
                temp_list = [
                    f"The macromolecules coupled with the small molecule {metId} are all closed, the generation rate is 0, and the gapfilling is performed:"
                ]
                new_way_rxn = self.gapfilling(
                    model, model_info, small_demand_id, general_library
                )
                need_add_rxn, gap_info = self.add_gapfilling_rxn(
                    model, model_info, check_model, general_library, new_way_rxn, obj
                )
                model_info["bio_modified"].extend(need_add_rxn)
                temp_list.extend(gap_info)
                if len(need_add_rxn) == 0:
                    temp_list.extend(
                        [
                            f"In biomass generation, {metId} must be coupled with macromolecules, the gapfilling cannot be performed, and there is an error! ! ! Please correct and try again"
                        ]
                    )
                    write_biomass(
                        model_control_info, temp_list, "Check_small_molecules"
                    )
                    return -1
                else:
                    final_flux = model.slim_optimize()
                    temp_list.extend(
                        [
                            f"The synthesis rate after gapfilling is: {round(final_flux,5)}"
                        ]
                    )
                    write_biomass(
                        model_control_info, temp_list, "Check_small_molecules"
                    )
                    return need_add_rxn, close_coupling_rxn
            else:
                temp_list = [
                    f"The macromolecules coupled with the small molecule {metId} are all closed, and the synthesis rate is {round(final_flux2,5)}, if it is not 0, no gap will be performed"
                ]
                write_biomass(model_control_info, temp_list, "Check_small_molecules")
                return need_add_rxn, close_coupling_rxn

    def small_fix(self, model_info, model, initial_rxn_id, smallmolecule, general_library, model_control_info, check_model, obj):
        """"""
        right_small_met, small_need_add_rxn, close_coupling_rxn = [], [], []
        for metId in smallmolecule:
            print(metId,'xxxxxxxxxx')
            with model:
                small_demand_id = add_demand(model_info, model, str(metId), independent=True)
                model.objective = small_demand_id
                # write_flux_file(model_info, model, small_demand_id, metId, self.cfg)
                small_info = self.smallmolecule_iscoupling(model, model_info, small_demand_id, initial_rxn_id, metId, general_library, model_control_info, obj, check_model)
                if small_info == -1 :
                    continue
                elif small_info == 0: 
                    right_small_met.append(metId)
                else:
                    small_need_add_rxn.extend(small_info[0])
                    close_coupling_rxn.extend(small_info[1])
            if close_coupling_rxn :
                model_info['biomass']['small_macro_modify'] = "yes"
                model_control_info['check_biomass_production']["score"] = 0
                for r in close_coupling_rxn:
                    model.reactions.get_by_id(r).bounds = (0,0)
            if small_need_add_rxn:
                model_control_info['check_biomass_production']["score"] = 0
                model_info['biomass']['small_macro_modify'] = "yes" 
                self.add_gap_rxn(model, small_need_add_rxn, general_library, check_model)
        temp = [f"The synthesis of small molecule metabolites {right_small_met} were not coupled to macromolecules, and the synthetic pathway were correct"]
        write_biomass(model_control_info, temp, "Check_small_molecules")
        model.reactions.get_by_id(initial_rxn_id).bounds = (0,1000)

    def Normalized_macromolecule(
        self,
        model,
        model_info,
        metId,
        initialMacromoleculeResult,
        macro_demand_id,
        initial_rxn_id,
        normalized_macro_time,
        model_control_info,
        general_library,
        obj,
        check_model,
    ):
        """"""
        (
            close_coupling_rxn,
            initial_macromolecule_rxn,
            macromoleculeRxnId,
            need_add_rxn,
            close_coupling_rxn,
        ) = ([], "", "", [], [])
        sink_rxn = model_info["sink_rxn"]
        mets = model.metabolites.get_by_id(metId)
        print(model.slim_optimize())
        print(macro_demand_id)
        print(model.reactions.get_by_id(macro_demand_id))
        if model.slim_optimize() <= 1e-6:
            temp = [
                f"error : Macromolecules {metId} cannot be synthesized! Unable to normalize"
            ]
            write_biomass(model_control_info, temp, "Check_macromolecules")
            return 0
        while initialMacromoleculeResult:
            status = 0  # Exit when no biomass or macromolecule coupling is detected in this iteration.
            if model.slim_optimize() <= 1e-6 or model.slim_optimize() == "nan":
                temp = [
                    f"error : After the gap addition reaction, Macromolecules {metId} cannot be synthesized! Unable to normalize"
                ]
                write_biomass(model_control_info, temp, "Check_macromolecules")
                return 0
            print("test:......................", model.slim_optimize())

            if model.slim_optimize() > 0:
                pfba_solution = pfba(model)
                need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-6]

            else:
                need_fluxes = {}
            for r in list(
                need_fluxes.keys()
            ):  # Resolve macromolecule synthesis that depends on the biomass reaction.
                rxns = model.reactions.get_by_id(r)
                all_mets = [m.id for m in rxns.metabolites]
                met_formula = list(rxns.metabolites.keys())[0].formula
                if r == initial_rxn_id or unresolved_boundary_source(rxns, need_fluxes[r]):
                    if (
                        r != macro_demand_id and r not in sink_rxn
                    ):  # Check for coupling through biomass or macromolecular exchange reactions.
                        print(r, sink_rxn, "xxxxxx")
                        met_formula = self.get_formula(
                            model, all_mets[0], model_info, check_model, obj
                        )
                        if met_formula != 0:
                            continue
                        status += 1
                        model.reactions.get_by_id(r).bounds = (0, 0)
                        close_coupling_rxn.append(r)
                        temp_Macromolecule = self.calculateResult(
                            model,
                            macro_demand_id,
                            model_info,
                            sink_rxn,
                            check_model,
                            obj,
                        )  # Save the current production rate when coupling is present.
                        temp = [
                            f"{macro_demand_id} composition is coupled with {r}:{rxns}"
                        ]
                        print("temp_Macromolecule: ", temp_Macromolecule)
                        if type(temp_Macromolecule) == str:
                            pass
                        else:
                            temp.extend(
                                [
                                    f"After closing {r}, the relative molecular mass of {macro_demand_id} is {round(temp_Macromolecule,5)}"
                                ]
                            )
                            write_biomass(
                                model_control_info, temp, "Check_macromolecules"
                            )
                        if (
                            temp_Macromolecule == 0
                        ):  # Complete a missing pathway if the target can no longer be synthesized.
                            temp = [
                                f"All macromolecules are closed, the relative molecular mass of {macro_demand_id} is 0, and the gap is performed:"
                            ]
                            new_way_rxn = self.gapfilling(
                                model, model_info, macro_demand_id, general_library
                            )
                            need_add_rxn, gap_info = self.add_gapfilling_rxn(
                                model,
                                model_info,
                                check_model,
                                general_library,
                                new_way_rxn,
                                obj,
                            )
                            model_info["bio_modified"].extend(need_add_rxn)
                            temp.extend(gap_info)
                            write_biomass(
                                model_control_info, temp, "Check_macromolecules"
                            )
                            if len(need_add_rxn) == 0:
                                temp = [
                                    f"In biomass generation, {metId} must be coupled with macromolecules, the gapfilling cannot be performed, there is an error, please correct"
                                ]
                                write_biomass(
                                    model_control_info, temp, "Check_macromolecules"
                                )
                                return 0
                            break
                        else:  # Repeat until synthesis no longer depends on other macromolecule reactions.
                            initialMacromoleculeResult = temp_Macromolecule
            if status == 0:
                break
        if model.slim_optimize() > 0:
            pfba_solution = pfba(model)
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-6]
        else:
            need_fluxes = {}
        for (
            r,
            v,
        ) in (
            need_fluxes.items()
        ):  # Read the macromolecule synthesis reaction and its stoichiometric coefficients.
            rxns = model.reactions.get_by_id(r)
            products_mets = [m.id for m in rxns.products]
            reactants_mets = [m.id for m in rxns.reactants]
            coefficient_dict = {}
            if (mets.id in products_mets and v > 0) or (
                mets.id in reactants_mets and v < 0
            ):
                macromoleculeRxnId = rxns.id
                model_info["bio_modified"].append(macromoleculeRxnId)
                model_info["revision_num"] += 1
                initial_macromolecule_rxn = model.reactions.get_by_id(
                    macromoleculeRxnId
                )
                for (
                    j
                ) in rxns.metabolites:  # Identify the macromolecule synthesis reaction.

                    if (
                        str(j) == mets.id
                    ):  # Preserve the coefficient of the target macromolecule.
                        continue
                    coefficient_dict[j] = rxns.get_coefficient(j) * (
                        1000 / initialMacromoleculeResult
                    ) - rxns.get_coefficient(
                        j
                    )  # Apply the additional coefficient required by the scaling factor.
                rxns.add_metabolites(coefficient_dict)
                break
        if type(initial_macromolecule_rxn) == str:
            temp = [
                f"error : After removing the coupled macromolecules, Macromolecules {metId} cannot be synthesized! Unable to normalize"
            ]
            write_biomass(model_control_info, temp, "Check_macromolecules")
            return 0
        initial_macromolecule_rxn, initial_macromolecule_rxn_id = keep_decimal_places(
            initial_macromolecule_rxn
        )
        if normalized_macro_time > 1:
            temp = [
                f"Round {normalized_macro_time}: {metId} Equation: {initial_macromolecule_rxn}"
            ]
            if status != 0:
                temp.extend(
                    [
                        f"Round {normalized_macro_time}: {metId} coupling macromolecule {','.join(close_coupling_rxn)}, and {','.join(close_coupling_rxn)} formula is not clear!!!, unable to calculate relative to initial {metId} Molecular mass"
                    ]
                )
            else:
                temp.extend(
                    [
                        f"Round {normalized_macro_time}: {mets.id} relative molecular mass results: {round(initialMacromoleculeResult,5)}"
                    ]
                )
        else:
            temp = [f"Initial equation: {initial_macromolecule_rxn}"]
            if status != 0:
                temp.extend(
                    [
                        f"{metId} coupling macromolecule {','.join(close_coupling_rxn)}, and {','.join(close_coupling_rxn)} formula is not clear!!!, the relative molecular mass of the initial {metId} cannot be calculated"
                    ]
                )
            else:
                temp.extend(
                    [
                        f"Initial relative molecular mass: {round(initialMacromoleculeResult,5)}"
                    ]
                )
        write_biomass(model_control_info, temp, "Check_macromolecules")
        return macromoleculeRxnId, need_add_rxn, close_coupling_rxn

    def recover_macromolecule_rxn(
        self,
        model,
        model_info,
        macromoleculeRxn,
        need_add_rxn,
        general_library,
        need_close_coupling_rxn,
        initial_rxn_id,
        model_control_info,
        check_model,
    ):
        """"""
        for k, v in macromoleculeRxn.items():
            model_control_info["check_biomass_production"]["score"] = 0
            model_info["biomass"]["small_macro_modify"] = "yes"
            model.reactions.get_by_id(k).reaction = v
        if need_add_rxn:
            model_control_info["check_biomass_production"]["score"] = 0
            model_info["biomass"]["small_macro_modify"] = "yes"
            self.add_gap_rxn(model, need_add_rxn, general_library, check_model)
        for r in need_close_coupling_rxn:
            model_control_info["check_biomass_production"]["score"] = 0
            model_info["biomass"]["small_macro_modify"] = "yes"
            model.reactions.get_by_id(r).bounds = (0, 0)
        model.reactions.get_by_id(initial_rxn_id).bounds = (0, 1000)

    def macro_fix(
        self,
        model_info,
        model,
        initial_rxn_id,
        macromolecule,
        general_library,
        model_control_info,
        check_model,
        obj,
    ):
        """"""
        gapfilling_rxn, need_close_coupling_rxn, macromoleculeRxn = [], [], {}
        sink_rxn = model_info["sink_rxn"]
        for metId in macromolecule:
            temp = [f"{metId}:"]
            write_biomass(model_control_info, temp, "Check_macromolecules")
            normalized_macro_time, while_time = 0, 0
            with model:
                macro_demand_id = add_demand(model_info, model, str(metId))
                model.objective = macro_demand_id
                print(model.slim_optimize(), macro_demand_id, "xxx")
                try:
                    model.reactions.get_by_id(macro_demand_id).bounds = (
                        check_model.reactions.get_by_id(macro_demand_id).bounds
                    )  # Restore macromolecule reaction bounds.
                except:
                    print(
                        f"{model.id} does not contain {macro_demand_id}; no reaction restoration is needed"
                    )
                print(f"{metId}_objective_value: ", model.slim_optimize())

                initialMacromoleculeResult = self.calculateResult(
                    model, macro_demand_id, model_info, sink_rxn, check_model, obj
                )
                print(model.slim_optimize(), macro_demand_id, "xxx")
                if initialMacromoleculeResult == 0:
                    continue
                if type(initialMacromoleculeResult) != str and (
                    990 < initialMacromoleculeResult < 1010
                ):
                    continue
                else:
                    normalized_macro_time += 1
                    print(model.slim_optimize(), macro_demand_id, "xxx")
                    macromoleculeInfo = self.Normalized_macromolecule(
                        model,
                        model_info,
                        metId,
                        initialMacromoleculeResult,
                        macro_demand_id,
                        initial_rxn_id,
                        normalized_macro_time,
                        model_control_info,
                        general_library,
                        obj,
                        check_model,
                    )
                    if macromoleculeInfo == -1:
                        return -1
                    elif macromoleculeInfo == 0:
                        continue
                    else:
                        macromoleculeRxnId, need_add_rxn, close_coupling_rxn = (
                            macromoleculeInfo
                        )
                    finalMacromoleculeResult = self.calculateResult(
                        model, macro_demand_id, model_info, sink_rxn, check_model, obj
                    )
                    gapfilling_rxn.extend(need_add_rxn)
                    need_close_coupling_rxn.extend(close_coupling_rxn)
                    while (
                        type(finalMacromoleculeResult) == str
                        or finalMacromoleculeResult < 990
                        or finalMacromoleculeResult > 1010
                    ):
                        normalized_macro_time += 1
                        macromoleculeInfo = self.Normalized_macromolecule(
                            model,
                            model_info,
                            metId,
                            finalMacromoleculeResult,
                            macro_demand_id,
                            initial_rxn_id,
                            normalized_macro_time,
                            model_control_info,
                            general_library,
                            obj,
                            check_model,
                        )
                        if macromoleculeInfo == -1:
                            return -1
                        elif macromoleculeInfo == 0:  # Exit the inner loop only.
                            while_time += 1
                            break
                        else:
                            macromoleculeRxnId, need_add_rxn, close_coupling_rxn = (
                                macromoleculeInfo
                            )
                        finalMacromoleculeResult = self.calculateResult(
                            model,
                            macro_demand_id,
                            model_info,
                            sink_rxn,
                            check_model,
                            obj,
                        )
                        gapfilling_rxn.extend(need_add_rxn)
                        need_close_coupling_rxn.extend(close_coupling_rxn)
                    if while_time == 1:
                        continue
                    finalRxn = model.reactions.get_by_id(macromoleculeRxnId)
                    finalRxn, finalRxn_id = keep_decimal_places(finalRxn)
                    final_rxn_id = model.reactions.get_by_id(
                        macromoleculeRxnId
                    ).build_reaction_string()
                    macromoleculeRxn[macromoleculeRxnId] = final_rxn_id
                    temp = [f"Corrected equation: {finalRxn}"]
                    temp.extend(
                        [
                            f"Corrected relative molecular mass: {round(finalMacromoleculeResult,5)}"
                        ]
                    )
                    write_biomass(model_control_info, temp, "Check_macromolecules")
            if need_close_coupling_rxn:
                for rxnId in need_close_coupling_rxn:
                    model.reactions.get_by_id(rxnId).bounds = (0, 0)
            if gapfilling_rxn:
                self.add_gap_rxn(model, gapfilling_rxn, general_library, check_model)
        self.recover_macromolecule_rxn(
            model,
            model_info,
            macromoleculeRxn,
            gapfilling_rxn,
            general_library,
            need_close_coupling_rxn,
            initial_rxn_id,
            model_control_info,
            check_model,
        )

    def modify_nadh(self, model_info, model, need_fluxes, grate_delete_rxn):
        """"""
        for ids, v in need_fluxes.items():
            for rxn in model_info["reactions"]:
                if ids == rxn["id"]:
                    rxns = model.reactions.get_by_id(ids)
                    rxn_bounds = rxns.bounds
                    if len(set(NADH_NAME) & set(rxn["reactants_mets"])) != 0 and v < 0:
                        rxns.bounds = (0, 1000)
                        rxn["bounds"] = [0, 1000]
                        rules2 = Rules2(ids, self.cfg, v)
                        rules2.get_all_rules(model_info, model)
                        if "rule_conflict" in rxn.keys():
                            rxns.bounds = rxn_bounds
                            rxn["bounds"] = rxn_bounds
                        else:
                            rxn["rules"]["modify_nadh_rxn"] = "true"
                            grate_delete_rxn.append(ids)
                            rxn["nadhs_modify"] = "true"
                            print(
                                ids,
                                "Constraining the implicated reaction set",
                                model.reactions.get_by_id(ids).build_reaction_string(
                                    use_metabolite_names=True
                                ),
                                model.reactions.get_by_id(ids).bounds,
                            )
                            return
                    if len(set(NADH_NAME) & set(rxn["products_mets"])) != 0 and v > 0:
                        rxns.bounds = (-1000, 0)
                        rxn["bounds"] = [-1000, 0]
                        rules2 = Rules2(ids, self.cfg, v)
                        rules2.get_all_rules(model_info, model)
                        if "rule_conflict" in rxn.keys():
                            rxns.bounds = rxn_bounds
                            rxn["bounds"] = rxn_bounds
                        else:
                            rxn["rules"]["modify_nadh_rxn"] = "true"
                            grate_delete_rxn.append(ids)
                            rxn["nadhs_modify"] = "true"
                            print(
                                ids,
                                "Constraining the implicated reaction set",
                                model.reactions.get_by_id(ids).build_reaction_string(
                                    use_metabolite_names=True
                                ),
                                model.reactions.get_by_id(ids).bounds,
                            )
                            return

    def close_max_value_nadh(
        self,
        model_info,
        model_control_info,
        model,
        check_model,
        obj,
        nadh_id,
        general_library,
    ):
        """
        Turn off the maximum value in each result
        """
        (
            grate_delete_rxn,
            count,
            fluxs,
            count2,
            now_rely_rxn,
            need_gap_rxn,
            temp_rxnId,
        ) = ([], 0, 0, 0, {}, [], "")
        if model.slim_optimize() <= 1e-6:
            print(model, "________NADH_yes")
        while model.slim_optimize() > 1e-6:
            count = count + 1

            if fluxs != round(model.slim_optimize(), 3):
                fluxs = round(model.slim_optimize(), 3)
                count2 = 0
            else:
                count2 += 1
            print(
                count,
                count2,
                fluxs,
                round(model.slim_optimize(), 3),
                sep="................",
            )
            if count2 > 25:
                print(
                    "Redox correction did not converge; check additional nutrient inputs or erroneous reactions"
                )
                write_flux_file(model_info, model, nadh_id, "biomass_nadh", self.cfg)
                return

            if model.slim_optimize() > 0:
                pfba_solution = pfba(model)
                need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-6]
            else:
                need_fluxes = {}
            if count2 > 5:
                self.modify_nadh(model_info, model, need_fluxes, grate_delete_rxn)
            fba_rxn_dic, flux_fic = {}, {}
            rely_rxn = get_rely_rxn(
                model, model_info, check_model, nadh_id, need_fluxes
            )
            print("rely_rxn:", rely_rxn)
            now_rely_rxn = {}
            for ids, v in need_fluxes.items():
                for rxn in model_info["reactions"]:
                    if ids == rxn["id"]:

                        if ids in rely_rxn:
                            now_rely_rxn.update({ids: rxn["net_penalty_points"]})
                        if ids == nadh_id or ids in rely_rxn:
                            continue
                        fba_rxn_dic[ids] = rxn["net_penalty_points"]
                        flux_fic[ids] = v
            if count2 == 15 and now_rely_rxn:
                print("now_rely_rxn:", now_rely_rxn)
                for rxnId, penalty_points in now_rely_rxn.items():
                    if penalty_points == max(now_rely_rxn.values()):
                        model.reactions.get_by_id(rxnId).bounds = (0, 0)
                        temp_rxnId = rxnId
                        print(
                            "Temporarily constraining growth-dependent reaction:",
                            temp_rxnId,
                            ".............",
                        )
                        break
                with model:
                    initial_rxn_id = model_info["initial_rxn"]["biomass_rxn_id"]
                    model.objective = initial_rxn_id
                    print(model.slim_optimize(), "...............")
                    big_model = self.get_big_model(model_info, model, general_library)
                    no_synthesized_mets = self.check_bio_component_is_zero(
                        model_info, model, model_control_info, initial_rxn_id
                    )
                    new_way_rxn = self.get_new_way_rxn(
                        model_info, big_model, no_synthesized_mets, model_control_info
                    )
                    need_add_rxn, gap_info = self.add_gapfilling_rxn(
                        model,
                        model_info,
                        check_model,
                        general_library,
                        new_way_rxn,
                        obj,
                    )
                    print("need_add_rxn: ", need_add_rxn)
                    print(model.slim_optimize(), "...............")
                    if need_add_rxn:
                        for rxnId in need_add_rxn:
                            bounds = model.reactions.get_by_id(rxnId).bounds
                            model.reactions.get_by_id(rxnId).bounds = (0, 0)
                            if (
                                model.slim_optimize() < 1e-6
                                or model.slim_optimize() == "nan"
                            ):
                                need_gap_rxn.append(rxnId)
                                model.reactions.get_by_id(rxnId).bounds = bounds
                if need_gap_rxn:
                    print("need_gap_rxn: ", need_gap_rxn)
                    self.add_gap_rxn(model, need_gap_rxn, general_library, check_model)
                    model_info["need_close_rely_rxn"] = temp_rxnId
                    rules = Rules(self.cfg)
                    run_rules(rules, model_info, model)
                if not need_add_rxn and temp_rxnId:
                    print(
                        "Restoring growth-dependent reaction: ",
                        temp_rxnId,
                        ".....................",
                    )
                    model.reactions.get_by_id(temp_rxnId).bounds = (
                        check_model.reactions.get_by_id(temp_rxnId).bounds
                    )
            for ids, penalty_points in fba_rxn_dic.items():

                if penalty_points == max(
                    fba_rxn_dic.values()
                ):  # Constrain all reactions tied for the highest score in the current route.
                    if ids != nadh_id and fba_rxn_dic[ids] != 0:
                        check_rule(model_info, model, ids, flux_fic[ids], self.cfg)

                        grate_delete_rxn.append(ids)
                        print(
                            ids,
                            " ",
                            model.reactions.get_by_id(ids).build_reaction_string(
                                use_metabolite_names=True
                            ),
                            model.reactions.get_by_id(ids).bounds,
                            fba_rxn_dic[ids],
                        )
            print("close :", model.slim_optimize(), "............", grate_delete_rxn)
        return grate_delete_rxn

    def check_nadh(
        self, model_info, model, check_model, model_control_info, obj, general_library
    ):
        """"""
        nadh_name, num = [
            "NADH",
            "NADPH",
            "FADH2",
            "FMNH2",
            "Q8H2",
            "MQL8",
            "DMMQL8",
        ], 1e-6
        reas = [NADH_NAME, NADPH, FADH2, FMNH2, Q8H2_NAME, MQL8, DMMQL8]
        pro1s = [NAD_NAME, NADP, FAD, FMN, Q8_NAME, MQN8, DMMQ8]
        pro2 = H_NAME
        for i, name in enumerate(nadh_name):
            nadh_id = get_other_nadh_rxn(model, reas[i], pro1s[i], pro2, name)
            if not nadh_id:
                continue
            model.reactions.get_by_id(nadh_id).bounds = (0, 1000)
            set_model_objective(model_info, model, nadh_id)
            close_autotrophic_or_c_source(model_info, model, "carbon_source")
            grate_delete_rxn = self.close_max_value_nadh(
                model_info,
                model_control_info,
                model,
                check_model,
                obj,
                nadh_id,
                general_library,
            )
            if grate_delete_rxn:
                all_need_fluxes, infinite_rxn = add_back_reactions_net(
                    model_info,
                    model,
                    check_model,
                    grate_delete_rxn,
                    "nadhs",
                    num,
                    self.cfg,
                )

    def Normalized_biomass(self, model, initial_rxn_id, final_biomass):
        """"""
        coefficient_dict = {}
        if final_biomass == 0:
            return -1
        rxns = model.reactions.get_by_id(initial_rxn_id)
        for j in rxns.metabolites:
            if (j in rxns.products) and (
                "biomass" in str(j).lower()
                or str(j).lower() in "biomass"
                or j.name in BIOMASS
            ):  # Preserve the coefficient of the biomass product metabolite.
                continue
            coefficient_dict[j] = rxns.get_coefficient(j) * (
                1000 / final_biomass
            ) - rxns.get_coefficient(j)
        rxns.add_metabolites(coefficient_dict)
        print(rxns.reaction)
        return 0

    def Normalized_sustaining_met(self, model, model_info):
        """Restore the recorded GAM using a signed, complete ATP-hydrolysis delta.

        An absent/ambiguous original ADP coefficient is not a measured zero GAM.
        In that case retain the current equation and record why restoration was
        skipped. Never change just part of the ATP hydrolysis reaction.
        """
        initial = model_info["initial_rxn"]
        ini_adp_coefficient = initial.get("adp_coefficient")
        initial_rxn_id = initial["biomass_rxn_id"]
        rxns = model.reactions.get_by_id(initial_rxn_id)
        state = model_info.setdefault("gam_restoration", {})
        if ini_adp_coefficient is None or not np.isfinite(ini_adp_coefficient) or ini_adp_coefficient < 0:
            state.update(status="skipped", reason="Original ADP coefficient missing or invalid")
            return
        groups = [(ATP_NAME, -1), (H2O_NAME, -1), (ADP_NAME, 1), (H_NAME, 1), (PI_NAME, 1)]
        selected = []
        for names, sign in groups:
            candidates = [met for met in rxns.metabolites if met.name in names]
            if len(candidates) != 1:
                state.update(status="skipped", reason="Missing or ambiguous ATP-hydrolysis metabolite")
                return
            selected.append((candidates[0], sign))
        if len({met.compartment for met, _ in selected}) != 1:
            state.update(status="skipped", reason="ATP-hydrolysis metabolites span multiple compartments")
            return
        fin_adp_coefficient = rxns.get_coefficient(selected[2][0].id)
        fin_coefficient = float(ini_adp_coefficient) - fin_adp_coefficient
        if not np.isfinite(fin_coefficient):
            state.update(status="skipped", reason="Non-finite current ADP coefficient")
            return
        rxns.add_metabolites({met: sign * fin_coefficient for met, sign in selected})
        state.update(status="restored", adp_target=float(ini_adp_coefficient), hydrolysis_delta=fin_coefficient)
        print(rxns.reaction)

    def biomass_fix(
        self, model_info, model, initial_rxn_id, model_control_info, check_model, obj
    ):
        """"""
        with model:
            middle_bio_value = model.slim_optimize()
            print(model.slim_optimize(), "xxxxxx")
            sink_rxn = model_info["sink_rxn"]

            print(model.slim_optimize(), "xxxxxxxx")

            modify_nor_bio = final_biomass = self.calculateResult(
                model, initial_rxn_id, model_info, sink_rxn, check_model, obj
            )
            modify_nor_bio_eq = model.reactions.get_by_id(
                initial_rxn_id
            ).build_reaction_string(use_metabolite_names=False)
            model_info["biomass"]["middle_bio_value"] = middle_bio_value
            model_info["biomass"]["modify_nor_bio_eq"] = modify_nor_bio_eq
            model_info["biomass"]["modify_nor_bio"] = modify_nor_bio
            if type(final_biomass) == str:
                model_control_info["check_biomass_production"]["score"] = 0
                print(model.slim_optimize(), "sssssssssssssssssssss")
                temp = [
                    f"The final biomass is still coupled with macromolecules {final_biomass}"
                ]
                write_biomass(model_control_info, temp, "Normalization_of_biomass")
                model_info["bio_coupling"] = "0"
                return -1
            model_info["ini_weight(g)"] = round(final_biomass / 1000, 3)
            if type(final_biomass) != str and (990 < final_biomass < 1010):
                temp = [
                    f"The biomass is 1g after treatment with small molecules and macromolecules"
                ]
            else:
                model_control_info["check_biomass_production"]["score"] = 0
                temp = [
                    f"When verifying the biomass weight of 1 g. The biomass was {round(final_biomass/1000,3)}g and did not 1g."
                ]
            write_biomass(model_control_info, temp, "Normalization_of_biomass")
            while final_biomass < 990 or final_biomass > 1010:
                get_final = self.Normalized_biomass(
                    model, initial_rxn_id, final_biomass
                )
                print("get_final:.............", get_final)
                if get_final == -1:
                    model_control_info["check_biomass_production"]["score"] = 0
                    temp = [f"The final flux final is 0, error"]
                    write_biomass(model_control_info, temp, "Normalization_of_biomass")
                    return -1
                final_biomass = self.calculateResult(
                    model, initial_rxn_id, model_info, sink_rxn, check_model, obj
                )
                print("final_biomass:.................", final_biomass)
            final_biomass_equation = model.reactions.get_by_id(initial_rxn_id)
            final_biomass_equation, final_biomass_equation_ids = keep_decimal_places(
                final_biomass_equation
            )
            final_biomass_equation_id = model.reactions.get_by_id(
                initial_rxn_id
            ).build_reaction_string()
        print("Recording final biomass output", final_biomass_equation_id)
        model_info["biomass"]["final_biomass"] = final_biomass
        model_info["biomass"]["final_biomass_equation"] = final_biomass_equation
        model.reactions.get_by_id(initial_rxn_id).reaction = final_biomass_equation_id

    def recover_atpm_bounds(self, model_info, model, check_model):
        """Undo only the temporary maintenance relaxation, including all aliases."""
        restore_bounds = model_info.get('atpm_restore_bounds', {})
        for rxn in model_info['reactions']:
            if 'atpms_modify' in rxn.keys() and rxn['atpms_modify'] == "true":
                bounds = restore_bounds.get(rxn['id'])
                if bounds is None:  # Compatibility with previously saved model_info.
                    bounds = check_model.reactions.get_by_id(rxn['id']).bounds
                model.reactions.get_by_id(rxn['id']).bounds = tuple(bounds)
                restore_bounds.pop(rxn['id'], None)
                rxn['atpms_modify'] = "false"

    def convert_nan_to_null(self, model):
        """"""
        for met in model.metabolites:

            if (
                pd.isna(met.formula)
                or met.formula == "nan"
                or met.formula is None
                or not met.formula
                or met.formula == "null"
            ):
                print(met, met.formula)
                met.formula = ""

    def get_exchange(self, model, model_control_info):
        if model.slim_optimize() > 0:
            pfba_solution = pfba(model)
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-6]
        else:
            need_fluxes = {}

        temp = [f"Model calculation conditions:"]
        for r, v in need_fluxes.items():
            rxn = model.reactions.get_by_id(r)
            all_mets = [m.id for m in rxn.metabolites]
            if len(all_mets) == 1:
                mets = model.metabolites.get_by_id(all_mets[0])
                temp.extend(
                    [
                        f"{r}: {rxn.build_reaction_string(use_metabolite_names=True)}    {mets.formula}  {round(v,5)}  {rxn.bounds}"
                    ]
                )
        write_biomass(model_control_info, temp, "model_calculation_conditions")

    def get_result(
        self, model_info, model, sustain_result, model_control_info, initial_rxn_id
    ):
        """"""
        final_bio_value = model.slim_optimize()
        initial_biomass = model_info["biomass"]["initial_biomass"]
        if type(initial_biomass) == str:
            sustain_result.extend(
                [
                    f"Biomass coupling macromolecule {initial_biomass}, and the formula of {initial_biomass} is unclear!!!, the relative molecular mass of the initial biomass cannot be calculated"
                ]
            )
        else:
            sustain_result.extend(
                [
                    f"Initial biomass relative molecular mass: {round(initial_biomass,5)}mg"
                ]
            )
        if "small_macro_modify" in model_info["biomass"].keys():
            modify_nor_bio = model_info["biomass"]["modify_nor_bio"]
            modify_nor_bio_eq = model_info["biomass"]["modify_nor_bio_eq"]
            middle_bio_value = model_info["biomass"]["middle_bio_value"]
            if type(modify_nor_bio) == str:
                sustain_result.extend(
                    [
                        f"Modifying the biomass equation after the reaction of small molecules and macromolecules: {modify_nor_bio_eq}"
                    ]
                )
                sustain_result.extend(
                    [
                        f"The biomass is coupled to the macromolecule {modify_nor_bio}, and the formula of {modify_nor_bio} is not clear!!!, it is impossible to calculate the relative molecular mass of the biomass after modifying the reaction of small molecules and macromolecules"
                    ]
                )
            else:
                sustain_result.extend(
                    [
                        f"Modifying the biomass equation after the reaction of small molecules and macromolecules: {modify_nor_bio_eq}"
                    ]
                )
                sustain_result.extend(
                    [
                        f"Modify the relative molecular mass results of biomass after the reaction of small molecules and macromolecules: {round(modify_nor_bio,5)}"
                    ]
                )
            sustain_result.extend(
                [
                    f"After modifying the reaction of small molecules and macromolecules, the growth rate of biomass is: {round(middle_bio_value,5)}"
                ]
            )
        final_biomass_equation = model_info["biomass"]["final_biomass_equation"]
        final_biomass = model_info["biomass"]["final_biomass"]
        sustain_result.extend(
            [
                f"After modifying the biomass, the final biomass equation: {final_biomass_equation}"
            ]
        )
        sustain_result.extend(
            [
                f"After modifying the biomass, the final biomass relative molecular mass result: {round(final_biomass,5)}"
            ]
        )
        sustain_result.extend(
            [
                f"The final relative molecular mass of biomass after modification is {round(final_biomass,5)}mg, and the final biomass growth rate is {round(final_bio_value,5)}"
            ]
        )
        write_biomass(model_control_info, sustain_result, "Normalization_of_biomass")
        model_info["bio_modified"].extend([model_info["initial_rxn"]["biomass_rxn_id"]])
        model_info["revision_num"] += 1
        fins = get_c_mol(model)
        model_info["c_bio"] = fins
        model_control_info["check_biomass_production"][
            "Final_biomass"
        ] = f"Initial growth rate: {model_info['biomass_original']}; after setting the carbon source, the value becomes: {model_info['biomass_by_c_source']}; after setting special boundaries, the value becomes: {model_info['biomass_by_special_boundary']}; after quality control, final biomass growth rate: {round(final_bio_value,3)}; single carbon source {fins}.\nThe equations modified during biomass normalization correction are placed in the 'Model Revision' list, and you can make changes and edits by clicking the 'Edit' button."
        model_info["control_analysis"].append(
            (
                model_info["biomass_original"],
                model_info["biomass_by_c_source"],
                round(final_bio_value, 3),
            )
        )

    def modify_bios(self, model_info, model, model_control_info):
        """"""
        temp = []
        bio_modified = list(set(model_info["bio_modified"]))
        for rxnId in bio_modified:

            rxns = model.reactions.get_by_id(rxnId)
            equation, equation_id = keep_decimal_places(rxns)
            modify_rxn = {
                "reaction_id": rxnId,
                "equation(ID)": equation_id,
                "equation(Name)": equation,
                "l_b": rxns.lower_bound,
                "u_b": rxns.upper_bound,
            }
            temp.append(modify_rxn)
        model_control_info["check_biomass_production"]["model_revision"].extend(temp)

    def recover_special_boundary_rxn(self, model_info, model, check_model):
        """"""
        for rxn in model_info["reactions"]:
            if "special_boundary" in rxn.keys():
                model.reactions.get_by_id(rxn["id"]).bounds = (
                    check_model.reactions.get_by_id(rxn["id"]).bounds
                )

    def get_final_fluxes(
        self, model_info, model, check_model, model_control_info, initial_rxn_id
    ):
        """"""
        atpm_id = get_atpm_rxn(model_info, model, "ATP", "ADP")
        set_model_objective(model_info, model, atpm_id)
        final_atp = round(model.slim_optimize(), 3)

        if "ADD_ATPM" in model.reactions:
            model.reactions.remove("ADD_ATPM")
        nadh_id = add_nadh_rxn(model)
        set_model_objective(model_info, model, nadh_id)
        final_nadh = round(model.slim_optimize(), 3)

        model.reactions.remove("ADD_NADH")
        final_netld_flux = get_final_net_fluxes(model_info, model, check_model)
        if is_autotrophic(model_info):
            final_yield_flux = ""
        else:
            final_yield_flux = get_final_yield_fluxes(model_info, model)
        model.objective = initial_rxn_id
        final_biomass_flux = round(model.slim_optimize(), 3)
        print("final_biomass_flux:", final_biomass_flux)
        print(model_control_info["final_flux"])
        model_control_info["final_flux"] = {
            "final_biomass": final_biomass_flux,
            "final_nadh": final_nadh,
            "final_atp": final_atp,
            "final_netld": final_netld_flux,
            "final_yield": final_yield_flux,
        }

    def biomass_control(
        self, model_info, model, check_model, model_control_info, obj, model_check_info
    ):
        """"""
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), ":biomass start")
        model.solver = "cplex"
        model_control_info["check_biomass_production"] = {}
        model_control_info["check_biomass_production"]["score"] = 1
        model_control_info["check_biomass_production"]["Biomass_information"] = []

        model_control_info["check_biomass_production"][
            "Check_synthesizability_of_biomass_components"
        ] = []
        model_control_info["check_biomass_production"][
            "Gapfilling_of_biomass_components"
        ] = []
        model_control_info["check_biomass_production"]["Check_small_molecules"] = []
        model_control_info["check_biomass_production"]["Check_macromolecules"] = []
        model_control_info["check_biomass_production"]["Normalization_of_biomass"] = []
        model_control_info["check_biomass_production"]["Final_biomass"] = ""
        model_control_info["check_biomass_production"]["model_revision"] = []

        (
            model_info["biomass"],
            model_info["bio_modified"],
            model_info["no_gap_rxn_input"],
        ) = ({}, [], [])

        if not is_bio_exp_correct(model_info):
            model_control_info["check_biomass_production"][
                "Biomass_information"
            ].append(
                "Could not find the biomass equation, please check and try again!!!"
            )
            model_control_info["check_biomass_production"]["score"] = 0
            model_info["bio_norxn"] = 0
            return 0
        initial_rxn_id = model_info["initial_rxn"]["biomass_rxn_id"]
        model.reactions.get_by_id(initial_rxn_id).bounds = (0, 1000)
        self.get_biomass_info(model_info, model_control_info)
        self.recover_bio_objective(model_info, model)
        print("Biomass output during curation:", initial_rxn_id, model.slim_optimize())
        if rule_policy_enabled(self.cfg):
            apply_native_rules(model, model_info)
            self.cfg.big_general_library = ""
        self.add_sink_rxn(model, model_info, check_model)

        general_library = self.get_general_library(model)
        model_info['gapfill_rejected_redox_shortcuts'] = self.rejected_gap_redox_shortcuts
        print(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            ": Candidate library loaded; checking biomass precursors",
        )
        no_synthesized_mets = self.check_bio_component_is_zero(
            model_info, model, model_control_info, initial_rxn_id
        )
        if (
            self.preprocess_initial_zero_bio(
                model_info,
                model,
                model_control_info,
                general_library,
                no_synthesized_mets,
                obj,
                check_model,
            )
            == -1
        ):
            allocate_auxiliary_supply(model, model_info, model_control_info)
            fins = get_c_mol(model)
            model_info["c_bio"] = fins
            write_flux_file(model_info, model, initial_rxn_id, "biomasss", self.cfg)
            model_control_info["check_biomass_production"][
                "Final_biomass"
            ] = f"Initial growth rate: {model_info['biomass_original']}; after setting the carbon source, the value becomes: {model_info['biomass_by_c_source']}; after setting special boundaries, the value becomes: {model_info['biomass_by_special_boundary']}; after quality control, final biomass growth rate: {round(model.slim_optimize(), 3)}; single carbon source {fins}.\nThe equations modified during biomass normalization correction are placed in the 'Model Revision' list, and you can make changes and edits by clicking the 'Edit' button."
            model_info["control_analysis"].append(
                (
                    model_info["biomass_original"],
                    model_info["biomass_by_c_source"],
                    round(model.slim_optimize(), 3),
                )
            )
            self.recover_atpm_bounds(model_info, model, check_model)
            self.modify_bios(model_info, model, model_control_info)

            return 0

        use_ATPS4rpp_or_ATPS4(model)
        print(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            ": Checking biomass-associated maintenance",
        )
        sustain_result = self.is_missing_sustaining_energy(
            model_info, model, model_control_info, initial_rxn_id
        )
        macromolecule = self.get_macromolecule(
            model, initial_rxn_id, model_control_info
        )
        smallmolecule = self.get_smallmolecule(
            model_info, model, initial_rxn_id, macromolecule, model_control_info
        )

        self.get_initial_biomass(
            model_info, model, initial_rxn_id, model_control_info, check_model, obj
        )
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), ":small_fix start")
        self.small_fix(
            model_info,
            model,
            initial_rxn_id,
            smallmolecule,
            general_library,
            model_control_info,
            check_model,
            obj,
        )

        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), ":macro_fix start")
        self.macro_fix(
            model_info,
            model,
            initial_rxn_id,
            macromolecule,
            general_library,
            model_control_info,
            check_model,
            obj,
        )
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), ":nadh_fix start")

        print(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), ":biomass_fix start"
        )
        if (
            self.biomass_fix(
                model_info, model, initial_rxn_id, model_control_info, check_model, obj
            )
            == -1
        ):
            allocate_auxiliary_supply(model, model_info, model_control_info)
            fins = get_c_mol(model)
            model_info["c_bio"] = fins
            write_flux_file(model_info, model, initial_rxn_id, "biomasss", self.cfg)
            model_control_info["check_biomass_production"][
                "Final_biomass"
            ] = f"Initial growth rate: {model_info['biomass_original']}; after setting the carbon source, the value becomes: {model_info['biomass_by_c_source']}; after setting special boundaries, the value becomes: {model_info['biomass_by_special_boundary']}; after quality control, final biomass growth rate: {round(model.slim_optimize(), 3)}; single carbon source {fins}.\nThe equations modified during biomass normalization correction are placed in the 'Model Revision' list, and you can make changes and edits by clicking the 'Edit' button."
            model_info["control_analysis"].append(
                (
                    model_info["biomass_original"],
                    model_info["biomass_by_c_source"],
                    round(model.slim_optimize(), 3),
                )
            )
            self.recover_atpm_bounds(model_info, model, check_model)
            self.modify_bios(model_info, model, model_control_info)

            return 0
        self.recover_atpm_bounds(model_info, model, check_model)

        write_flux_file(model_info, model, initial_rxn_id, "biomasss", self.cfg)
        self.Normalized_sustaining_met(model, model_info)
        allocate_auxiliary_supply(model, model_info, model_control_info)
        self.get_result(
            model_info, model, sustain_result, model_control_info, initial_rxn_id
        )

        print("Biomass output before maintenance normalization:", model.slim_optimize())
        write_flux_file(model_info, model, initial_rxn_id, "biomasss2", self.cfg)
        self.modify_bios(model_info, model, model_control_info)
