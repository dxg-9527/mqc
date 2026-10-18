from mqc.config import Config
from mqc.utils import *
from mqc.control.rules_control import Rules
from mqc.control.rules_control2 import Rules2
from mqc.control.biomass_control import Biomasses


class Nadhs:

    def __init__(self, cfg: Config):
        """"""
        self.cfg = cfg

    def old_get_model_initial_flux(self, model_info, model, model_control_info):
        """"""
        with model:
            model_info["initial_flux"] = {}
            model_info["initial_flux"]["initial_biomass"] = model_info["initial_rxn"]
            model_control_info["initial_flux"]["initial_biomass"] = round(
                model_info["initial_rxn"]["biomass_rxn_flux"], 3
            )
            nadh_id = add_nadh_rxn(model)
            model.reactions.get_by_id(nadh_id).bounds = (0, 1000)
            set_model_objective(model_info, model, nadh_id)
            model_info["initial_flux"]["initial_nadh"] = model_info["now_obj_rxn"]
            model_control_info["initial_flux"]["initial_nadh"] = round(
                model_info["now_obj_rxn"]["rxn_flux"], 3
            )

            atpm_id = get_atpm_rxn(model_info, model, "ATP", "ADP")
            model.reactions.get_by_id(atpm_id).bounds = (0, 1000)
            set_model_objective(model_info, model, atpm_id)
            model_info["initial_flux"]["initial_atp"] = model_info["now_obj_rxn"]
            model_control_info["initial_flux"]["initial_atp"] = round(
                model_info["now_obj_rxn"]["rxn_flux"], 3
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
        name,
        temp,
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
                write_flux_file(model_info, model, nadh_id, "nadh", self.cfg)
                model_control_info["check_reducing_equivalents_production"][name][
                    "Optimal_rate_in_initial_model"
                ] = ""
                model_control_info["check_reducing_equivalents_production"][name][
                    "Details"
                ] = "error: Infinite cycle of reducing equivalents, there may be additional substrate inputs or erroneous reactions in the model"
                model_control_info["check_reducing_equivalents_production"][name][
                    "model_revision"
                ] = []
                model_control_info["check_reducing_equivalents_production"][name][
                    "related_reaction"
                ] = []
                model_control_info["check_reducing_equivalents_production"][name][
                    "Optimal_rate_in_revised_model"
                ] = ""

                if (
                    "Check_synthesizability_of_biomass_components"
                    in model_control_info["check_biomass_production"]
                ):
                    del model_control_info["check_biomass_production"][
                        "Check_synthesizability_of_biomass_components"
                    ]
                if (
                    "Gapfilling_of_biomass_components"
                    in model_control_info["check_biomass_production"]
                ):
                    del model_control_info["check_biomass_production"][
                        "Gapfilling_of_biomass_components"
                    ]
                for rxnId in grate_delete_rxn:
                    model.reactions.get_by_id(rxnId).bounds = (
                        check_model.reactions.get_by_id(rxnId).bounds
                    )
                raise RuntimeError("error: nadh infinite loop")

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
                    biomasses = Biomasses(self.cfg)
                    general_library = biomasses.get_general_library(model)

                    model_control_info["check_biomass_production"][
                        "Check_synthesizability_of_biomass_components"
                    ] = []
                    model_control_info["check_biomass_production"][
                        "Gapfilling_of_biomass_components"
                    ] = []
                    big_model = biomasses.get_big_model(
                        model_info, model, general_library
                    )
                    no_synthesized_mets = biomasses.check_bio_component_is_zero(
                        model_info, model, model_control_info, initial_rxn_id
                    )
                    new_way_rxn = biomasses.get_new_way_rxn(
                        model_info, big_model, no_synthesized_mets, model_control_info
                    )
                    need_add_rxn, gap_info = biomasses.add_gapfilling_rxn(
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
                    biomasses.add_gap_rxn(
                        model, need_gap_rxn, general_library, check_model
                    )
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

    def get_final_fluxes(self, model_info, model, check_model, model_control_info):
        """"""
        with model:
            for rxn in model_info["reactions"]:
                if "nadhs_modify" not in rxn.keys():
                    rxns = model.reactions.get_by_id(rxn["id"])
                    if rxns.id != "ADD_NADH":
                        rxns.bounds = check_model.reactions.get_by_id(rxn["id"]).bounds
            close_autotrophic_or_c_source(model_info, model, "carbon_source")
            set_c_source_supply(model_info, model, "nadhs", check_model)
            model_control_info["final_nadh_flux"] = model.slim_optimize()
            print(model.slim_optimize())
            print(model.optimize().objective_value)
            nadh_id = model_info["now_obj_rxn"]["rxn_id"]

            model.reactions.remove("ADD_NADH")
            model.objective = model_info["initial_rxn"]["biomass_rxn_id"]
            model_control_info["final_biomass_flux"] = model.slim_optimize()

    def nadh_control(
        self, model_info, model, check_model, model_control_info, model_check_info, obj
    ):
        """"""
        temp, model_info["modified_rxn_id"] = {}, {}

        model_control_info["check_reducing_equivalents_production"]["score"] = 1
        num = 1e-6
        nadh_name, temp = [
            "NADH",
            "NADPH",
            "FADH2",
            "FMNH2",
            "Q8H2",
            "MQL8",
            "DMMQL8",
        ], {}
        reas = [NADH_NAME, NADPH, FADH2, FMNH2, Q8H2_NAME, MQL8, DMMQL8]
        pro1s = [NAD_NAME, NADP, FAD, FMN, Q8_NAME, MQN8, DMMQ8]
        pro2 = H_NAME
        for i, name in enumerate(nadh_name):
            model_control_info["check_reducing_equivalents_production"][name] = {}
            nadh_id = get_other_nadh_rxn(model, reas[i], pro1s[i], pro2, name)
            if not nadh_id:
                model_control_info["check_reducing_equivalents_production"][name][
                    "score"
                ] = 2
                model_control_info["check_reducing_equivalents_production"][name][
                    "Details"
                ] = f"The reducing equivalent {name} is not present in the model."
                model_info["control_analysis_final"].append("")
                continue

            model.reactions.get_by_id(nadh_id).bounds = (0, 1000)
            set_model_objective(model_info, model, nadh_id)

            close_autotrophic_or_c_source(model_info, model, "carbon_source")

            if pd.isna(model.slim_optimize()) or model.slim_optimize() <= 1e-6:
                model_control_info["check_reducing_equivalents_production"][name][
                    "score"
                ] = 1
            else:
                model_control_info["check_reducing_equivalents_production"][name][
                    "score"
                ] = 0
                model_control_info["check_reducing_equivalents_production"]["score"] = 0

            grate_delete_rxn = self.close_max_value_nadh(
                model_info,
                model_control_info,
                model,
                check_model,
                obj,
                nadh_id,
                name,
                temp,
            )

            all_need_fluxes, infinite_rxn = add_back_reactions_net(
                model_info, model, check_model, grate_delete_rxn, "nadhs", num, self.cfg
            )
            print("infinite_rxn", infinite_rxn)
            temp["model_revision"] = get_net_infinite_info(model_info, model, "nadhs")
            model_revision_id = [rea["reaction_id"] for rea in temp["model_revision"]]
            temp["related_reaction"] = get_other_net_infinite_info(
                model_info, model, all_need_fluxes
            )
            with model:
                set_c_source_supply(model_info, model, "nadhs", check_model)
                set_auto_source_supply(model_info, model, check_model, "initial")
                set_c_source_supply(model_info, model, "nadhs", check_model)
                final1 = round(model.slim_optimize(), 3)
                finals = get_c_mol(model)
                print(finals, round(model.slim_optimize(), 3), ".............")

            model_info["control_analysis_final"].append(finals)
            if (
                model_control_info["check_reducing_equivalents_production"][name][
                    "score"
                ]
                == 0
            ):
                temp["score"] = 0
                model_control_info["check_reducing_equivalents_production"][name][
                    "Optimal_rate_in_initial_model"
                ] = f"Optimal {name} production rate was {model_info['control_analysis_initial'][i]} mmol/gDCW/h when no carbon sourse was supplied."
                model_control_info["check_reducing_equivalents_production"][name][
                    "Details"
                ] = f"Identifying erroneous reaction set {model_revision_id} was crucial in assressing the erroneous net generation of {name}.\n In the 'Model revision' section, the equation {model_revision_id} and the modifications made to it were presented. The pathway responsible for the erroneous net generation of {name} could be visualized by clicking on 'Visualisation'.\n In the 'Related reactions' section, reactions associated with the pathways contributing to the erroneous net generation of {name} were listed, and modifications to these reactions can be performed by clicking 'Edit'."
                model_control_info["check_reducing_equivalents_production"][name][
                    "model_revision"
                ] = temp["model_revision"]
                model_control_info["check_reducing_equivalents_production"][name][
                    "related_reaction"
                ] = temp["related_reaction"]
                model_control_info["check_reducing_equivalents_production"][name][
                    "Optimal_rate_in_revised_model"
                ] = f"Optimal {name} production rate was {final1} mmol/gDCW/h when carbon source in the list of “Carbon source supply ” was supplied."

            else:
                model_control_info["check_reducing_equivalents_production"][name][
                    "Summary"
                ] = f"Pass. The reducing equivalent {name} did not generate infinitely. "
                model_control_info["check_reducing_equivalents_production"][name][
                    "Details"
                ] = f"Reducing equivalent {name} is not generated when the carbon source supply is turned off, namely there is no erroneous net generation of reducing equivalents in the model."
                model_control_info["check_reducing_equivalents_production"][name][
                    "Optimal_rate"
                ] = f"Optimal {name} production rate was {final1} mmol/gDCW/h when carbon source set in the list of “Carbon source supply ” was supplied."

            if f"ADD_{name}" in model.reactions:
                model.remove_reactions(f"ADD_{name}")
