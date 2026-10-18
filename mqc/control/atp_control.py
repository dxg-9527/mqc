from mqc.config import Config
from mqc.utils import *
from mqc.control.rules_control import Rules
from mqc.control.rules_control2 import Rules2
from mqc.control.biomass_control import Biomasses


class Atps:

    def __init__(self, cfg: Config):
        """"""
        self.cfg = cfg

    def modify_nadh(self, model_info, model, need_fluxes, grate_delete_rxn):
        """"""
        if model.slim_optimize() > 0:
            pfba_solution = pfba(model)
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-6]
        else:
            need_fluxes = {}
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
                            print(
                                ids,
                                " ",
                                " ",
                                v,
                                " ",
                                model.reactions.get_by_id(ids).build_reaction_string(
                                    use_metabolite_names=True
                                ),
                            )
                            rxn["atps_modify"] = "true"
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
                            print(
                                ids,
                                " ",
                                " ",
                                v,
                                " ",
                                model.reactions.get_by_id(ids).build_reaction_string(
                                    use_metabolite_names=True
                                ),
                            )
                            rxn["atps_modify"] = "true"
                            return

    def get_exchange(self, model_info, need_fluxes):
        for ids, v in need_fluxes.items():
            for rxn in model_info["reactions"]:
                if ids == rxn["id"] and rxn["id"] in model_info["exchange_rxns"]:
                    if (
                        len(set(H_NAME) & set(rxn["reactants_mets"])) != 0 and v < 0
                    ) or (len(set(H_NAME) & set(rxn["products_mets"])) != 0 and v > 0):
                        return ids, v

    def get_transport(self, model_info, model, need_fluxes, h_close):
        for ids, v in need_fluxes.items():
            for rxn in model_info["reactions"]:
                if ids == rxn["id"] and rxn["id"] in model_info["transport_rxns"]:
                    rxns = model.reactions.get_by_id(ids)
                    reactants_mets = [m.id for m in rxns.reactants]
                    products_mets = [m.id for m in rxns.products]
                    print("transport:", h_close)
                    if (len(set(h_close) & set(reactants_mets)) != 0 and v < 0) or (
                        len(set(h_close) & set(products_mets)) != 0 and v > 0
                    ):
                        return ids, v

    def get_other(self, model_info, model, need_fluxes, h_close):
        for ids, v in need_fluxes.items():
            for rxn in model_info["reactions"]:
                if (
                    ids == rxn["id"]
                    and rxn["id"] not in model_info["exchange_rxns"]
                    and rxn["id"] not in model_info["transport_rxns"]
                ):
                    if "rule_conflict" in rxn.keys():
                        print(
                            "Proton adjustment conflicts with reaction rules; retaining reaction:",
                            ids,
                            rxn["bounds"],
                        )
                        break
                    rxns = model.reactions.get_by_id(ids)
                    reactants_mets = [m.id for m in rxns.reactants]
                    products_mets = [m.id for m in rxns.products]
                    if (len(set(h_close) & set(reactants_mets)) != 0 and v < 0) or (
                        len(set(h_close) & set(products_mets)) != 0 and v > 0
                    ):
                        return ids, v

    def set_bounds(self, ids, v, grate_delete_rxn, model_info, model):
        for rxn in model_info["reactions"]:
            if ids == rxn["id"]:
                rxns = model.reactions.get_by_id(ids)
                print(
                    ids,
                    " ",
                    " ",
                    v,
                    " ",
                    model.reactions.get_by_id(ids).build_reaction_string(
                        use_metabolite_names=True
                    ),
                )
                if v < 0:
                    rxns.bounds = (0, 1000)
                    rxn["bounds"] = [0, 1000]
                    rxn["rules"]["modify_pmf_rxn"] = "true"
                    rxn["atps_modify"] = "true"
                    grate_delete_rxn.append(ids)
                if v > 0:
                    rxns.bounds = (-1000, 0)
                    rxn["bounds"] = [-1000, 0]
                    rxn["rules"]["modify_pmf_rxn"] = "true"
                    rxn["atps_modify"] = "true"
                    grate_delete_rxn.append(ids)

    def set_bounds2(self, ids, v, grate_delete_rxn, model_info, model):
        for rxn in model_info["reactions"]:
            if ids == rxn["id"]:
                rxns = model.reactions.get_by_id(ids)
                rxn_bounds = rxns.bounds
                if v < 0:
                    rxns.bounds = (0, 1000)
                    rxn["bounds"] = [0, 1000]
                    rules2 = Rules2(ids, self.cfg, v)
                    rules2.get_all_rules(model_info, model)
                    if "rule_conflict" in rxn.keys():
                        rxns.bounds = rxn_bounds
                        rxn["bounds"] = rxn_bounds
                    else:
                        rxn["rules"]["modify_pmf_rxn"] = "true"
                        grate_delete_rxn.append(ids)
                        print(
                            ids,
                            " ",
                            " ",
                            v,
                            " ",
                            model.reactions.get_by_id(ids).build_reaction_string(
                                use_metabolite_names=True
                            ),
                        )
                        rxn["atps_modify"] = "true"
                if v > 0:
                    rxns.bounds = (-1000, 0)
                    rxn["bounds"] = [-1000, 0]
                    rules2 = Rules2(ids, self.cfg, v)
                    rules2.get_all_rules(model_info, model)
                    if "rule_conflict" in rxn.keys():
                        rxns.bounds = rxn_bounds
                        rxn["bounds"] = rxn_bounds
                    else:
                        rxn["rules"]["modify_pmf_rxn"] = "true"
                        grate_delete_rxn.append(ids)
                        print(
                            ids,
                            " ",
                            " ",
                            v,
                            " ",
                            model.reactions.get_by_id(ids).build_reaction_string(
                                use_metabolite_names=True
                            ),
                        )
                        rxn["atps_modify"] = "true"

    def modify_pmf(self, model_info, model, need_fluxes, grate_delete_rxn, h_close):
        """"""
        if model.slim_optimize() > 0:
            pfba_solution = pfba(model)
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-6]
        else:
            need_fluxes = {}
        print("h_close:", h_close)
        print("Adjusting extracellular proton flux")
        pmf_info = self.get_exchange(model_info, need_fluxes)
        if pmf_info:
            print("pmf_info", pmf_info)
            exchange_rxn_id, v = pmf_info
            self.set_bounds(exchange_rxn_id, v, grate_delete_rxn, model_info, model)
            return
        pmf_info = self.get_transport(model_info, model, need_fluxes, h_close)
        if pmf_info:
            transport_rxn_id, v = pmf_info
            self.set_bounds(transport_rxn_id, v, grate_delete_rxn, model_info, model)
            return
        pmf_info = self.get_other(model_info, model, need_fluxes, h_close)
        if pmf_info:
            other_rxn_id, v = pmf_info
            self.set_bounds2(other_rxn_id, v, grate_delete_rxn, model_info, model)
            print(
                "Examining proton flux in an otherwise valid reaction", other_rxn_id, v
            )
            return

    def modify_atp(self, model_info, model, need_fluxes, grate_delete_rxn):
        """"""
        if model.slim_optimize() > 0:
            pfba_solution = pfba(model)
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-6]
        else:
            need_fluxes = {}
        for ids, v in need_fluxes.items():
            for rxn in model_info["reactions"]:
                if ids == rxn["id"]:
                    if (
                        len(set(ATP_NAME) & set(rxn["reactants_mets"])) != 0
                        and len(set(ATP_NAME) & set(rxn["products_mets"])) != 0
                    ):
                        break
                    rxns = model.reactions.get_by_id(ids)
                    rxn_bounds = rxns.bounds
                    if (
                        len(set(ATP_NAME) & set(rxn["reactants_mets"])) != 0
                        and v < 0
                        and "right_chain_rxn" not in rxn["rules"].keys()
                    ):
                        rxns.bounds = (0, 1000)
                        rxn["bounds"] = [0, 1000]
                        rules2 = Rules2(ids, self.cfg, v)
                        rules2.get_all_rules(model_info, model)
                        if "rule_conflict" in rxn.keys():
                            rxns.bounds = rxn_bounds
                            rxn["bounds"] = rxn_bounds
                        else:
                            rxn["rules"]["modify_atp_rxn"] = "true"
                            grate_delete_rxn.append(ids)
                            print(
                                ids,
                                " ",
                                v,
                                " ",
                                model.reactions.get_by_id(ids).build_reaction_string(
                                    use_metabolite_names=True
                                ),
                            )
                            rxn["atps_modify"] = "true"
                            return
                    if (
                        len(set(ATP_NAME) & set(rxn["products_mets"])) != 0
                        and v > 0
                        and "right_chain_rxn" not in rxn["rules"].keys()
                    ):
                        rxns.bounds = (-1000, 0)
                        rxn["bounds"] = [-1000, 0]
                        rules2 = Rules2(ids, self.cfg, v)
                        rules2.get_all_rules(model_info, model)
                        if "rule_conflict" in rxn.keys():
                            rxns.bounds = rxn_bounds
                            rxn["bounds"] = rxn_bounds
                        else:
                            rxn["rules"]["modify_atp_rxn"] = "true"
                            grate_delete_rxn.append(ids)
                            print(
                                ids,
                                " ",
                                " ",
                                v,
                                " ",
                                model.reactions.get_by_id(ids).build_reaction_string(
                                    use_metabolite_names=True
                                ),
                            )
                            rxn["atps_modify"] = "true"
                            return

    def close_max_value_atp(
        self,
        model_info,
        model_control_info,
        model,
        h_close,
        check_model,
        atpm_id,
        obj,
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
            temp_rxnId,
            need_gap_rxn,
        ) = ([], 0, 0, 0, {}, "", [])
        if model.slim_optimize() <= 1e-6:
            print(model, "________ATP_yes")
        while model.slim_optimize() > 1e-6:
            count = count + 1
            print(count, count2, atpm_id, sep="................")
            if fluxs != round(model.slim_optimize(), 3):
                fluxs = round(model.slim_optimize(), 3)
                count2 = 0
            else:
                count2 += 1
            if count2 > 35:
                print(
                    "Energy correction did not converge; check additional nutrient inputs or erroneous reactions"
                )
                write_flux_file(model_info, model, atpm_id, "atp", self.cfg)
                model_control_info["check_energy_production"][name][
                    "Optimal_rate_in_initial_model"
                ] = ""
                model_control_info["check_energy_production"][name][
                    "Details"
                ] = "error: Infinite energy cycle, there may be additional substrate inputs or erroneous reactions in the model"
                model_control_info["check_energy_production"][name][
                    "model_revision"
                ] = []
                model_control_info["check_energy_production"][name][
                    "related_reaction"
                ] = []
                model_control_info["check_energy_production"][name][
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
                raise RuntimeError("error: atp infinite loop")

            if count2 > 3:
                self.modify_nadh(model_info, model, need_fluxes, grate_delete_rxn)
            if count2 > 6:
                self.modify_pmf(
                    model_info, model, need_fluxes, grate_delete_rxn, h_close
                )
            if count2 > 9:
                self.modify_atp(model_info, model, need_fluxes, grate_delete_rxn)
            fba_rxn_dic, flux_fic = {}, {}

            if model.slim_optimize() > 0:
                pfba_solution = pfba(model)
                need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-6]
            else:
                need_fluxes = {}
            rely_rxn = get_rely_rxn(
                model, model_info, check_model, atpm_id, need_fluxes
            )
            print("......................................")
            print("rely_rxn: ", rely_rxn)
            now_rely_rxn = {}
            for ids, v in need_fluxes.items():
                for rxn in model_info["reactions"]:
                    if ids == rxn["id"]:

                        if ids in rely_rxn:
                            now_rely_rxn.update({ids: rxn["net_penalty_points"]})
                        if ids == atpm_id or ids in rely_rxn:
                            continue
                        fba_rxn_dic[ids] = rxn["net_penalty_points"]
                        flux_fic[ids] = v

            if count2 == 12 and now_rely_rxn:

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
                        print(model.slim_optimize(), "ATP production rate:")
                        break

                with model:
                    set_c_source_supply(model_info, model, "atps", check_model)
                    set_auto_source_supply(model_info, model, check_model, "atps")
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
                    if ids != atpm_id and ids not in ATPM and fba_rxn_dic[ids] != 0:
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
        model_info["temp_final_flux"] = {}

        with model:
            set_c_source_supply(model_info, model, "atps", check_model)
            set_auto_source_supply(model_info, model, check_model, "atps")
            model_info["temp_final_flux"]["final_atp"] = model.slim_optimize()

            nadh_id = model_info["initial_flux"]["initial_nadh"]["rxn_id"]
            set_model_objective(model_info, model, nadh_id)
            model_info["temp_final_flux"]["final_nadh"] = model.slim_optimize()

            if "ADD_ATPM" in model.reactions:
                model.reactions.remove("ADD_ATPM")
            model.reactions.remove("ADD_NADH")
            model.objective = model_info["initial_rxn"]["biomass_rxn_id"]
            model_info["temp_final_flux"]["final_biomass"] = model.slim_optimize()

    def atp_control2(
        self, model_info, model, check_model, model_control_info, model_check_info, obj
    ):
        """"""
        xtp = ["A", "C", "G", "U", "I"]
        model_control_info["check_energy_production"]["score"] = 1
        num, temp = 1e-6, {}
        atpm_id = get_atpm_rxn(model_info, model, "ATP", "ADP")
        set_model_objective(model_info, model, atpm_id)
        close_autotrophic_or_c_source(model_info, model, "carbon_source")

        if model.slim_optimize() <= 1e-6:
            model_control_info["check_energy_production"][
                "special_correction:"
            ] = "Energy is not generated  when the carbon source supply is turned off, namely there is no erroneous net generation of energy in the model"
        else:
            model_control_info["check_energy_production"]["score"] = 0

        h_close = model_info["h_close"]
        grate_delete_rxn = self.close_max_value_atp(
            model_info, model_control_info, model, h_close, check_model, obj
        )
        all_need_fluxes, infinite_rxn = add_back_reactions_net(
            model_info, model, check_model, grate_delete_rxn, "atps", num, self.cfg
        )
        temp["model_revision"] = get_net_infinite_info(model_info, model, "atps")
        temp["related_reaction"] = get_other_net_infinite_info(
            model_info, model, all_need_fluxes
        )

        if model_control_info["check_energy_production"]["score"] == 0:
            temp["score"] = 0
            model_control_info["check_energy_production"][
                "description1"
            ] = "energy correction response"
            model_control_info["check_energy_production"]["model_revision"] = temp[
                "model_revision"
            ]
            print(model_control_info["check_energy_production"])
            model_control_info["check_energy_production"][
                "description2"
            ] = "Energy-related pathways available for revision"
            model_control_info["check_energy_production"]["related_reaction"] = temp[
                "related_reaction"
            ]
            for k in xtp:
                model_check_info["check_energy_production"]["initial_flux"] = (
                    model_control_info["initial_flux"][f"initial_{k}tp"]
                )
            model_check_info["check_energy_production"]["Summary"] = "energy not passed"
        else:
            model_control_info["check_energy_production"][
                "Summary"
            ] = "energy quality control passed"
            model_check_info["check_energy_production"][
                "Summary"
            ] = "energy quality control passed"
        if "ADD_ATPM" in model.reactions:
            model.reactions.remove("ADD_ATPM")

    def atp_control(
        self, model_info, model, check_model, model_control_info, model_check_info, obj
    ):
        """"""
        xtp = ["A", "C", "G", "U", "I"]
        model_control_info["check_energy_production"]["score"] = 1
        num, temp, ini_num = 1e-6, {}, 0
        for k in xtp:
            ini_num += 1
            name = f"{k}TP"
            model_control_info["check_energy_production"][name] = {}
            model_check_info["check_energy_production"][name] = {}
            atpm_id = get_atpm_rxn(model_info, model, name, f"{k}DP")
            if not atpm_id:
                model_control_info["check_energy_production"][name]["score"] = 2
                model_control_info["check_energy_production"][name][
                    "Details"
                ] = f"The energy {k}TP is not present in the model."
                model_info["control_analysis_final"].append("")
                continue
            model.reactions.get_by_id(atpm_id).bounds = (0, 1000)
            set_model_objective(model_info, model, atpm_id)

            close_autotrophic_or_c_source(model_info, model, "carbon_source")

            if pd.isna(model.slim_optimize()) or model.slim_optimize() <= 1e-6:
                model_control_info["check_energy_production"][name]["score"] = 1

            else:
                model_control_info["check_energy_production"][name]["score"] = 0
                model_control_info["check_energy_production"]["score"] = 0

            h_close = model_info["h_close"]
            grate_delete_rxn = self.close_max_value_atp(
                model_info,
                model_control_info,
                model,
                h_close,
                check_model,
                atpm_id,
                obj,
                name,
                temp,
            )
            all_need_fluxes, infinite_rxn = add_back_reactions_net(
                model_info, model, check_model, grate_delete_rxn, "atps", num, self.cfg
            )
            temp["model_revision"] = get_net_infinite_info(model_info, model, "atps")
            model_revision_id = [rea["reaction_id"] for rea in temp["model_revision"]]
            temp["related_reaction"] = get_other_net_infinite_info(
                model_info, model, all_need_fluxes
            )

            with model:
                set_c_source_supply(model_info, model, "initial", check_model)
                set_auto_source_supply(model_info, model, check_model, "initial")
                set_c_source_supply(model_info, model, "initial", check_model)
                final1 = round(model.slim_optimize(), 3)
                finals = get_c_mol(model)

            model_info["control_analysis_final"].append(finals)
            if model_control_info["check_energy_production"][name]["score"] == 0:
                temp["score"] = 0
                model_control_info["check_energy_production"][name][
                    "Optimal_rate_in_initial_model"
                ] = f"Optimal {k}TP production rate was {model_info['control_analysis_initial'][6+ini_num]} mmol/gDCW/h when no carbon sourse was supplied."
                model_control_info["check_energy_production"][name][
                    "Details"
                ] = f"Identifying erroneous reaction set {model_revision_id} was crucial in assressing the erroneous net generation of {k}TP.\n In the 'Model revision' section, the equation {model_revision_id} and the modifications made to it were presented. The pathway responsible for the erroneous net generation of {k}TP could be visualized by clicking on 'Visualisation'.\n In the 'Related reactions' section, reactions associated with the pathways contributing to the erroneous net generation of {k}TP were listed, and modifications to these reactions can be performed by clicking 'Edit'."
                model_control_info["check_energy_production"][name][
                    "model_revision"
                ] = temp["model_revision"]
                model_control_info["check_energy_production"][name][
                    "related_reaction"
                ] = temp["related_reaction"]
                model_control_info["check_energy_production"][name][
                    "Optimal_rate_in_revised_model"
                ] = f"Optimal {k}TP production rate was {final1} mmol/gDCW/h when carbon source in the list of “Carbon source supply ” was supplied."

            else:
                model_control_info["check_energy_production"][name][
                    "Summary"
                ] = f"Pass. The energy {k}TP did not generate infinitely. "
                model_control_info["check_energy_production"][name][
                    "Details"
                ] = f"Energy {k}TP is not generated when the carbon source supply is turned off, namely there is no erroneous net generation of energy in the model."
                model_control_info["check_energy_production"][name][
                    "Optimal_rate"
                ] = f"Optimal {k}TP production rate was {final1} mmol/gDCW/h when carbon source set in the list of “Carbon source supply ” was supplied."

            if f"ADD_{k}TPM" in model.reactions:
                model.remove_reactions([f"ADD_{k}TPM"])
        control_analysis_initial = model_info["control_analysis_initial"]
        control_analysis_final = model_info["control_analysis_final"]
        for i in range(len(control_analysis_initial)):
            model_info["control_analysis"].append(
                (control_analysis_initial[i], control_analysis_final[i])
            )
