import re

from mqc.config import Config
from mqc.utils import *
from mqc.defaults import *
from mqc.control.biomass_control import Biomasses
from mqc.control.rules_control import Rules


class Nets:
    """
    Get net substance output information.

    """

    def __init__(self, cfg: Config):
        """
        Define net substance output.

        """
        self.cfg = cfg
        self.nets = []

    def find_net_generation(self, model_info, model, model_control_info):
        """
        Determining whether there is a net formation of substances
        """
        model_control_info["check_metabolite_production"]["score"] = 1
        r_met, net_generation, only_pro_mets, met_num, c_met_num = [], [], [], 0, 0
        initial_net_flux, temps = [], {}
        for rxn in model_info["reactions"]:
            if rxn["id"] in model_info["exchange_rxns"] and rxn["bounds"][0] < 0:
                r_met += rxn["all_mets"]
            if len(rxn["reactants_mets"]) == 0 and len(rxn["products_mets"]) != 0:
                only_pro_mets += rxn["products_mets"]
        for met in model.metabolites:
            met_num += 1
            if met.name not in (FREE_METS + only_pro_mets):
                if re.search("C[A-Z]", str(met.formula)) or re.search(
                    r"C[\d]", str(met.formula)
                ):
                    c_met_num += 1
                    with model:
                        if met.name in ACCOA and met.compartment in C_COMPARTMENT:
                            objectiveId = add_accoa(model, str(met.id))
                        else:
                            objectiveId = add_demand(model_info, model, met.id)
                        if not objectiveId:
                            continue
                        model.objective = objectiveId
                        if model.slim_optimize() > 1e-5 and met.name not in r_met:
                            net_generation.append(met.id)

                            temps[f"{met.id}"] = round(model.slim_optimize(), 3)
        model_info["net_generation"] = net_generation
        if len(net_generation) != 0:
            model_control_info["check_metabolite_production"]["score"] = 0
        model_info["total_number_of_metabolites"] = met_num
        model_info["total_number_of_carbon_containing_metabolites"] = c_met_num

        model_info["total_net_matter"] = len(net_generation)
        model_info["net_matter"] = ", ".join(net_generation)

        model_info["initial_flux"]["initial_net"] = temps

    def close_max_value_net(
        self,
        model_info,
        model_control_info,
        model,
        object_rxn_id,
        netId,
        check_model,
        obj,
        temp,
    ):
        """
        Turn off the maximum value in each result
        """
        (
            grate_delete_rxn,
            fluxs,
            objects,
            count2,
            now_rely_rxn,
            temp_rxnId,
            need_gap_rxn,
        ) = ([], 0, "", 0, {}, "", [])
        count = 0
        if model.slim_optimize() <= 1e-5:

            pass
        while model.slim_optimize() > 1e-5:
            count += 1
            print(count, " --- ", object_rxn_id)
            if objects == object_rxn_id and fluxs == round(model.slim_optimize(), 3):
                count2 += 1
            else:
                fluxs = round(model.slim_optimize(), 3)
                objects = object_rxn_id
                count2 = 0
            if count2 > 7:
                print(
                    "Unexpected production persists; check additional nutrient inputs"
                )
                write_flux_file(model_info, model, object_rxn_id, "net", self.cfg)
                model_control_info["check_metabolite_production"][
                    "Summary"
                ] = f"Calculating the optimal production rate of {temp['total_number_of_carbon_containing_metabolites']} carbon-containing metabolites, respectively when no carbon source was supplied. The erroneous net generation of {len(model_info['net_generation'])} metabolites is as follows:\n {model_info['net_matter']}"
                model_control_info["check_metabolite_production"][
                    "Details"
                ] = "error: There may be additional substrate inputs in the model"
                model_control_info["check_metabolite_production"]["model_revision"] = (
                    temp["model_revision"]
                )
                model_control_info["check_metabolite_production"][
                    "related_reaction"
                ] = temp["related_reaction"]

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
                raise RuntimeError("error: net infinite loop")

            fba_rxn_dic, flux_fic = {}, {}
            if model.slim_optimize() > 0:
                pfba_solution = pfba(model)
                need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes) > 1e-6]
            else:
                need_fluxes = {}

            rely_rxn = get_rely_rxn(
                model, model_info, check_model, object_rxn_id, need_fluxes
            )
            print("......................................")
            print("rely_rxn: ", rely_rxn)
            now_rely_rxn = {}
            for ids, v in need_fluxes.items():
                for rxn in model_info["reactions"]:
                    if ids == rxn["id"]:

                        if ids in rely_rxn:
                            now_rely_rxn.update({ids: rxn["net_penalty_points"]})
                        if ids == object_rxn_id or ids in rely_rxn:
                            continue
                        fba_rxn_dic[ids] = rxn["net_penalty_points"]
                        flux_fic[ids] = v

            if count2 == 3 and now_rely_rxn:

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
                    if ids != object_rxn_id and fba_rxn_dic[ids] != 0:
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
            print("close :", model.slim_optimize())
        return grate_delete_rxn

    def old_get_associated_substances(
        self,
        modelInfo,
        model,
        temp,
        infinite_rxn,
        all_net,
        check_model,
        model_control_info,
    ):
        """"""
        temp_net, modified_rxn = [], []

        for netId in all_net:
            with model:
                print(
                    "ABAH",
                    model.reactions.get_by_id("ABAH").bounds,
                    "...........................................",
                )
                print(
                    "EX_h_e",
                    model.reactions.get_by_id("EX_h_e").bounds,
                    "...........................................",
                )
                met = model.metabolites.get_by_id(netId)
                if met.name in ACCOA and met.compartment in C_COMPARTMENT:
                    object_rxn_id = add_accoa(model, str(met.id))
                else:
                    object_rxn_id = add_demand(modelInfo, model, netId)
                model.objective = object_rxn_id
                if model.slim_optimize() <= 1e-6:
                    temp_net.append(netId)

        temp["correlation_between_net_substances"][f"{','.join(infinite_rxn)}"] = (
            ",".join(temp_net)
        )
        return temp_net

    def get_associated_substances(self, modelInfo, model, temp, all_net, check_model):
        """"""
        temp_net = []
        for rxnId in modelInfo["modified_rxn_id"].keys():
            model.reactions.get_by_id(rxnId).bounds = check_model.reactions.get_by_id(
                rxnId
            ).bounds
        for rxnId, bounds in modelInfo["modified_rxn_id"].items():
            model.reactions.get_by_id(rxnId).bounds = bounds
            for netId in all_net:
                with model:
                    met = model.metabolites.get_by_id(netId)
                    if met.name in ACCOA and met.compartment in C_COMPARTMENT:
                        object_rxn_id = add_accoa(model, str(met.id))
                    else:
                        object_rxn_id = add_demand(modelInfo, model, netId)
                    if not object_rxn_id:
                        continue
                    model.objective = object_rxn_id
                    if model.slim_optimize() <= 1e-6:
                        temp_net.append(netId)

            if temp_net:
                temp["correlation_between_net_substances"][rxnId] = ",".join(temp_net)

                temp_net = []

    def net_control(
        self, model_info, model, check_model, model_control_info, model_check_info, obj
    ):
        """
        nets correction process
        """
        all_net, num, temp, all_need_fluxes = [], 1e-6, {}, {}
        close_autotrophic_or_c_source(model_info, model, "carbon_source")

        self.find_net_generation(model_info, model, model_control_info)
        print(model_info["total_number_of_metabolites"])
        temp["total_number_of_metabolites"] = model_info["total_number_of_metabolites"]
        temp["total_number_of_carbon_containing_metabolites"] = model_info[
            "total_number_of_carbon_containing_metabolites"
        ]
        temp["total_net_matter"] = model_info["total_net_matter"]
        temp["net_matter"] = model_info["net_matter"]
        all_net.extend(model_info["net_generation"])
        temp["correlation_between_net_substances"] = {}
        temp["substance_approach"], temp["model_revision"], temp["related_reaction"] = (
            [],
            [],
            [],
        )
        for netId in model_info["net_generation"]:
            with model:
                met = model.metabolites.get_by_id(netId)
                if met.name in ACCOA and met.compartment in C_COMPARTMENT:
                    object_rxn_id = add_accoa(model, str(met.id))
                else:
                    object_rxn_id = add_demand(model_info, model, netId)
                if not object_rxn_id:
                    continue
                model.objective = object_rxn_id
                grate_delete_rxn = self.close_max_value_net(
                    model_info,
                    model_control_info,
                    model,
                    object_rxn_id,
                    netId,
                    check_model,
                    obj,
                    temp,
                )
                if len(grate_delete_rxn) == 0:
                    continue
                need_flux, infinite_rxn = add_back_reactions_net(
                    model_info,
                    model,
                    check_model,
                    grate_delete_rxn,
                    "nets",
                    num,
                    self.cfg,
                )
                all_need_fluxes.update(need_flux)
                temp_model_revision = get_net_infinite_info(model_info, model, "nets")
                temp["model_revision"].extend(temp_model_revision)
                model_revision_id = [
                    rea["reaction_id"] for rea in temp["model_revision"]
                ]
                temp_related_reaction = get_other_net_infinite_info(
                    model_info, model, all_need_fluxes
                )
                temp["related_reaction"].extend(temp_related_reaction)

            boundary_restoration(model_info, model, "nets")

        self.get_associated_substances(model_info, model, temp, all_net, check_model)

        if model_control_info["check_metabolite_production"]["score"] == 0:

            model_control_info["check_metabolite_production"][
                "Summary"
            ] = f"Calculating the optimal production rate of {temp['total_number_of_carbon_containing_metabolites']} carbon-containing metabolites, respectively when no carbon source was supplied. The erroneous net generation of {len(model_info['net_generation'])} metabolites is as follows:\n {model_info['net_matter']}"
            model_control_info["check_metabolite_production"][
                "Details"
            ] = f"Identifying erroneous reaction set {model_revision_id} was crucial in assressing the erroneous net generation of these metabolites.\n In the 'Model revision' section, the equation set {model_revision_id} and the modifications made to it were presented. The pathway responsible for the erroneous net generation of these metabiltes could be visualized by clicking on 'Visualisation'.\nIn the 'Related reactions' section, reactions associated with the pathways contributing to the erroneous net generation of metabolites were listed, and modifications to these reactions can be performed by clicking 'Edit'."
            model_control_info["check_metabolite_production"]["model_revision"] = temp[
                "model_revision"
            ]
            model_control_info["check_metabolite_production"]["related_reaction"] = (
                temp["related_reaction"]
            )
            model_info["control_analysis"].append("error")

        else:
            model_control_info["check_metabolite_production"][
                "Summary"
            ] = f"There are a total of {temp['total_number_of_carbon_containing_metabolites']} carbon-containing metabolites, and none of them can product without carbon substrate supply. Therefore, this model does not exhibit errors in net production of metabolites."
            model_info["control_analysis"].append("yes")
