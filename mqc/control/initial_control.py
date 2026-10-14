from cobra.util.solver import linear_reaction_coefficients
import json
import multiprocessing
from multiprocessing import Pipe, Process
import time

from mqc.config import Config
from mqc.control.rules_control import Rules
from mqc.utils import *


class InitialPreprocess:
    """ """

    def __init__(self, cfg: Config) -> None:
        """ """
        self.cfg = cfg

    def get_exchange(self, model, model_info, check_model):
        """"""
        temp, exchange_info = [], {}

        for rxn in model.reactions:
            all_mets = [m.id for m in rxn.metabolites]
            if len(all_mets) == 1:
                check_rxns = check_model.reactions.get_by_id(rxn.id)
                mets = model.metabolites.get_by_id(all_mets[0])
                exchange_info = {
                    "Reaction_id": rxn.id,
                    "Equation(ID)": rxn.build_reaction_string(
                        use_metabolite_names=True
                    ),
                    "Formula": mets.formula,
                    "Original": str(check_rxns.bounds),
                    "L_b": rxn.lower_bound,
                    "U_b": rxn.upper_bound,
                }
                temp.append(exchange_info)
        return temp

    def find_net_initial_flux(self, model_info, model, model_control_info):
        """
        Determining whether there is a net formation of substances
        """
        model_control_info["check_metabolite_production"]["score"] = 1
        r_met, net_generation, only_pro_mets, met_num = [], [], [], 0
        initial_net_flux, temp = [], {}
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
                            initial_net_flux.extend(
                                [f"{met.id} : {round(model.slim_optimize(),3)}"]
                            )
                            temp[f"{met.id}"] = round(model.slim_optimize(), 3)
        model_info["net_generation"] = net_generation
        if len(net_generation) != 0:
            model_control_info["check_metabolite_production"]["score"] = 0
        model_info["total_number_of_metabolites:"] = met_num
        model_info["total_net_matter:"] = len(net_generation)
        model_info["net_matter:"] = ", ".join(net_generation)

    def get_reducing_power_initial_flux(self, model_info, model, model_control_info):
        """"""
        nadh_name, nadh_flux = [
            "NADH",
            "NADPH",
            "FADH2",
            "FMNH2",
            "Q8H2",
            "MQL8",
            "DMMQL8",
        ], []
        reas = [NADH_NAME, NADPH, FADH2, FMNH2, Q8H2_NAME, MQL8, DMMQL8]
        pro1s = [NAD_NAME, NADP, FAD, FMN, Q8_NAME, MQN8, DMMQ8]
        pro2 = H_NAME
        for i, name in enumerate(nadh_name):
            nadh_id = get_other_nadh_rxn(model, reas[i], pro1s[i], pro2, name)
            if not nadh_id:
                continue
            model.reactions.get_by_id(nadh_id).bounds = (0, 1000)
            set_model_objective(model_info, model, nadh_id)
            initial_nadh = f"initial_{name}".lower()
            model_info["initial_flux"][initial_nadh] = model_info["now_obj_rxn"]

            nadh_flux.extend(
                [f"{name} : {round(model_info['now_obj_rxn']['rxn_flux'],3)}"]
            )

    def get_energy_initial_flux(self, model_info, model, model_control_info):
        """"""
        xtp, xtp_flux = ["A", "C", "G", "U", "I"], []
        for k in xtp:
            atpm_id = get_atpm_rxn(model_info, model, f"{k}TP", f"{k}DP")
            if not atpm_id:
                continue
            model.reactions.get_by_id(atpm_id).bounds = (0, 1000)
            set_model_objective(model_info, model, atpm_id)
            model_info["initial_flux"][f"initial_{k}tp"] = model_info["now_obj_rxn"]

            xtp_flux.extend(
                [f"{k}TP : {round(model_info['now_obj_rxn']['rxn_flux'],3)}"]
            )

    def check_yield_formula(self, model, metId):
        """"""
        calculate_met = ["C", "H", "O", "N", "P", "S"]
        met = model.metabolites.get_by_id(metId)
        if met.elements and (not set(met.elements.keys()).issubset(set(calculate_met))):
            return 0
        if (
            pd.isna(met.formula)
            or "R" in met.formula
            or "X" in met.formula
            or met.formula == "nan"
            or not met.formula
            or met.formula == "null"
        ):
            return 0
        return 1

    def find_yield_generation_subprocess(self, model_info, model, metabolites, conn):
        """Get all carbonaceous species that exceed the maximum theoretical yield"""
        yield_generation, c_num = [], 0
        initial_yield_flux, temp, substance_flux = [], {}, []
        temps = []

        pid = multiprocessing.current_process().pid

        for met in metabolites:
            if met.id in model_info["need_set_carbon_mets"]:
                continue
            if re.search("C[A-Z]", str(met.formula)) or re.search(
                r"C[\d]", str(met.formula)
            ):
                c_num += 1
                if self.check_yield_formula(model, met.id) == 0:

                    continue
                try:
                    with model:
                        if met.name in ACCOA and met.compartment in C_COMPARTMENT:
                            objectiveId = add_accoa(model, str(met.id))
                        else:
                            objectiveId = add_demand(model_info, model, met.id)
                        model.objective = objectiveId

                        max_rate, temp = max_reduced_degree_rate(
                            model, met.id, model_info
                        )
                        model_info["yield_note"].extend(temp)
                        if (
                            max_rate == 0 or "R" in met.formula or "X" in met.formula
                        ):  # Skip zero product limits and formulas containing elements outside CHONPS.
                            continue

                        model_opt = model.slim_optimize()

                        if round(model_opt, 2) > round(max_rate, 2):

                            yield_generation.append(met.id)
                            initial_yield_flux.extend(
                                [f"{met.id} : {round(model_opt,3)}"]
                            )
                            temp[f"{met.id}"] = round(model_opt, 3)
                            approach_info = get_net_approach(model, met.id, max_rate)
                            temps.append(approach_info)

                except:
                    pass

        conn.send(
            {
                "yield_generation": yield_generation,
                "c_num": c_num,
                "initial_yield_flux": initial_yield_flux,
                "met_flux": temp,
                "initial_substance_approach": temps,
            }
        )
        conn.close()

    def find_yield_initial_flux(self, model_info, model, model_control_info):
        """ """
        print("Running in parallel")

        metabolites = model.metabolites
        print("metabolites num: ", len(metabolites))
        metabolites_num = len(metabolites)
        PROCESS_NUM = 10

        average = metabolites_num // 10
        remains = metabolites_num % 10

        processes = []

        parent_connections = []

        for i in range(PROCESS_NUM):
            if i == PROCESS_NUM - 1:
                tmp_metabolites = metabolites[i * average :]
            else:
                tmp_metabolites = metabolites[i * average : (i + 1) * average]

            parent_conn, child_conn = Pipe()
            parent_connections.append(parent_conn)

            process = Process(
                target=self.find_yield_generation_subprocess,
                args=(
                    model_info,
                    model,
                    tmp_metabolites,
                    child_conn,
                ),
            )
            processes.append(process)

        for process in processes:
            process.start()

        for process in processes:
            process.join()

        yield_generation = []
        c_num = 0
        initial_yield_flux = []
        temp_met_flux = {}
        initial_substance_approach = []
        substance_flux = []
        for parent_connection in parent_connections:
            response = parent_connection.recv()
            yield_generation += response["yield_generation"]
            c_num += response["c_num"]
            initial_yield_flux += response["initial_yield_flux"]
            temp_met_flux.update(response["met_flux"])
            initial_substance_approach += response["initial_substance_approach"]

        model_info["yield_generation"] = yield_generation

        if len(yield_generation) != 0:
            model_control_info["check_metabolite_yield"]["score"] = 0
        model_info["total_carbon_content"] = c_num
        model_info["yield_exceeds_the_theoretical_maximum_number"] = len(
            yield_generation
        )
        model_info["yield_net_matter"] = ", ".join(yield_generation)

        model_control_info["check_metabolite_yield"][
            "description0"
        ] = "Initial information on all substances exceeding the maximum theoretical yield"
        model_control_info["check_metabolite_yield"][
            "initial_substance_approach"
        ] = initial_substance_approach

    def get_model_initial_flux(
        self, model_info, model, model_control_info, model_check_info
    ):
        """"""
        with model:
            model_info["initial_flux"] = {}
            model_info["initial_flux"]["initial_biomass"] = model_info["initial_rxn"]

            self.get_reducing_power_initial_flux(model_info, model, model_control_info)
            self.get_energy_initial_flux(model_info, model, model_control_info)

            model_control_info["check_metabolite_yield"]["score"] = 1

            if is_autotrophic(model_info):
                model_control_info["check_metabolite_yield"][
                    "autotrophic"
                ] = "This is an autotrophic strain, no need to calculate the yield"
                model_check_info["check_metabolite_yield"][
                    "autotrophic"
                ] = "This is an autotrophic strain, no need to calculate the yield"

                return 0
            else:
                get_c_source_supply(model_info, model)
                self.find_yield_initial_flux(model_info, model, model_control_info)

    def get_initial_nadh_atp(self, model_info, model):
        """"""
        nadh_name = ["NADH", "NADPH", "FADH2", "FMNH2", "Q8H2", "MQL8", "DMMQL8"]
        reas = [NADH_NAME, NADPH, FADH2, FMNH2, Q8H2_NAME, MQL8, DMMQL8]
        pro1s = [NAD_NAME, NADP, FAD, FMN, Q8_NAME, MQN8, DMMQ8]
        pro2 = H_NAME
        for i, name in enumerate(nadh_name):
            with model:
                nadh_id = get_other_nadh_rxn(model, reas[i], pro1s[i], pro2, name)
                if not nadh_id:
                    model_info["control_analysis_initial"].append("")
                    continue
                model.reactions.get_by_id(nadh_id).bounds = (0, 1000)
                set_model_objective(model_info, model, nadh_id)

                print(model.slim_optimize(), "................")

                initial = round(model.slim_optimize(), 3)

                model_info["control_analysis_initial"].append(initial)

        xtp = ["A", "C", "G", "U", "I"]
        for k in xtp:
            name = f"{k}TP"
            with model:
                atpm_id = get_atpm_rxn(model_info, model, name, f"{k}DP")
                if not atpm_id:
                    model_info["control_analysis_initial"].append("")
                    continue
                model.reactions.get_by_id(atpm_id).bounds = (0, 1000)
                set_model_objective(model_info, model, atpm_id)
                initial2 = round(model.slim_optimize(), 3)

                model_info["control_analysis_initial"].append(initial2)

    def initial_control(
        self, model_info, model, check_model, model_control_info, model_check_info
    ):
        """"""

        model_info["control_analysis"].append(model.id)

        initial_rxn_id = model_info["initial_rxn"]["biomass_rxn_id"]
        has_biomass = bool(initial_rxn_id)
        if has_biomass:
            model.objective = initial_rxn_id
            write_flux_file2(
                model_info, model, initial_rxn_id, "ini_biomasss", self.cfg
            )
            get_biomass_adp_coefficient(model, model_info)

        del_model_boundary(model_info, model, check_model, model_check_info)

        model_check_info["boundary_information"]["biomass_boundary"] = model_info[
            "biomass_boundary"
        ]
        model_control_info["boundary_information"]["carbon_source_supply"] = ", ".join(
            model_info["carbon_source_supply"]
        )
        if has_biomass:
            print(
                "Initial biomass output:",
                model.slim_optimize(),
                "........................................111.......................",
            )

        get_carbon_proportion(model_info, model, check_model, model_control_info)
        set_c_source_supply(model_info, model, "initial", check_model)

        model_info["biomass_by_c_source"] = (
            round(model.slim_optimize(), 3) if has_biomass else 0
        )
        print("Carbon source:", model_info["carbon_source_boundary"])
        if has_biomass:
            print(
                "Biomass output after carbon-source configuration:",
                model.slim_optimize(),
                ".......................................222........................",
            )
            write_flux_file2(
                model_info, model, initial_rxn_id, "set_carbon_biomasss", self.cfg
            )

        self.get_initial_nadh_atp(model_info, model)

        model_info["initial_flux"] = {}
        set_special_boundary_rxn(model_info, model, check_model)

        set_bio_boundary(model_info, model)
        set_no2_no3_bounds(model_info, model)
        if has_biomass:
            non_carbon_limit = get_not_carbon_source_limiting_metabolites(
                model_info, model, check_model
            )
        else:
            non_carbon_limit = (
                "Not evaluated because no biomass reaction was identified."
            )
        model_control_info["boundary_information"][
            "non_carbon_source_limiting_substrate"
        ] = non_carbon_limit

        model_check_info["boundary_information"]["special_boundary_reaction"] = ""
        model_check_info["boundary_information"]["min_carbon_source_set"] = []
        model_info["biomass_by_special_boundary"] = (
            round(model.slim_optimize(), 3) if has_biomass else 0
        )
        if has_biomass:
            print(
                "Biomass output after special-boundary configuration:",
                model.slim_optimize(),
                ".......................................333........................",
            )

        set_auto_source_supply(model_info, model, check_model, "initial")

        model_control_info["boundary_information"]["carbon_source_boundary"] = ",".join(
            model_info["carbon_source_boundary"]
        )
        exchange_rxns = self.get_exchange(model, model_info, check_model)
        model_check_info["boundary_information"][
            "exchange_reaction_boundary"
        ] = exchange_rxns
        model_control_info["boundary_information"]["reaction_with_special_boundary"] = (
            ",".join(model_info["special_boundary_reaction"])
        )
        model_control_info["boundary_information"][
            "exchange_reaction_boundary"
        ] = exchange_rxns

        print(
            model.slim_optimize(),
            "...............................................................",
        )
        rules = Rules(self.cfg)
        run_rules(rules, model_info, model)

        add_respiratory_chain(model_info, model, check_model)
