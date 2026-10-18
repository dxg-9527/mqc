from cobra.util.solver import linear_reaction_coefficients
import json
import multiprocessing
from multiprocessing import Pipe, Process
import matplotlib.pyplot as plt
import time
import numpy as np
from mqc.config import Config
from mqc.utils import *


class Check:
    """ """

    def __init__(self, cfg: Config) -> None:
        """ """
        self.cfg = cfg
        self.check_data = []

    def calculateResult(self, model, initial_rxn_id):
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
                    return all_mets[0]
                result -= v * rxns.get_coefficient(all_mets[0]) * relative_molecular_mass(model, all_mets[0])
        try:
            final = abs(result / model.slim_optimize())
        except ZeroDivisionError:
            final = 0
        return final

    def find_net_initial_flux(self, model_info, model, model_control_info):
        """
        Determining whether there is a net formation of substances
        """
        r_met, net_generation, only_pro_mets, met_num = [], [], [], 0
        initial_net_flux, temp = [], []
        for rxn in model_info["reactions"]:
            if rxn["id"] in model_info["exchange_rxns"] and rxn["bounds"][0] < 0:
                r_met += rxn["all_mets"]
            if len(rxn["reactants_mets"]) == 0 and len(rxn["products_mets"]) != 0:
                only_pro_mets += rxn["products_mets"]
        for met in model.metabolites:
            if met.name not in (FREE_METS + only_pro_mets):
                if re.search("C[A-Z]", str(met.formula)) or re.search(
                    r"C[\d]", str(met.formula)
                ):
                    met_num += 1
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

                            approach_info = get_net_approach2(model, met.id)
                            temp.append(approach_info)
        model_info["net_generation"] = net_generation

        if len(net_generation) != 0:
            self.check_data.append("0")
            model_control_info["check_metabolite_production"]["score"] = 0
            model_control_info["check_metabolite_production"][
                "Summary"
            ] = "Not passed. Erroneous net generation of some metabolites."
            model_control_info["check_metabolite_production"][
                "Details"
            ] = f"Calculating the optimal production rate of {met_num} carbon-containing metabolites, respectively when no carbon source was supplied. The erroneous net generation of {len(net_generation)} metabolites is as follows: \n{', '.join(net_generation)}"
            model_control_info["check_metabolite_production"][
                "net_dec"
            ] = "List of net generation metabolites and optimal rate without carbon source supplement"
            model_control_info["check_metabolite_production"]["Net_info"] = temp
            net_num = [len(net_generation), met_num - len(net_generation)]
            lables = ["error", "right"]
            plt.pie(net_num, labels=lables, autopct="%1.1f%%")

            plt.axis("equal")
            plt.savefig(self.cfg.net_png, dpi=300)
            plt.clf()
        else:
            self.check_data.append("1")
            model_control_info["check_metabolite_production"]["score"] = 1
            model_control_info["check_metabolite_production"][
                "Summary"
            ] = "Pass. All the metabolites did not generate form nothing."
            model_control_info["check_metabolite_production"][
                "Details"
            ] = f"Calculating the optimal production rate of {met_num} carbon-containing metabolites, respectively when no carbon source was supplied. All metabolites did not generate when no carbon source was supplied."

    def get_reducing_power_initial_flux(self, model_info, model, model_control_info):
        """"""
        flag, exist_nadh, nadh_flux = 0, [], []
        nadh_name, temp = [
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
            exist_nadh.append(name)
            model.reactions.get_by_id(nadh_id).bounds = (0, 1000)
            set_model_objective(model_info, model, nadh_id)
            nadh_flux.append(round(model.slim_optimize(), 3))
            if pd.isna(model.slim_optimize()) or model.slim_optimize() <= 1e-6:
                temp.extend([f"{name} - does not generate infinitely"])
            else:
                flag += 1
                temp.extend([f"{name} - generate infinitely, not passed"])

        if flag > 1:
            self.check_data.append("0")
            exist_nadh_dict = dict(zip(exist_nadh, nadh_flux))
            model_control_info["check_reducing_equivalents_production"]["score"] = 0
            model_control_info["check_reducing_equivalents_production"][
                "Summary"
            ] = "Not passed. Erroneous net generation of reducing equivalents."
            model_control_info["check_reducing_equivalents_production"][
                "Details"
            ] = f"Calculating the optimal production rate of {', '.join(exist_nadh)} respectively when no carbon source was supplied. The erroneous reducing equivalent generation is as follows: \n{[(k,v) for k,v in exist_nadh_dict.items() if v > 1e-6]}"
            plt.bar(exist_nadh, nadh_flux)
            plt.xlabel("Reducing equivalent substances")
            plt.ylabel("Optimal production rate")
            plt.savefig(self.cfg.nadh_png, dpi=300)
            plt.clf()
        else:
            self.check_data.append("1")
            model_control_info["check_reducing_equivalents_production"]["score"] = 1
            model_control_info["check_reducing_equivalents_production"][
                "Summary"
            ] = "Pass. All the reducing equivalent substances did not generate infinitely."
            model_control_info["check_reducing_equivalents_production"][
                "Details"
            ] = f"Calculating the optimal production rate of {', '.join(exist_nadh)} respectively when no carbon source is supplied. All the reducing equivalent substances did not generate when no carbon source was supplied"
        return self.cfg.nadh_png

    def get_energy_initial_flux(self, model_info, model, model_control_info):
        """"""
        xtp, temp = ["A", "C", "G", "U", "I"], []
        flag, exist_atp, atp_flux = 0, [], []
        for k in xtp:
            atpm_id = get_atpm_rxn(model_info, model, f"{k}TP", f"{k}DP")
            if not atpm_id:
                continue
            exist_atp.append(f"{k}TP")
            model.reactions.get_by_id(atpm_id).bounds = (0, 1000)
            set_model_objective(model_info, model, atpm_id)
            atp_flux.append(round(model.slim_optimize(), 3))
            if pd.isna(model.slim_optimize()) or model.slim_optimize() <= 1e-6:
                temp.extend([f"{k}TP - does not generate infinitely"])
            else:
                flag += 1
                temp.extend([f"{k}TP - generate infinitely, not passed"])

        if flag > 1:
            self.check_data.append("0")
            exist_atp_dict = dict(zip(exist_atp, atp_flux))
            model_control_info["check_energy_production"]["score"] = 0
            model_control_info["check_energy_production"][
                "Summary"
            ] = "Not passed. Erroneous net generation of energy."
            model_control_info["check_energy_production"][
                "Details"
            ] = f"Calculating the optimal production rate of {', '.join(exist_atp)} respectively when no carbon source was supplied. The erroneous energy generation is as follows: \n{[(k,v) for k,v in exist_atp_dict.items() if v > 1e-6]}"
            plt.bar(exist_atp, atp_flux)
            plt.xlabel("Energy substances")
            plt.ylabel("Optimal production rate")
            plt.savefig(self.cfg.atp_png, dpi=300)
            plt.clf()
        else:
            self.check_data.append("1")
            model_control_info["check_energy_production"]["score"] = 1
            model_control_info["check_energy_production"][
                "Summary"
            ] = "Pass. All the energy substances did not generate infinitely."
            model_control_info["check_energy_production"][
                "Details"
            ] = f"Calculating the optimal production rate of {', '.join(exist_atp)} respectively when no carbon source was supplied. All the energy substances did not generate when no carbon source was supplied."
        return self.cfg.atp_png

    def check_yield_formula(self, model, metId):
        """"""
        calculate_met = ["C", "H", "O", "N", "P", "S"]
        met = model.metabolites.get_by_id(metId)
        if not set(met.elements.keys()).issubset(set(calculate_met)):
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
        initial_yield_flux, temp, temps, substance_flux = [], [], [], []
        temp_degree, netName, maxRate, objRate = [], [], [], []
        pid = multiprocessing.current_process().pid
        print(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            f"pid {pid} model.metabolites num:",
            len(metabolites),
        )

        for met in metabolites:
            if met.id in model_info["need_set_carbon_mets"]:
                continue
            if re.search("C[A-Z]", str(met.formula)) or re.search(
                r"C[\d]", str(met.formula)
            ):
                c_num += 1
                if self.check_yield_formula(model, met.id) == 0:
                    print(met, ": Formula unresolved; skipping calculation")
                    continue
                try:
                    with model:
                        if met.name in ACCOA and met.compartment in C_COMPARTMENT:
                            objectiveId = add_accoa(model, str(met.id))
                        else:
                            objectiveId = add_demand(model_info, model, met.id)
                        if not objectiveId:
                            continue
                        model.objective = objectiveId

                        max_rate, temps = max_reduced_degree_rate(
                            model, met.id, model_info
                        )
                        temp_degree.extend(temps)
                        if (
                            max_rate == 0
                        ):  # Skip zero product limits and formulas containing elements outside CHONPS.
                            continue

                        model_opt = model.slim_optimize()

                        if round(model_opt, 2) > round(max_rate, 2):

                            yield_generation.append(met.id)
                            initial_yield_flux.extend(
                                [f"{met.id} : {round(model_opt,3)}"]
                            )

                            approach_info = get_net_approach(model, met.id, max_rate)
                            temp.append(approach_info)
                            netName.append(met.id)
                            maxRate.append(round(max_rate, 2))
                            objRate.append(round(model_opt, 2))

                except:
                    pass

        conn.send(
            {
                "yield_generation": yield_generation,
                "c_num": c_num,
                "initial_yield_flux": initial_yield_flux,
                "met_info": temp,
                "temp_degree": temp_degree,
                "netName": netName,
                "maxRate": maxRate,
                "objRate": objRate,
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
        temp_met_info = []
        temp_degree = []
        substance_flux, netName, maxRate, objRate = [], [], [], []
        for parent_connection in parent_connections:
            response = parent_connection.recv()
            yield_generation += response["yield_generation"]
            c_num += response["c_num"]
            initial_yield_flux += response["initial_yield_flux"]
            temp_met_info += response["met_info"]
            temp_degree += response["temp_degree"]
            netName += response["netName"]
            maxRate += response["maxRate"]
            objRate += response["objRate"]

        model_control_info["check_metabolite_yield"]["notes"] = ", ".join(temp_degree)
        model_info["yield_generation"] = yield_generation
        netName = [netName[k] for k in range(10) if k < len(netName)]
        new_maxRate = [maxRate[k] for k in range(10) if k < len(maxRate)]
        new_objRate = [objRate[k] for k in range(10) if k < len(objRate)]

        if len(yield_generation) != 0:
            self.check_data.append("0")
            model_control_info["check_metabolite_yield"]["score"] = 0
            model_control_info["check_metabolite_yield"][
                "Summary"
            ] = "Not passed. The calculated optimal yield of some metabolites exceeded the theoretical maximum yield."
            model_control_info["check_metabolite_yield"][
                "Details"
            ] = f"Calculating the optimal production rate of {c_num} carbon-containing metabolites, respectively when carbon source was supplied. The optimal yield of {len(yield_generation)} metabolites exceeded their theoretical maximum yield as follows:\n{', '.join(initial_yield_flux)}"
            model_control_info["check_metabolite_yield"][
                "yield_dec"
            ] = "List of metabolites whose calculated optimal yield exceeded their theoretical maximum yield"
            model_control_info["check_metabolite_yield"]["yield_info"] = temp_met_info
            yield_num = [len(yield_generation), c_num - len(yield_generation)]
            lables = ["error", "right"]
            plt.pie(yield_num, labels=lables, autopct="%1.1f%%")

            plt.axis("equal")
            plt.savefig(self.cfg.yield_png, dpi=300)
            plt.clf()
            width = 0.4  # Set the bar width.
            x = np.arange(len(netName))  # Set the x coordinates.
            x_middle = x + width / 2  # Center each pair of bars on its category.

            plt.bar(x_middle, new_maxRate, width=width, label="maximum_yield")
            plt.bar(x_middle + width, new_objRate, width=width, label="objective_value")

            plt.xticks(x_middle + width / 2, netName)  # Set tick positions and labels.
            plt.xticks(rotation=30, ha="center")
            plt.xticks(fontsize=7)
            plt.xlabel("Substance")
            plt.ylabel("Rate")

            plt.legend()  # Show the legend.
            plt.savefig(self.cfg.yield2_png, dpi=300)
            plt.clf()
        else:
            self.check_data.append("1")
            model_control_info["check_metabolite_yield"]["score"] = 1
            model_control_info["check_metabolite_yield"][
                "Summary"
            ] = "Pass. The calculated optimal yield of all metabolites did not exceed the theoretical maximum yield. "
            model_control_info["check_metabolite_yield"][
                "Details"
            ] = f"Calculating the optimal production rate of {c_num} carbon-containing metabolites, respectively when carbon source was supplied. The optimal yield of  all metabolites did not exceed their theoretical maximum."

        return {
            "yield_met_id": yield_generation,
            "max_rate": maxRate,
            "obj_rate": objRate,
        }

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

    def check_control(
        self, model_info, model, check_model, model_control_info, model_check_info
    ):
        """"""

        self.check_data.append(f"{model.id}")
        initial_rxn_id = model_info["initial_rxn"]["biomass_rxn_id"]
        if initial_rxn_id:
            model.objective = initial_rxn_id
        del_model_boundary(model_info, model, check_model, model_control_info)

        set_special_boundary_rxn(model_info, model, check_model)
        set_bio_boundary(model_info, model)

        close_autotrophic_or_c_source(model_info, model, "carbon_source")

        reducing_pic_path = self.get_reducing_power_initial_flux(
            model_info, model, model_check_info
        )
        energy_pic_path = self.get_energy_initial_flux(
            model_info, model, model_check_info
        )
        self.find_net_initial_flux(model_info, model, model_check_info)

        get_carbon_proportion(model_info, model, check_model, model_control_info)
        set_c_source_supply(model_info, model, "yields", check_model)
        set_auto_source_supply(model_info, model, check_model, "yields")

        yield_data = {}
        if is_autotrophic(model_info):
            model_check_info["check_metabolite_yield"][
                "autotrophic"
            ] = "This is an autotrophic strain, no need to calculate the yield"
        else:
            get_c_source_supply(model_info, model)
            yield_data = self.find_yield_initial_flux(
                model_info, model, model_check_info
            )

        if initial_rxn_id:
            model.objective = initial_rxn_id
            temp = [
                f"{model_info['initial_rxn']['biomass_rxn_id']} : {model_info['initial_rxn']['biomass_rxn_exp']}"
            ]
            temp.append(model_info["initial_rxn"]["check_biomass_composition"])
            model_check_info["check_biomass_production"]["biomass_info"] = ", ".join(
                temp
            )
            ini = round(model_info["initial_rxn"]["biomass_rxn_flux"], 3)
            if pd.isna(model.slim_optimize()) or ini <= 1e-6:
                model_check_info["check_biomass_production"][
                    "Summary"
                ] = f"Not passed.The biomass growth rate was 0 h-1. It is necessary to fill the gaps for the biomass components that cannot be generated"
                model_check_info["check_biomass_production"]["score"] = 0
                self.check_data.append("0")
            else:
                biomass_result = self.calculateResult(model, initial_rxn_id)
                if type(biomass_result) != str and (990 < biomass_result < 1010):
                    model_check_info["check_biomass_production"][
                        "Summary"
                    ] = f"Pass. The biomass composition was 1g."
                    model_check_info["check_biomass_production"]["score"] = 1
                    self.check_data.append("1")
                else:
                    model_check_info["check_biomass_production"][
                        "Summary"
                    ] = f"Not passed.The biomass composition was {biomass_result}mg, not 1g. It is necessary to perform quality control for biomass to ensure the mass of 1g"
                    model_check_info["check_biomass_production"]["score"] = 0
                    self.check_data.append("0")
        else:
            model_check_info["check_biomass_production"][
                "Summary"
            ] = "No biomass reaction was identified; this check did not pass."
            model_check_info["check_biomass_production"]["biomass_info"] = ""
            model_check_info["check_biomass_production"]["score"] = 0
            self.check_data.append("0")
        pear_data = {
            "reduce": reducing_pic_path,
            "energy": energy_pic_path,
            "yield": yield_data,
        }
        return self.check_data, pear_data
