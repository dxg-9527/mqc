"""Detection and correction workflows used by the MQC web service."""

import argparse
import os
import time
from pathlib import Path

from mqc.config import Config
from mqc.utils import (
    change_model_control_info,
    convert_list_to_string,
    convert_nan_to_null,
    get_model_info,
    get_new_model,
    read_model,
    reduct_modelinfo_modify,
    return_restoration,
    run_rules,
    write_check_model,
    write_final_model,
    write_model_info,
    write_result,
    write_result2,
    write_result3,
    write_result4,
)
import cobra
from mqc.control.preprocessing_control import Preprocess
from mqc.control.model_control import ModelPreprocess
from mqc.control.initial_control import InitialPreprocess
from mqc.control.check_control import Check
from mqc.control.rules_control import Rules
from mqc.control.nadh_control import Nadhs
from mqc.control.atp_control import Atps
from mqc.control.net_control import Nets
from mqc.control.yield_control import Yields
from mqc.control.biomass_control import Biomasses
from mqc.control.postgap_control import run_postgap_control, PostGapValidationError
from mqc.control.electron_stoichiometry import validate_final_templates
from mqc.control.gap_rule_policy import validate_final_rule_bounds


class ModelCheck:
    """
    Obtain the total output information of the model quality control.

    """

    def __init__(self, model_file: str, output_dir: str):
        """
        define output dictionary.

        """

        self.model_file = model_file
        self.output_dir = self.create_outputdir(output_dir)
        self.cfg = Config(self.output_dir)
        self.model_control_info = {}
        self.model_check_info = {}
        self.model_check_info["boundary_information"] = {}
        self.all_data = []
        self.model_check_info["check_reducing_equivalents_production"] = {}
        self.model_check_info["check_energy_production"] = {}
        self.model_check_info["check_metabolite_production"] = {}
        self.model_check_info["check_metabolite_yield"] = {}
        self.model_check_info["check_biomass_production"] = {}

        self.model_control_info["boundary_information"] = {}
        self.model_control_info["check_reducing_equivalents_production"] = {}
        self.model_control_info["check_energy_production"] = {}
        self.model_control_info["check_metabolite_production"] = {}
        self.model_control_info["check_metabolite_yield"] = {}
        self.model_control_info["check_biomass_production"] = {}

    def create_outputdir(self, output_dir):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        return output_dir

    def model_check(self):
        """"""
        t1 = time.time()
        controler = Preprocess(self.model_file, self.cfg)
        if not controler.model:
            model_control_info = change_model_control_info(
                self.model_file, self.model_control_info
            )
            model_control_infos = write_result2(model_control_info, self.cfg)
            t2 = time.time()
            print("Total time: ", t2 - t1)
            return model_control_infos
        model_pre_process = ModelPreprocess(self.cfg)
        model_pre_process.get_model_info(controler)
        checks = Check(self.cfg)
        all_data = checks.check_control(
            model_pre_process.model_info,
            controler.model,
            controler.check_model,
            self.model_control_info,
            self.model_check_info,
        )

        del self.model_check_info["boundary_information"]
        model_control_infos = write_result(self.model_check_info, self.cfg)
        t2 = time.time()
        print("Total time: ", t2 - t1)
        return model_control_infos, all_data

    def model_check2(self):
        """
        The overall process of model quality control.

        """
        t1 = time.time()
        model_control_infos, final_model = "", ""
        self.model_check_info["boundary_information"] = {}
        controler = Preprocess(self.model_file, self.cfg)
        if not controler.model:
            model_control_info = change_model_control_info(
                self.model_file, self.model_control_info
            )
            model_control_infos = write_result2(model_control_info, self.cfg)
            final_model = self.model_file
            t2 = time.time()
            print("Total time: ", t2 - t1)
            return model_control_infos

        model_pre_process = ModelPreprocess(self.cfg)
        model_pre_process.get_model_info(controler)
        model_pre_process.model_info["model_file"] = self.model_file

        try:
            initial_pre_process = InitialPreprocess(self.cfg)
            initial_pre_process.initial_control(
                model_pre_process.model_info,
                controler.model,
                controler.check_model,
                self.model_control_info,
                self.model_check_info,
            )

            final_model = write_final_model(controler.model, self.cfg, self.model_file)

            controler.model = read_model(final_model)

            nadhs = Nadhs(self.cfg)
            nadhs.nadh_control(
                model_pre_process.model_info,
                controler.model,
                controler.check_model,
                self.model_control_info,
                self.model_check_info,
                controler,
            )
            final_model = write_final_model(controler.model, self.cfg, self.model_file)
            controler.model = read_model(final_model)

            atps = Atps(self.cfg)
            atps.atp_control(
                model_pre_process.model_info,
                controler.model,
                controler.check_model,
                self.model_control_info,
                self.model_check_info,
                controler,
            )
            final_model = write_final_model(controler.model, self.cfg, self.model_file)
            controler.model = read_model(final_model)

            nets = Nets(self.cfg)
            nets.net_control(
                model_pre_process.model_info,
                controler.model,
                controler.check_model,
                self.model_control_info,
                self.model_check_info,
                controler,
            )
            final_model = write_final_model(controler.model, self.cfg, self.model_file)
            controler.model = read_model(final_model)

            yields = Yields(self.cfg)
            yield_one = yields.yield_control(
                model_pre_process.model_info,
                controler.model,
                controler.check_model,
                self.model_control_info,
                controler,
                self.model_check_info,
            )
            if yield_one == 1:
                return_restoration(
                    model_pre_process.model_info,
                    controler.model,
                    controler.check_model,
                    Biomasses(self.cfg),
                )
                yields.yield_control(
                    model_pre_process.model_info,
                    controler.model,
                    controler.check_model,
                    self.model_control_info,
                    controler,
                    self.model_check_info,
                )
            final_model = write_final_model(controler.model, self.cfg, self.model_file)
            controler.model = read_model(final_model)

            biomasses = Biomasses(self.cfg)
            pre_gapfill_ids = {r.id for r in controler.model.reactions}
            biomasses.biomass_control(
                model_pre_process.model_info,
                controler.model,
                controler.check_model,
                self.model_control_info,
                controler,
                self.model_check_info,
            )
            run_postgap_control(model_pre_process.model_info, controler.model,
                                self.model_control_info, self.cfg, pre_gapfill_ids)
            validate_final_templates(controler.model, model_pre_process.model_info,
                                     self.model_control_info)
            validate_final_rule_bounds(controler.model, model_pre_process.model_info,
                                       self.model_control_info)

            convert_nan_to_null(controler.model)
            convert_nan_to_null(controler.check_model)

            convert_list_to_string(self.model_control_info, self.model_check_info)

            final_model = write_final_model(
                controler.model, self.cfg, self.model_file,
                biomass_rxn_id=model_pre_process.model_info['initial_rxn']['biomass_rxn_id'] or None)
            check_model = write_check_model(
                controler.check_model, self.cfg, self.model_file
            )
        except PostGapValidationError:
            write_result2(self.model_control_info, self.cfg)
            write_result3(self.model_control_info, self.cfg)
            write_model_info(model_pre_process.model_info, self.cfg)
            raise
        except RuntimeError as e:
            final_model = write_final_model(controler.model, self.cfg, self.model_file)
            check_model = write_check_model(
                controler.check_model, self.cfg, self.model_file
            )
            model_control_infos = write_result2(self.model_control_info, self.cfg)
            print(repr(e), ".............")

        except Exception as e:
            model_control_infos = write_result2(self.model_control_info, self.cfg)
            final_model = write_final_model(controler.model, self.cfg, self.model_file)
            check_model = write_check_model(
                controler.check_model, self.cfg, self.model_file
            )
            print(repr(e), ".............")
            raise

        write_model_info(model_pre_process.model_info, self.cfg)

        model_control_infos = write_result3(self.model_control_info, self.cfg)

        t2 = time.time()
        print(
            controler.model.slim_optimize(),
            "...............................................................",
        )

        print("Total time: ", t2 - t1)
        return model_control_infos, final_model, check_model


class ModelControl:
    """
    Obtain the total output information of the model quality control.

    """

    def __init__(
        self,
        model_file: str,
        output_dir: str,
        types: str,
        rxns_job: list,
        check_model: str,
    ):
        """
        define output dictionary.

        """

        self.output_dir = self.create_outputdir(output_dir)
        self.cfg = Config(self.output_dir)
        self.model_file = model_file
        self.types = types
        self.rxns_job = rxns_job
        self.check_model = check_model
        self.model_check_info = {}
        self.model_check_info["boundary_information"] = {}
        self.model_check_info["check_reducing_equivalents_production"] = {}
        self.model_check_info["check_energy_production"] = {}
        self.model_check_info["check_metabolite_production"] = {}
        self.model_check_info["check_metabolite_yield"] = {}
        self.model_check_info["check_biomass_production"] = {}

        self.model_control_info = {}
        self.model_control_info["boundary_information"] = {}
        self.model_control_info["check_reducing_equivalents_production"] = {}
        self.model_control_info["check_energy_production"] = {}
        self.model_control_info["check_metabolite_production"] = {}
        self.model_control_info["check_metabolite_yield"] = {}
        self.model_control_info["check_biomass_production"] = {}

        self.temp_model_control_info = {}
        self.model_control_info["check_biomass_production"][
            "Check_synthesizability_of_biomass_components"
        ] = []
        self.model_control_info["check_biomass_production"][
            "Gapfilling_of_biomass_components"
        ] = []

    def create_outputdir(self, output_dir):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        return output_dir

    def get_initial_check_model(self, controler, check_model):
        """"""
        if check_model.endswith(".json"):
            controler.check_model = cobra.io.load_json_model(check_model)
        elif check_model.endswith(".yaml"):
            controler.check_model = cobra.io.load_yaml_model(check_model)
        elif check_model.endswith(".mat"):
            controler.check_model = cobra.io.load_matlab_model(check_model)
        else:
            controler.check_model = cobra.io.read_sbml_model(check_model)

    def model_control(self):
        """
        The overall process of model quality control.

        """
        t1 = time.time()

        model_path = get_new_model(self.model_file, self.rxns_job, self.cfg)
        model_info = get_model_info(self.cfg)
        reduct_modelinfo_modify(model_info)
        controler = Preprocess(model_path, self.cfg)
        self.get_initial_check_model(controler, self.check_model)

        print(
            controler.model.slim_optimize(),
            "...............................................................",
        )
        try:
            initial_pre_process = InitialPreprocess(self.cfg)
            final_model = write_final_model(controler.model, self.cfg, self.model_file)
            controler.model = cobra.io.read_sbml_model(final_model)
            rules = Rules(self.cfg)
            run_rules(rules, model_info, controler.model)

            if self.types == "boundary_information":
                initial_pre_process.initial_control(
                    model_info,
                    controler.model,
                    controler.check_model,
                    self.model_control_info,
                    self.model_check_info,
                )

            if self.types == "check_reducing_equivalents_production":
                nadhs = Nadhs(self.cfg)
                nadhs.nadh_control(
                    model_info,
                    controler.model,
                    controler.check_model,
                    self.model_control_info,
                    self.model_check_info,
                    controler,
                )

            if self.types == "check_energy_production":
                atps = Atps(self.cfg)
                atps.atp_control(
                    model_info,
                    controler.model,
                    controler.check_model,
                    self.model_control_info,
                    self.model_check_info,
                    controler,
                )

                write_model_info(model_info, self.cfg)

            if self.types == "check_metabolite_production":
                nets = Nets(self.cfg)
                nets.net_control(
                    model_info,
                    controler.model,
                    controler.check_model,
                    self.model_control_info,
                    self.model_check_info,
                    controler,
                )

            if self.types == "check_metabolite_yield":
                yields = Yields(self.cfg)
                yields.yield_control(
                    model_info,
                    controler.model,
                    controler.check_model,
                    self.model_control_info,
                    controler,
                    self.model_check_info,
                )

            if self.types == "check_biomass_production":
                biomasses = Biomasses(self.cfg)
                biomasses.biomass_control(
                    model_info,
                    controler.model,
                    controler.check_model,
                    self.model_control_info,
                    controler,
                    self.model_check_info,
                )
                convert_list_to_string(self.model_control_info, self.model_check_info)

            convert_nan_to_null(controler.model)
            final_model = write_final_model(controler.model, self.cfg, self.model_file)
        except RuntimeError as e:
            final_model = ""
            print(repr(e), ".............")
        except Exception as e:
            print(repr(e), ".............")
            raise
        self.temp_model_control_info[self.types] = self.model_control_info[self.types]
        model_control_infos = write_result4(self.temp_model_control_info, self.cfg)
        modelInfo = write_model_info(model_info, self.cfg)
        t2 = time.time()
        print("Total time: ", t2 - t1)
        return model_control_infos, final_model, modelInfo


def run(model_file, output_dir, mode="correct"):
    """Detect quantitative defects and optionally correct the model.

    ``check`` writes the detection report and returns its path and module flags.
    ``correct`` first runs detection, then returns the correction report path,
    revised model path and reference model path, following the web workflow.
    """
    model_file = Path(model_file).expanduser().resolve()
    if not model_file.is_file():
        raise FileNotFoundError("Model file does not exist: {}".format(model_file))
    if mode not in {"check", "correct"}:
        raise ValueError("Mode must be 'check' or 'correct'.")
    output_dir = Path(output_dir).expanduser().resolve()
    workflow = ModelCheck(str(model_file), str(output_dir))
    detected = workflow.model_check()
    if not isinstance(detected, tuple):
        raise RuntimeError("Model detection failed; inspect {}".format(detected))
    if mode == "check":
        return detected
    corrected = workflow.model_check2()
    if not isinstance(corrected, tuple):
        raise RuntimeError("Model correction failed; inspect {}".format(corrected))
    return corrected


def main(argv=None):
    """Run MQC from a terminal."""
    parser = argparse.ArgumentParser(description="MQC metabolic model quality control")
    parser.add_argument(
        "model", type=Path, help="Input SBML, JSON, YAML or MATLAB model"
    )
    parser.add_argument(
        "output", type=Path, help="Directory for reports and model outputs"
    )
    parser.add_argument(
        "--mode",
        choices=["check", "correct"],
        default="correct",
        help="Run detection only or detection followed by correction",
    )
    args = parser.parse_args(argv)
    try:
        result = run(args.model, args.output, args.mode)
    except (FileNotFoundError, ValueError, RuntimeError) as error:
        parser.exit(1, "MQC failed: {}\n".format(error))
    print("MQC completed. Outputs:")
    for value in result:
        print(value)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
