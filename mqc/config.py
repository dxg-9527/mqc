"""Portable paths for packaged reference data and per-run outputs."""

from pathlib import Path

ROOT = Path(__file__).resolve().parent


class Config:
    """Resolve reference files relative to the installed package."""

    def __init__(self, output_dir):
        output = Path(output_dir).expanduser().resolve()
        output.mkdir(parents=True, exist_ok=True)
        self.output_dir = str(output)
        self.ROOT = ROOT
        self.flux_path = str(output / "flux.txt")
        self.flux_path2 = output / "flux2.txt"
        self.flux_path3 = output / "flux3.txt"
        self.flux_path4 = output / "flux.tsv"
        self.result_file = output / "check_result.json"
        self.result_file2 = output / "fail.txt"
        self.result_file3 = output / "control_result.json"
        self.result_file4 = output / "new_result.json"
        for extension in ["xml", "json", "yaml", "mat"]:
            setattr(
                self,
                "big_general_library_" + extension,
                str(output / ("big_model." + extension)),
            )
        self.big_general_library = ""
        self.biomass_rule_consistency = True
        # Independent integration trial only; the accepted package is unchanged.
        self.native_reaction_coverage = True
        self.repair_single_electron_templates = True
        self.model_info = str(output / "model_info.json")
        data = ROOT / "summary"
        self.general_library_path = str(data / "general_library.xml.gz")
        self.general_library_path_meta = str(data / "merged_meta_model.xml.gz")
        self.metacyc_rules = data / "MetacycRule.xlsx"
        self.modelseed_met = data / "modelseed_met.xlsx"
        self.meta_met = data / "meta_met.xlsx"
        self.kegg_met = data / "kegg_met.xlsx"
        self.bigg_metacyc = data / "bigg_metacyc.xlsx"
        self.modelseed_rxn = data / "modelseed_reactions.xlsx"
        self.virtual_met = data / "virtual_met.xlsx"
        for stem in ["nadh", "atp", "net", "yield", "yield2", "yield_pear"]:
            setattr(self, stem + "_png", str(output / (stem + ".png")))
