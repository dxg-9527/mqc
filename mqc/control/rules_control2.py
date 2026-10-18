import re
import pandas as pd
import fnmatch
from mqc.defaults import *
from mqc.config import Config
from mqc.control.rules_control import (
    ammonia_assimilation_supported, glutamate_synthesis_direction,
    respiratory_membership_supported,
)


class Rules2:
    """
    Get all rule information
    """

    def find_sugar_identity_conflicts(self, model_info, model):
        from mqc.control.pathway_rules import sugar_identity_conflict, SUGAR_IDENTITY_RULE
        if abs(self.flux) <= 1e-9:
            return
        if sugar_identity_conflict(model.reactions.get_by_id(self.rxnId)):
            for record in model_info['reactions']:
                if record['id'] == self.rxnId:
                    record['rule_conflict'] = 'true'
                    record['rules2'][SUGAR_IDENTITY_RULE] = 'true'

    def find_biotin_synthesis_direction(self, model_info, model):
        from mqc.control.biotin_direction import synthesis_direction, BIOTIN_RULE
        if not getattr(self.cfg, 'biomass_rule_consistency', True):
            return
        reaction = model.reactions.get_by_id(self.rxnId)
        direction = synthesis_direction(reaction)
        if direction is not None and self.flux * direction < -1e-9:
            for record in model_info['reactions']:
                if record['id'] == self.rxnId:
                    record['rule_conflict'] = 'true'
                    record['rules2'][BIOTIN_RULE] = 'true'

    def find_deacetylase_direction(self, model_info, model):
        from mqc.control.deacetylase_direction import deacetylase_direction, ACETATE_RULE
        if not getattr(self.cfg, 'biomass_rule_consistency', True):
            return
        direction = deacetylase_direction(model.reactions.get_by_id(self.rxnId))
        if direction is not None and self.flux * direction < -1e-9:
            for record in model_info['reactions']:
                if record['id'] == self.rxnId:
                    record['rule_conflict'] = 'true'
                    record['rules2'][ACETATE_RULE] = 'true'

    def __init__(self, rxnId: str, cfg: Config, v: str):
        """"""
        self.rxnId = rxnId
        self.cfg = cfg
        self.flux = v

    def find_desulfurase_direction(self, model_info, model):
        from mqc.control.desulfurase_direction import desulfuration_direction, DESULFURASE_RULE
        if not getattr(self.cfg, 'biomass_rule_consistency', True):
            return
        direction = desulfuration_direction(model.reactions.get_by_id(self.rxnId))
        if direction is not None and self.flux * direction < -1e-9:
            for record in model_info['reactions']:
                if record['id'] == self.rxnId:
                    record['rule_conflict'] = 'true'
                    record['rules2'][DESULFURASE_RULE] = 'true'

    def initial_rule(self, model_info):
        """"""
        for rxn in model_info["reactions"]:
            if "rule_conflict" in rxn.keys():
                del rxn["rule_conflict"]
            rxn["rules2"] = {}

    def find_incomplete_quinone_methylation(self, model_info, model):
        from mqc.control.quinone_methylation import incomplete_methylation, METHYL_REDOX_RULE
        if not getattr(self.cfg, 'biomass_rule_consistency', True) or abs(self.flux) <= 1e-9:
            return
        if incomplete_methylation(model.reactions.get_by_id(self.rxnId)):
            for record in model_info['reactions']:
                if record['id'] == self.rxnId:
                    record['rule_conflict'] = 'true'
                    record['rules2'][METHYL_REDOX_RULE] = 'true'

    def find_o2s(self, model_info):
        """
        Check external uptake, with the same scope as Rules.find_o2s.

        Intracellular production is not uptake from the environment. Keep the
        original species list and flux signs; only restore the exchange guard.
        """
        if self.rxnId not in model_info["exchange_rxns"]:
            return
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if (
                    len(set(rxn["reactants_mets"]) & set(O2S + SO3_NAME)) != 0
                    and self.flux < 0
                ) or (
                    len(set(rxn["products_mets"]) & set(O2S + SO3_NAME)) != 0
                    and self.flux > 0
                ):
                    rxn["rule_conflict"] = "true"
                    rxn["rules2"][
                        "Superoxide anion cannot be taken up by cells directly"
                    ] = "true"

    def find_O2_inorganic_rxn(self, model_info):
        """
        find inorganic reactions involving oxygen
        """
        O2_inorganic_rxn = []
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if (
                    rxn["id"] not in model_info["transport_rxns"]
                    and rxn["id"] not in model_info["exchange_rxns"]
                ):
                    n, s = 0, 0
                    for met in rxn["all_mets"]:
                        s += 1
                        if met in O2_METS:
                            n += 1
                    if s == n:
                        O2_inorganic_rxn.append(rxn["id"])
        return O2_inorganic_rxn

    def find_superoxide_rxn(self, model_info):
        """
        superoxide, reaction of peroxide to produce oxygen

        Notes
        -----
        SPODMpp: 2.0 h_p + 2.0 o2s_p --> h2o2_p + o2_p
        SOD_m: 2.0 sox_m --> h2o2_m + o2_m
        GTHP_CAT: 2.0 gthrd_c + 3.0 h2o2_c --> gthox_c + 4.0 h2o_c + o2_c   Reaction of peroxide and reduced glutathione to produce oxygen

        """
        superoxide_rxn = []
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if (
                    len(set(O2S) & set(rxn["reactants_mets"])) != 0
                    or len(set(SOX) & set(rxn["reactants_mets"])) != 0
                ):
                    if (
                        len(set(H2O2_NAME) & set(rxn["products_mets"])) != 0
                        and len(set(O2_NAME) & set(rxn["products_mets"])) != 0
                    ):
                        superoxide_rxn.append(rxn["id"])
                if (
                    len(set(O2S) & set(rxn["products_mets"])) != 0
                    or len(set(SOX) & set(rxn["products_mets"])) != 0
                ):
                    if (
                        len(set(H2O2_NAME) & set(rxn["reactants_mets"])) != 0
                        and len(set(O2_NAME) & set(rxn["reactants_mets"])) != 0
                    ):
                        superoxide_rxn.append(rxn["id"])
                if (
                    len(set(GTHRD) & set(rxn["reactants_mets"])) != 0
                    or len(set(H2O2_NAME) & set(rxn["reactants_mets"])) != 0
                ):
                    if (
                        len(set(GTHOX) & set(rxn["products_mets"])) != 0
                        and len(set(O2_NAME) & set(rxn["products_mets"])) != 0
                    ):
                        superoxide_rxn.append(rxn["id"])
                if (
                    len(set(GTHRD) & set(rxn["products_mets"])) != 0
                    or len(set(H2O2_NAME) & set(rxn["products_mets"])) != 0
                ):
                    if (
                        len(set(GTHOX) & set(rxn["reactants_mets"])) != 0
                        and len(set(O2_NAME) & set(rxn["reactants_mets"])) != 0
                    ):
                        superoxide_rxn.append(rxn["id"])
        return superoxide_rxn

    def find_photosynthetic_rxn(self, model_info):
        """
        What are the photosynthetic bacteria, the reaction of photosynthesis to produce oxygen
        """
        photosynthetic_rxn = []
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if re.search(r"photosystem", rxn["name"].lower()):
                    photosynthetic_rxn.append(rxn["id"])
                if (
                    len(set(PHOTON_NAME) & set(rxn["reactants_mets"])) != 0
                    and len(set(O2_NAME) & set(rxn["products_mets"])) != 0
                ):
                    photosynthetic_rxn.append(rxn["id"])
                if (
                    len(set(PHOTON_NAME) & set(rxn["products_mets"])) != 0
                    and len(set(O2_NAME) & set(rxn["reactants_mets"])) != 0
                ):
                    photosynthetic_rxn.append(rxn["id"])
        return photosynthetic_rxn

    def find_O2_rxn(self, model_info):
        """
        Find the reaction that obeys the oxygen rule
        """
        O2_rxn = []
        O2_inorganic_rxn = self.find_O2_inorganic_rxn(model_info)
        superoxide_rxn = self.find_superoxide_rxn(model_info)
        photosynthetic_rxn = self.find_photosynthetic_rxn(model_info)
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if (
                    rxn["id"] not in model_info["exchange_rxns"]
                    and rxn["id"] not in model_info["transport_rxns"]
                    and rxn["id"] not in O2_inorganic_rxn
                    and rxn["id"] not in superoxide_rxn
                    and rxn["id"] not in photosynthetic_rxn
                ):
                    if len(set(O2_NAME) & set(rxn["reactants_mets"])) != 0:
                        if self.flux < 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"][
                                "The direction of the reaction containing oxygen is the direction of oxygen consumption"
                            ] = "true"
                            O2_rxn.append(rxn["id"])
                    if len(set(O2_NAME) & set(rxn["products_mets"])) != 0:
                        if self.flux > 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"][
                                "The direction of the reaction containing oxygen is the direction of oxygen consumption"
                            ] = "true"
                            O2_rxn.append(rxn["id"])
        return O2_rxn

    def find_nh4_inorganic_rxn(self, model_info):
        """
        find the inorganic reaction of nh3 and nh4
        """
        nh4_inorganic_rxn = []
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if (
                    rxn["id"] not in model_info["exchange_rxns"]
                    and rxn["id"] not in model_info["transport_rxns"]
                ):
                    if (
                        len(set(NH3_NAME) & set(rxn["reactants_mets"])) != 0
                        and len(set(NH4_NAME) & set(rxn["products_mets"])) != 0
                    ):
                        nh4_inorganic_rxn.append(rxn["id"])
                    if (
                        len(set(NH3_NAME) & set(rxn["products_mets"])) != 0
                        and len(set(NH4_NAME) & set(rxn["reactants_mets"])) != 0
                    ):
                        nh4_inorganic_rxn.append(rxn["id"])
        return nh4_inorganic_rxn

    def find_nh3_nh4_rxn(self, model_info, O2_rxn):
        """
        Except nh3, nh4 and ATP, nadh, nadph, chor, prpp, udpacblfuc reactions can fix ammonia, others are reactions to generate ammonia
        """
        nh4_inorganic_rxn = self.find_nh4_inorganic_rxn(model_info)
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if (
                    rxn["id"] not in model_info["exchange_rxns"]
                    and rxn["id"] not in model_info["transport_rxns"]
                    and rxn["id"] not in nh4_inorganic_rxn
                    and rxn["id"] not in O2_rxn
                ):
                    if (
                        len(set(NH4_NAME) & set(rxn["reactants_mets"])) != 0
                        or len(set(NH3_NAME) & set(rxn["reactants_mets"])) != 0
                    ):
                        if not ammonia_assimilation_supported(rxn["reactants_mets"], rxn["products_mets"]):
                            if self.flux > 0:
                                rxn["rule_conflict"] = "true"
                                rxn["rules2"][
                                    "Most reactions producing NH3/NH4+ were assumed to consume NH3/NH4+ when it reacted with ATP, 2-oxoglutarate, chorismate, 5-phospho-alpha-D-ribose-1-diphosphate, and UDP-N-acetyl-beta-L-fucosamine"
                                ] = "true"
                    if (
                        len(set(NH4_NAME) & set(rxn["products_mets"])) != 0
                        or len(set(NH3_NAME) & set(rxn["products_mets"])) != 0
                    ):
                        if not ammonia_assimilation_supported(rxn["products_mets"], rxn["reactants_mets"]):
                            if self.flux < 0:
                                rxn["rule_conflict"] = "true"
                                rxn["rules2"][
                                    "Most reactions producing NH3/NH4+ were assumed to consume NH3/NH4+ when it reacted with ATP, 2-oxoglutarate, chorismate, 5-phospho-alpha-D-ribose-1-diphosphate, and UDP-N-acetyl-beta-L-fucosamine"
                                ] = "true"

    def find_CO2_rxn_all(self, model_info):
        """
        Find all CO2 reactions
        """
        CO2_rxn_all = []
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if (
                    rxn["id"] not in model_info["exchange_rxns"]
                    and rxn["id"] not in model_info["transport_rxns"]
                ):
                    if (rxn["bounds"][0] < 0 and rxn["bounds"][1] > 0) and len(
                        set(CO2_NAME) & set(rxn["all_mets"])
                    ) != 0:
                        CO2_rxn_all.append(rxn["id"])
                    if (rxn["bounds"][0] == 0 and rxn["bounds"][1] > 0) and len(
                        set(CO2_NAME) & set(rxn["reactants_mets"])
                    ) != 0:
                        CO2_rxn_all.append(rxn["id"])
                    if (rxn["bounds"][0] < 0 and rxn["bounds"][1] == 0) and len(
                        set(CO2_NAME) & set(rxn["products_mets"])
                    ) != 0:
                        CO2_rxn_all.append(rxn["id"])
        return CO2_rxn_all

    def find_co2_inorganic_rxn(self, model_info):
        """
        find the inorganic reaction involving co2
        """
        co2_inorganic_rxn = []
        CO2_rxn_all = self.find_CO2_rxn_all(model_info)
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if rxn["id"] in CO2_rxn_all and set(rxn["all_mets"]).issubset(
                    CO2_INORGANIC
                ):
                    co2_inorganic_rxn.append(rxn["id"])
        return co2_inorganic_rxn

    def find_co2_atp_pep_rxn(self, model_info):
        """
        co2 and high energy compound atp pep can fix carbon
        """
        CO2_rxn_all = self.find_CO2_rxn_all(model_info)
        co2_atp_pep_rxn = []
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if rxn["id"] in CO2_rxn_all:
                    if len(set(CO2_NAME) & set(rxn["reactants_mets"])) != 0:
                        if (
                            len(set(ATP_NAME) & set(rxn["reactants_mets"])) != 0
                            or len(set(PEP) & set(rxn["reactants_mets"])) != 0
                        ):
                            co2_atp_pep_rxn.append(rxn["id"])

                    if len(set(CO2_NAME) & set(rxn["products_mets"])) != 0:
                        if (
                            len(set(ATP_NAME) & set(rxn["products_mets"])) != 0
                            or len(set(PEP) & set(rxn["products_mets"])) != 0
                        ):
                            co2_atp_pep_rxn.append(rxn["id"])
        return co2_atp_pep_rxn

    def find_natural_co2_rxn(self, model_info):
        """
        6 natural carbon fixation reactions
        """
        natural_co2_rxn = []
        CO2_rxn_all = self.find_CO2_rxn_all(model_info)
        co2_inorganic_rxn = self.find_co2_inorganic_rxn(model_info)

        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if rxn["id"] in CO2_rxn_all and rxn["id"] not in co2_inorganic_rxn:
                    if (
                        len(set(CO2_NAME) & set(rxn["reactants_mets"])) != 0
                        and len(set(RB15BP) & set(rxn["reactants_mets"])) != 0
                        and len(set(PG3) & set(rxn["products_mets"])) != 0
                    ):
                        natural_co2_rxn.append(rxn["id"])
                    if (
                        len(set(CO2_NAME) & set(rxn["products_mets"])) != 0
                        and len(set(RB15BP) & set(rxn["products_mets"])) != 0
                        and len(set(PG3) & set(rxn["reactants_mets"])) != 0
                    ):
                        natural_co2_rxn.append(rxn["id"])

                    if (
                        len(set(CO2_NAME) & set(rxn["reactants_mets"])) != 0
                        and len(set(FOR) & set(rxn["products_mets"])) != 0
                    ):
                        natural_co2_rxn.append(rxn["id"])
                    if (
                        len(set(CO2_NAME) & set(rxn["products_mets"])) != 0
                        and len(set(FOR) & set(rxn["reactants_mets"])) != 0
                    ):
                        natural_co2_rxn.append(rxn["id"])

                    if (
                        len(set(CO2_NAME) & set(rxn["reactants_mets"])) != 0
                        and len(set(FOR) & set(rxn["products_mets"])) != 0
                    ):
                        natural_co2_rxn.append(rxn["id"])
                    if (
                        len(set(CO2_NAME) & set(rxn["products_mets"])) != 0
                        and len(set(FOR) & set(rxn["reactants_mets"])) != 0
                    ):
                        natural_co2_rxn.append(rxn["id"])

                    natural_co2_rxn.append("CODH_ACS")
                    natural_co2_rxn.append("MNXR96830")

                    if set(CO2_TCA_I).issubset(rxn["all_mets"]):
                        natural_co2_rxn.append(rxn["id"])
                    if set(CO2_TCA_II).issubset(rxn["all_mets"]):
                        natural_co2_rxn.append(rxn["id"])

                    if set(CO2_ACCOA_I).issubset(rxn["all_mets"]):
                        natural_co2_rxn.append(rxn["id"])
                    if set(CO2_ACCOA_II).issubset(rxn["all_mets"]):
                        natural_co2_rxn.append(rxn["id"])
        return natural_co2_rxn

    def find_CO2_rxn(self, model_info):
        """
        Find reactions that do not comply with the CO2 rule
        """
        CO2_rxn_all = self.find_CO2_rxn_all(model_info)
        co2_inorganic_rxn = self.find_co2_inorganic_rxn(model_info)
        co2_atp_pep_rxn = self.find_co2_atp_pep_rxn(model_info)
        natural_co2_rxn = self.find_natural_co2_rxn(model_info)
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if rxn["id"] in list(
                    set(CO2_rxn_all)
                    - set(co2_inorganic_rxn)
                    - set(co2_atp_pep_rxn)
                    - set(natural_co2_rxn)
                ):
                    if len(set(CO2_NAME) & set(rxn["reactants_mets"])) != 0:
                        if self.flux > 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"][
                                "Most reactions cannot proceed in the direction of CO2 fixation, except for naturally carbon fixation reactions and and reactions utilizing CO2 and high-energy substrates such as phosphoenolpyruvate (PEP) and ATP"
                            ] = "true"
                    if len(set(CO2_NAME) & set(rxn["products_mets"])) != 0:
                        if self.flux < 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"][
                                "Most reactions cannot proceed in the direction of CO2 fixation, except for naturally carbon fixation reactions and and reactions utilizing CO2 and high-energy substrates such as phosphoenolpyruvate (PEP) and ATP"
                            ] = "true"

    def find_ATP_synthase_rxn(self, model_info):
        """
        Find out the list of ATP synthase reactions
        """
        ATP_synthase_rxn = []
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if sorted(rxn["all_mets"]) == sorted(ATP_SYNTHASE):
                    ATP_synthase_rxn.append(rxn["id"])
        return ATP_synthase_rxn

    def find_ATP_rxn(self, model_info):
        """
        Find reactions that do not comply with ATP rules
        """
        ATP_synthase_rxn = self.find_ATP_synthase_rxn(model_info)
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if (
                    rxn["id"] not in model_info["exchange_rxns"]
                    and rxn["id"] not in model_info["transport_rxns"]
                    and rxn["id"] not in ATP_synthase_rxn
                ):
                    if len(set(XTP) & set(rxn["reactants_mets"])) != 0:
                        if len(set(XPI) & set(rxn["products_mets"])) != 0:
                            if (
                                len(set(rxn["products_mets"]) & set(YLCOA)) == 0
                                and self.flux < 0
                            ):
                                rxn["rule_conflict"] = "true"
                                rxn["rules2"][
                                    "Most reactions with ATP and phosphate/diphosphate/triphosphate are in the direction of ATP consumption ATP synthase reaction"
                                ] = "true"
                    if len(set(XTP) & set(rxn["products_mets"])) != 0:
                        if len(set(XPI) & set(rxn["reactants_mets"])) != 0:
                            if (
                                len(set(rxn["reactants_mets"]) & set(YLCOA)) == 0
                                and self.flux > 0
                            ):
                                rxn["rule_conflict"] = "true"
                                rxn["rules2"][
                                    "Most reactions with ATP and phosphate/diphosphate/triphosphate are in the direction of ATP consumption ATP synthase reaction"
                                ] = "true"

    def find_proton_rxn(self, model_info, model):
        """
        find the proton transport reaction
        """
        h_close = []
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                h_coefficient, hm = [], []
                rxns = model.reactions.get_by_id(rxn["id"])
                h_met = [met.id for met in rxns.metabolites if met.name in H_NAME]
                if (
                    len(rxn["all_mets"]) == 6
                    and set(rxn["all_mets"]).issubset(ATP_SYNTHASE)
                    and len(h_met) == 2
                ) or (
                    len(rxn["all_mets"]) == 7
                    and set(rxn["all_mets"]).issubset(ATP_SYNTHASE)
                ):
                    for met in h_met:
                        h_coefficient.append(abs(rxns.get_coefficient(met)))
                    for mets in rxns.metabolites:
                        if mets.name in ATP_NAME:
                            atp_coefficient = abs(rxns.get_coefficient(mets.id))
                            break
                    try:
                        if any((i / atp_coefficient) > 2 for i in h_coefficient):
                            pass
                    except:
                        print("find_proton_rxn:..............", rxn)

                    if any((i / atp_coefficient) > 2 for i in h_coefficient):
                        print(rxn["rxn_exp_name"], rxn["rxn_exp_id"])
                        if len(set(ATP_NAME) & set(rxn["reactants_mets"])) != 0:
                            if self.flux > 0:
                                rxn["rule_conflict"] = "true"
                                rxn["rules2"][
                                    "The direction of the reaction for proton-driven generation of ATP is incorrect"
                                ] = "true"
                            else:
                                rxn["rules"]["right_chain_rxn"] = "true"
                                rxn["rules2"]["right_chain_rxn"] = "true"
                            h_met.reverse()
                            hm = h_met

                        if len(set(ATP_NAME) & set(rxn["products_mets"])) != 0:
                            if self.flux < 0:
                                rxn["rule_conflict"] = "true"
                                rxn["rules2"][
                                    "The direction of the reaction for proton-driven generation of ATP is incorrect"
                                ] = "true"
                            else:
                                rxn["rules"]["right_chain_rxn"] = "true"
                                rxn["rules2"]["right_chain_rxn"] = "true"
                            hm = h_met
                        h_close.append(hm[0])
                        print("rules2-h_close", h_close)
                        for now_rxn in model.reactions:
                            rxn_met_c = [met.id for met in now_rxn.metabolites]
                            r_mets_id = [met.id for met in now_rxn.reactants]
                            pro_mets_id = [met.id for met in now_rxn.products]
                            if len(set(hm) & set(rxn_met_c)) == 2:
                                if now_rxn.id in model_info["transport_rxns"]:
                                    if (
                                        hm[0] in r_mets_id
                                        and hm[1] in pro_mets_id
                                        and (
                                            now_rxn.lower_bound < 0
                                            and now_rxn.upper_bound == 0
                                        )
                                    ):
                                        rxn["rule_conflict"] = "true"
                                        rxn["rules2"][
                                            "The direction of the reaction for proton-driven generation of ATP is incorrect"
                                        ] = "true"
                                    elif (
                                        hm[0] in pro_mets_id
                                        and hm[1] in r_mets_id
                                        and (
                                            now_rxn.lower_bound == 0
                                            and now_rxn.upper_bound > 0
                                        )
                                    ):
                                        rxn["rule_conflict"] = "true"
                                        rxn["rules2"][
                                            "The direction of the reaction for proton-driven generation of ATP is incorrect"
                                        ] = "true"

    def find_Sugar_hydrolysis_rxn(self, model_info, model):
        """
        Find reactions that don't follow the rules for polysaccharides and glycogen
        """
        sugar = []
        for sugar_met in SUGAR:
            for met in model.metabolites:
                if sugar_met in met.name:
                    sugar.append(met.name)

        sugar = list(set(sugar) - set(EXCEPT_SUGAR))

        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if (
                    rxn["id"] not in model_info["exchange_rxns"]
                    and rxn["id"] not in model_info["transport_rxns"]
                ):
                    if (
                        len(set(sugar) & set(rxn["reactants_mets"])) == 1
                        and len(set(H2O_NAME) & set(rxn["reactants_mets"])) != 0
                    ):
                        if self.flux < 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"][
                                "Polysaccharide hydrolysis reaction is the direction of hydrolysis"
                            ] = "true"
                    if (
                        len(set(sugar) & set(rxn["products_mets"])) == 1
                        and len(set(H2O_NAME) & set(rxn["products_mets"])) != 0
                    ):
                        if self.flux > 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"][
                                "Polysaccharide hydrolysis reaction is the direction of hydrolysis"
                            ] = "true"

    def find_PTS_transport(self, model_info):
        """
        Find reactions that do not correspond to the PTS pathway for material transport. PTS pathway for material transport, irreversible
        """
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if "PTS" in rxn["name"] and (
                    rxn["bounds"][0] < 0 and rxn["bounds"][1] > 0
                ):
                    rxn["rule_conflict"] = "true"
                    rxn["rules2"][
                        "Phosphotransferase reactions (PTS) are irreversible"
                    ] = "true"

    def find_ppi_h2o(self, model_info):
        """
        Hydrolysis of polyphosphoric acid to generate pi ppi
        """
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if rxn["id"] not in model_info["exchange_rxns"]:
                    if len(set(H2O_NAME) & set(rxn["reactants_mets"])) != 0:
                        for r_m in rxn["reactants_mets"]:
                            if (
                                re.search(
                                    r"phosphate|ppi|triphosphate|diphosphate|bisphosphate",
                                    r_m.lower(),
                                )
                                and self.flux < 0
                            ):
                                if (
                                    len(set(ALL_PI_NAME) & set(rxn["products_mets"]))
                                    != 0
                                ):
                                    rxn["rule_conflict"] = "true"
                                    rxn["rules2"][
                                        "The hydrolysis reaction of polyphosphoric acid is in the direction of hydrolysis"
                                    ] = "true"
                    if len(set(H2O_NAME) & set(rxn["products_mets"])) != 0:
                        for p_m in rxn["products_mets"]:
                            if (
                                re.search(
                                    r"phosphate|ppi|triphosphate|diphosphate|bisphosphate",
                                    p_m.lower(),
                                )
                                and self.flux > 0
                            ):
                                if (
                                    len(set(ALL_PI_NAME) & set(rxn["reactants_mets"]))
                                    != 0
                                ):
                                    rxn["rule_conflict"] = "true"
                                    rxn["rules2"][
                                        "The hydrolysis reaction of polyphosphoric acid is in the direction of hydrolysis"
                                    ] = "true"

    def find_acyl_h2o(self, model_info):
        """
        Acyl-CoA is hydrolyzed to produce CoA
        """
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if rxn["id"] not in model_info["exchange_rxns"]:
                    if (
                        len(set(H2O_NAME) & set(rxn["reactants_mets"])) != 0
                        and len(set(YLCOA) & set(rxn["reactants_mets"])) != 0
                        and len(rxn["reactants_mets"]) == 2
                        and len(set(COA) & set(rxn["products_mets"])) != 0
                    ):
                        if self.flux < 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"][
                                "The hydrolysis of acyl-CoA to form CoA is in the direction of hydrolysis"
                            ] = "true"
                    if (
                        len(set(H2O_NAME) & set(rxn["products_mets"])) != 0
                        and len(set(YLCOA) & set(rxn["products_mets"])) != 0
                        and len(rxn["products_mets"]) == 2
                        and len(set(COA) & set(rxn["reactants_mets"])) != 0
                    ):
                        if self.flux > 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"][
                                "The hydrolysis of acyl-CoA to form CoA is in the direction of hydrolysis"
                            ] = "true"

    def find_ac_acid(self, model_info):
        """
        In addition to the reaction of acetic acid and high-energy substances (atp, acyl-CoA) substances, the reaction containing acetic acid is the direction of generating acetic acid
        """
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if (
                    rxn["id"] not in model_info["exchange_rxns"]
                    and rxn["id"] not in model_info["transport_rxns"]
                ):
                    if len(set(AC) & set(rxn["products_mets"])) != 0:
                        if (
                            len(set(ATP_NAME) & set(rxn["products_mets"])) == 0
                            and len(set(rxn["products_mets"]) & set(YLCOA)) == 0
                        ):
                            if self.flux < 0:
                                rxn["rule_conflict"] = "true"
                                rxn["rules2"][
                                    "The reaction containing acetate is in the direction of producing acetate except for reactions with another high energy metabolite such as ATP, Acyl-CoA"
                                ] = "true"
                    if len(set(AC) & set(rxn["reactants_mets"])) != 0:
                        if (
                            len(set(ATP_NAME) & set(rxn["reactants_mets"])) == 0
                            and len(set(rxn["reactants_mets"]) & set(YLCOA)) == 0
                        ):
                            if self.flux > 0:
                                rxn["rule_conflict"] = "true"
                                rxn["rules2"][
                                    "The reaction containing acetate is in the direction of producing acetate except for reactions with another high energy metabolite such as ATP, Acyl-CoA"
                                ] = "true"

    def find_h2o2_rxn(self, model_info, model):
        """
        When hydrogen peroxide reacts with reducing substances (nadh, ndph, Reduced glutathione, Reduced thioredoxin, Ferrocytochrome c-553), it is irreversible
        """
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if rxn["id"] not in model_info["exchange_rxns"]:
                    rxns = model.reactions.get_by_id(rxn["id"])
                    reactants_mets_name = [met.name for met in rxns.reactants]
                    products_mets_name = [met.name for met in rxns.products]
                    if (
                        len(set(H2O2_NAME) & set(reactants_mets_name)) != 0
                        and len(rxn["reactants_mets"]) == 2
                    ):
                        if (
                            len(set(NADH_NAME) & set(rxn["reactants_mets"])) != 0
                            or len(set(NADPH) & set(rxn["reactants_mets"])) != 0
                            or len(set(H2O2_REDUCED) & set(reactants_mets_name)) != 0
                        ):
                            if self.flux < 0:
                                rxn["rule_conflict"] = "true"
                                rxn["rules2"][
                                    "When hydrogen peroxide reacts with reducing substances such as NADH, NADPH, reduced glutathione, reduced thioredoxin, and ferrocytochrome c-553, the reaction is irreversible"
                                ] = "true"
                    if (
                        len(set(H2O2_NAME) & set(products_mets_name)) != 0
                        and len(rxn["products_mets"]) == 2
                    ):
                        if (
                            len(set(NADH_NAME) & set(rxn["products_mets"])) != 0
                            or len(set(NADPH) & set(rxn["products_mets"])) != 0
                            or len(set(H2O2_REDUCED) & set(products_mets_name)) != 0
                        ):
                            if self.flux > 0:
                                rxn["rule_conflict"] = "true"
                                rxn["rules2"][
                                    "When hydrogen peroxide reacts with reducing substances such as NADH, NADPH, reduced glutathione, reduced thioredoxin, and ferrocytochrome c-553, the reaction is irreversible"
                                ] = "true"

    def find_fe3_fe2_rxn(self, model_info):
        """
        Fe3 reacts with nadh and fadh2 to Fe2, which is irreversible
        """
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if rxn["id"] not in model_info["exchange_rxns"]:
                    if (
                        len(set(FE3) & set(rxn["reactants_mets"])) != 0
                        and len(set(FE2) & set(rxn["products_mets"])) != 0
                    ):
                        if (
                            len(set(NADH_NAME) & set(rxn["reactants_mets"])) != 0
                            or len(set(FADH2) & set(rxn["reactants_mets"])) != 0
                        ):
                            if self.flux < 0:
                                rxn["rule_conflict"] = "true"
                                rxn["rules2"][
                                    "The reduction of Fe3 to Fe2 by nadh and fadh2 is irreversible"
                                ] = "true"

                    if (
                        len(set(FE3) & set(rxn["products_mets"])) != 0
                        and len(set(FE2) & set(rxn["reactants_mets"])) != 0
                    ):
                        if (
                            len(set(NADH_NAME) & set(rxn["products_mets"])) != 0
                            or len(set(FADH2) & set(rxn["products_mets"])) != 0
                        ):
                            if self.flux > 0:
                                rxn["rule_conflict"] = "true"
                                rxn["rules2"][
                                    "The reduction of Fe3 to Fe2 by nadh and fadh2 is irreversible"
                                ] = "true"

    def find_ate_met(self, model):
        """
        Acquire acid metabolites
        """
        ate_met = []
        for met in model.metabolites:
            if (
                met.name.endswith("ate")
                or met.name.endswith("Acid")
                or met.name.endswith("acid")
            ):
                ate_met.append(met.name)
        return ate_met

    def find_aldehyde_met(self, model):
        """
        Access to aldehyde metabolites
        """
        aldehyde_met = []
        for met in model.metabolites:
            if (
                met.name.endswith("aldehyde")
                or met.name.endswith("anal")
                or re.search(r"erythrose", met.name.lower())
            ):
                aldehyde_met.append(met.name)
        return aldehyde_met

    def find_aldehyde_ate_rxn(self, model_info, model):
        """
        Acyl-CoA is hydrolyzed to produce CoA
        """
        ate_met = self.find_ate_met(model)
        aldehyde_met = self.find_aldehyde_met(model)
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if rxn["id"] not in model_info["exchange_rxns"]:
                    if len(set(NAD_NADP_FDXO) & set(rxn["reactants_mets"])) != 0:
                        if (
                            len(set(rxn["reactants_mets"]) & set(aldehyde_met)) != 0
                            and len(set(H2O_NAME) & set(rxn["reactants_mets"])) != 0
                            and len(set(rxn["products_mets"]) & set(ate_met)) != 0
                            and len(set(O2_NAME) & set(rxn["products_mets"])) == 0
                        ):
                            if self.flux < 0:
                                rxn["rule_conflict"] = "true"
                                rxn["rules2"][
                                    "The oxidation of aldehydes to acids by NA(D)P and ferredoxin is irreversible"
                                ] = "true"
                    if len(set(NAD_NADP_FDXO) & set(rxn["products_mets"])) != 0:
                        if (
                            len(set(rxn["products_mets"]) & set(aldehyde_met)) != 0
                            and len(set(H2O_NAME) & set(rxn["products_mets"])) != 0
                            and len(set(rxn["reactants_mets"]) & set(ate_met)) != 0
                            and len(set(O2_NAME) & set(rxn["reactants_mets"])) == 0
                        ):
                            if self.flux > 0:
                                rxn["rule_conflict"] = "true"
                                rxn["rules2"][
                                    "The oxidation of aldehydes to acids by NA(D)P and ferredoxin is irreversible"
                                ] = "true"

    def find_aldehyde_ate_rxns(self, model_info, model):
        """
        Hydrolysis of aldehydes to acids and alcohols
        """
        anol_met = []
        for met in model.metabolites:
            if met.name.endswith("anol"):
                anol_met.append(met.name)
        ate_met = self.find_ate_met(model)
        aldehyde_met = self.find_aldehyde_met(model)

        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if rxn["id"] not in model_info["exchange_rxns"]:
                    if (
                        len(set(H2O_NAME) & set(rxn["reactants_mets"])) != 0
                        and len(set(rxn["reactants_mets"]) & set(aldehyde_met)) != 0
                        and len(set(rxn["products_mets"]) & set(ate_met)) != 0
                        and len(set(rxn["products_mets"]) & set(anol_met)) != 0
                    ):
                        if self.flux < 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"][
                                "The hydrolysis of aldehydes to acids and alcohols is an irreversible reaction."
                            ] = "true"

                    if (
                        len(set(H2O_NAME) & set(rxn["products_mets"])) != 0
                        and len(set(rxn["products_mets"]) & set(aldehyde_met)) != 0
                        and len(set(rxn["reactants_mets"]) & set(ate_met)) != 0
                        and len(set(rxn["reactants_mets"]) & set(anol_met)) != 0
                    ):
                        if self.flux > 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"][
                                "The hydrolysis of aldehydes to acids and alcohols is an irreversible reaction."
                            ] = "true"

    def find_sugar_pi_rxn(self, model_info):
        """
        Phosphorylated sugar hydrolysis reaction, irreversible
        """
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if rxn["id"] not in model_info["exchange_rxns"]:
                    if (
                        len(set(H2O_NAME) & set(rxn["reactants_mets"])) != 0
                        and len(set(SUGAR_PHOSPHATE) & set(rxn["reactants_mets"])) != 0
                        and len(rxn["reactants_mets"]) == 2
                    ):
                        if self.flux < 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"][
                                "The hydrolysis of phosphorylated sugars is an irreversible reaction"
                            ] = "true"

                    if (
                        len(set(H2O_NAME) & set(rxn["products_mets"])) != 0
                        and len(set(SUGAR_PHOSPHATE) & set(rxn["products_mets"])) != 0
                        and len(rxn["products_mets"]) == 2
                    ):
                        if self.flux > 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"][
                                "The hydrolysis of phosphorylated sugars is an irreversible reaction"
                            ] = "true"

    def find_glutamate_synthesis(self, model_info):
        """Identify reactions that do not follow the rules for glutamate synthesis"""
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                direction = glutamate_synthesis_direction(rxn)
                if direction is not None and self.flux * direction < -1e-9:
                    rxn["rule_conflict"] = "true"
                    rxn["rules2"]["glutamate_synthesis_rxn"] = "true"

    def find_quino_rxn(self, model_info, model):
        """Find reactions that do not follow the rules for reactions of quinones"""
        quino_mets = []
        for met in model.metabolites:
            if re.search(
                r"ubiquinone|ubiquinol|menaquinone|menaquinol|plastoquinol|plastoquinone|2-demethylmenaquinol|2-demethylmenaquinone|flavin adenine dinucleotide|ubiquinone-8|qh2|ubiquinol-\d|mql7|menaquinone \d|menaquinol|2-demethylmenaquinone 8|2-demethylmenaquinol 8|2-demethylmenaquinol8|2dmmq7|fad",
                met.name.lower(),
            ):
                quino_mets.append(met.name)
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if (
                    len(set(rxn["reactants_mets"]) & set(quino_mets)) != 0
                    and len(set(NADH_NAME) & set(rxn["reactants_mets"])) != 0
                    and len(set(NAD_NAME) & set(rxn["products_mets"])) != 0
                ):
                    if self.flux < 0:
                        rxn["rule_conflict"] = "true"
                        rxn["rules2"][
                            "The reduction of ubiquinones via NADH is an irreversible process."
                        ] = "true"
                if (
                    len(set(rxn["products_mets"]) & set(quino_mets)) != 0
                    and len(set(NADH_NAME) & set(rxn["products_mets"])) != 0
                    and len(set(NAD_NAME) & set(rxn["reactants_mets"])) != 0
                ):
                    if self.flux > 0:
                        rxn["rule_conflict"] = "true"
                        rxn["rules2"][
                            "The reduction of ubiquinones via NADH is an irreversible process."
                        ] = "true"

    def find_respiratory_chain_rxn(self, model_info, model):
        """Get Respiratory Chain Reaction"""
        for rxn in model_info["reactions"]:
            # self.flux belongs only to self.rxnId, not every proton carrier.
            if rxn["id"] != self.rxnId:
                continue
            rxns = model.reactions.get_by_id(rxn["id"])
            if (
                respiratory_membership_supported(rxn["all_mets"])
                and len(set(rxn["reactants_mets"]) & set(H_NAME)) != 0
                and len(set(rxn["products_mets"]) & set(H_NAME)) != 0
            ):
                if not (
                    len(set(rxn["reactants_mets"]) & set(NADH_NAME)) != 0
                    and len(set(rxn["products_mets"]) & set(NADPH)) != 0
                ):
                    if any(
                        met.compartment in C_COMPARTMENT
                        for met in rxns.reactants
                        if met.name in H_NAME
                    ):
                        if self.flux < 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"]["respiratory_chian_rxn"] = "true"
                        else:
                            rxn["rules"]["right_chain_rxn"] = "true"
                            rxn["rules2"]["right_chain_rxn"] = "true"
                    if any(
                        met.compartment in C_COMPARTMENT
                        for met in rxns.products
                        if met.name in H_NAME
                    ):
                        if self.flux > 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"]["respiratory_chian_rxn"] = "true"
                        else:
                            rxn["rules"]["right_chain_rxn"] = "true"
                            rxn["rules2"]["right_chain_rxn"] = "true"

    def find_respiratory_chain_rxn2(self, model_info):
        """
        Get Respiratory Chain Reaction
        """
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"]:
                if len(rxn["all_mets"]) == 7 and set(rxn["all_mets"]).issubset(
                    ATP_SYNTHASE
                ):
                    if len(set(ATP_NAME) & set(rxn["reactants_mets"])) != 0:
                        if self.flux > 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"]["respiratory_chian_rxn"] = "true"
                        else:
                            rxn["rules"]["right_chain_rxn"] = "true"
                            rxn["rules2"]["right_chain_rxn"] = "true"
                    if len(set(ATP_NAME) & set(rxn["products_mets"])) != 0:
                        if self.flux < 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"]["respiratory_chian_rxn"] = "true"
                        else:
                            rxn["rules"]["right_chain_rxn"] = "true"
                            rxn["rules2"]["right_chain_rxn"] = "true"

    def get_metacyc_rxn(self, rxn, l_list, r_list, lower_dict, upper_dict, dfId, metName, met):
        """Compatibility wrapper using the actual cycle flux direction."""
        from mqc.control.metacyc_direction import reference_direction, violates
        allowed, _, _ = reference_direction(rxn, lower_dict[dfId], upper_dict[dfId],
                                             ' + '.join(l_list) + ' --> ' + ' + '.join(r_list))
        return int(violates(rxn, {'allowed': allowed}, self.flux))

    def find_metacyc_bounds_rxn(self, model_info):
        """Only a forbidden current pFBA direction sets rule_conflict."""
        from mqc.control.metacyc_direction import apply_rule
        apply_rule(model_info, self.cfg.metacyc_rules, self.rxnId, self.flux)

    def find_energy_reducing_power_rxn(self, model_info, model):
        """Find reactions that do not follow the rules for reactions of energy or reducing power"""
        energy_reducing_power_mets = []
        for met in model.metabolites:
            if re.search(r"hydrogen|phosphonate", met.name.lower()):
                energy_reducing_power_mets.append(met.name)
        for rxn in model_info["reactions"]:
            if self.rxnId == rxn["id"] and rxn["id"] in model_info["exchange_rxns"]:
                if (
                    len(set(rxn["reactants_mets"]) & set(energy_reducing_power_mets))
                    != 0
                ):
                    if self.flux < 0:
                        rxn["rule_conflict"] = "true"
                        rxn["rules2"]["energy_reducing_power"] = "true"
                if (
                    len(set(rxn["products_mets"]) & set(energy_reducing_power_mets))
                    != 0
                ):
                    if self.flux > 0:
                        rxn["rule_conflict"] = "true"
                        rxn["rules2"]["energy_reducing_power"] = "true"

    def find_modelseed_error_bound_rxn(self, model_info):
        """"""
        df = pd.read_excel(self.cfg.modelseed_rxn)
        df_name = df.iloc[:, 0]
        df_lower_bound = df.iloc[:, 3]
        df_upper_bound = df.iloc[:, 4]
        lower_dict = dict(zip(df_name, df_lower_bound))
        upper_dict = dict(zip(df_name, df_upper_bound))
        for rxn in model_info["reactions"]:
            if rxn["id"] != self.rxnId:
                continue
            for dfId in list(df_name):
                if rxn["id"] == dfId:
                    if lower_dict[dfId] == -1000 and upper_dict[dfId] == 0:
                        if self.flux > 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"][
                                "the boundary reversibility of modelseed does not match metecyc"
                            ] = "true"
                    if lower_dict[dfId] == 0 and upper_dict[dfId] == 1000:
                        if self.flux < 0:
                            rxn["rule_conflict"] = "true"
                            rxn["rules2"][
                                "the boundary reversibility of modelseed does not match metecyc"
                            ] = "true"

    def get_all_rules(self, model_info, model):
        """"""
        self.initial_rule(model_info)
        self.find_o2s(model_info)
        O2_rxn = self.find_O2_rxn(model_info)
        self.find_nh3_nh4_rxn(model_info, O2_rxn)
        self.find_CO2_rxn(model_info)
        self.find_ATP_rxn(model_info)
        self.find_proton_rxn(model_info, model)
        self.find_Sugar_hydrolysis_rxn(model_info, model)
        self.find_PTS_transport(model_info)
        self.find_ppi_h2o(model_info)
        self.find_acyl_h2o(model_info)
        self.find_ac_acid(model_info)
        self.find_h2o2_rxn(model_info, model)
        self.find_fe3_fe2_rxn(model_info)
        self.find_aldehyde_ate_rxn(model_info, model)
        self.find_aldehyde_ate_rxns(model_info, model)
        self.find_sugar_pi_rxn(model_info)
        self.find_glutamate_synthesis(model_info)
        self.find_quino_rxn(model_info, model)
        self.find_respiratory_chain_rxn(model_info, model)
        self.find_respiratory_chain_rxn2(model_info)
        self.find_metacyc_bounds_rxn(model_info)
        self.find_energy_reducing_power_rxn(model_info, model)
        self.find_modelseed_error_bound_rxn(model_info)
        self.find_biotin_synthesis_direction(model_info, model)
        self.find_deacetylase_direction(model_info, model)
        self.find_desulfurase_direction(model_info, model)
        self.find_incomplete_quinone_methylation(model_info, model)
        self.find_sugar_identity_conflicts(model_info, model)
