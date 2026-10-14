from drawing import ROOT, MSO_SHAPE
from evidence import initial, flux


def panel_b3(d):
    p = d.group("b3_membrane")
    d.membrane(p, 3.78, 7.10, 3.60, 0.37, gap=(5.61, 6.03))
    d.rect(p, 5.61, 6.96, 0.42, 0.68, "blue", 1.2, "pale_blue", True)
    d.rect(p, 5.73, 7.61, 0.19, 0.28, "blue", 1.0, "pale_blue", True)
    d.shape(p, MSO_SHAPE.OVAL, 5.30, 7.82, 1.05, 0.56, "blue", 1.2, "pale_blue")
    d.text(p, 5.33, 7.98, 0.99, 0.24, "ATPS4r", 15.5, "blue", bold=True, italic=True)
    d.text(p, 6.22, 6.74, 1.38, 0.28, "extracellular", 13.5, "muted")
    d.text(p, 3.72, 7.55, 1.36, 0.28, "cytoplasm", 13.5, "muted")
    q = d.group("b3_coupling")
    d.text(q, 3.64, 6.37, 1.04, 0.25, "EX_h_e", 15.5, "orange", bold=True, italic=True)
    d.line(
        q,
        4.12,
        6.66,
        5.71,
        6.66,
        color="orange",
        width=2.3,
        name="EX_h_e proton uptake (blocked)",
    )
    d.cross(q, 4.72, 6.66, 0.077)
    d.text(q, 4.04, 6.17, 1.56, 0.28, "uptake blocked", 13.5, "stop")
    d.text(q, 5.38, 6.35, 0.88, 0.26, "4 H⁺(e)", 15.5, "blue", bold=True)
    d.path(
        q,
        [
            ("M", 5.82, 6.67),
            ("L", 5.82, 7.37),
            ("C", 5.82, 7.67, 6.23, 7.68, 6.68, 7.68),
        ],
        width=2.15,
        name="Proton flow through ATP synthase",
    )
    d.text(q, 6.65, 7.53, 0.95, 0.24, "3 H⁺(c)", 15.5, "blue")
    # Two small input and output nodes connected around the catalytic head.
    d.path(
        q, [("M", 4.65, 8.03), ("C", 4.85, 8.03, 5.03, 8.01, 5.31, 8.07)], width=1.25
    )
    d.path(
        q, [("M", 4.65, 8.64), ("C", 4.90, 8.64, 5.00, 8.28, 5.31, 8.18)], width=1.25
    )
    d.path(
        q, [("M", 6.34, 8.07), ("C", 6.64, 8.00, 6.75, 8.03, 6.97, 8.03)], width=1.25
    )
    d.path(
        q, [("M", 6.34, 8.18), ("C", 6.58, 8.24, 6.69, 8.64, 6.97, 8.64)], width=1.25
    )
    for xx, yy, col in [
        (4.60, 8.03, "blue"),
        (4.60, 8.64, "blue"),
        (7.02, 8.03, "energy"),
        (7.02, 8.64, "blue"),
    ]:
        d.circle(q, xx, yy, 0.052, col)
    d.text(q, 3.84, 7.89, 0.63, 0.26, "ADP", 15.5, align="right")
    d.text(q, 3.98, 8.50, 0.49, 0.26, "Pi", 15.5, align="right")
    d.text(q, 7.13, 7.89, 0.57, 0.26, "ATP", 15.5, "ink", bold=True, align="left")
    d.text(q, 7.13, 8.50, 0.57, 0.26, "H₂O", 15.5, align="left")
    d.text(
        q,
        3.38,
        9.08,
        4.23,
        0.29,
        "Proton uptake supports NTP production",
        14.5,
        bold=True,
    )
    d.text(
        q,
        3.38,
        9.39,
        4.23,
        0.27,
        "ATPS4r flux: " + flux("EX_h_e", "ATPS4r"),
        14.5,
        bold=True,
    )
    for name in ["b3_membrane", "b3_coupling"]:
        d.shift_group(name, -0.20)


def panel_b4(d):
    p = d.group("b4_madh")
    d.text(p, 7.94, 6.13, 2.08, 0.49, "Methylmalonate\nsemialdehyde", 15.5, bold=True)
    d.text(p, 10.29, 6.51, 2.08, 0.31, "methylmalonate", 15.5, bold=True)
    d.line(
        p,
        8.65,
        6.93,
        11.50,
        6.93,
        color="orange",
        width=2.35,
        name="MADH faulty reaction (closed)",
    )
    d.circle(p, 8.60, 6.93, 0.062)
    d.circle(p, 11.57, 6.93, 0.062)
    d.text(p, 9.44, 6.63, 1.03, 0.27, "MADH", 15.5, "orange", bold=True, italic=True)
    d.text(
        p,
        9.36,
        7.04,
        1.18,
        0.24,
        "v = " + flux("MADH", "MADH", signed=True),
        14.5,
        "orange",
        bold=True,
    )
    d.cross(p, 10.76, 6.93, 0.08)
    d.text(p, 11.63, 7.01, 0.77, 0.25, "closed", 13.5, "stop")
    d.path(
        p, [("M", 8.90, 7.56), ("C", 9.20, 7.56, 9.27, 7.20, 9.48, 6.96)], width=1.25
    )
    d.circle(p, 8.84, 7.56, 0.039, "redox")
    d.text(p, 8.10, 7.43, 0.66, 0.25, "NAD⁺", 13.5, "muted", align="right")
    d.path(
        p,
        [("M", 10.56, 6.96), ("C", 10.42, 7.31, 10.42, 7.62, 10.63, 7.68)],
        width=1.25,
    )
    d.circle(p, 10.69, 7.69, 0.039, "redox")
    d.text(p, 10.81, 7.57, 0.65, 0.24, "NADH", 13.5, "muted", align="left")
    d.path(
        p, [("M", 10.89, 6.96), ("C", 10.96, 7.24, 11.28, 7.32, 11.38, 7.43)], width=1.1
    )
    d.circle(p, 11.42, 7.47, 0.036)
    d.text(p, 11.55, 7.34, 0.59, 0.25, "CO₂", 13.5, "muted", align="left")
    d.rect(p, 8.16, 7.95, 3.02, 0.31, fill="pale_orange", rounding=True)
    d.text(p, 8.18, 7.95, 2.98, 0.31, "Mass / charge imbalance", 13.5, "orange")
    q = d.group("b4_products")
    d.line(
        q,
        11.42,
        7.53,
        11.42,
        7.80,
        width=1.6,
        dash=True,
        arrow=True,
        name="CO2:co2_c>downstream",
    )
    d.path(
        q,
        [("M", 11.28, 8.25), ("C", 11.03, 8.27, 10.81, 8.33, 10.62, 8.37)],
        width=1.4,
        dash=True,
        arrow=False,
    )
    d.text(q, 11.08, 7.72, 1.24, 0.66, "CO₂ pool\n+ NH₄⁺", 13.5, "muted")
    d.circle(q, 10.62, 8.37, 0.025, "blue", "blue", width=0.6)
    d.path(q, [("M", 10.62, 8.37), ("C", 10.06, 8.45, 9.50, 8.60, 9.11, 8.94)], width=2)
    # ALPHNH is supported by model stoichiometry, not by an archived pFBA trace
    # for the same urea objective. Keep that distinction visible in the artwork.
    d.path(
        q,
        [("M", 10.62, 8.37), ("C", 11.15, 8.42, 11.57, 8.62, 11.71, 8.94)],
        width=2,
        dash=True,
        name="ALPHNH:model-derived separate-objective branch",
    )
    d.text(
        q,
        8.01,
        8.34,
        1.43,
        0.47,
        "UREA\n(reverse)",
        14.5,
        "blue",
        bold=True,
        italic=True,
    )
    d.text(
        q,
        11.48,
        8.24,
        0.92,
        0.47,
        "ALPHNH\n(reverse)",
        13.5,
        "blue",
        bold=True,
        italic=True,
    )
    d.circle(q, 9.06, 9.00, 0.06, "green")
    d.circle(q, 11.73, 9.00, 0.06, "purple")
    d.text(q, 8.56, 9.10, 1.00, 0.25, "Urea", 15.5, "ink", bold=True)
    d.text(q, 10.99, 9.10, 1.48, 0.25, "Allophanate", 15.5, "ink", bold=True)
    d.text(q, 8.59, 9.39, 0.94, 0.27, initial("urea_c"), 15.5, bold=True)
    d.text(q, 11.26, 9.39, 0.94, 0.27, initial("allphn_c"), 15.5, bold=True)
    for name in ["b4_madh", "b4_products"]:
        d.shift_group(name, -0.20)


def panel_c(d):
    p = d.group("c_data_driven_validation")
    pic = p.add_picture(
        str(ROOT / "Fig4_panel_c" / "output" / "Fig4_panel_c_post_QC_validation.png"),
        int(12.65 * 914400),
        int(0.41 * 914400),
        width=int(3.20 * 914400),
    )
    pic.name = "Post-MQC validation — regenerated from archived CSV data"
    # Put the panel identifier above the opaque plot image, not underneath it.
    d.text(p, 12.68, 0.25, 0.27, 0.32, "c", 20, bold=True)


def legend(d):
    p = d.group("legend")
    d.line(p, 0.35, 9.81, 15.85, 9.81, color="rule", width=0.75, arrow=False)
    y = 10.11
    d.line(p, 0.45, y, 0.80, y, width=1.8)
    d.text(p, 0.89, 9.96, 1.76, 0.32, "Active pFBA route", 13.5, align="left")
    d.line(p, 2.81, y, 3.14, y, color="orange", width=2.2)
    d.text(p, 3.23, 9.96, 1.93, 0.32, "Problematic reaction", 13.5, align="left")
    d.cross(p, 5.43, y, 0.064)
    d.text(p, 5.58, 9.96, 1.16, 0.32, "MQC stop", 13.5, align="left")
    d.circle(p, 7.01, y, 0.052)
    d.text(p, 7.16, 9.96, 1.14, 0.32, "Metabolite", 13.5, align="left")
    d.membrane(p, 8.57, 10.01, 0.34, 0.20, density=0.085)
    d.text(p, 9.03, 9.96, 1.21, 0.32, "Membrane", 13.5, align="left")
    d.line(p, 10.47, y, 10.80, y, width=1.5, dash=True)
    d.text(p, 10.91, 9.96, 2.07, 0.32, "Model-derived link", 13.5, align="left")
    d.text(
        p,
        13.04,
        9.94,
        2.62,
        0.46,
        "Flux: mmol gDCW⁻¹ h⁻¹\nc: cytoplasm; e: extracellular",
        13.5,
        "muted",
    )
