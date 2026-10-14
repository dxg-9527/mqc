from drawing import MSO_SHAPE
from evidence import initial, flux


def frame(d):
    p = d.group("00_titles_and_panel_rules")
    for x in [3.08, 12.52]:
        d.line(p, x, 0.29, x, 9.50, color="rule", width=0.75, arrow=False)
    d.line(p, 3.27, 5.15, 12.34, 5.15, color="rule", width=0.7, arrow=False)
    d.text(p, 0.35, 0.25, 0.27, 0.32, "a", 20, bold=True)
    d.text(p, 0.66, 0.25, 2.17, 0.50, "No-carbon\nphenotype screen", 17, bold=True)
    d.text(p, 3.28, 0.25, 0.27, 0.32, "b", 20, bold=True)
    d.text(
        p,
        3.65,
        0.25,
        8.22,
        0.30,
        "Active pathways and targeted constraints",
        17,
        bold=True,
    )
    d.text(p, 3.27, 0.66, 4.47, 0.40, "(1) FDR-dependent redox cycle", 16, bold=True)
    d.text(p, 7.93, 0.66, 4.47, 0.40, "(2) ABAH-dependent redox cycle", 16, bold=True)
    d.text(
        p, 3.27, 5.30, 4.47, 0.34, "(3) Proton-supported ATP synthesis", 16, bold=True
    )
    d.text(
        p,
        7.93,
        5.30,
        4.47,
        0.34,
        "(4) Unsupported metabolite production",
        16,
        bold=True,
    )


def panel_a(d):
    p = d.group("a_condition_model")
    d.shape(p, MSO_SHAPE.HEXAGON, 0.53, 2.41, 0.44, 0.42, "blue", 1.1, "pale_blue")
    d.text(p, 0.38, 2.89, 0.76, 0.47, "carbon\nsource", 15.5)
    d.line(p, 1.02, 2.62, 1.51, 2.62, width=1.8)
    d.cross(p, 1.25, 2.62, 0.087)
    pts = [
        (1.81, 2.34),
        (2.13, 2.13),
        (2.48, 2.32),
        (2.69, 2.74),
        (2.40, 3.04),
        (2.02, 3.12),
        (1.69, 2.83),
        (2.17, 2.66),
    ]
    for i, j in [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 6),
        (6, 0),
        (0, 7),
        (1, 7),
        (2, 7),
        (3, 7),
        (4, 7),
        (5, 7),
        (6, 7),
        (2, 4),
    ]:
        d.line(p, *pts[i], *pts[j], width=0.85, arrow=False)
    for x, y in pts:
        d.circle(p, x, y, 0.045, width=1)
    d.text(p, 1.64, 3.23, 1.05, 0.29, "iCN718", 17, bold=True)
    d.text(p, 0.39, 3.65, 2.55, 0.22, "Carbon-source uptake closed", 13.5, "muted")
    d.text(p, 0.51, 4.00, 2.30, 0.32, "Objective screening", 15.5, "blue")
    d.line(p, 1.66, 4.35, 1.66, 4.68, color="blue", width=1.25)
    d.shift_group("a_condition_model", -1.18)
    q = d.group("a_objectives")
    d.text(q, 0.40, 3.72, 2.53, 0.40, "Positive objectives", 16, bold=True)
    d.text(q, 0.40, 4.30, 2.53, 0.30, "Reducing equivalents", 15.5, "redox", bold=True)
    d.circle(q, 0.67, 4.99, 0.060, "redox")
    d.text(q, 0.83, 4.82, 0.73, 0.30, "NADH", 15.5, "ink", bold=True, align="left")
    d.text(q, 1.81, 4.82, 1.1, 0.32, initial("NADH"), 15.5, bold=True, align="left")
    d.circle(q, 0.67, 5.96, 0.06, "energy")
    d.text(q, 0.83, 5.79, 0.73, 0.30, "NTPs", 15.5, "ink", bold=True, align="left")
    d.text(
        q,
        1.53,
        5.79,
        1.41,
        0.32,
        initial("ATP") + " each",
        15.5,
        bold=True,
        align="right",
    )
    d.text(
        q,
        0.50,
        6.17,
        2.42,
        0.45,
        "ATP / CTP / GTP\nUTP / ITP",
        13.5,
        "muted",
        align="left",
    )
    d.circle(q, 0.75, 7.05, 0.055, "green")
    d.text(q, 0.45, 7.20, 0.74, 0.27, "Urea", 15.5, "ink", bold=True)
    d.text(q, 0.45, 7.50, 0.74, 0.24, initial("urea_c"), 15.5, bold=True)
    d.circle(q, 2.13, 7.05, 0.055, "purple")
    d.text(q, 1.45, 7.20, 1.38, 0.27, "Allophanate", 15.5, "ink", bold=True)
    d.text(q, 1.74, 7.50, 0.77, 0.24, initial("allphn_c"), 15.5, bold=True)
    d.path(
        q,
        [("M", 0.75, 7.91), ("L", 0.75, 8.04), ("L", 2.13, 8.04), ("L", 2.13, 7.91)],
        color="rule",
        width=0.8,
        arrow=False,
    )
    d.line(q, 1.44, 8.04, 1.44, 8.29, color="rule", width=0.8, arrow=False)
    d.text(q, 0.45, 8.39, 2.43, 0.29, "Trace active pathways", 15.5, "blue", bold=True)
    d.line(q, 1.03, 8.85, 3.24, 8.85, width=1.65)
    d.text(q, 0.45, 9.12, 2.43, 0.23, "mmol gDCW⁻¹ h⁻¹", 13.5, "muted", italic=True)


def panel_b1(d):
    p = d.group("b1_folate_cycle")
    cx = 5.51
    top = 2.54
    bot = 3.40
    d.path(
        p,
        [
            ("M", cx + 0.065, top),
            ("C", 6.83, top + 0.04, 6.83, bot - 0.06, cx + 0.09, bot),
        ],
        width=2,
        name="MTHFR3 forward",
    )
    d.path(
        p,
        [
            ("M", cx - 0.065, bot),
            ("C", 4.19, bot - 0.04, 4.19, top + 0.06, cx - 0.09, top),
        ],
        width=2,
        name="MTHFR2 reverse",
    )
    for yy in [top, bot]:
        d.circle(p, cx, yy, 0.058)
    d.text(p, 4.16, 2.16, 2.70, 0.30, "5,10-methylene-THF", 15.5, bold=True)
    d.text(p, 4.46, 3.46, 2.10, 0.30, "5-methyl-THF", 15.5, bold=True)
    d.text(p, 6.52, 2.69, 1.18, 0.31, "MTHFR3", 15, "blue", bold=True, italic=True)
    d.text(p, 6.52, 3.03, 1.22, 0.45, "NADPH\n→ NADP⁺", 13.5, "muted")
    d.text(
        p,
        3.30,
        2.58,
        1.18,
        0.46,
        "MTHFR2\n(reverse)",
        15,
        "blue",
        bold=True,
        italic=True,
    )
    d.text(p, 3.28, 3.05, 1.30, 0.43, "NAD⁺\n→ NADH", 13.5, "muted")
    q = d.group("b1_rubredoxin_cycle")
    top = 4.03
    bot = 4.91
    d.path(
        q,
        [
            ("M", cx + 0.065, top),
            ("C", 6.80, top + 0.04, 6.80, bot - 0.04, cx + 0.09, bot),
        ],
        color="orange",
        width=2.3,
        name="FDR reverse (closed)",
    )
    d.path(
        q,
        [
            ("M", cx - 0.065, bot),
            ("C", 4.22, bot - 0.04, 4.22, top + 0.04, cx - 0.09, top),
        ],
        width=2,
        name="RDXR return route",
    )
    for yy in [top, bot]:
        d.circle(q, cx, yy, 0.058)
    d.text(q, 4.21, 3.72, 2.60, 0.24, "Reduced rubredoxin", 14.5, bold=True)
    d.text(q, 4.16, 4.99, 2.70, 0.30, "Oxidized rubredoxin", 15.5, bold=True)
    d.text(
        q,
        6.56,
        4.00,
        1.16,
        0.46,
        "FDR\n(reverse)",
        15,
        "orange",
        bold=True,
        italic=True,
    )
    d.text(q, 6.48, 4.75, 1.26, 0.45, "NADP⁺ + H⁺\n→ NADPH", 13.5, "muted")
    d.text(q, 3.31, 4.07, 1.14, 0.31, "RDXR", 15, "blue", bold=True, italic=True)
    d.text(q, 3.29, 4.39, 1.29, 0.45, "NADH →\nNAD⁺ + H⁺", 13.5, "muted")
    d.cross(q, 6.475, 4.48, 0.077)
    d.text(q, 6.65, 4.48, 1.0, 0.27, "closed", 13.5, "stop")
    d.text(
        q,
        3.40,
        5.27,
        4.20,
        0.24,
        "Net: H⁺ + NAD⁺ → NADH; output " + f'{float(flux("FDR", "ADD_NADH")):g}',
        13.5,
        bold=True,
    )
    d.shift_group("b1_folate_cycle", -1.10)
    d.shift_group("b1_rubredoxin_cycle", -0.55)


def panel_b2(d):
    p = d.group("b2_outer_cycle")
    d.rect(p, 8.40, 2.21, 3.55, 0.43, fill="pale_green", rounding=True)
    d.text(p, 8.43, 2.21, 3.49, 0.43, "2 ATP + NAD⁺ →\n2 ADP + 2 Pi + NADH", 13.5)
    lx, rx = 8.97, 11.20
    yt, yb = 3.20, 4.42
    d.path(
        p,
        [
            ("M", rx - 0.08, yt - 0.04),
            ("C", 10.71, 3.09, 9.49, 3.09, lx + 0.11, yt - 0.025),
        ],
        width=2,
        name="OMMBLHX3:2ommbl_c>2omhmbl_c",
    )
    d.path(
        p,
        [
            ("M", lx - 0.07, yt + 0.07),
            ("C", 8.41, 3.50, 8.43, 4.14, lx - 0.09, yb - 0.07),
        ],
        width=2,
        name="DMQMT:2omhmbl_c>q8h2_c",
    )
    d.path(
        p,
        [
            ("M", rx + 0.07, yb - 0.08),
            ("C", 11.70, 4.12, 11.69, 3.51, rx + 0.085, yt + 0.08),
        ],
        color="orange",
        width=2.35,
        name="ABAH:q8_c>2ommbl_c",
    )
    d.path(
        p,
        [("M", lx + 0.05, yb + 0.02), ("C", 9.45, 4.66, 9.65, 4.73, 9.80, 4.73)],
        width=2,
        name="CYTBD:q8h2_c>protein",
    )
    d.path(
        p,
        [("M", 10.49, 4.73), ("C", 10.70, 4.68, 10.91, 4.58, rx - 0.06, yb + 0.02)],
        width=2,
        name="CYTBD:protein>q8_c",
    )
    for xx, yy in [(lx, yt), (rx, yt), (lx, yb), (rx, yb)]:
        d.circle(p, xx, yy, 0.072, fill="pale_blue", width=1.2)
    d.text(p, 8.12, 2.87, 1.42, 0.29, "2OMHMBL", 15.5, bold=True)
    d.text(p, 10.83, 2.87, 1.44, 0.29, "2OMMBL", 15.5, bold=True)
    d.text(p, 9.31, 2.77, 1.71, 0.26, "OMMBLHX3", 14.5, "blue", bold=True, italic=True)
    d.text(p, 8.75, 3.57, 0.96, 0.29, "DMQMT", 13.5, "blue", bold=True, italic=True)
    d.text(
        p,
        11.61,
        3.23,
        0.79,
        0.47,
        "ABAH\n(reverse)",
        13.5,
        "orange",
        bold=True,
        italic=True,
    )
    d.cross(p, 11.55, 3.97, 0.075)
    d.text(p, 11.63, 3.84, 0.73, 0.28, "closed", 13.5, "stop")
    d.text(p, 8.12, 4.26, 0.74, 0.29, "Q₈H₂", 15.5, bold=True)
    d.text(p, 11.36, 4.26, 0.66, 0.26, "Q₈", 15.5, bold=True)
    q = d.group("b2_inner_cycle")
    d.path(
        q,
        [("M", 10.12, 4.13), ("C", 9.50, 4.12, 9.50, 3.46, 10.12, 3.45)],
        color="blue",
        width=1.4,
        name="DMQMT:amet_c>ahcys_c",
    )
    d.path(
        q,
        [("M", 10.20, 3.45), ("C", 10.82, 3.46, 10.82, 4.12, 10.20, 4.13)],
        color="orange",
        width=1.4,
        name="ABAH:ahcys_c>amet_c",
    )
    d.circle(q, 10.16, 3.45, 0.060, "energy", "pale_yellow")
    d.circle(q, 10.16, 4.13, 0.060, "energy", "pale_yellow")
    d.text(q, 9.38, 3.12, 1.56, 0.24, "AdoHcy + H⁺", 13.5, bold=True)
    d.text(q, 9.59, 4.20, 1.14, 0.27, "AdoMet", 14.5, bold=True)
    d.text(q, 9.52, 3.57, 1.28, 0.42, "Cofactor\nrecycling", 13.5, "muted", italic=True)
    m = d.group("b2_membrane")
    d.membrane(m, 8.68, 4.70, 3.18, 0.24, gap=(9.80, 10.49))
    d.rect(m, 9.80, 4.58, 0.69, 0.49, "blue", 1.15, "pale_blue", True)
    d.text(m, 9.70, 4.70, 0.89, 0.29, "CYTBD", 13.5, "blue", bold=True, italic=True)
    d.text(m, 8.15, 4.58, 0.47, 0.27, "c", 13.5, "muted")
    d.text(m, 8.15, 4.96, 0.47, 0.27, "e", 13.5, "muted")
    d.text(m, 8.63, 5.08, 3.69, 0.24, "½ O₂ + 2 H⁺(c) → H₂O + 2 H⁺(e)", 13.5, "muted")
    d.text(
        m,
        8.00,
        5.32,
        4.31,
        0.24,
        "NADH output: " + flux("ABAH", "ADD_NADH"),
        13.5,
        bold=True,
    )
    for name in ["b2_outer_cycle", "b2_inner_cycle", "b2_membrane"]:
        d.shift_group(name, -0.83)
