#!/usr/bin/env python3
"""Rebuild all standalone MQC publication panels from packaged CSV data."""

from __future__ import annotations

import fig2_panels
import fig3_panels
import fig5_panels
import fig6_panels
import fig4_panels
import compose_full_figures


def main() -> None:
    fig2_panels.main()
    fig3_panels.main()
    fig4_panels.main()
    fig5_panels.main()
    fig6_panels.main()
    compose_full_figures.main()
    print("Complete: wrote 41 standalone panels and five combined figures")


if __name__ == "__main__":
    main()
