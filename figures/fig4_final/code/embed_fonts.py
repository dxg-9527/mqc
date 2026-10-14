"""Embed complete, freely embeddable Liberation Sans faces as Office EOT parts.

Format references:
https://www.w3.org/submissions/EOT/ (EOT version 0x00020001)
https://learn.microsoft.com/sv-se/openspecs/office_standards/ms-oe376/1663dabc-5d98-463f-889e-bcd9b77c3d34
The source font's license and fsType are preserved; no subsetting or modification.
"""

import struct, zipfile
from lxml import etree as ET
from drawing import FONTDIR, ROOT


def to_eot(path):
    b = path.read_bytes()
    u16 = lambda o: struct.unpack_from(">H", b, o)[0]
    u32 = lambda o: struct.unpack_from(">I", b, o)[0]
    tables = {
        b[12 + i * 16 : 16 + i * 16].decode("ascii"): (
            u32(20 + i * 16),
            u32(24 + i * 16),
        )
        for i in range(u16(4))
    }
    os2 = tables["OS/2"][0]
    head = tables["head"][0]
    nt = tables["name"][0]
    assert u16(os2 + 8) == 0, "Only installable embedding licensed fonts supported here"
    names = {}
    for i in range(u16(nt + 2)):
        plat, enc, lang, nid, nlen, noff = struct.unpack_from(">6H", b, nt + 6 + i * 12)
        if plat == 3 and lang == 0x409:
            raw = b[nt + u16(nt + 4) + noff : nt + u16(nt + 4) + noff + nlen]
            names[nid] = raw.decode("utf-16-be").encode("utf-16-le")
    hdr = bytearray(struct.pack("<4I", 0, len(b), 0x00020001, 0))
    hdr.extend(b[os2 + 32 : os2 + 42])
    hdr.extend(
        struct.pack("<BBIHH", 1, u16(os2 + 62) & 1, u16(os2 + 4), u16(os2 + 8), 0x504C)
    )
    hdr.extend(struct.pack("<4I", *[u32(os2 + 42 + i * 4) for i in range(4)]))
    hdr.extend(struct.pack("<2I", u32(os2 + 78), u32(os2 + 82)))
    hdr.extend(struct.pack("<5I", u32(head + 8), 0, 0, 0, 0))
    assert len(hdr) == 80
    for nid in [1, 2, 5, 4]:
        v = names[nid]
        hdr.extend(struct.pack("<HH", 0, len(v)))
        hdr.extend(v)
    hdr.extend(struct.pack("<HH", 0, 0))
    struct.pack_into("<I", hdr, 0, len(hdr) + len(b))
    return bytes(hdr) + b


def embed(pptx):
    P = "http://schemas.openxmlformats.org/presentationml/2006/main"
    R = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
    PKG = "http://schemas.openxmlformats.org/package/2006/relationships"
    CT = "http://schemas.openxmlformats.org/package/2006/content-types"
    with zipfile.ZipFile(pptx) as z:
        files = {n: z.read(n) for n in z.namelist()}
    # Office templates include localized fallback-font names. The artwork uses
    # embedded Liberation Sans; keep unused theme metadata in English as well.
    for name in list(files):
        if name.startswith("ppt/theme/") and name.endswith(".xml"):
            theme = ET.fromstring(files[name])
            for element in theme.iter():
                typeface = element.get("typeface", "")
                if typeface and not typeface.isascii():
                    element.set("typeface", "Arial")
            files[name] = ET.tostring(theme, encoding="UTF-8", xml_declaration=True)
    pres = ET.fromstring(files["ppt/presentation.xml"])
    old = pres.find("{" + P + "}embeddedFontLst")
    if old is not None:
        pres.remove(old)
    pres.set("embedTrueTypeFonts", "1")
    pres.set("saveSubsetFonts", "0")
    lst = ET.Element("{" + P + "}embeddedFontLst")
    ent = ET.SubElement(lst, "{" + P + "}embeddedFont")
    ET.SubElement(ent, "{" + P + "}font", typeface="Liberation Sans")
    rels = ET.fromstring(files["ppt/_rels/presentation.xml.rels"])
    for r in list(rels):
        if r.get("Id", "").startswith("rIdLiberation"):
            rels.remove(r)
    for idx, (tag, variant) in enumerate(
        [
            ("regular", "Regular"),
            ("bold", "Bold"),
            ("italic", "Italic"),
            ("boldItalic", "BoldItalic"),
        ],
        1,
    ):
        rid = "rIdLiberation" + str(idx)
        target = f"fonts/LiberationSans-{variant}.fntdata"
        ET.SubElement(ent, "{" + P + "}" + tag, {"{" + R + "}id": rid})
        ET.SubElement(
            rels, "{" + PKG + "}Relationship", Id=rid, Type=R + "/font", Target=target
        )
        files["ppt/" + target] = to_eot(FONTDIR / f"LiberationSans-{variant}.ttf")
    before = pres.find("{" + P + "}defaultTextStyle")
    if before is None:
        pres.append(lst)
    else:
        pres.insert(list(pres).index(before), lst)
    types = ET.fromstring(files["[Content_Types].xml"])
    if not any(e.get("Extension") == "fntdata" for e in types):
        ET.SubElement(
            types,
            "{" + CT + "}Default",
            Extension="fntdata",
            ContentType="application/x-fontdata",
        )
    for name, xml in [
        ("ppt/presentation.xml", pres),
        ("ppt/_rels/presentation.xml.rels", rels),
        ("[Content_Types].xml", types),
    ]:
        files[name] = ET.tostring(
            xml, encoding="UTF-8", xml_declaration=True, standalone=True
        )
    with zipfile.ZipFile(pptx, "w", zipfile.ZIP_DEFLATED) as z:
        for n, data in files.items():
            z.writestr(n, data)
    print("Embedded all four Liberation Sans faces:", pptx)


if __name__ == "__main__":
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument("pptx", nargs="?", type=Path, default=ROOT / "Fig4_final.pptx")
    embed(parser.parse_args().pptx)
