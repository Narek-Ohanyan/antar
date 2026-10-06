"""Copyright and citation: one statement, the same everywhere."""
import html
import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
NAME = "ANTAR — Assessment of Niche, Treeline & Analogue Refugia"
URL = "https://antar.narekohanyan.com"
REPO = "https://github.com/Narek-Ohanyan/antar"


def text_of(path):
    return html.unescape(re.sub(r"<[^>]+>", "", (ROOT / path).read_text()))


def test_licence_is_all_rights_reserved_not_mit():
    lic = (ROOT / "LICENSE").read_text()
    assert "Copyright (c) 2026 Narek Ohanyan. All rights reserved." in lic
    assert "Permission is hereby granted" not in lic
    assert "MIT" not in (ROOT / "pyproject.toml").read_text().split("[project.optional-dependencies]")[0]
    assert "MIT licensed" not in (ROOT / "README.md").read_text()


def test_citation_file_names_the_work_exactly():
    c = yaml.safe_load((ROOT / "CITATION.cff").read_text())
    assert c["title"] == NAME and c["preferred-citation"]["title"] == NAME
    assert len(c["authors"]) == 1
    assert (c["authors"][0]["family-names"], c["authors"][0]["given-names"]) == ("Ohanyan", "Narek")
    assert c["authors"][0]["website"] == "https://www.narekohanyan.com"
    assert c["preferred-citation"]["year"] == 2026 and c["repository-code"] == REPO and c["url"] == URL and c["preferred-citation"]["url"] == URL
    assert "2.0.0" in c["version"]


def test_footer_and_readme_carry_the_same_reference():
    footer = text_of("ui/index.html")
    readme = (ROOT / "README.md").read_text()
    ref = f"Ohanyan, N. (2026). {NAME} (Version 2.0.0-alpha) [Computer software and web interface]. {URL} (source code: {REPO})"
    assert ref in re.sub(r"\s+", " ", footer)
    assert ref.replace("*", "") in re.sub(r"\s+", " ", readme.replace("\n> ", " ").replace("*", ""))
    assert "© 2026 Narek Ohanyan. All rights reserved." in footer
    assert "@misc{ohanyan2026antar" in footer and "@misc{ohanyan2026antar" in readme
    assert "Treeline \\& Analogue Refugia" in footer and "Treeline \\& Analogue Refugia" in readme      # the ampersand is escaped for LaTeX


def test_no_doi_or_other_identifier_is_invented():
    for p in ("CITATION.cff", "README.md", "ui/index.html"):
        assert "doi.org" not in (ROOT / p).read_text().lower(), p
