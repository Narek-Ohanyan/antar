"""The acknowledgments page: the supplied names, and the FORACCA facts as published by the sources it cites."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
APP = (ROOT / "ui" / "app.js").read_text()
ACK = APP[APP.index("function renderAck"):APP.index("/* ---------- ", APP.index("function renderAck"))]
TEXT = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", ACK.replace("&amp;", "&")))


def test_every_person_and_institution_is_named_as_supplied():
    for s in ("Franziska Zilker", "Tobias Kühnhanss", "Dr. Michael James McCarthy", "Dynamic Macroecology group", "PD Dr. Marco Pütz",
              "Dr. Dominik Braunschweiger", "Alen Amirkhanian, Director of the AUA Acopian Center for the Environment",
              "Swiss Federal Institute for Forest, Snow and Landscape Research WSL"):
        assert s in TEXT, s
    assert "his supervisor" in TEXT and TEXT.count("The author") >= 3            # third person throughout


def test_foracca_section_has_the_published_facts_and_links_to_its_sources():
    for s in ("Forest Restoration and Climate Change in Armenia (FORACCA)", "Swiss Agency for Development and Cooperation (SDC)",
              "10 years, 2023–2033, CHF 10 million", "Main phase 2025–2028", "Forest Alliance", "Shen NGO",
              "Food and Agriculture Organization of the United Nations (FAO)",
              "Advance scientific understanding of Armenia’s capacity to address climate change and sustainably manage its forests.",
              "Promote climate-smart practices in rural areas.",
              "Ensure evidence-based policymaking for climate adaptation and efficient forest management."):
        assert s in TEXT, s
    assert "https://www.wsl.ch/en/projects/foracca/" in ACK and "https://armenpress.am/en/article/1126549" in ACK
    assert "Recreation" not in ACK
    assert 'id="foracca"' in ACK and "data-scroll" in ACK


def test_external_links_open_safely_and_announce_it():
    for m in re.finditer(r"<a [^>]*target=\"_blank\"[^>]*>", ACK):
        assert 'rel="noopener"' in m.group(0)
    assert ACK.count("(opens in a new tab)") >= 1
