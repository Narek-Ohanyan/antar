"""The author page: the bio as supplied, the portrait, the link, and the route."""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
APP = (ROOT / "ui" / "app.js").read_text()
SITE = "https://www.narekohanyan.com"


def test_route_footer_link_and_assets_exist():
    assert "author: renderAuthor" in APP and "function renderAuthor" in APP
    assert 'href="#/author"' in (ROOT / "ui" / "index.html").read_text()
    img = ROOT / "ui" / "assets" / "author.jpg"
    assert img.exists() and img.stat().st_size < 400_000                      # a web-sized portrait, not the 2000 px original


def test_bio_is_the_supplied_text_and_links_to_the_personal_site():
    body = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", APP.replace("&amp;", "&")))
    for sentence in (
        "Narek Ohanyan is a young climate leader from Armenia.",
        "He is a climate & environmental researcher at the AUA Acopian Center for the Environment.",
        "He has also served as a Guest Scientist at the Swiss Federal Institute for Forest, Snow and Landscape Research WSL.",
        "His current research focuses on modeling forest climate resilience.",
        "he authored the book The Overshoot: Life After the 1.5°C Limit.",
        "he has managed 50+ projects and secured tens of thousands of dollars in grants, empowering 10k+ Youth for Climate Action.",
        "Narek serves as the Lead Organizer for LCOY Armenia (2025-2026).",
        "He advises the Ministry of Environment and UNICEF through the Youth Climate Council,",
        "previously served as a UNFCCC COP27 Party Delegate and UN Youth and Children High-Level Climate Champion.",
        "Narek holds a BS in Computer Science (’26) from the American University of Armenia.",
    ):
        assert sentence in body, sentence
    assert f'href="{SITE}"' in APP and 'rel="noopener"' in APP


def test_portrait_has_a_text_alternative_and_declared_size():
    imgs = re.findall(r"<img[^>]*assets/author\.jpg[^>]*>", APP)
    assert len(imgs) == 2
    assert all('alt="Portrait of Narek Ohanyan"' in i and 'width="720"' in i and 'height="720"' in i for i in imgs)
