"""The deployable bundle (scripts/build_site.py) is complete and consistent, and deploy/.htaccess does what it says on a real Apache.

The Apache tests start the system httpd on loopback ports with a throwaway certificate (mod_ssl, mod_rewrite, mod_headers, mod_deflate); they are skipped, not
passed, when httpd, openssl or its modules are missing."""
import http.client
import json
import shutil
import socket
import ssl
import subprocess
import sys
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import build_site  # noqa: E402

SUB = "antar.narekohanyan.com"


@pytest.fixture(scope="module")
def bundle(tmp_path_factory):
    out = build_site.assemble(tmp_path_factory.mktemp("bundle") / "site")
    return out


def test_the_bundle_is_complete_and_consistent(bundle):
    assert build_site.check_bundle(bundle) == []


def test_python_and_caches_stay_behind_and_the_server_files_are_in(bundle):
    names = {p.name for p in bundle.rglob("*")}
    assert not any(n.endswith((".py", ".pyc")) or n == "__pycache__" for n in names)
    assert ".htaccess" in names and "robots.txt" in names and "hero.mp4" in names and "credits.json" in names
    assert (bundle / "robots.txt").read_text().startswith("User-agent: *")


def test_a_broken_bundle_is_reported_not_passed(bundle, tmp_path):
    broken = tmp_path / "broken"
    shutil.copytree(bundle, broken)
    (broken / "app.js").unlink()
    (broken / "assets" / "photos" / "aragats-1600.jpg").unlink()
    (broken / "stray.py").write_text("x")
    text = " | ".join(build_site.check_bundle(broken))
    assert "app.js" in text and "aragats-1600.jpg" in text and "stray.py" in text
    (bundle / "index.html").read_text()
    stale = tmp_path / "stale"
    shutil.copytree(bundle, stale)
    idx = (stale / "index.html").read_text().replace("styles.css?v=", "styles.css?v=zzzz", 1)
    (stale / "index.html").write_text(idx)
    assert any("stamp" in p for p in build_site.check_bundle(stale))


def test_the_manifest_in_the_bundle_lists_the_photographs_and_grids(bundle):
    m = json.loads((bundle / "data" / "manifest.json").read_text())
    assert len(m["photos"]) == 4 and set(m["grids"]) >= {"validation"}


# ------------------------------------------------------------------ real Apache
MODS = Path("/usr/libexec/apache2")
NEEDED = ["mpm_prefork", "authz_core", "unixd", "mime", "dir", "alias", "rewrite", "headers", "deflate", "ssl", "socache_shmcb", "log_config", "filter", "env"]
HTTPD = shutil.which("httpd") or ("/usr/sbin/httpd" if Path("/usr/sbin/httpd").exists() else None)
can_apache = bool(HTTPD and shutil.which("openssl") and all((MODS / f"mod_{m}.so").exists() for m in NEEDED))


def free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


@pytest.fixture(scope="module")
def apache(bundle, tmp_path_factory):
    if not can_apache:
        pytest.skip("no system Apache with mod_ssl / openssl")
    d = tmp_path_factory.mktemp("httpd")
    key, crt = d / "k.pem", d / "c.pem"
    subprocess.run(["openssl", "req", "-x509", "-newkey", "rsa:2048", "-nodes", "-keyout", str(key), "-out", str(crt), "-days", "1", "-subj", "/CN=localhost"],
                   check=True, capture_output=True)
    hp, sp = free_port(), free_port()
    loads = "\n".join(f"LoadModule {m}_module {MODS / ('mod_' + m + '.so')}" for m in NEEDED)
    conf = d / "httpd.conf"
    conf.write_text(f"""ServerRoot "{d}"
DefaultRuntimeDir "{d}"
Mutex file:{d} default
PidFile "{d}/httpd.pid"
ServerName localhost
Listen 127.0.0.1:{hp}
Listen 127.0.0.1:{sp}
{loads}
ErrorLog "{d}/error.log"
LogLevel warn
DocumentRoot "{bundle}"
TypesConfig /private/etc/apache2/mime.types
<Directory "{bundle}">
  AllowOverride All
  Require all granted
</Directory>
<VirtualHost 127.0.0.1:{sp}>
  SSLEngine on
  SSLCertificateFile "{crt}"
  SSLCertificateKeyFile "{key}"
</VirtualHost>
""")
    check = subprocess.run([HTTPD, "-f", str(conf), "-t"], capture_output=True, text=True)
    if check.returncode:
        pytest.skip("this Apache cannot load the test configuration: " + check.stderr.strip()[:200])
    proc = subprocess.Popen([HTTPD, "-f", str(conf), "-D", "FOREGROUND"], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, start_new_session=True)   # own session: Apache signals its whole process group when it stops
    for _ in range(100):
        try:
            socket.create_connection(("127.0.0.1", sp), timeout=0.2).close()
            break
        except OSError:
            time.sleep(0.1)
    else:
        proc.terminate()
        pytest.skip("Apache did not start: " + (proc.stderr.read().decode()[:300] if proc.stderr else ""))
    yield {"http": hp, "https": sp}
    proc.terminate()
    proc.wait(timeout=10)


def get(apache, path, host=SUB, tls=True, headers=None):
    ctx = ssl._create_unverified_context()
    conn = (http.client.HTTPSConnection("127.0.0.1", apache["https"], context=ctx, timeout=10) if tls else http.client.HTTPConnection("127.0.0.1", apache["http"], timeout=10))
    conn.request("GET", path, headers={"Host": host, "Accept-Encoding": "identity", **(headers or {})})
    r = conn.getresponse()
    body = r.read()
    out = (r.status, {k.lower(): v for k, v in r.getheaders()}, body)
    conn.close()
    return out


def test_plain_http_goes_to_https_on_the_subdomain(apache):
    s, h, _ = get(apache, "/data/manifest.json?v=1", tls=False)
    assert s == 301 and h["location"] == f"https://{SUB}/data/manifest.json?v=1"          # the query string is kept


def test_the_same_folder_under_the_main_domain_redirects_to_the_subdomain(apache):
    for tls in (True, False):
        s, h, _ = get(apache, "/app.js?v=1", host="narekohanyan.com", tls=tls)
        assert s == 301 and h["location"] == f"https://{SUB}/app.js?v=1"
    s, h, _ = get(apache, "/", host="www.narekohanyan.com")
    assert s == 301 and h["location"] == f"https://{SUB}/"


def test_the_page_is_served_and_always_revalidated(apache):
    s, h, body = get(apache, "/")
    assert s == 200 and b"ANTAR" in body and h["cache-control"] == "no-cache"
    s, h, _ = get(apache, "/index.html")
    assert s == 301 and h["location"].startswith(f"https://{SUB}") and h["location"].endswith("/") and "index.html" not in h["location"]    # Apache makes it absolute


def test_scripts_data_and_images_get_their_cache_rules(apache):
    assert get(apache, "/app.js?v=x")[1]["cache-control"] == "public, max-age=3600, must-revalidate"
    assert get(apache, "/data/manifest.json?v=x")[1]["cache-control"] == "public, max-age=3600, must-revalidate"
    assert get(apache, "/assets/photos/aragats-800.jpg")[1]["cache-control"] == "public, max-age=604800"
    assert get(apache, "/assets/fonts/crimson-pro-var.woff2")[1]["content-type"] == "font/woff2"
    assert get(apache, "/assets/map/borders.geojson")[1]["content-type"].startswith("application/geo+json")


def test_text_is_compressed_and_video_supports_ranges(apache):
    s, h, _ = get(apache, "/app.js", headers={"Accept-Encoding": "gzip"})
    assert s == 200 and h.get("content-encoding") == "gzip"
    s, h, body = get(apache, "/assets/hero.mp4", headers={"Range": "bytes=0-99"})
    assert s == 206 and len(body) == 100


def test_folders_are_not_listed_and_missing_files_are_404(apache):
    s, _, body = get(apache, "/assets/")
    assert s in (403, 404) and b"Index of" not in body and b"hero.mp4" not in body
    assert get(apache, "/nothing-here.html")[0] == 404


def test_the_map_rasters_reach_the_browser_byte_for_byte_and_untransformed(apache, bundle):
    import hashlib
    for name in ("elevation", "region", "forest"):
        # even a browser that announces WebP support must get the original bytes: they are data, and a converted image would change the values
        s, h, body = get(apache, f"/assets/map/{name}.bin?v=x", headers={"Accept": "image/avif,image/webp,image/apng,image/svg+xml,image/*,*/*;q=0.8"})
        assert s == 200 and h["content-type"] == "application/octet-stream" and "no-transform" in h["cache-control"]
        assert hashlib.sha256(body).hexdigest() == hashlib.sha256((bundle / "assets" / "map" / f"{name}.bin").read_bytes()).hexdigest()
