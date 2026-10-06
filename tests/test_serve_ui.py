"""scripts/serve_ui.py: Range requests (Safari needs them for video) and revalidation instead of stale caching."""
import http.client
import sys
import threading
from http.server import ThreadingHTTPServer
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import serve_ui  # noqa: E402


@pytest.fixture(scope="module")
def server():
    srv = ThreadingHTTPServer(("127.0.0.1", 0), serve_ui.Handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield srv.server_address[1]
    srv.shutdown()


def get(port, path, headers=None):
    c = http.client.HTTPConnection("127.0.0.1", port)
    c.request("GET", path, headers=headers or {})
    r = c.getresponse()
    return r.status, dict(r.getheaders()), r.read()


def test_scripts_are_revalidated_not_reused(server):
    status, h, _ = get(server, "/app.js")
    assert status == 200 and h["Cache-Control"] == "no-cache"


def test_range_request_returns_exactly_the_requested_bytes(server):
    full_status, full_h, full = get(server, "/assets/hero.mp4")
    status, h, body = get(server, "/assets/hero.mp4", {"Range": "bytes=100-199"})
    assert status == 206 and body == full[100:200]
    assert h["Content-Range"] == f"bytes 100-199/{len(full)}"
    status, _, tail = get(server, "/assets/hero.mp4", {"Range": "bytes=-50"})
    assert status == 206 and tail == full[-50:]
    status, _, open_end = get(server, "/assets/hero.mp4", {"Range": f"bytes={len(full) - 10}-"})
    assert status == 206 and open_end == full[-10:]
    assert get(server, "/assets/hero.mp4", {"Range": f"bytes={len(full) + 5}-"})[0] == 416
