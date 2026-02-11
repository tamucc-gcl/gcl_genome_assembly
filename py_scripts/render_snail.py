from pathlib import Path
from playwright.sync_api import sync_playwright

HTML = Path("snail.html").resolve()
OUT  = Path("snail.png").resolve()

with sync_playwright() as p:
    browser = p.chromium.launch(args=["--no-sandbox", "--disable-setuid-sandbox"])
    page = browser.new_page(viewport={"width": 1200, "height": 1200}, device_scale_factor=2)

    page.goto(HTML.as_uri(), wait_until="networkidle")
    page.wait_for_timeout(1000)  # let D3 finish

    el = page.query_selector("#assembly_stats")
    if el is None:
        raise SystemExit("ERROR: could not find #assembly_stats in snail.html")

    el.screenshot(path=str(OUT))
    browser.close()

print(f"Wrote {OUT}")
