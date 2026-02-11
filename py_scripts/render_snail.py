#!/usr/bin/env python3

import argparse
from pathlib import Path
from playwright.sync_api import sync_playwright


def main():
    parser = argparse.ArgumentParser(
        description="Render an assembly-stats snail HTML plot to a static PNG using Playwright."
    )
    parser.add_argument(
        "--html",
        required=True,
        help="Input HTML file (e.g. snail.html)"
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output PNG filename (e.g. snail.png)"
    )
    parser.add_argument(
        "--width",
        type=int,
        default=1200,
        help="Viewport width in pixels (default: 1200)"
    )
    parser.add_argument(
        "--height",
        type=int,
        default=1200,
        help="Viewport height in pixels (default: 1200)"
    )
    parser.add_argument(
        "--scale",
        type=int,
        default=2,
        help="Device scale factor for higher resolution (default: 2)"
    )
    parser.add_argument(
        "--selector",
        default="#assembly_stats",
        help="CSS selector to screenshot (default: #assembly_stats)"
    )
    parser.add_argument(
        "--wait-ms",
        type=int,
        default=1000,
        help="Milliseconds to wait after load (default: 1000)"
    )

    args = parser.parse_args()

    html_path = Path(args.html).resolve()
    out_path = Path(args.out).resolve()

    if not html_path.exists():
        raise SystemExit(f"ERROR: input HTML file does not exist: {html_path}")

    with sync_playwright() as p:
        browser = p.chromium.launch(args=["--no-sandbox", "--disable-setuid-sandbox"])
        page = browser.new_page(
            viewport={"width": args.width, "height": args.height},
            device_scale_factor=args.scale
        )

        page.goto(html_path.as_uri(), wait_until="networkidle")
        page.wait_for_timeout(args.wait_ms)

        el = page.query_selector(args.selector)
        if el is None:
            raise SystemExit(f"ERROR: could not find selector '{args.selector}' in {html_path}")

        el.screenshot(path=str(out_path))
        browser.close()

    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
