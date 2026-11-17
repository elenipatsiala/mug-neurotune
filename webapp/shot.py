import sys
from playwright.sync_api import sync_playwright
with sync_playwright() as p:
    b = p.chromium.launch(executable_path="/usr/bin/google-chrome", args=["--no-sandbox"])
    pg = b.new_page(viewport={"width": 1280, "height": 900})
    errs = []
    pg.on("pageerror", lambda e: errs.append(str(e)))
    pg.on("console", lambda m: m.type == "error" and errs.append(m.text))
    pg.goto("file://" + __import__("os").path.abspath("index.html"))
    pg.click("#demoBtn")
    pg.wait_for_timeout(800)
    pg.screenshot(path=sys.argv[1], full_page=True)
    print("errors:", [e for e in errs if "fonts" not in e])
    print(pg.inner_text("#verdict"))
    b.close()
