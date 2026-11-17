import json, pathlib
h = pathlib.Path("template.html").read_text()
h = h.replace("/*CORE*/", pathlib.Path("core.js").read_text())
h = h.replace("/*GENESETS*/", pathlib.Path("genesets.json").read_text())
h = h.replace("/*DEMO_COUNTS*/", json.dumps(pathlib.Path("demo_counts.csv").read_text()))
h = h.replace("/*DEMO_META*/", json.dumps(pathlib.Path("demo_meta.csv").read_text()))
pathlib.Path("index.html").write_text(h)
print("built", len(h)//1024, "KB")
