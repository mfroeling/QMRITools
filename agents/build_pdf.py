"""Builds a document PDF: Markdown -> HTML (pandoc, styled by manual.css) -> PDF (Chrome, headless).

Run from the repository root:
    python agents\\build_pdf.py             the BIDS manual (QMRITools/BIDS Example/BIDS-config.md)
    python agents\\build_pdf.py Other.md    any other Markdown file in that folder
Needs pandoc on the PATH and Google Chrome.  Output: <name>.pdf in the same folder (the HTML and CSS it is made from are removed again), with
the document title (its first heading) in the page footer.  The Markdown file is the source, never edit the PDF.
"""

import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DOCS = os.path.join(os.path.dirname(HERE), "QMRITools", "BIDS Example")
SOURCE = os.path.join(DOCS, sys.argv[1] if len(sys.argv) > 1 else "BIDS-config.md")
NAME = os.path.splitext(os.path.basename(SOURCE))[0]
CSS = os.path.join(HERE, "manual.css")
OUT = DOCS
HTML = os.path.join(OUT, NAME + ".html")
PDF = os.path.join(OUT, NAME + ".pdf")
BROWSERS = (r"C:\Program Files\Google\Chrome\Application\chrome.exe",
            r"C:\Program Files (x86)\Google\Chrome\Application\chrome.exe",
            os.path.expandvars(r"%LOCALAPPDATA%\Google\Chrome\Application\chrome.exe"))


def slug(heading: str) -> str:
    """The anchor GitHub, VS Code and pandoc (gfm_auto_identifiers) all make from a heading."""
    text = heading.replace("`", "").strip().lower()
    return re.sub(r"[^\w\- ]", "", text).replace(" ", "-")


def update_toc(path: str) -> None:
    """Write the table of contents (headings of levels 2 and 3) between <!-- TOC --> and <!-- /TOC --> in the Markdown
    file itself, so that it is right in the editor preview and in the PDF.  Nothing happens if there are no markers."""
    with open(path, encoding="utf-8", newline="") as f:
        text = f.read().replace("\r\n", "\n")
    if "<!-- TOC -->" not in text:
        return
    entries, fenced = [], False
    for line in text.split("\n"):
        if line.startswith("```"):
            fenced = not fenced
        match = None if fenced else re.match(r"^(#{2,3}) (.+)$", line)
        if match and match.group(2).strip() != "Contents":
            title = match.group(2).strip()
            entries.append("%s- [%s](#%s)" % ("  " * (len(match.group(1)) - 2), title.replace("`", ""), slug(title)))
    head, rest = text.split("<!-- TOC -->\n", 1)
    tail = rest.split("<!-- /TOC -->", 1)[1]
    new = head + "<!-- TOC -->\n" + "\n".join(entries) + "\n<!-- /TOC -->" + tail
    if new != text:
        with open(path, "w", encoding="utf-8", newline="") as f:
            f.write(new)


def main() -> int:
    browser = next((p for p in BROWSERS if os.path.exists(p)), None)
    if browser is None:
        print("Google Chrome not found")
        return 1
    if not os.path.exists(SOURCE):
        print("not found:", SOURCE)
        return 1
    os.makedirs(OUT, exist_ok=True)
    update_toc(SOURCE)
    with open(SOURCE, encoding="utf-8") as f:                    # the title is the first heading: it goes in the footer
        title = next((line[2:].strip() for line in f if line.startswith("# ")), NAME)
    with open(CSS, encoding="utf-8") as f:
        css = f.read().replace("Muscle Stimulator Software: User Manual", title.replace('"', ""))
    css_path = os.path.join(OUT, NAME + ".css")
    with open(css_path, "w", encoding="utf-8") as f:
        f.write(css)
    subprocess.run(["pandoc", SOURCE, "-f", "markdown+gfm_auto_identifiers", "-s", "--embed-resources", "--syntax-highlighting=none", "--resource-path", DOCS, "--css", css_path,
                    "--metadata", "pagetitle=" + title, "-o", HTML], check=True)

    with open(HTML, encoding="utf-8") as f:
        html = f.read()
    # the table of contents is the list under the heading "Contents" in the Markdown: give it the TOC look (manual.css)
    html = re.sub(r'<h2 id="contents">Contents</h2>\s*<!-- TOC -->\s*(<ul>.*?</ul>)\s*(?=<h2)',
                  lambda m: '<nav id="TOC"><h2 class="toc-title">Contents</h2>' + m.group(1).replace(" <!-- /TOC -->", "") + "</nav>", html, count=1, flags=re.S)
    html = html.replace("\\|", "|")                               # a pipe escaped for GitHub tables, inside code
    html = re.sub(r"(<td[^>]*>)(\d{4})-(\d{2})-(\d{2})(</td>)",     # dates must not wrap at their hyphens
                  lambda m: m.group(1) + "‑".join(m.group(2, 3, 4)) + m.group(5), html)
    with open(HTML, "w", encoding="utf-8") as f:
        f.write(html)

    url = "file:///" + HTML.replace("\\", "/")
    subprocess.run([browser, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
                    "--generate-pdf-document-outline", "--print-to-pdf=" + PDF, url],
                   check=True, stderr=subprocess.DEVNULL)
    print("written", PDF)
    for temporary in (HTML, css_path):                           # only the PDF is kept in docs/build
        os.remove(temporary)
    return 0


if __name__ == "__main__":
    sys.exit(main())
