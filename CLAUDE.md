# Claude Instructions

Read `AGENTS.md` first.

`AGENTS.md` is the canonical repository instruction file for coding
agents in this repo. The package design standards, interface rules, and
validation contract are defined there and in the linked files under
`markdown_notes/design/`; the compact live roadmap is
`markdown_notes/rdistnn_roadmap.md`, and completed plans are archived
under `markdown_notes/archive/`.

Before changing code, especially in exported functions or shared
helpers, follow `AGENTS.md` and its required read order.

For mathematical derivations, symbolic checks, or finite combinatorial
verification, Mathematica / Wolfram Language is available through
`wolframscript` and should be used when helpful; record nontrivial
conclusions in code, tests, or design notes rather than relying on an
ephemeral session.

Local CLI tooling on the parent Windows workspace was verified on
2026-05-06: `rg`, `rga`, and `fd` for search/discovery; `jq`, `yq`, and
`duckdb` for structured data; `pdfinfo`, `pdftotext`, `pdftoppm`,
`pdffonts`, `qpdf`, `gswin64c`, `magick`, `tesseract`, `ocrmypdf`, and
`exiftool` for PDF/OCR/image work; and `Rscript`, `pandoc`, `latexmk`,
`lualatex`, `gh`, `gcc`, `g++`, `make`, `cmake`, and `wolframscript` for
package/document/build work. Use `magick`, not Windows `convert.exe`.
Tesseract has `eng` and `osd`; use
`ocrmypdf --skip-text --deskew --rotate-pages` for scanned PDFs and
`pdftotext -layout` for searchable PDFs. Code-quality tools now
available include `actionlint`, `shellcheck`, `shfmt`, `pre-commit`, and
`ruff`. `pip`, `pip3`, `pre-commit`, `ruff`, and `exiftool` are exposed
through shims in `C:\Users\jakob\.local\bin`.
