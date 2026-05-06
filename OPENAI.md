# OpenAI / Codex Instructions

Read `AGENTS.md` first.

`AGENTS.md` is the canonical repository instruction file for coding agents in
this repo. The design-note system for `RdistNN` flows from that file and from
its linked markdown notes under `markdown_notes/design/`.

When working on user-facing APIs, validation, error messages, interfaces,
documentation, or tests, follow the read order and rules in `AGENTS.md` before
making changes.

Local CLI tooling on the parent Windows workspace was verified on 2026-05-06:
`rg`, `rga`, and `fd` for search/discovery; `jq`, `yq`, and `duckdb` for
structured data; `pdfinfo`, `pdftotext`, `pdftoppm`, `pdffonts`, `qpdf`,
`gswin64c`, `magick`, `tesseract`, and `ocrmypdf` for PDF/OCR/image work; and
`Rscript`, `pandoc`, `latexmk`, `lualatex`, `gh`, `gcc`, `g++`, `make`,
`cmake`, and `wolframscript` for package/document/build work. Use `magick`,
not Windows `convert.exe`. Tesseract has `eng` and `osd`; use
`ocrmypdf --skip-text --deskew --rotate-pages` for scanned PDFs and
`pdftotext -layout` for searchable PDFs.
