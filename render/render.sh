#!/usr/bin/env bash
# Renders a handout from Markdown to one self-contained HTML file in GitHub style.
# Images and mermaid.js are fetched once, at render time, and embedded. Mermaid diagrams are
# drawn in the browser by the embedded mermaid.js, so the file works offline.
# Usage: render.sh IN.md "Browser tab title" OUT.html
here=$(dirname "$0")
pandoc "$1" -f gfm -t html5 --standalone \
  --lua-filter="$here/mermaid.lua" --include-in-header="$here/header.html" \
  --highlight-style=pygments --embed-resources --resource-path="$(dirname "$1")" --metadata pagetitle="$2" -o "$3"
