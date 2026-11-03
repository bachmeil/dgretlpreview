#!/bin/bash
# Convert the guide's markdown files to html with pandoc.
# Usage: guide/make-html.sh [output directory]      (default: guide/html)
# Every .md file under guide/ becomes a .html file in the same place under the output directory,
# with pandoc/style.css next to index.html and the images copied across.
# Needs: pandoc (3.x). Syntax highlighting uses pandoc's own, with the D definition.
set -e
cd "$(dirname "$0")"
OUT="${1:-html}"
mkdir -p "$OUT"
cp pandoc/style.css "$OUT/style.css"

find . -name '*.md' -not -path './build/*' -not -path "./$OUT/*" | sort | while read -r f; do
    f="${f#./}"
    dest="$OUT/${f%.md}.html"
    mkdir -p "$(dirname "$dest")"
    # the path from this page back to the top of the guide, for the style sheet and the links in the template
    depth=$(grep -o '/' <<< "$f" | wc -l)
    root=""
    for ((i = 0; i < depth; i++)); do root="../$root"; done
    pandoc "$f" --from gfm --to html5 --standalone --template pandoc/template.html \
        --lua-filter pandoc/links.lua --highlight-style tango -V root="$root" -o "$dest"
done

# images the pages refer to
find . \( -name '*.png' -o -name '*.svg' -o -name '*.jpg' \) -not -path './build/*' -not -path "./$OUT/*" | while read -r f; do
    f="${f#./}"
    mkdir -p "$OUT/$(dirname "$f")"
    cp "$f" "$OUT/$f"
done
echo "wrote $OUT/index.html and $(find "$OUT" -name '*.html' | wc -l) pages"
