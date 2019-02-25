#!/bin/bash
for file in $(find . -type f -name '*.md') ; do
  pandoc ${file} -f markdown+lists_without_preceding_blankline -t html5 \
  -s -o ${file%.md}.html --metadata pagetitle=${file%.md} --mathjax
done

# pandoc options
# -lists_without_preceding_blankline: allow a list to occur right after a paragraph, with no intervening blank space
# -s, --standalone
# --mathjax (or --mathml)
#   - Only works with pandoc markdown format (-f markdown), not GitHub-flavored markdown (-f gfm)
#   - mathml only renders properly in Firefox and Safari
