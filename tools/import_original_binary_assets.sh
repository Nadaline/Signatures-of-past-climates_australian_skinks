#!/usr/bin/env bash
# Add original binary assets from the private research repository without
# overwriting the chapter figures already published from the R7 thesis.
set -euo pipefail

SOURCE_URL="https://github.com/Nadaline/Lampropholis-delicata-genomic-history-staging.git"
PUBLIC_URL="https://github.com/Nadaline/Signatures-of-past-climates_australian_skinks.git"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

git clone --depth 1 --branch docs/reader-guided-reproducibility "$SOURCE_URL" "$WORK/source"
git clone --depth 1 "$PUBLIC_URL" "$WORK/public"

rsync -a --ignore-existing \
  --include='*/' --include='*.png' --include='*.pdf' \
  --exclude='*' "$WORK/source/" "$WORK/public/"

cd "$WORK/public"
git add assets/ design/ results/ || true
if git diff --cached --quiet; then
  echo "All original binary assets are already present."
  exit 0
fi

git commit -m "Add original figure and guide exports"
git push origin main
echo "Original missing PNG and PDF assets have been added."
