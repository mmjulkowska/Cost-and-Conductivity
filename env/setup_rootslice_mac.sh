#!/usr/bin/env bash
# Install RootSlice (Jagdeep Sidhu's Apple Silicon build) and the GRANAR R package on macOS.
# RootSlice is third-party software and is NOT in this repository. Get RootSlice_for_Ruth.zip from Box (see data/README.md),
# then run:   bash env/setup_rootslice_mac.sh /path/to/RootSlice_for_Ruth.zip
# Installs to $ROOTSLICE_HOME (default ~/root_models/RootSlice). Needs R (https://cran.r-project.org).
set -e
ZIP="${1:?give the path to RootSlice_for_Ruth.zip}"; RS="${ROOTSLICE_HOME:-$HOME/root_models/RootSlice}"
[ "$(uname -m)" = "arm64" ] || echo "WARNING: this RootSlice build runs on Apple Silicon only"
mkdir -p "$(dirname "$RS")"; tmp="$(mktemp -d)"; unzip -q "$ZIP" -d "$tmp"
src="$(dirname "$(find "$tmp" -path '*/app/RootSlice' | head -1)")/.."; rm -rf "$RS"; cp -R "$src" "$RS"; rm -rf "$tmp"
(cd "$RS" && xattr -dr com.apple.quarantine . 2>/dev/null || true; chmod +x run.sh app/RootSlice 2>/dev/null || true)
echo "RootSlice installed at $RS"; ls "$RS"
Rscript -e 'if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", repos="https://cloud.r-project.org");
  if (!"granar" %in% rownames(installed.packages())) remotes::install_github("granar/granar", upgrade="never");
  cat("granar", as.character(packageVersion("granar")), "(analysis used 0.0.1.2)\n")'
