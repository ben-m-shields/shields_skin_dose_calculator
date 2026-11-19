#!/usr/bin/env bash
# Builds the Windows portable zip on Linux. Does the same as build-windows.ps1 but without running anything Windows:
# R's installer is unpacked with innoextract, and the R and Python packages are the Windows builds from CRAN and PyPI.
# The result can't be tested here, so try it on a Windows PC before giving it to anyone.
#
# Needs R 4.5 (to work out the package list), pip, npm and innoextract. Set INNOEXTRACT if innoextract isn't on PATH.
set -euo pipefail

R_VERSION=4.5.2
PYTHON_VERSION=3.11.9
R_PACKAGES="'ggplot2','dplyr','lubridate','DT','tidyr','plotly','shiny','reticulate','gdata','bslib','shinyMatrix','shinybusy','officer','flextable','ragg'"
# pandas's dependencies are pinned too, to the versions the app was tested with
PYTHON_PACKAGES="pydicom==2.4.0 rdsr_navigator==0.5.0 pandas==2.2.3 numpy==2.2.6 python-dateutil==2.9.0.post0 six==1.17.0 pytz==2026.3.post1 tzdata==2026.3"
APP_FILES="app.R skin_dose_map.R report_docx_skin_dose.R read_RDSR_irradiation_events.py read_RDSR_acc_dose_data.py environment.RData www"
INNOEXTRACT=${INNOEXTRACT:-innoextract}

desktop=$(cd "$(dirname "$0")" && pwd)
appsrc="$(dirname "$desktop")/app"
bundle="$desktop/bundle"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
rm -rf "$bundle"
mkdir -p "$bundle"

echo "== R $R_VERSION"
curl -fsSL -o "$tmp/R-setup.exe" "https://cloud.r-project.org/bin/windows/base/old/$R_VERSION/R-$R_VERSION-win.exe"
"$INNOEXTRACT" --silent --output-dir "$tmp/R" "$tmp/R-setup.exe"
mv "$tmp/R/app" "$bundle/R"
rm -rf "$bundle/R/doc" "$bundle/R/tests"

echo "== R packages (Windows binaries)"
mkdir -p "$tmp/rpkgs"
# work out every package needed, skipping ones R already ships with, and download their Windows builds
Rscript --vanilla -e "
  contrib <- 'https://cloud.r-project.org/bin/windows/contrib/4.5'
  ap <- available.packages(contriburl = contrib, type = 'win.binary')
  pkgs <- c($R_PACKAGES)
  need <- unique(c(pkgs, unlist(tools::package_dependencies(pkgs, db = ap, which = c('Depends', 'Imports'), recursive = TRUE))))
  need <- setdiff(need, list.files('$bundle/R/library'))
  got <- download.packages(need, '$tmp/rpkgs', contriburl = contrib, type = 'win.binary', quiet = TRUE)
  missing <- setdiff(need, got[, 1])
  if (length(missing)) stop('Could not download: ', paste(missing, collapse = ', '))
  cat('Downloaded', nrow(got), 'packages\n')"
for z in "$tmp"/rpkgs/*.zip; do unzip -q -o "$z" -d "$bundle/R/library"; done
for p in "$bundle"/R/library/*/; do rm -rf "$p/doc" "$p/tests"; done

echo "== Python $PYTHON_VERSION"
# the python.org NuGet package is a complete, relocatable Python (unlike a venv, which hard-codes its location)
curl -fsSL -o "$tmp/python.zip" "https://www.nuget.org/api/v2/package/python/$PYTHON_VERSION"
unzip -q "$tmp/python.zip" "tools/*" -d "$tmp/python"
mv "$tmp/python/tools" "$bundle/python"
python3 -m pip install --quiet --disable-pip-version-check --no-compile --target "$bundle/python/Lib/site-packages" \
  --platform win_amd64 --python-version "${PYTHON_VERSION%.*}" --implementation cp --only-binary=:all: $PYTHON_PACKAGES

# Python runs inside R's process, and Windows looks for DLLs next to the program that started the process (R), not
# next to Python. PCs without the Visual C++ runtime installed then can't load numpy, so put Python's copy beside R.
cp "$bundle"/python/vcruntime140*.dll "$bundle/R/bin/x64/"

echo "== App files"
mkdir -p "$bundle/shiny"
for f in $APP_FILES; do cp -r "$appsrc/$f" "$bundle/shiny/"; done

echo "== Zip"
cd "$desktop"
[ -d node_modules ] || npm ci
npm run dist
ls -lh "$desktop"/dist/*.zip
