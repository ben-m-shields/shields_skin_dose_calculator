# Builds windows\bundle\ (a private copy of R with the app's packages, Python with the RDSR packages, and the shiny app)
# ready for electron-builder to package, then builds the zip into windows\dist\. Run it in PowerShell on Windows.
$ErrorActionPreference = 'Stop'
$ProgressPreference = 'SilentlyContinue'   # Invoke-WebRequest is very slow with the progress bar on

$RVersion = '4.5.2'
$PythonVersion = '3.11.9'
$RPackages = "'ggplot2','dplyr','lubridate','DT','tidyr','plotly','shiny','reticulate','gdata','bslib','shinyMatrix','shinybusy','officer','flextable','ragg'"
# pandas's dependencies are pinned too, to the versions the app was tested with
$PythonPackages = @('pydicom==2.4.0', 'rdsr_navigator==0.5.0', 'pandas==2.2.3', 'numpy==2.2.6', 'python-dateutil==2.9.0.post0', 'six==1.17.0', 'pytz==2026.3.post1', 'tzdata==2026.3')
$AppFiles = @('app.R', 'skin_dose_map.R', 'report_docx_skin_dose.R', 'read_RDSR_irradiation_events.py',
              'read_RDSR_acc_dose_data.py', 'environment.RData', 'www')

$appsrc = Join-Path (Split-Path -Parent $PSScriptRoot) 'app'
$bundle = Join-Path $PSScriptRoot 'bundle'
$tmp = Join-Path ([IO.Path]::GetTempPath()) 'skin-dose-build'
Remove-Item -Recurse -Force $bundle, $tmp -ErrorAction SilentlyContinue
New-Item -ItemType Directory $bundle, $tmp | Out-Null

# --- R -------------------------------------------------------------------------------------------
Write-Host "Installing R $RVersion"
Invoke-WebRequest "https://cloud.r-project.org/bin/windows/base/old/$RVersion/R-$RVersion-win.exe" -OutFile "$tmp\R-setup.exe"
$p = Start-Process "$tmp\R-setup.exe" -ArgumentList '/VERYSILENT', '/SUPPRESSMSGBOXES', '/NORESTART', '/CURRENTUSER', '/NOICONS', "/DIR=$bundle\R" -Wait -PassThru
if ($p.ExitCode -ne 0) { throw "R installer failed with exit code $($p.ExitCode)" }
Remove-Item -Recurse -Force "$bundle\R\doc", "$bundle\R\tests" -ErrorAction SilentlyContinue

Write-Host 'Installing R packages'
$rscript = "$bundle\R\bin\Rscript.exe"
# install into R's own library so the bundle doesn't depend on anything outside its folder
& $rscript --vanilla -e "install.packages(c($RPackages), lib = .Library, repos = 'https://cloud.r-project.org', type = 'binary')"
if ($LASTEXITCODE -ne 0) { throw 'Installing R packages failed' }
# install.packages only warns when a package fails, so check they all load
& $rscript --vanilla -e "for (p in c($RPackages)) library(p, character.only = TRUE); cat('All R packages load\n')"
if ($LASTEXITCODE -ne 0) { throw 'An R package failed to load' }
Get-ChildItem "$bundle\R\library" -Directory | ForEach-Object {
  Remove-Item -Recurse -Force "$($_.FullName)\doc", "$($_.FullName)\tests" -ErrorAction SilentlyContinue
}

# --- Python --------------------------------------------------------------------------------------
# The python.org NuGet package is a complete, relocatable Python (unlike a venv, which hard-codes its location)
Write-Host "Installing Python $PythonVersion"
Invoke-WebRequest "https://www.nuget.org/api/v2/package/python/$PythonVersion" -OutFile "$tmp\python.zip"
Expand-Archive "$tmp\python.zip" "$tmp\python-nuget"
Move-Item "$tmp\python-nuget\tools" "$bundle\python"
$py = "$bundle\python\python.exe"
& $py -m ensurepip --default-pip
& $py -m pip install --no-warn-script-location --disable-pip-version-check @PythonPackages
if ($LASTEXITCODE -ne 0) { throw 'Installing Python packages failed' }

# Python runs inside R's process, and Windows looks for DLLs next to the program that started the process (R), not
# next to Python. PCs without the Visual C++ runtime installed then can't load numpy, so put Python's copy beside R.
Copy-Item "$bundle\python\vcruntime140*.dll" "$bundle\R\bin\x64\"

# --- The shiny app -------------------------------------------------------------------------------
Write-Host 'Copying the app'
New-Item -ItemType Directory "$bundle\shiny" | Out-Null
foreach ($f in $AppFiles) { Copy-Item -Recurse (Join-Path $appsrc $f) "$bundle\shiny\" }

# --- Check the bundled R can load the Python packages ----------------------------------------------
Write-Host 'Checking R can load the Python packages'
$check = "library(reticulate); use_python('$($py -replace '\\', '/')', required = TRUE); " +
         "for (m in c('pydicom', 'rdsr_navigator', 'pandas')) import(m); cat('Python packages load\n')"
$env:PYTHONNOUSERSITE = '1'
& $rscript --vanilla -e $check
if ($LASTEXITCODE -ne 0) { throw 'The bundled R could not load the bundled Python packages' }
Remove-Item -Recurse -Force $tmp

# --- Zip -----------------------------------------------------------------------------------------
Write-Host 'Building the zip'
Push-Location $PSScriptRoot
if (-not (Test-Path node_modules)) { npm ci; if ($LASTEXITCODE -ne 0) { throw 'npm ci failed' } }
npm run dist
$ok = $LASTEXITCODE
Pop-Location
if ($ok -ne 0) { throw 'Building the zip failed' }
Get-ChildItem "$PSScriptRoot\dist\*.zip"
