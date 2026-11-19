# Shields skin dose calculator

Estimates peak skin dose for fluoroscopically guided procedures from the RDSR (radiation dose structured report)
that the X-ray system produces.

- `app/` - the source code for the Shiny app
- `windows/` - the Windows version: download it from Releases, and the scripts that build it

## Running it on Windows (no install needed)

1. Download `Shields-Skin-Dose-Calculator-<version>-win-x64.zip` from the
   [latest release](https://github.com/ben-m-shields/shields_skin_dose_calculator/releases/latest).
2. Right-click the zip, choose **Extract All**, and extract it somewhere like your Documents folder.
3. Open the extracted folder and double-click **Shields Skin Dose Calculator.exe**. Right-click it and choose
   **Send to > Desktop (create shortcut)** if you want a desktop shortcut.

The first start can take up to a minute. Closing the window closes the calculator.

If Windows shows "Windows protected your PC", click **More info** then **Run anyway**. This appears because the
app isn't code signed.

Everything runs on your own computer - RDSR files and patient details are never uploaded anywhere. To remove it,
delete the folder.

If the calculator won't start, the error message shows where the log file is (in your temp folder,
`shields-skin-dose-calculator.log`). Please include it when reporting a problem.

## Running from the source code

Install R 4.5 and RStudio. Open `app/app.R` in RStudio and click 'Run App'. The first run installs the R packages
and creates a Python environment in `app/python_venv/`, which needs Python 3 installed.

Or with Docker, from the top of the repo:

```
docker build -t shields-skin-dose app
docker run -p 3838:3838 shields-skin-dose
```

then open http://localhost:3838.
