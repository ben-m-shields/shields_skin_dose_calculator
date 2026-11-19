# Windows version

**Download:** `Shields-Skin-Dose-Calculator-<version>-win-x64.zip` from the
[latest release](https://github.com/ben-m-shields/shields_skin_dose_calculator/releases/latest). The zip is too
big to keep in the repo itself (GitHub's limit is 100 MB per file), so each version is attached to a release.
See the [main README](../README.md) for how to run it.

The zip contains a small Electron app, a private copy of R 4.5.2 with the app's packages, Python 3.11.9 with the
RDSR packages, and the app from `../app`. `main.js` starts the bundled R on a free port, shows a loading screen
until the app is ready, shows the app in its own window, and stops R when the window is closed.

## Building the zip

On Linux, run `./build-windows-on-linux.sh`. It needs R 4.5, pip, npm and `innoextract`
(`sudo apt install innoextract`). It downloads the Windows builds of R, the R packages and Python, and writes the
zip to `dist/`. Nothing Windows can run on Linux, so test the zip on a Windows PC before releasing it - Windows
Sandbox is a good test, as it doesn't have the extra runtimes most PCs pick up from other software.

On Windows, run `.\build-windows.ps1` in PowerShell. It needs Node.js. It installs R and Python into `bundle/`,
checks the bundled R can load the Python packages, and writes the zip to `dist/`.

Set the version in `package.json` before building - it goes in the zip's name.

To try the Electron window using your own R and `app/python_venv/`, run `npm install` then `npm start` here.

## Releasing a new version

Build the zip, test it on Windows, then create a release on GitHub (Releases > Draft a new release), tag it with
the version (e.g. `v1.0.1`) and attach the zip.
