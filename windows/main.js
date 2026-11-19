const { app, BrowserWindow, dialog, shell, Menu } = require('electron');
const { spawn, execFileSync } = require('child_process');
const fs = require('fs');
const http = require('http');
const net = require('net');
const os = require('os');
const path = require('path');

// The packaged app carries its own R, Python and copy of the shiny app in resources/.
// Running `npm start` from the repo uses the system R and app/python_venv instead.
const res = app.isPackaged ? process.resourcesPath : null;
const rscript = res ? path.join(res, 'R', 'bin', 'Rscript.exe') : 'Rscript';
const python = res ? path.join(res, 'python', 'python.exe')
  : path.join(__dirname, '..', 'app', 'python_venv', process.platform === 'win32' ? 'Scripts/python.exe' : 'bin/python');
const appDir = res ? path.join(res, 'shiny') : path.join(__dirname, '..', 'app');

const logFile = path.join(os.tmpdir(), 'shields-skin-dose-calculator.log');
let rProcess = null;
let win = null;
let quitting = false;
let ready = false;

function log(line) {
  fs.appendFileSync(logFile, line + '\n');
}

function lastLogLines(n) {
  try {
    return fs.readFileSync(logFile, 'utf8').trim().split('\n').slice(-n).join('\n');
  } catch {
    return '';
  }
}

function freePort() {
  return new Promise((resolve, reject) => {
    const srv = net.createServer();
    srv.on('error', reject);
    srv.listen(0, '127.0.0.1', () => {
      const port = srv.address().port;
      srv.close(() => resolve(port));
    });
  });
}

function startR(port) {
  const env = { ...process.env, SKIN_DOSE_PYTHON: python, RETICULATE_PYTHON: python, PYTHONNOUSERSITE: '1' };
  // Stop anything the user already has installed (their own R, R packages or Python) from being picked up instead of ours
  for (const k of ['R_HOME', 'R_LIBS', 'R_LIBS_SITE', 'PYTHONHOME', 'PYTHONPATH']) delete env[k];
  if (res) env.R_LIBS_USER = path.join(appDir, 'no-user-library'); // R ignores library paths that don't exist

  const code = `shiny::runApp('${appDir.replace(/\\/g, '/')}', port = ${port}, host = '127.0.0.1', launch.browser = FALSE)`;
  rProcess = spawn(rscript, ['--vanilla', '-e', code], { cwd: appDir, env, windowsHide: true });
  rProcess.stdout.on('data', (d) => fs.appendFileSync(logFile, d));
  rProcess.stderr.on('data', (d) => fs.appendFileSync(logFile, d));
  rProcess.on('error', (err) => {
    log(`Could not start R: ${err.message}`);
    rProcess = null;
  });
  rProcess.on('exit', (code) => {
    log(`R exited with code ${code}`);
    rProcess = null;
    // before the app is ready, waitForApp notices R has gone and reports it
    if (ready && !quitting) {
      dialog.showErrorBox('Shields Skin Dose Calculator',
        `The calculator stopped unexpectedly.\n\nDetails were saved to:\n${logFile}\n\n${lastLogLines(15)}`);
      app.quit();
    }
  });
}

function stopR() {
  if (!rProcess) return;
  // Rscript starts R as a child process, so kill the whole tree on Windows
  if (process.platform === 'win32') {
    try { execFileSync('taskkill', ['/pid', String(rProcess.pid), '/T', '/F'], { windowsHide: true }); } catch {}
  } else {
    rProcess.kill();
  }
  rProcess = null;
}

function waitForApp(url, timeoutMs) {
  const end = Date.now() + timeoutMs;
  return new Promise((resolve, reject) => {
    const attempt = () => {
      if (!rProcess) return reject(new Error('R stopped before the app was ready'));
      http.get(url, (r) => {
        r.resume();
        if (r.statusCode === 200) resolve();
        else retry();
      }).on('error', retry);
    };
    const retry = () => (Date.now() > end ? reject(new Error('Timed out waiting for the app to start')) : setTimeout(attempt, 500));
    attempt();
  });
}

async function createWindow() {
  fs.writeFileSync(logFile, `Shields Skin Dose Calculator ${app.getVersion()} started ${new Date().toISOString()}\n`);
  Menu.setApplicationMenu(null);
  win = new BrowserWindow({ width: 1500, height: 950, title: 'Shields Skin Dose Calculator', show: false });
  win.once('ready-to-show', () => win.show());
  win.on('closed', () => { win = null; });
  win.loadFile(path.join(__dirname, 'loading.html'));

  const port = await freePort();
  const url = `http://127.0.0.1:${port}/`;
  startR(port);

  // Shiny's download buttons open in a new tab, so turn those into downloads (Electron asks where to save).
  // Any other link goes to the user's normal browser.
  win.webContents.setWindowOpenHandler(({ url: target }) => {
    if (target.startsWith(url)) win.webContents.downloadURL(target);
    else shell.openExternal(target);
    return { action: 'deny' };
  });
  // Keep the app's own page titled as the app rather than the shiny page title
  win.on('page-title-updated', (e) => e.preventDefault());

  try {
    await waitForApp(url, 180000);
  } catch (err) {
    log(err.message);
    if (!quitting) {
      quitting = true;
      dialog.showErrorBox('Shields Skin Dose Calculator',
        `The calculator could not start.\n\nDetails were saved to:\n${logFile}\n\n${lastLogLines(15)}`);
      stopR();
      app.quit();
    }
    return;
  }
  ready = true;
  log(`Shiny app ready at ${url}`);
  if (win) win.loadURL(url);
}

if (!app.requestSingleInstanceLock()) {
  app.quit();
} else {
  app.on('second-instance', () => {
    if (win) {
      if (win.isMinimized()) win.restore();
      win.focus();
    }
  });
  app.whenReady().then(createWindow);
  app.on('window-all-closed', () => app.quit());
  app.on('before-quit', () => {
    quitting = true;
    stopR();
  });
}
