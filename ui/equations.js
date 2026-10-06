/* Symbolic state equations for the editor's equations window.
 *
 * Runs bondgraph.py (embedded by embed-python.js as self.BONDGRAPH_PY) with
 * SymPy in the browser through Pyodide, in a worker so a slow solve never
 * freezes the editor, and typesets its LaTeX with KaTeX. Both load from a CDN
 * the first time the window opens (about 11 MB, then cached).
 */
(function (root) {
  "use strict";

  // 0.29 is the newest Pyodide that runs in a classic worker; later ones need a
  // module worker, which Chrome won't start from a page opened as file://.
  const PYODIDE = "https://cdn.jsdelivr.net/pyodide/v0.29.3/full/";
  const KATEX = "https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.19.0/";
  const KATEX_SRI = {
    js: "sha384-QFFtAGzvvj+bfgCGxXJlNZZR1nXEZgvG8tDLCCY1F19xl20WlfTYgguB4VcNdxYk",
    css: "sha384-3rdsX6e5mueWyoweR9NIVmtEsUkokpBT/0ALqKKIBMr9j4qhHkaIkAcGgsE6uVlp",
  };
  // A solve running longer than this is abandoned (worker restarted) when a newer request arrives.
  const STUCK_MS = 8000;

  let katexReady = null;
  function loadKatex() {
    if (!katexReady) {
      katexReady = new Promise((resolve, reject) => {
        const css = document.createElement("link");
        Object.assign(css, { rel: "stylesheet", href: KATEX + "katex.min.css", integrity: KATEX_SRI.css, crossOrigin: "anonymous" });
        const js = document.createElement("script");
        Object.assign(js, { src: KATEX + "katex.min.js", integrity: KATEX_SRI.js, crossOrigin: "anonymous" });
        js.onload = () => resolve(root.katex);
        js.onerror = () => { katexReady = null; reject(new Error("couldn't load KaTeX")); };
        document.head.append(css, js);
      });
    }
    return katexReady;
  }

  // Body of the worker, run from a Blob URL (file:// pages can't start a worker from a file).
  function workerMain() {
    let booting = null, py = null;
    const post = (m) => self.postMessage(m);
    async function boot(cfg) {
      post({ type: "status", text: "Loading Python…" });
      importScripts(cfg.pyodide + "pyodide.js");
      py = await self.loadPyodide({ indexURL: cfg.pyodide });
      post({ type: "status", text: "Loading SymPy…" });
      await py.loadPackage("sympy", { messageCallback: () => {} });
      py.FS.writeFile("/home/pyodide/bondgraph.py", cfg.source);
      py.runPython("import sys; sys.path.insert(0, '/home/pyodide'); import bondgraph");
    }
    self.onmessage = async ({ data }) => {
      if (data.type === "init") {
        booting = boot(data);
        booting.catch((err) => post({ type: "failed", error: String((err && err.message) || err) }));
        return;
      }
      try { await booting; } catch (_) { return; }
      post({ type: "started", id: data.id });
      try {
        py.globals.set("bond_text", data.bond);
        py.globals.set("integral_json", JSON.stringify(data.integral));
        const out = py.runPython(
          "import json, bondgraph as bg\n" +
          "_m = bg.parse_text(bond_text)\n" +
          "json.dumps(bg.latex_report(_m, bg.storage_causality(_m, json.loads(integral_json))))");
        post({ type: "result", id: data.id, result: JSON.parse(out) });
      } catch (err) {
        // A PythonError carries the traceback; its last line is "ValueError: message".
        const last = String((err && err.message) || err).trim().split("\n").pop();
        post({ type: "result", id: data.id, error: last.replace(/^\w*(Error|Exception): /, "") });
      }
    };
  }

  // Latest-wins solver: while one solve runs, only the newest request waits,
  // and `done` hears back only about the newest.
  function createSolver(onStatus) {
    let worker = null, running = null, queued = null, latest = 0, done = null;
    function start() {
      const url = URL.createObjectURL(new Blob([`(${workerMain})()`], { type: "text/javascript" }));
      worker = new Worker(url);
      URL.revokeObjectURL(url);
      worker.onmessage = ({ data }) => {
        if (data.type === "status") onStatus(data.text);
        else if (data.type === "started") { if (running) running.since = Date.now(); onStatus("Solving…"); }
        else if (data.type === "failed") fail("Couldn't load SymPy (" + data.error + "). The equations window needs an internet connection the first time it opens.");
        else if (data.type === "result") {
          running = null;
          if (data.id === latest && done) done(data.error ? { error: data.error } : { result: data.result });
          if (queued) { const q = queued; queued = null; send(q); }
        }
      };
      worker.onerror = (e) => { e.preventDefault(); fail("The solver stopped: " + (e.message || "unknown error")); };
      worker.postMessage({ type: "init", pyodide: PYODIDE, source: root.BONDGRAPH_PY });
    }
    function fail(msg) {
      stop();
      if (done) done({ error: msg, fatal: true });
    }
    function stop() {
      if (worker) worker.terminate();
      worker = null; running = null; queued = null;
    }
    function send(req) {
      running = { id: req.id, since: Infinity }; // the clock starts once Python is loaded
      worker.postMessage(Object.assign({ type: "solve" }, req));
    }
    return {
      solve(bond, integral, callback) {
        done = callback;
        const req = { id: ++latest, bond, integral };
        if (!root.BONDGRAPH_PY) { callback({ error: "bondgraph-py.js is missing. Run: node ui/embed-python.js", fatal: true }); return; }
        if (running && Date.now() - running.since > STUCK_MS) stop();
        if (!worker) start();
        if (running) queued = req; else send(req);
      },
      stop,
    };
  }

  root.BondEquations = { loadKatex, createSolver };
})(typeof self !== "undefined" ? self : this);
