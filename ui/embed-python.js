// Copies ../bondgraph.py into bondgraph-py.js, so the equations window can run
// it in the browser: a page opened from file:// can't fetch the .py file itself.
// Run after changing bondgraph.py (the engine tests check it's up to date):
//   node ui/embed-python.js
const fs = require("node:fs");
const path = require("node:path");

const OUT = path.join(__dirname, "bondgraph-py.js");

function render() {
  const src = fs.readFileSync(path.join(__dirname, "..", "bondgraph.py"), "utf8");
  return "// Generated from ../bondgraph.py by embed-python.js. Don't edit; rerun: node ui/embed-python.js\n" +
    "self.BONDGRAPH_PY = " + JSON.stringify(src) + ";\n";
}

if (require.main === module) {
  fs.writeFileSync(OUT, render());
  console.log("Wrote " + path.relative(process.cwd(), OUT));
}
module.exports = { render, OUT };
