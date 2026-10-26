// Copies the examples that have a "#@title" line (examples/*.bond) into
// examples.js, for the editor's "Load an example" menu: a page opened from
// file:// can't read them itself. Run after changing one (the engine tests
// check it's up to date):
//   node ui/embed-examples.js
const fs = require("node:fs");
const path = require("node:path");

const DIR = path.join(__dirname, "..", "examples");
const OUT = path.join(__dirname, "examples.js");

function render() {
  const examples = fs.readdirSync(DIR).filter((f) => f.endsWith(".bond")).sort().flatMap((file) => {
    const text = fs.readFileSync(path.join(DIR, file), "utf8");
    const title = /^#@title\s+(.+)$/m.exec(text);
    return title ? [{ file, title: title[1].trim(), text }] : [];
  }).sort((a, b) => a.title.localeCompare(b.title));
  return "// Generated from ../examples/*.bond by embed-examples.js. Don't edit; rerun: node ui/embed-examples.js\n" +
    "self.BOND_EXAMPLES = " + JSON.stringify(examples, null, 1) + ";\n";
}

if (require.main === module) {
  fs.writeFileSync(OUT, render());
  console.log("Wrote " + path.relative(process.cwd(), OUT));
}
module.exports = { render, OUT };
