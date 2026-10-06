// Typesets every formula of the engine descriptions with KaTeX in strict mode; run by tests/test_method_page.py.
import fs from 'fs';
import { createRequire } from 'module';
const root = process.argv[2];
const require = createRequire(import.meta.url);
const katex = require(root + '/ui/assets/katex/katex.min.js');
const methodology = JSON.parse(fs.readFileSync(root + '/ui/data/methodology.json', 'utf8'));
let n = 0, bad = 0;
for (const e of methodology.engines) {
  for (const t of e.technical) {
    for (const m of t.matchAll(/\$([^$]+)\$/g)) {
      n++;
      try { katex.renderToString(m[1], { throwOnError: true, displayMode: false, strict: 'error', trust: false }); }
      catch (err) { bad++; console.log('ERR', e.id, err.message.slice(0, 200), '\n   ', m[1].slice(0, 160)); }
    }
  }
}
console.log(n + ' formulas, ' + bad + ' errors');
process.exit(bad || !n ? 1 : 0);
