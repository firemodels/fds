// Exercise availability filtering without a browser or network.
const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const source = fs.readFileSync(__dirname + '/reader.js', 'utf8');
const start = source.indexOf('async function refreshLanguages()');
const end = source.indexOf('\nrefreshLanguages();', start);
async function optionsFor(languages, reject = false) {
  const select = {options: ['current edition'], replaceChildren(...items) { this.options = items; }};
  const context = {
    URL, location: {href: 'http://localhost/en/mass-transport.html#section', pathname: '/en/mass-transport.html'},
    document: {body: {dataset: {languageManifest: '../languages.json'}}, documentElement: {lang: 'en'},
      querySelector: () => select, createElement: () => ({})},
    fetch: async () => { if (reject) throw new Error('unavailable'); return {ok: true, json: async () => ({languages})}; }
  };
  vm.createContext(context);
  vm.runInContext(source.slice(start, end), context);
  await context.refreshLanguages();
  return select.options;
}
async function manualsLink(href, available, spanish = false) {
  const header = {link: null, querySelector() { return this.link; }, prepend(link) { this.link = link; }};
  const context = {URL, spanish, location: {href},
    document: {body: {dataset: {languageManifest: '../languages.json'}},
      querySelector: () => header, createElement: () => ({})},
    fetch: async () => ({ok: available, json: async () => ({kind: 'fds-manual-catalog'})})};
  vm.createContext(context);
  const start = source.indexOf('async function addManualsNavigation()');
  const end = source.indexOf('\naddManualsNavigation();', start);
  vm.runInContext(source.slice(start, end), context);
  await context.addManualsNavigation();
  await context.addManualsNavigation();
  return header.link;
}
(async () => {
  const en = {code: 'en', label: 'English', pages: ['mass-transport.html']};
  const es = {code: 'es', label: 'Español', pages: ['overview.html']};
  let options = await optionsFor([en, es]);
  assert.equal(options.length, 1);
  assert.equal(options[0].value, 'http://localhost/en/mass-transport.html');
  es.pages.push('mass-transport.html');
  options = await optionsFor([en, es]);
  assert.equal(options.length, 2);
  assert.equal(options[1].value, 'http://localhost/es/mass-transport.html');
  assert.equal(options[0].selected, true);
  assert.deepEqual(await optionsFor([], true), ['current edition']);
  const link = await manualsLink('http://localhost/docs/user-guide/html/en/chapter-1.html', true);
  assert.equal(link.href, 'http://localhost/docs/index.html');
  assert.equal(link.textContent, 'All manuals');
  assert.equal((await manualsLink('http://localhost/docs/tech-guide/html/es/overview.html', true, true)).textContent, 'Todos los manuales');
  assert.equal(await manualsLink('http://localhost/en/overview.html', false), null);
  console.log('Language menu checks passed: page coverage, relative URLs, current selection, and unavailable registry.');
})().catch(error => { console.error(error); process.exitCode = 1; });
