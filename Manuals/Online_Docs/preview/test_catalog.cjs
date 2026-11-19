// Test localization and edition routing with a small in-memory DOM substitute.
const fs = require("node:fs");
const vm = require("node:vm");
const assert = require("node:assert/strict");
const source = fs.readFileSync(__dirname + "/catalog.js", "utf8");
function node() {
  return {children: [], textContent: "", append(...items) { this.children.push(...items); },
    replaceChildren(...items) { this.children = items; },
    setAttribute(key, value) { this[key] = value; }};
}
const elements = new Map();
const document = {documentElement: {}, createElement: node,
  querySelector(selector) {
    if (!elements.has(selector)) elements.set(selector, node());
    return elements.get(selector);
  }};
const context = {document};
vm.createContext(context);
vm.runInContext(source.slice(0, source.indexOf("// Discover editions")), context);
const manuals = [
  {key: "tech", title: "Technical Reference Guide", title_es: "Guía de referencia técnica",
   subtitle: "Volume 1", subtitle_es: "Volumen 1",
   editions: [{code: "en", url: "tech/en/index.html"}, {code: "es", url: "tech/es/index.html"}]},
  {key: "user", title: "User Guide", title_es: "Guía del usuario",
   subtitle: "Fire Dynamics Simulator", subtitle_es: "Simulador de dinámica de incendios",
   editions: [{code: "en", url: "user/en/index.html"}]}
];
context.renderCatalog(manuals, "es");
assert.equal(document.documentElement.lang, "es");
assert.equal(elements.get("h1").textContent, "Manuales de FDS en línea");
let cards = elements.get("#manuals").children.flatMap(section => section.children[1].children);
assert.equal(cards[0].children[1].children[0].textContent, "Guía de referencia técnica");
assert.equal(cards[0].children[1].children[0].href, "tech/es/index.html");
assert.equal(cards[1].children[1].children[0].href, "user/en/index.html");
assert.equal(cards[1].children[2].textContent, "Disponible en inglés");
context.renderCatalog(manuals, "en");
cards = elements.get("#manuals").children.flatMap(section => section.children[1].children);
assert.equal(cards[0].children[1].children[0].href, "tech/en/index.html");
assert.equal(cards[1].children[1].children[0].textContent, "User Guide");
assert.equal(cards[1].children.length, 2);
console.log("Catalog localization and edition routing checks passed.");
