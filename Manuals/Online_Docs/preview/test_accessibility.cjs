// Deterministic state tests, not a substitute for browser/assistive technology.
const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const source = fs.readFileSync(__dirname + '/reader.js', 'utf8');
const events = {};
const viewport = {matches: true, addEventListener(type, fn) { this.change = fn; }};
function element() {
  const classes = new Set();
  return {attributes: {}, events: {},
    classList: {toggle(name, on) { on ? classes.add(name) : classes.delete(name); }, contains: name => classes.has(name)},
    setAttribute(name, value) { this.attributes[name] = value; },
    getAttribute(name) { return this.attributes[name]; },
    removeAttribute(name) { delete this.attributes[name]; },
    addEventListener(type, fn) { this.events[type] = fn; },
    focus() { context.document.activeElement = this; }};
}
const menu = element(), sidebar = element(), inside = element();
sidebar.contains = node => node === inside;
const wide = Object.assign(element(), {scrollWidth: 700, clientWidth: 300});
const fits = Object.assign(element(), {scrollWidth: 250, clientWidth: 300});
const context = {menu, sidebar, spanish: false, dialog: {open: false},
  window: {matchMedia: () => viewport, addEventListener() {}},
  ResizeObserver: class { observe() {} },
  document: {activeElement: null, fonts: {ready: Promise.resolve()},
    addEventListener(type, fn) { events[type] = fn; }, querySelectorAll: () => [wide, fits]}};
vm.createContext(context);
vm.runInContext(source.slice(source.indexOf('const mobile ='), source.indexOf('const sections =')), context);
assert.equal(sidebar.inert, true);
menu.events.click();
assert.equal(sidebar.inert, false);
assert.equal(menu.attributes['aria-expanded'], 'true');
context.dialog.open = true;
events.keydown({key: 'Escape'});
assert.equal(sidebar.inert, false, 'Search modal owns Escape first');
context.dialog.open = false;
inside.focus();
events.keydown({key: 'Escape'});
assert.equal(sidebar.inert, true);
assert.equal(context.document.activeElement, menu);
viewport.matches = false;
viewport.change();
assert.equal(sidebar.inert, false, 'Desktop navigation stays usable');
inside.focus();
viewport.matches = true;
viewport.change();
assert.equal(sidebar.inert, true);
assert.equal(context.document.activeElement, menu);
assert.equal(wide.attributes.tabindex, '0');
assert.equal(fits.attributes.tabindex, undefined);
wide.clientWidth = 800;
context.updateScrollRegions();
assert.equal(wide.attributes.tabindex, undefined);
assert.equal(wide.attributes['aria-label'], undefined);
console.log('Accessibility state checks passed: mobile visibility, Escape, resize focus, and overflow tab stops.');
