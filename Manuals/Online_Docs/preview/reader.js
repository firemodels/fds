"use strict";
const spanish = document.documentElement.lang === "es";
document.querySelector(".language-picker select")?.addEventListener("change", event => {
  const target = new URL(event.target.value, location.href);
  if (event.target.selectedOptions[0].dataset.keepSection !== "false") target.hash = location.hash;
  location.assign(target.href);
});
async function refreshLanguages() {
  const manifest = document.body.dataset.languageManifest;
  const select = document.querySelector(".language-picker select");
  if (!manifest || !select) return;
  try {
    const response = await fetch(manifest, {cache: "no-store"});
    if (!response.ok) return;
    const {languages} = await response.json();
    const page = location.pathname.split("/").pop() || "index.html";
    const options = languages.filter(edition => edition.pages.includes(page)).map(edition => {
      const option = document.createElement("option");
      option.textContent = edition.label;
      option.value = new URL(`${edition.code}/${page}`, new URL(manifest, location.href)).href;
      option.selected = edition.code === document.documentElement.lang;
      return option;
    });
    if (options.some(option => option.selected)) select.replaceChildren(...options);
  } catch { /* Keep the working, current-edition option if the registry is unavailable. */ }
}
refreshLanguages();
const dialog = document.querySelector("#search-dialog");
const input = document.querySelector("#search-input");
const results = document.querySelector("#search-results");
const status = document.querySelector("#search-status");
let searchIndex;
let searchRequest;
async function loadIndex() {
  if (!searchRequest) searchRequest = fetch("search-index.json").then(response => {
    if (!response.ok) throw new Error("Search index unavailable");
    return response.json();
  });
  searchIndex = await searchRequest;
}
async function openSearch() {
  dialog.showModal();
  input.focus();
  try { await loadIndex(); renderSearch(); }
  catch { status.textContent = spanish ? "No se pudo cargar la búsqueda. Recargue la página." : "Search could not load. Please reload the page."; }
}
document.querySelector(".search-trigger").addEventListener("click", openSearch);
document.addEventListener("keydown", event => {
  if ((event.metaKey || event.ctrlKey) && event.key.toLowerCase() === "k") {
    event.preventDefault();
    if (dialog.open) dialog.close(); else openSearch();
  }
});
dialog.addEventListener("click", event => { if (event.target === dialog) {
  const box = dialog.getBoundingClientRect();
  if (event.clientX < box.left || event.clientX > box.right || event.clientY < box.top || event.clientY > box.bottom) dialog.close();
}});
function renderSearch() {
  results.replaceChildren();
  const query = input.value.trim().toLowerCase();
  if (!query) { status.textContent = document.body.dataset.languageManifest ? (spanish ? "Buscar en este manual." : "Search this manual.") : (spanish ? "Buscar en el capítulo traducido." : "Search the chapters included in this preview."); return; }
  if (!searchIndex) return;
  const terms = query.split(/\s+/);
  const matches = searchIndex.map(item => {
    const text = (item.title + " " + item.text).toLowerCase();
    return {item, score: terms.every(term => text.includes(term)) ?
      1 + (item.title.toLowerCase().includes(query) ? 10 : 0) + (text.includes(query) ? 3 : 0) : 0};
  }).filter(result => result.score).sort((a, b) => b.score - a.score);
  status.textContent = spanish ? `${matches.length} ${matches.length === 1 ? "sección encontrada" : "secciones encontradas"}` : `${matches.length} matching section${matches.length === 1 ? "" : "s"}`;
  for (const {item} of matches.slice(0, 25)) {
    const link = document.createElement("a");
    link.className = "search-result";
    link.href = item.url;
    const chapter = document.createElement("small"); chapter.textContent = item.chapter;
    const title = document.createElement("strong"); title.textContent = item.title;
    const snippet = document.createElement("p");
    const at = Math.max(0, item.text.toLowerCase().indexOf(terms[0]) - 70);
    snippet.textContent = (at ? "…" : "") + item.text.slice(at, at + 240) + "…";
    link.append(chapter, title, snippet);
    link.addEventListener("click", () => dialog.close());
    results.append(link);
  }
}
input.addEventListener("input", renderSearch);
const menu = document.querySelector(".menu-toggle");
const sidebar = document.querySelector(".sidebar");
const mobile = window.matchMedia("(max-width: 720px)");
function setMenu(open, restoreFocus = false) {
  sidebar.classList.toggle("open", open);
  menu.setAttribute("aria-expanded", String(open));
  sidebar.inert = mobile.matches && !open;
  if (restoreFocus) menu.focus();
}
menu.addEventListener("click", () => setMenu(!sidebar.classList.contains("open")));
mobile.addEventListener("change", () => {
  const wasInside = sidebar.contains(document.activeElement);
  setMenu(false, mobile.matches && wasInside);
});
setMenu(false);
document.addEventListener("keydown", event => {
  // Escape first belongs to the modal search dialog.
  if (event.key === "Escape" && !dialog.open && mobile.matches && sidebar.classList.contains("open")) {
    setMenu(false, true);
  }
});
function updateScrollRegions() {
  for (const region of document.querySelectorAll(".table-scroll")) {
    const overflowing = region.scrollWidth > region.clientWidth + 1;
    if (overflowing) {
      region.setAttribute("tabindex", "0");
      region.setAttribute("role", "region");
      region.setAttribute("aria-label", spanish ? "Ecuación o tabla desplazable" : "Scrollable equation or table");
    } else {
      region.removeAttribute("tabindex");
      region.removeAttribute("role");
      region.removeAttribute("aria-label");
    }
  }
}
const resizeObserver = new ResizeObserver(updateScrollRegions);
for (const region of document.querySelectorAll(".table-scroll")) resizeObserver.observe(region);
window.addEventListener("load", updateScrollRegions);
document.fonts.ready.then(updateScrollRegions);
updateScrollRegions();
const sections = [...document.querySelectorAll("article h2[id], article h3[id]")];
const tocLinks = [...document.querySelectorAll(".page-toc nav a")];
let queued = false;
function highlightSection() {
  const current = [...sections].reverse().find(element => element.getBoundingClientRect().top <= 170) || sections[0];
  for (const link of tocLinks) {
    const active = current && link.hash === "#" + current.id;
    link.classList.toggle("current", !!active);
    if (active) link.setAttribute("aria-current", "location"); else link.removeAttribute("aria-current");
  }
  queued = false;
}
document.addEventListener("scroll", () => { if (!queued) { queued = true; requestAnimationFrame(highlightSection); } }, {passive:true});
highlightSection();
