"use strict";
const catalogLabels = {
  en: {
    title: "FDS online manuals", eyebrow: "DOCUMENTATION LIBRARY",
    intro: "Choose a manual to start reading.", language: "Language",
    technical: "TECHNICAL DOCUMENTATION", skip: "Skip to content",
    empty: "No manuals have completed building yet.",
    fallback: "Available in English",
    technicalGroup: "NIST SP 1018 · Technical Reference Guide",
    userGroup: "NIST SP 1019 · User Guide",
    footer: "These manuals are generated from the FDS LaTeX sources. Each manual includes search, chapter navigation, and access to its PDF edition."
  },
  es: {
    title: "Manuales de FDS en línea", eyebrow: "BIBLIOTECA DE DOCUMENTACIÓN",
    intro: "Seleccione un manual para comenzar a leer.", language: "Idioma",
    technical: "DOCUMENTACIÓN TÉCNICA", skip: "Saltar al contenido",
    empty: "Todavía no hay manuales disponibles.",
    fallback: "Disponible en inglés",
    technicalGroup: "NIST SP 1018 · Guía de referencia técnica",
    userGroup: "NIST SP 1019 · Guía del usuario",
    footer: "Estos manuales se generan a partir de los archivos LaTeX de FDS. Cada manual incluye búsqueda, navegación por capítulos y acceso a su edición PDF."
  }
};

function catalogEdition(manual, language) {
  return manual.editions.find(e => e.code === language) ||
    manual.editions.find(e => e.code === "en") || manual.editions[0];
}

function renderCatalog(manuals, language) {
  const labels = catalogLabels[language] || catalogLabels.en;
  const spanish = language === "es";
  document.documentElement.lang = language;
  document.title = labels.title;
  document.querySelector("h1").textContent = labels.title;
  document.querySelector(".eyebrow").textContent = labels.eyebrow;
  document.querySelector(".intro").textContent = labels.intro;
  document.querySelector(".language-picker>span").textContent = labels.language;
  document.querySelector("#catalog-language").setAttribute("aria-label", labels.language);
  document.querySelector(".brand small").textContent = labels.technical;
  document.querySelector(".skip-link").textContent = labels.skip;
  document.querySelector("footer").textContent = labels.footer;
  const sections = [];
  for (const [group, heading] of [["technical", labels.technicalGroup], ["user", labels.userGroup]]) {
    const members = manuals.filter(m => (m.key === "user") === (group === "user"));
    if (!members.length) continue;
    const section = document.createElement("section");
    section.className = "publication";
    const groupTitle = document.createElement("h2");
    groupTitle.textContent = heading;
    const grid = document.createElement("div");
    grid.className = "manual-grid";
    const cards = members.map(manual => {
    const edition = catalogEdition(manual, language);
    const card = document.createElement("section");
    card.className = "manual-card";
    const volume = document.createElement("p");
    volume.className = "volume";
    volume.textContent = spanish ? manual.subtitle_es : manual.subtitle;
    const title = document.createElement("h3");
    const link = document.createElement("a");
    link.href = edition.url;
    link.textContent = spanish ? manual.title_es : manual.title;
    title.append(link);
    card.append(volume, title);
    if (edition.code !== language) {
      const note = document.createElement("p");
      note.className = "edition-note";
      note.textContent = labels.fallback;
      card.append(note);
    }
    return card;
    });
    grid.append(...cards);
    section.append(groupTitle, grid);
    sections.push(section);
  }
  document.querySelector("#manuals").replaceChildren(...sections);
  document.querySelector("#empty").textContent = manuals.length ? "" : labels.empty;
}

// Discover editions at load time, including builds started before the catalog.
(async function refreshManuals() {
  try {
    const response = await fetch("manuals.json", {cache: "no-store"});
    if (!response.ok) return;
    const registry = await response.json();
    if (registry.kind !== "fds-manual-catalog") return;
    const found = await Promise.all(registry.locations.map(async manual => {
      const editions = await Promise.all(Object.entries(registry.languages).map(async ([code, label]) => {
        const folder = manual.output + "/html/" + code + "/";
        try {
          const response = await fetch(folder + "build-manifest.json", {cache: "no-store"});
          if (!response.ok) return null;
          const proof = await response.json();
          if (!proof.complete || proof.language !== code || (proof.manual || "tech") !== manual.key || !proof.pages?.length) return null;
          const index = await fetch(folder + "index.html", {method: "HEAD", cache: "no-store"});
          return index.ok ? {code, label, url: folder + "index.html"} : null;
        } catch { return null; }
      }));
      const available = editions.filter(Boolean);
      return available.length ? {...manual, editions: available} : null;
    }));
    const manuals = found.filter(Boolean);
    const available = new Set(manuals.flatMap(m => m.editions.map(e => e.code)));
    const picker = document.querySelector("#catalog-language");
    const requested = new URL(location.href).searchParams.get("language") || picker.value || "en";
    picker.replaceChildren(...Object.entries(registry.languages).filter(([code]) => available.has(code)).map(([code, label]) => {
      const option = document.createElement("option");
      option.value = code;
      option.textContent = label;
      return option;
    }));
    const language = available.has(requested) ? requested : (picker.value || "en");
    picker.value = language;
    renderCatalog(manuals, language);
    picker.addEventListener("change", () => {
      renderCatalog(manuals, picker.value);
      const url = new URL(location.href);
      url.searchParams.set("language", picker.value);
      history.replaceState(null, "", url);
    });
  } catch { /* Keep the statically generated library on fetch failure. */ }
})();
