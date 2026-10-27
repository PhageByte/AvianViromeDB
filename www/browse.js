"use strict";

(() => {
  const byId = (id) => document.getElementById(id);
  const size = 50;
  let catalog = null;
  let records = [];
  let filtered = [];
  let page = 1;
  const controls = ["query", "quality", "source", "min-length", "max-length", "coverage"];
  const missing = (value) => value === null || value === "" || value === "NA";
  const display = (value, index) => missing(value) ? "Unknown" :
    index >= 8 && String(value).startsWith("[") ? "Unassigned" :
    index === 2 ? Number(value).toLocaleString("en-US") : String(value);

  function render() {
    const pages = Math.max(1, Math.ceil(filtered.length / size));
    page = Math.min(Math.max(1, page), pages);
    const start = (page - 1) * size;
    const body = byId("results");
    body.replaceChildren();
    const fragment = document.createDocumentFragment();
    for (const row of filtered.slice(start, start + size)) {
      const tr = document.createElement("tr");
      row.forEach((value, index) => {
        const td = document.createElement("td");
        td.textContent = display(value, index);
        tr.append(td);
      });
      fragment.append(tr);
    }
    if (!filtered.length) {
      const tr = document.createElement("tr");
      const td = document.createElement("td");
      td.colSpan = catalog.fields.length;
      td.textContent = "No records match these filters. Adjust a filter or reset your selection.";
      tr.append(td);
      fragment.append(tr);
    }
    body.append(fragment);
    byId("result-count").textContent = `${filtered.length.toLocaleString("en-US")} of ${records.length.toLocaleString("en-US")} records`;
    byId("page-summary").textContent = filtered.length ? `Showing ${start + 1}–${Math.min(start + size, filtered.length)} · Page ${page} of ${pages}` : "0 records · 50 records per page";
    byId("first").disabled = byId("previous").disabled = page === 1;
    byId("next").disabled = byId("last").disabled = page === pages;
    byId("page-number").max = String(pages);
    byId("page-number").value = String(page);
    byId("page-number").disabled = byId("go").disabled = !filtered.length;
    byId("download").disabled = !filtered.length;
  }

  function applyFilters() {
    if (!catalog) return;
    const query = byId("query").value.trim().toLowerCase();
    const quality = byId("quality").value;
    const source = byId("source").value;
    const coverage = byId("coverage").value;
    const min = byId("min-length").value === "" ? 0 : Number(byId("min-length").value);
    const max = byId("max-length").value === "" ? Infinity : Number(byId("max-length").value);
    const valid = min >= 0 && max >= 0 && min <= max && byId("min-length").validity.valid && byId("max-length").validity.valid;
    byId("error").hidden = valid;
    byId("error").textContent = valid ? "" : "Enter whole-number length limits of zero or greater, with the minimum at or below the maximum.";
    filtered = valid ? records.filter((row) =>
      (!query || row[0].toLowerCase().includes(query) || row[1].toLowerCase().includes(query)) &&
      (!quality || row[3] === quality) && row[2] >= min && row[2] <= max &&
      (!source || (row[7] || "").split(";").includes(source)) && (!coverage || row[6] === coverage)) : [];
    page = 1;
    render();
  }

  controls.forEach((id) => byId(id).addEventListener("input", applyFilters));
  byId("filter-form").addEventListener("submit", (event) => event.preventDefault());
  byId("filter-form").addEventListener("reset", () => setTimeout(applyFilters, 0));
  for (const [id, move] of Object.entries({first: () => 1, previous: () => page - 1, next: () => page + 1, last: () => Math.ceil(filtered.length / size)})) {
    byId(id).addEventListener("click", () => { page = move(); render(); });
  }
  byId("page-form").addEventListener("submit", (event) => {
    event.preventDefault();
    const requested = Number(byId("page-number").value);
    if (Number.isInteger(requested)) { page = requested; render(); }
  });
  byId("download").addEventListener("click", () => {
    const cell = (value) => missing(value) ? "NA" : String(value).replace(/[\t\r\n]/g, " ");
    const tsv = [catalog.fields.join("\t"), ...filtered.map((row) => row.map(cell).join("\t"))].join("\n") + "\n";
    const url = URL.createObjectURL(new Blob([tsv], {type: "text/tab-separated-values;charset=utf-8"}));
    const link = document.createElement("a");
    link.href = url;
    link.download = "AvianViromeDB-metadata-preview-filtered.tsv";
    document.body.append(link);
    link.click();
    link.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  });

  fetch("catalog-summary.json").then((response) => {
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    return response.json();
  }).then((data) => {
    if (!Array.isArray(data.records) || data.records.length !== data.row_count || data.fields.length !== 10) throw new Error("Invalid catalog");
    catalog = data;
    records = data.records;
    const sources = new Set(records.flatMap((row) => missing(row[7]) ? [] : row[7].split(";")));
    [...sources].sort().forEach((source) => {
      const option = document.createElement("option");
      option.value = source;
      option.textContent = source;
      byId("source").append(option);
    });
    byId("filter-fields").disabled = false;
    applyFilters();
  }).catch(() => {
    byId("result-count").textContent = "Catalog unavailable";
    byId("error").hidden = false;
    byId("error").textContent = "The preview could not be loaded. Reload the page to try again, or use the database downloads.";
  });
})();
