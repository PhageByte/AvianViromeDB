const $ = id => document.getElementById(id);
const fmt = value => value.toLocaleString('en-US');
const {Deck, OrthographicView, ScatterplotLayer, LineLayer, COORDINATE_SYSTEM} = deck;
const params = new URLSearchParams(location.search);
const palette = {vcx: [91, 143, 249], catalog: [184, 184, 184], refseq: [0, 158, 115], missing: [205, 211, 215]};
const common = {coordinateSystem: COORDINATE_SYSTEM.CARTESIAN, parameters: {depthCompare: 'always', depthWriteEnabled: false}};
let nodes, modules, moduleMap, members, positions, viewer, worker, viewState;
let nodeLayer, vcxLayer, backboneLayer, moduleLayer, fullLayer, selectionLayer, selectedLayer;
let selection = null, selected = new Set(), neighbors = new Set(), selectionResult = null;
let ready = false, request = 0, revision = 0, fullRequested = false, urlTimer, exportRequested = false;
let firstFrame, frameStarted = 0, fullFrames = 0, fullTime = 0;
let color = 'vcx', edgeMode = 'selection', searchRows = [], matches = [], pendingSelection = params.get('sel');
const allIds = [], vcxIds = [];

function categorical(id) {
  const hue = ((id * 137.508) % 360) / 60;
  return [5, 3, 1].map(n => Math.round(255 * (.46 - .30 * Math.max(-1, Math.min((n + hue) % 6 - 3, 9 - (n + hue) % 6, 1)))));
}

function nodeColor(i) {
  let rgb;
  if (color === 'vcx') rgb = nodes.module[i] === 44 ? palette.vcx : nodes.set[i] ? palette.refseq : palette.catalog;
  else if (color === 'set') rgb = nodes.set[i] ? palette.refseq : palette.catalog;
  else rgb = nodes[color][i] < 0 ? palette.missing : categorical(nodes[color][i]);
  return [...rgb, selection && !selected.has(i) && !neighbors.has(i) ? 65 : 230];
}

function binaryPoints(ids) {
  const value = new Float32Array(ids.length * 2);
  ids.forEach((i, j) => value.set(positions.subarray(i * 2, i * 2 + 2), j * 2));
  return {length: ids.length, attributes: {getPosition: {value, size: 2}}};
}

function lineLayer(id, source, target, props) {
  return new LineLayer({...common, id, pickable: false, widthUnits: 'pixels',
    data: {length: source.length / 2, attributes: {
      getSourcePosition: {value: source, size: 2}, getTargetPosition: {value: target, size: 2}}}, ...props});
}

function tooltip(text, info) {
  $('tooltip').hidden = !text;
  if (!text) return;
  $('tooltip').textContent = text;
  $('tooltip').style.left = `${Math.max(4, Math.min(info.x + 12, $('viewport').clientWidth - 270))}px`;
  $('tooltip').style.top = `${Math.max(4, Math.min(info.y + 12, $('viewport').clientHeight - 65))}px`;
}

function nodeHover(info, ids) {
  const i = info.index < 0 ? -1 : ids ? ids[info.index] : info.index;
  tooltip(i < 0 ? '' : `${nodes.vOTU[i] || nodes.label[i]} · ${nodes.label[i]}\n${moduleName(nodes.module[i])}`, info);
}

function moduleName(id) { return id === 44 ? 'VCX · Module 44' : id < 0 ? 'No assigned module' : `Module ${id}`; }
function countNote(id) {
  const note = `Edge counts use the ${fmt(nodes.edges)} edges drawn in Fig. 6b.`;
  return id === 44 ? `${note} Table S9 counts VCX edges in the vConTACT3 network that defined the modules: 25,590 to catalog vOTUs in other modules and 156 to RefSeq genomes.` : note;
}

function createLayers() {
  nodeLayer = new ScatterplotLayer({...common, id: 'nodes',
    data: {length: nodes.n, attributes: {getPosition: {value: positions, size: 2}}},
    radiusUnits: 'pixels', getRadius: (_, {index}) => selected.has(index) ? 3.2 : neighbors.has(index) ? 2.5 : 1.5,
    getFillColor: (_, {index}) => nodeColor(index), pickable: true,
    onHover: info => nodeHover(info), onClick: info => ready && info.index >= 0 && select({kind: 'node', id: info.index})});
  vcxLayer = new ScatterplotLayer({...common, id: 'vcx-outline', data: binaryPoints(vcxIds),
    radiusUnits: 'pixels', getRadius: 2.5, stroked: true, lineWidthUnits: 'pixels', getLineWidth: .65,
    getFillColor: (_, {index}) => nodeColor(vcxIds[index]), getLineColor: [29, 53, 87, 230], pickable: true,
    onHover: info => nodeHover(info, vcxIds), onClick: info => ready && info.index >= 0 && select({kind: 'node', id: vcxIds[info.index]})});
  const source = new Float32Array(modules.backbone.length * 2), target = new Float32Array(source.length);
  modules.backbone.forEach((edge, i) => {
    const a = moduleMap.get(edge.m1), b = moduleMap.get(edge.m2);
    source.set([a.x, a.y], i * 2); target.set([b.x, b.y], i * 2);
  });
  backboneLayer = lineLayer('backbone', source, target, {visible: edgeMode === 'backbone', getColor: [26, 75, 107, 75],
    getWidth: (_, {index}) => Math.min(5, .4 + Math.sqrt(modules.backbone[index].edges) / 35)});
  moduleLayer = new ScatterplotLayer({...common, id: 'module-centres', visible: edgeMode === 'backbone',
    data: {length: modules.modules.length, attributes: {getPosition: {
      value: new Float32Array(modules.modules.flatMap(m => [m.x, m.y])), size: 2}}},
    radiusUnits: 'pixels', getRadius: 6, getFillColor: [255, 255, 255, 210], stroked: true,
    lineWidthUnits: 'pixels', getLineWidth: 1.5, getLineColor: [13, 107, 109], pickable: true,
    onHover: info => tooltip(info.index < 0 ? '' : `${moduleName(modules.modules[info.index].module_id)}\nClick to select module`, info),
    onClick: info => {
      if (ready && info.index >= 0) {
        $('edges').value = 'selection'; setEdgeMode();
        select({kind: 'mod', id: modules.modules[info.index].module_id}, true);
      }
    }});
}

function refreshNodes() {
  // clone preserves the layer's data identity and GPU position buffer; only these accessors are invalidated.
  nodeLayer = nodeLayer.clone({updateTriggers: {getFillColor: [color, revision], getRadius: revision}});
  vcxLayer = vcxLayer.clone({visible: color === 'vcx', updateTriggers: {getFillColor: [color, revision]}});
  const ids = [...selected];
  selectedLayer = ids.length ? new ScatterplotLayer({...common, id: 'selected', data: binaryPoints(ids),
    radiusUnits: 'pixels', getRadius: 3.6, getFillColor: (_, {index}) => nodeColor(ids[index]),
    stroked: true, lineWidthUnits: 'pixels', getLineWidth: 1.2, getLineColor: [29, 53, 87], pickable: true,
    onHover: info => nodeHover(info, ids), onClick: info => info.index >= 0 && select({kind: 'node', id: ids[info.index]})}) : null;
}

function render() {
  viewer.setProps({layers: [fullLayer, backboneLayer, selectionLayer, nodeLayer, vcxLayer, selectedLayer, moduleLayer].filter(Boolean)});
}

function fit(ids) {
  let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
  for (const i of ids) {
    minX = Math.min(minX, nodes.x[i]); maxX = Math.max(maxX, nodes.x[i]);
    minY = Math.min(minY, nodes.y[i]); maxY = Math.max(maxY, nodes.y[i]);
  }
  const w = $('viewport').clientWidth, h = $('viewport').clientHeight;
  const mobile = matchMedia('(max-width: 720px)').matches, panel = $('inspector');
  const right = !panel.hidden && !mobile ? panel.offsetWidth + 28 : 0;
  const bottom = mobile ? (panel.hidden ? 0 : panel.offsetHeight + 16) + $('legend').offsetHeight + 16 : 0;
  const usableW = Math.max(100, w - right - 64), usableH = Math.max(100, h - bottom - 100);
  const scale = Math.min(usableW / Math.max(240, (maxX - minX) * 1.16), usableH / Math.max(240, (maxY - minY) * 1.16));
  viewState = {target: [(minX + maxX) / 2 + right / (2 * scale), (minY + maxY) / 2 - bottom / (2 * scale), 0],
    zoom: Math.min(3, Math.log2(scale)), minZoom: -8, maxZoom: 6};
  viewer?.setProps({viewState});
  scheduleURL();
}

function saveURL() {
  clearTimeout(urlTimer);
  const url = new URL(location.href);
  for (const [key, value] of Object.entries({x: viewState.target[0].toFixed(3), y: viewState.target[1].toFixed(3),
    z: viewState.zoom.toFixed(5), color, edges: edgeMode})) url.searchParams.set(key, value);
  if (selection) url.searchParams.set('sel', `${selection.kind}:${selection.kind === 'node' ? nodes.label[selection.id] : selection.id}`);
  else if (pendingSelection) url.searchParams.set('sel', pendingSelection);
  else url.searchParams.delete('sel');
  history.replaceState(null, '', url);
}

function scheduleURL() { clearTimeout(urlTimer); urlTimer = setTimeout(saveURL, 180); }

function element(tag, text, className) {
  const el = document.createElement(tag);
  if (text !== undefined) el.textContent = text;
  if (className) el.className = className;
  return el;
}

function selectionButton(text, value) {
  const button = element('button', text, 'text-button');
  button.onclick = () => select(value, true);
  return button;
}

function attributes(rows) {
  const dl = element('dl');
  for (const [key, value] of rows) {
    const dd = element('dd'); dd.append(value instanceof Node ? value : String(value || 'Unknown'));
    dl.append(element('dt', key), dd);
  }
  return dl;
}

function inspect() {
  const panel = $('details'); panel.replaceChildren();
  if (!selection) {
    panel.append(element('h2', 'Explore Fig. 6b'), element('p', 'Select a node to inspect its genome and neighbors, or choose an Infomap module.'),
      element('p', `${fmt(modules.modules.length)} modules · ${fmt(nodes.edges)} gene-sharing edges`, 'muted'),
      element('p', 'Search a vOTU, contig, RefSeq accession, genome name, family, “VCX” or “module 44”.'));
    return;
  }
  if (selection.kind === 'node') {
    const i = selection.id, m = nodes.module[i];
    panel.append(element('h2', nodes.vOTU[i] || nodes.label[i]), element('div', nodes.set[i] ? 'RefSeq genome' : 'Catalog sequence', 'badge'),
      attributes([['Label', nodes.label[i]], ['Name', nodes.name[i]], ['vOTU', nodes.vOTU[i]],
        ['Module', m < 0 ? 'Unassigned' : selectionButton(moduleName(m), {kind: 'mod', id: m})],
        ['Length', nodes.length[i] < 0 ? 'Unknown' : `${fmt(nodes.length[i])} bp`], ['CheckV', nodes.quality[i]],
        ['Lifestyle', nodes.lifestyle[i]], ['Bird sources', nodes.birds[i].split(';').join(', ')],
        ['Family', nodes.family_levels[nodes.family[i]]], ['Class', nodes.class_levels[nodes.class[i]]],
        ['Host phylum', nodes.host_phylum_levels[nodes.host_phylum[i]]]]));
  } else {
    const m = moduleMap.get(selection.id);
    panel.append(element('h2', moduleName(m.module_id)), element('div', `${fmt(m.members)} nodes selected`, 'badge'),
      attributes([['Members', fmt(m.members)], ['Catalog', fmt(m.catalog)], ['RefSeq', fmt(m.refseq)],
        ['Composition', m.composition], ['Internal edges', fmt(m.internal)]]));
  }
  if (selectionResult) {
    panel.append(element('p', `${fmt(selectionResult.source.length / 2)} incident edges · ${fmt(neighbors.size)} neighboring nodes`, 'muted'),
      element('p', `${fmt(selectionResult.internalCount)} within-module · ${fmt(selectionResult.externalCount)} between-module edges`, 'muted'));
    if (selection.kind === 'mod') panel.append(element('p', countNote(selection.id), 'muted'));
    panel.append(element('h3', 'Top neighboring modules'));
    const list = element('ul', undefined, 'neighbors');
    for (const [id, count] of selectionResult.shared.slice(0, 10)) {
      const row = element('li');
      row.append(id < 0 ? element('span', 'Unassigned') : selectionButton(moduleName(id), {kind: 'mod', id}), element('small', `${fmt(count)} edges`));
      list.append(row);
    }
    panel.append(list.childElementCount ? list : element('p', 'No edges to other modules.'));
  } else panel.append(element('p', 'Calculating incident edges…', 'muted'));
  if (selection.kind === 'mod') {
    const details = element('details'), list = element('div', undefined, 'member-list');
    details.append(element('summary', `Browse all ${fmt(selected.size)} members`));
    for (const i of selected) list.append(selectionButton(`${nodes.vOTU[i] || nodes.label[i]} · ${nodes.label[i]}`, {kind: 'node', id: i}));
    details.append(list); panel.append(details);
  }
}

function select(value, frame = false) {
  if (!ready) return;
  selection = value; selectionResult = null; selectionLayer = null; neighbors = new Set();
  selected = new Set(value.kind === 'mod' ? members.get(value.id) : [value.id]);
  $('viewport').dataset.selection = `${value.kind}:${value.id}`;
  $('viewport').dataset.selectedCount = selected.size;
  $('viewport').dataset.selectionReady = 'false';
  $('inspector').hidden = false; $('toggle').setAttribute('aria-expanded', 'true');
  $('suggestions').hidden = true; $('search').setAttribute('aria-expanded', 'false');
  $('progress').textContent = 'Calculating selection edges…';
  revision++; refreshNodes(); inspect(); render();
  if (frame) fit([...selected]);
  saveURL();
  worker.postMessage({type: 'select', selection, request: ++request});
}

function whole() {
  pendingSelection = null;
  selection = null; selected.clear(); neighbors.clear(); selectionResult = null; selectionLayer = null; request++;
  $('viewport').dataset.selection = ''; $('viewport').dataset.selectedCount = '0';
  $('progress').textContent = ready ? 'Selection ready' : 'Loading adjacency data…';
  revision++; refreshNodes(); inspect(); fit(allIds); render();
}

function legend() {
  $('legend-title').textContent = $('color').selectedOptions[0].textContent;
  const items = $('legend-items'); items.replaceChildren();
  const rows = color === 'vcx' ? [['VCX (266)', palette.vcx], ['Other catalog', palette.catalog], ['RefSeq', palette.refseq]] :
    color === 'set' ? [['Catalog', palette.catalog], ['RefSeq', palette.refseq]] :
    color === 'module' ? modules.modules.map(m => [moduleName(m.module_id), categorical(m.module_id), m.module_id]) :
    nodes[`${color}_levels`].map((name, i) => [name, categorical(i)]);
  if (!['vcx', 'set'].includes(color)) rows.push(['Unassigned', palette.missing]);
  for (const [name, rgb, id] of rows) {
    const item = id === undefined ? element('span') : selectionButton('', {kind: 'mod', id});
    const swatch = element('i', undefined, 'swatch'); swatch.style.background = `rgb(${rgb.join(',')})`;
    item.append(swatch, document.createTextNode(name)); items.append(item);
  }
  $('edge-legend').textContent = {selection: 'Edges: dark blue within module · pale blue between modules',
    backbone: 'Module centres: teal rings · Lines: ≥50 shared edges; width shows count', all: 'All 2,301,056 edges · faint gray', none: 'Edges hidden'}[edgeMode];
}

function setEdgeMode() {
  edgeMode = $('edges').value;
  backboneLayer = backboneLayer.clone({visible: edgeMode === 'backbone'});
  moduleLayer = moduleLayer.clone({visible: edgeMode === 'backbone'});
  if (selectionLayer) selectionLayer = selectionLayer.clone({visible: edgeMode === 'selection'});
  $('notice').textContent = edgeMode === 'all' ? 'All edges: 2.3 million lines; rendering may be slow. Use Selection or Backbone for navigation.' : '';
  if (edgeMode === 'all' && !fullRequested) {
    fullRequested = true; $('progress').textContent = 'Preparing all edges…'; worker.postMessage({type: 'all'});
  }
  if (fullLayer) fullLayer = fullLayer.clone({visible: edgeMode === 'all'});
  legend(); render(); scheduleURL();
}

function search() {
  const query = $('search').value.trim().toLowerCase();
  const id = /^(?:module\s+)?(\d+)$/.exec(query);
  matches = [];
  if (query === 'vcx' || (id && moduleMap.has(Number(id[1])))) {
    const mod = query === 'vcx' ? 44 : Number(id[1]);
    matches.push({value: {kind: 'mod', id: mod}, title: moduleName(mod), detail: `${moduleMap.get(mod).members} members`});
  } else if (query) {
    const exact = [], prefix = [];
    searchRows.forEach((fields, i) => {
      const row = {value: {kind: 'node', id: i}, title: nodes.vOTU[i] || nodes.label[i],
        detail: `${nodes.label[i]} · ${nodes.name[i]} · ${nodes.family_levels[nodes.family[i]] || 'Unassigned family'}`};
      if (fields.some(field => field === query)) exact.push(row);
      else if (fields.some(field => field.startsWith(query))) prefix.push(row);
    });
    matches = exact.concat(prefix).slice(0, 8);
  }
  const list = $('suggestions'); list.replaceChildren();
  for (const match of matches) {
    const li = element('li'), button = element('button', match.title); button.type = 'button';
    button.append(element('small', match.detail)); button.onclick = () => select(match.value, true);
    li.append(button); list.append(li);
  }
  if (query && !matches.length) list.append(element('li', 'No matches. Try a vOTU, family or module number.'));
  list.hidden = !query; $('search').setAttribute('aria-expanded', String(!!query));
}

function downloadPNG() {
  const canvas = $('network-canvas'), output = document.createElement('canvas');
  output.width = canvas.width; output.height = canvas.height;
  const ctx = output.getContext('2d'); ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, output.width, output.height);
  ctx.drawImage(canvas, 0, 0);
  output.toBlob(blob => {
    const url = URL.createObjectURL(blob), link = element('a');
    link.href = url; link.download = 'AvianViromeDB-network.png'; link.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }, 'image/png');
}

function workerMessage({data}) {
  if (data.type === 'error') { $('progress').textContent = `Adjacency data failed: ${data.message}`; return; }
  if (data.type === 'ready') {
    ready = true; $('viewport').dataset.csrReady = 'true';
    for (const id of ['vcx', 'search']) $(id).disabled = false;
    $('edges').querySelector('[value="all"]').disabled = false;
    $('progress').textContent = 'Selection ready';
    const saved = pendingSelection;
    pendingSelection = null;
    if (saved) {
      const split = saved.indexOf(':'), kind = saved.slice(0, split), key = saved.slice(split + 1);
      const id = kind === 'node' ? nodes.label.indexOf(key) : Number(key);
      if ((kind === 'node' && id >= 0) || (kind === 'mod' && moduleMap.has(id))) select({kind, id}, !['x', 'y', 'z'].every(k => params.has(k)));
      else $('notice').textContent = 'The linked selection was not found.';
    }
    if (params.get('edges') === 'all') { $('edges').value = 'all'; setEdgeMode(); }
  }
  if (data.type === 'selection' && data.request === request) {
    selectionResult = data; neighbors = new Set(data.neighbors);
    selectionLayer = lineLayer('selection-edges', data.source, data.target, {visible: edgeMode === 'selection',
      getColor: (_, {index}) => data.internal[index] ? [29, 83, 132, 100] : [91, 158, 197, 32], getWidth: .7});
    $('viewport').dataset.selectionReady = 'true'; $('viewport').dataset.incidentEdges = data.source.length / 2;
    $('progress').textContent = 'Selection ready'; revision++; refreshNodes(); inspect(); render();
  }
  if (data.type === 'all') {
    fullLayer = lineLayer('all-edges', data.source, data.target, {getColor: [93, 111, 122, 9], getWidth: .35, visible: edgeMode === 'all'});
    $('viewport').dataset.allEdges = data.source.length / 2;
    $('progress').textContent = 'All edges ready'; render();
  }
}

async function start() {
  [nodes, modules] = await Promise.all(['nodes', 'modules'].map(async name => {
    const response = await fetch(new URL(`data/${name}.json`, import.meta.url));
    if (!response.ok) throw new Error(`${name}.json: HTTP ${response.status}`);
    return response.json();
  }));
  moduleMap = new Map(modules.modules.map(m => [m.module_id, m]));
  members = new Map(modules.modules.map(m => [m.module_id, []]));
  positions = new Float32Array(nodes.n * 2);
  for (let i = 0; i < nodes.n; i++) {
    positions.set([nodes.x[i], nodes.y[i]], i * 2); allIds.push(i);
    members.get(nodes.module[i])?.push(i);
    if (nodes.module[i] === 44) vcxIds.push(i);
    searchRows.push([nodes.vOTU[i], nodes.label[i], nodes.name[i], nodes.family_levels[nodes.family[i]] || ''].map(v => v.toLowerCase()));
  }
  if ([...$('color').options].some(o => o.value === params.get('color'))) color = params.get('color');
  if (['selection', 'backbone', 'none'].includes(params.get('edges'))) edgeMode = params.get('edges');
  $('color').value = color; $('edges').value = edgeMode;
  fit(allIds);
  if (['x', 'y', 'z'].every(k => params.has(k) && Number.isFinite(Number(params.get(k))))) {
    viewState = {...viewState, target: [Number(params.get('x')), Number(params.get('y')), 0], zoom: Math.max(-8, Math.min(6, Number(params.get('z'))))};
  }
  createLayers(); refreshNodes();
  const firstDraw = new Promise(resolve => { firstFrame = resolve; });
  viewer = new Deck({canvas: 'network-canvas', views: new OrthographicView({flipY: false}), viewState,
    controller: {dragRotate: false}, useDevicePixels: Math.min(devicePixelRatio, 2),
    parameters: {clearColor: [1, 1, 1, 1]},
    onViewStateChange: ({viewState: next}) => { viewState = next; viewer.setProps({viewState}); scheduleURL(); },
    onError: error => { $('progress').textContent = `Rendering failed: ${error.message}`; console.error(error); },
    onBeforeRender: () => { frameStarted = performance.now(); },
    onAfterRender: () => {
      if (firstFrame) { firstFrame(); firstFrame = null; $('viewport').dataset.nodes = nodes.n; }
      if (exportRequested) { exportRequested = false; downloadPNG(); }
      if (edgeMode === 'all' && fullLayer && fullFrames < 5) {
        fullTime += performance.now() - frameStarted; fullFrames++;
        if (fullTime / fullFrames > 80) $('notice').textContent = 'All edges render slowly on this device. Choose Selection or Backbone for smoother navigation.';
      }
    }});
  render(); inspect(); legend();
  $('counts').textContent = `${fmt(nodes.n)} nodes · ${fmt(nodes.edges)} edges`;
  for (const id of ['whole', 'color', 'edges', 'copy', 'png']) $(id).disabled = false;
  $('whole').onclick = whole;
  $('vcx').onclick = () => { $('edges').value = 'selection'; setEdgeMode(); select({kind: 'mod', id: 44}, true); };
  $('color').onchange = () => { color = $('color').value; refreshNodes(); legend(); render(); scheduleURL(); };
  $('edges').onchange = setEdgeMode;
  $('search').oninput = search;
  $('search-form').onsubmit = event => { event.preventDefault(); if (matches[0]) select(matches[0].value, true); };
  $('search').onkeydown = event => {
    if (event.key === 'Escape') { $('suggestions').hidden = true; $('search').setAttribute('aria-expanded', 'false'); }
    if (event.key === 'ArrowDown') { event.preventDefault(); $('suggestions').querySelector('button')?.focus(); }
  };
  $('copy').onclick = async () => {
    saveURL();
    try { await navigator.clipboard.writeText(location.href); $('notice').textContent = 'Link copied'; }
    catch (error) { $('notice').textContent = `Copy failed: ${error.message}. Copy the address bar URL.`; }
  };
  $('png').onclick = () => { exportRequested = true; viewer.redraw('export'); };
  $('network-canvas').onpointerleave = () => { $('tooltip').hidden = true; };
  $('toggle').onclick = () => { $('inspector').hidden = !$('inspector').hidden; $('toggle').setAttribute('aria-expanded', String(!$('inspector').hidden)); };
  await firstDraw;
  $('progress').textContent = 'Nodes ready · Loading adjacency data…';
  worker = new Worker(new URL('worker.js', import.meta.url), {type: 'module'});
  worker.onmessage = workerMessage;
  worker.onerror = event => { $('progress').textContent = `Network worker failed: ${event.message}`; };
  worker.postMessage({type: 'init', positions, moduleOf: new Int32Array(nodes.module)});
}

start().catch(error => { $('progress').textContent = `Network could not load: ${error.message}`; console.error(error); });
