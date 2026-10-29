let positions, moduleOf, offsets, targets;
const members = new Map();

async function initialize(data) {
  positions = data.positions;
  moduleOf = data.moduleOf;
  for (let i = 0; i < moduleOf.length; i++) {
    if (!members.has(moduleOf[i])) members.set(moduleOf[i], []);
    members.get(moduleOf[i]).push(i);
  }
  const buffers = await Promise.all(['csr-offsets.u32', 'csr-targets.u16'].map(async name => {
    const response = await fetch(new URL(`data/${name}`, import.meta.url));
    if (!response.ok) throw new Error(`${name}: HTTP ${response.status}`);
    return response.arrayBuffer();
  }));
  offsets = new Uint32Array(buffers[0]);
  targets = new Uint16Array(buffers[1]);
  postMessage({type: 'ready'});
}

function edgeBuffers(count) {
  return {source: new Float32Array(count * 2), target: new Float32Array(count * 2)};
}

function putEdge(buffer, k, a, b) {
  buffer.source[k * 2] = positions[a * 2];
  buffer.source[k * 2 + 1] = positions[a * 2 + 1];
  buffer.target[k * 2] = positions[b * 2];
  buffer.target[k * 2 + 1] = positions[b * 2 + 1];
}

function select({selection, request}) {
  const ids = selection.kind === 'node' ? [selection.id] : members.get(selection.id);
  const selected = new Set(ids), neighbors = new Set(), shared = new Map();
  let count = 0;
  // Two degree-bounded passes allocate exactly the buffers transferred to the renderer.
  for (const a of ids) for (let p = offsets[a]; p < offsets[a + 1]; p++) {
    const b = targets[p];
    if (!selected.has(b) || a < b) count++;
  }
  const buffer = edgeBuffers(count), internal = new Uint8Array(count);
  let k = 0, internalCount = 0;
  for (const a of ids) for (let p = offsets[a]; p < offsets[a + 1]; p++) {
    const b = targets[p];
    if (selected.has(b) && a >= b) continue;
    putEdge(buffer, k, a, b);
    const same = moduleOf[a] >= 0 && moduleOf[a] === moduleOf[b];
    internal[k++] = Number(same);
    internalCount += Number(same);
    if (!selected.has(b)) neighbors.add(b);
    if (!same) shared.set(moduleOf[b], (shared.get(moduleOf[b]) || 0) + 1);
  }
  const result = {type: 'selection', request, ...buffer, internal, internalCount,
    externalCount: count - internalCount, neighbors: [...neighbors],
    shared: [...shared].sort((a, b) => b[1] - a[1] || a[0] - b[0])};
  postMessage(result, [buffer.source.buffer, buffer.target.buffer, internal.buffer]);
}

function allEdges() {
  const buffer = edgeBuffers(targets.length / 2);
  let k = 0;
  for (let a = 0; a < moduleOf.length; a++) for (let p = offsets[a]; p < offsets[a + 1]; p++) {
    const b = targets[p];
    if (a < b) putEdge(buffer, k++, a, b);
  }
  postMessage({type: 'all', ...buffer}, [buffer.source.buffer, buffer.target.buffer]);
}

self.onmessage = ({data}) => {
  if (data.type === 'init') initialize(data).catch(error => postMessage({type: 'error', message: error.message}));
  if (data.type === 'select') select(data);
  if (data.type === 'all') allEdges();
};
