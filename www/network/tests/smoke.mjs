// From www/: python3 -m http.server 8765 --bind 127.0.0.1
// Run: node network/tests/smoke.mjs (Playwright installed; optional CHROMIUM_PATH / PLAYWRIGHT_MODULE).
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {fileURLToPath} from 'node:url';
const {chromium} = await import(process.env.PLAYWRIGHT_MODULE || 'playwright');
const base = process.env.NETWORK_BASE_URL || 'http://127.0.0.1:8765/';
const browser = await chromium.launch({headless: true, executablePath: process.env.CHROMIUM_PATH,
  args: ['--no-sandbox', '--use-angle=swiftshader', '--enable-unsafe-swiftshader']});
const context = await browser.newContext({viewport: {width: 1440, height: 1000}, deviceScaleFactor: 1});
const page = await context.newPage();
const errors = [], requests = [];
page.on('console', message => { if (message.type() === 'error') errors.push(message.text()); });
page.on('pageerror', error => errors.push(error.message));
page.on('request', request => requests.push(request.url()));
page.setDefaultTimeout(120000);
const data = JSON.parse(await readFile(new URL('../data/nodes.json', import.meta.url)));
const offsetsBytes = await readFile(new URL('../data/csr-offsets.u32', import.meta.url));
const targetsBytes = await readFile(new URL('../data/csr-targets.u16', import.meta.url));
const offsets = new Uint32Array(offsetsBytes.buffer, offsetsBytes.byteOffset, offsetsBytes.byteLength / 4);
const targets = new Uint16Array(targetsBytes.buffer, targetsBytes.byteOffset, targetsBytes.byteLength / 2);

async function waitSelection(key, count) {
  await page.waitForFunction(({key, count}) => {
    const state = document.getElementById('viewport').dataset;
    return state.selection === key && Number(state.selectedCount) === count && state.selectionReady === 'true';
  }, {key, count});
}

async function explorerSmoke() {
  // Delay the CSR to verify that the actual canvas renders before adjacency is available.
  await page.route('**/csr-*.u*', async route => { await new Promise(resolve => setTimeout(resolve, 1000)); await route.continue(); });
  await page.goto(`${base}network.html`);
  await page.waitForFunction(() => document.getElementById('viewport').dataset.nodes === '17838');
  assert.equal(await page.locator('#counts').textContent(), '17,838 nodes · 2,301,056 edges');
  assert.equal(await page.locator('#vcx').isDisabled(), true);
  await page.waitForFunction(() => document.getElementById('viewport').dataset.csrReady === 'true');
  await page.screenshot({path: fileURLToPath(new URL('whole.png', import.meta.url))});
  await page.locator('#vcx').click();
  await waitSelection('mod:44', 266);
  let internal = 0, external = 0;
  const shared = new Map(), neighbors = new Set();
  for (let a = 0; a < data.n; a++) if (data.module[a] === 44) {
    for (let p = offsets[a]; p < offsets[a + 1]; p++) {
      const b = targets[p];
      if (data.module[b] === 44) internal++;
      else { external++; neighbors.add(b); shared.set(data.module[b], (shared.get(data.module[b]) || 0) + 1); }
    }
  }
  assert.equal(internal / 2, 14668);
  assert.equal(Number(await page.locator('#viewport').getAttribute('data-incident-edges')), internal / 2 + external);
  assert.match(await page.locator('#details').textContent(), new RegExp(`${neighbors.size.toLocaleString('en-US')} neighboring nodes`));
  assert.match(await page.locator('#details').textContent(), /Table S9 counts VCX edges .*25,590 to catalog vOTUs in other modules and 156 to RefSeq genomes/);
  const [topModule, topCount] = [...shared].sort((a, b) => b[1] - a[1])[0];
  assert.match(await page.locator('.neighbors li').first().textContent(), new RegExp(`Module ${topModule}.*${topCount.toLocaleString('en-US')} edges`));
  await page.screenshot({path: fileURLToPath(new URL('vcx.png', import.meta.url))});
  console.log(`PASS VCX: 266 selected; ${internal / 2} internal + ${external} external edges; top neighbor ${topModule} (${topCount} edges)`);
  await page.locator('#search').fill('vOTU000463'); await page.locator('#search').press('Enter');
  const node = data.vOTU.indexOf('vOTU000463'); await waitSelection(`node:${node}`, 1);
  assert.match(await page.locator('#details').textContent(), /vOTU000463/);
  assert.equal(Number(await page.locator('#viewport').getAttribute('data-incident-edges')), offsets[node + 1] - offsets[node]);
  assert.match(await page.locator('#details').textContent(), /280 neighboring nodes/);
  console.log(`PASS search vOTU000463: ${offsets[node + 1] - offsets[node]} incident edges; inspector correct`);
  await page.goto(`${base}network.html?sel=mod:44&color=vcx&edges=selection`);
  await waitSelection('mod:44', 266);
  await page.waitForFunction(() => new URL(location.href).searchParams.has('z'));
  const saved = page.url(); await page.goto(saved); await waitSelection('mod:44', 266);
  assert.equal(new URL(page.url()).searchParams.get('z'), new URL(saved).searchParams.get('z'));
  console.log('PASS URL: module 44, color, edges and camera restored');
}

async function modesAndExport() {
  // Pick an isolated VCX member using the public URL camera and the original, upward-positive y.
  const isolated = data.module.reduce((best, m, i) => m === 44 && (best < 0 || data.y[i] > data.y[best]) ? i : best, -1);
  const point = await page.evaluate(({x, y}) => {
    const camera = new URL(location.href).searchParams, scale = 2 ** Number(camera.get('z'));
    const box = document.getElementById('network-canvas').getBoundingClientRect();
    return {x: box.x + box.width / 2 + (x - Number(camera.get('x'))) * scale,
      y: box.y + box.height / 2 - (y - Number(camera.get('y'))) * scale};
  }, {x: data.x[isolated], y: data.y[isolated]});
  await page.mouse.move(point.x, point.y);
  await page.locator('#tooltip').waitFor({state: 'visible'});
  assert.match(await page.locator('#tooltip').textContent(), new RegExp(data.label[isolated]));
  await page.mouse.click(point.x, point.y); await waitSelection(`node:${isolated}`, 1);
  await page.locator('#details button').filter({hasText: 'Module 44'}).click(); await waitSelection('mod:44', 266);
  for (const color of ['set', 'module', 'family', 'host_phylum', 'vcx']) {
    await page.locator('#color').selectOption(color);
    assert.ok(await page.locator('#legend-items').textContent());
  }
  await page.locator('#edges').selectOption('backbone');
  await page.screenshot(); // Flush the newly visible picking layer before clicking its marker.
  assert.match(await page.locator('#edge-legend').textContent(), /50 shared edges/);
  const modules = JSON.parse(await readFile(new URL('../data/modules.json', import.meta.url)));
  const m = modules.modules.find(m => m.module_id === 44);
  const anchor = await page.evaluate(({x, y}) => {
    const camera = new URL(location.href).searchParams, scale = 2 ** Number(camera.get('z'));
    const box = document.getElementById('network-canvas').getBoundingClientRect();
    return {x: box.x + box.width / 2 + (x - Number(camera.get('x'))) * scale,
      y: box.y + box.height / 2 - (y - Number(camera.get('y'))) * scale};
  }, m);
  await page.mouse.move(anchor.x, anchor.y);
  await page.waitForFunction(() => document.getElementById('tooltip').textContent.includes('Click to select module'));
  await page.mouse.click(anchor.x, anchor.y);
  await page.waitForFunction(() => document.getElementById('edges').value === 'selection');
  await waitSelection('mod:44', 266);
  assert.equal(await page.locator('#edges').inputValue(), 'selection');
  console.log('PASS canvas hover, node picking, inspector module link and backbone module picking');
  await page.locator('#whole').click();
  await page.locator('#edges').selectOption('all');
  await page.waitForFunction(() => document.getElementById('viewport').dataset.allEdges === '2301056');
  await page.waitForTimeout(1000);
  assert.match(await page.locator('#notice').textContent(), /slow/);
  await page.locator('#edges').selectOption('none');
  await page.locator('#edges').selectOption('selection');
  const downloadPromise = page.waitForEvent('download'); await page.locator('#png').click();
  const download = await downloadPromise;
  const png = await readFile(await download.path());
  assert.equal(png.subarray(1, 4).toString(), 'PNG'); assert.ok(png.length > 5000);
  const pixel = await page.evaluate(async base64 => {
    const image = new Image(); image.src = `data:image/png;base64,${base64}`; await image.decode();
    const canvas = document.createElement('canvas'); canvas.width = image.width; canvas.height = image.height;
    const ctx = canvas.getContext('2d'); ctx.drawImage(image, 0, 0);
    return [...ctx.getImageData(0, 0, 1, 1).data];
  }, png.toString('base64'));
  assert.deepEqual(pixel, [255, 255, 255, 255]);
  console.log(`PASS color modes, all 2,301,056 edges and PNG export (${png.length} bytes)`);
  await page.setViewportSize({width: 390, height: 844});
  await page.locator('#vcx').click(); await waitSelection('mod:44', 266);
  const panel = await page.locator('#inspector').boundingBox();
  assert.ok(panel.y > 400 && panel.width === 390);
  await page.screenshot({path: fileURLToPath(new URL('mobile.png', import.meta.url))});
  await page.locator('#toggle').click(); assert.equal(await page.locator('#inspector').isVisible(), false);
  assert.equal(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), true);
  console.log('PASS mobile bottom sheet, collapse and viewport width');
}

try {
  await explorerSmoke(); await modesAndExport();
  assert.deepEqual(errors, []);
  assert.ok(requests.every(url => url.startsWith(base)), 'All runtime requests must use local assets');
  console.log('PASS console: 0 errors; runtime requests: local assets only');
  console.log('Screenshots: network/tests/vcx.png, network/tests/whole.png');
} finally { await browser.close(); }
