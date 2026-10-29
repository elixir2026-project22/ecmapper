// Load an Escher map (with a model) in Escher in a headless browser, save a screenshot and report the
// reactions and metabolites of the map that are not in the model.
// usage: ESCHER_JS=.../escher.min.js node render_escher.js map.json model.json|- out.png [width height] ["x y w h" (map units)]
// escher.min.js: e.g. escher/static/escher.min.js of the escher Python package; needs the npm package playwright.
const { chromium } = require('playwright');
const fs = require('fs');
(async () => {
  const [mapFile, modelFile, out, w = '2400', h = '1500', region] = process.argv.slice(2);
  const escher = fs.readFileSync(process.env.ESCHER_JS, 'utf8');
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: parseInt(w), height: parseInt(h) } });
  const messages = [];
  page.on('console', m => messages.push(m.text()));
  page.on('pageerror', e => messages.push('PAGEERROR ' + e.message));
  await page.setContent(`<html><body style="margin:0"><div id="map" style="width:${w}px;height:${h}px"></div></body></html>`);
  await page.addScriptTag({ content: escher });
  const mapData = fs.readFileSync(mapFile, 'utf8');
  const modelData = modelFile !== '-' ? fs.readFileSync(modelFile, 'utf8') : 'null';
  await page.evaluate(([mapData, modelData]) => {
    window.builder = escher.Builder(JSON.parse(mapData), JSON.parse(modelData), null, escher.libs.d3_select('#map'),
      { menu: 'zoom', fill_screen: false, never_ask_before_quit: true, enable_tooltips: false,
        identifiers_on_map: 'bigg_id' });
  }, [mapData, modelData]);
  await page.waitForFunction(() => window.builder && window.builder.map && window.builder.map.reactions, null, { timeout: 300000 });
  const info = await page.evaluate((region) => {
    const b = window.builder, m = b.map;
    if (region) {
      const [x, y, rw, rh] = region.split(' ').map(Number);
      const size = m.get_size(); const s = Math.min(size.width / rw, size.height / rh);
      m.zoomContainer.goTo(s, { x: -x * s, y: -y * s });
    } else {
      m.zoom_extent_canvas();
    }
    const model = b.cobra_model;
    const rs = Object.values(m.reactions);
    return { reactions: rs.length, nodes: Object.keys(m.nodes).length,
             reactions_not_in_model: model ? rs.filter(r => !(r.bigg_id in model.reactions)).length : null,
             metabolites_not_in_model: model ? Object.values(m.nodes).filter(n => n.node_type === 'metabolite' && !(n.bigg_id in model.metabolites)).length : null };
  }, region);
  await page.waitForTimeout(2000);
  await page.screenshot({ path: out });
  console.log(JSON.stringify(info));
  console.log(messages.slice(0, 10).join('\n'));
  await browser.close();
})();
