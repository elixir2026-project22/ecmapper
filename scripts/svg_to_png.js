// Render an SVG (e.g. a draft map) as PNG, optionally a part of it; needs the npm package playwright.
// usage: node svg_to_png.js in.svg out.png [width] [viewBox "x y w h"]
const { chromium } = require('playwright');
const fs = require('fs');
(async () => {
  const [inp, out, widthArg, vb] = process.argv.slice(2);
  let svg = fs.readFileSync(inp, 'utf8');
  const width = parseInt(widthArg || '2400');
  // size from the viewBox (given or the SVG's own)
  svg = svg.replace(/<svg([^>]*?)\swidth="[^"]*"/, '<svg$1').replace(/<svg([^>]*?)\sheight="[^"]*"/, '<svg$1');
  if (vb) svg = svg.replace(/<svg([^>]*?)\sviewBox="[^"]*"/, '<svg$1').replace('<svg', `<svg viewBox="${vb}"`);
  const m = svg.match(/viewBox="([^"]+)"/);
  const [, , w, h] = m[1].split(/[\s,]+/).map(Number);
  const height = Math.round(width * h / w);
  svg = svg.replace('<svg', `<svg width="${width}" height="${height}"`);
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width, height } });
  await page.setContent(`<html><body style="margin:0;background:white">${svg}</body></html>`);
  await page.screenshot({ path: out, fullPage: true });
  await browser.close();
})();
