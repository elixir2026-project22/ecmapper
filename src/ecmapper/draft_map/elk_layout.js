// reads an ELK JSON graph from stdin, writes the laid out graph to stdout
const ELK = require('elkjs/lib/elk.bundled.js');
let input = '';
process.stdin.on('data', d => (input += d));
process.stdin.on('end', async () => {
  const elk = new ELK();
  const graph = await elk.layout(JSON.parse(input));
  process.stdout.write(JSON.stringify(graph));
});
