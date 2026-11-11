/* Optional offline browser audit; npm dependency: playwright. */
const crypto = require('node:crypto');
const fs = require('node:fs'), path = require('node:path'), assert = require('node:assert/strict');
const {pathToFileURL} = require('node:url');
const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const root = path.resolve(process.argv[2] || 'tmp/multi-annotation-audit/results');
const checks = [], failures = [], images = path.join(root, 'browser');
fs.mkdirSync(images, {recursive: true});
function files(dir, suffix) {
  return fs.readdirSync(dir, {withFileTypes: true}).flatMap(e => e.isDirectory() && e.name !== 'browser' ? files(path.join(dir, e.name), suffix) : e.isFile() && e.name.endsWith(suffix) ? [path.join(dir, e.name)] : []);
}
async function unclipped(page, selector) {
  return page.locator(selector).evaluateAll(svgs => svgs.flatMap(svg => {
    const view = svg.viewBox.baseVal, inverse = svg.getScreenCTM()?.inverse(), bad = [];
    if (!inverse || !svg.getBoundingClientRect().height) return [];
    for (const text of svg.querySelectorAll('text')) {
      const b = text.getBBox(), m = inverse.multiply(text.getScreenCTM());
      for (const [x,y] of [[b.x,b.y],[b.x+b.width,b.y+b.height]]) {
        const p = new DOMPoint(x,y).matrixTransform(m);
        if (p.x < view.x-1 || p.y < view.y-1 || p.x > view.x+view.width+1 || p.y > view.y+view.height+1) bad.push(text.textContent);
      }
    }
    return bad;
  }));
}
(async () => {
  const browser = await chromium.launch({headless: true, executablePath: process.env.ANNOSTAT_BROWSER || undefined});
  try {
    for (const file of files(root, '.html')) for (const width of [1440,430]) {
      const page = await browser.newPage({viewport: {width, height: 1100}, locale: 'en-IN', acceptDownloads: true});
      const errors = [];
      page.on('pageerror', e => errors.push(e.message));
      page.on('dialog', async d => {errors.push('unexpected dialog'); await d.dismiss();});
      await page.route(/^https?:/, route => {errors.push('remote request '+route.request().url()); return route.abort();});
      const name = path.relative(root, file);
      try {
        await page.goto(pathToFileURL(file).href, {timeout: 60000});
        await page.evaluate(() => document.fonts.ready);
        const dataFile = path.join(path.dirname(file), 'comparison.json');
        const data = fs.existsSync(dataFile) ? JSON.parse(fs.readFileSync(dataFile)) : null;
        const summaryFile = path.join(path.dirname(file),'summary.json');
        if(fs.existsSync(summaryFile)) {
          const s = JSON.parse(fs.readFileSync(summaryFile));
          const cards = await page.locator('.metric').evaluateAll(nodes => Object.fromEntries(nodes.map(n => [n.querySelector('.metric-label').textContent,n.querySelector('.metric-value').textContent])));
          assert.equal(cards['Genome size'],s.genome_length.toLocaleString('en-US')+' bp');
          assert.equal(cards['Features'],s.total_features.toLocaleString('en-US'));
          assert.equal(cards['Genome GC'],s.quality_control.genome_gc_percent.toFixed(2)+'%');
          assert.equal(cards['Coding density'],s.quality_control.coding_density_percent.toFixed(2)+'%');
          assert.equal(cards['COG coverage'],s.cog_data_available?(100*s.cds_with_cog_count/s.cds_count).toFixed(2)+'%':'Not available');
        }
        const cohortFile = path.join(path.dirname(file),'cohort.json');
        if(fs.existsSync(cohortFile)) {
          const c = JSON.parse(fs.readFileSync(cohortFile));
          const rows = await page.locator('tbody tr').evaluateAll(nodes => nodes.map(n => [...n.querySelectorAll('td')].map(t => t.textContent)));
          assert.equal(rows.length,c.sample_count);
          c.samples.forEach((s,i) => {
            assert.equal(rows[i][0],s.sample);assert.equal(rows[i][2],s.genome_length_bp.toLocaleString('en-US'));
            assert.equal(rows[i][3],s.cds_count.toLocaleString('en-US'));assert.equal(rows[i][4],s.genome_gc_percent.toFixed(2));
          });
        }
        const gene = data?.gene_level_comparison;
        if (gene) {
          const info = () => page.evaluate(() => new DataTable('#gene-table').page.info());
          assert.equal((await info()).recordsTotal, gene.relationships.length);
          assert.equal(await page.locator('#gene-view').inputValue(), 'review');
          const review = gene.relationships.filter(r => r.status !== 'matched' || r.relationship !== 'exact' || r.protein_consequence !== 'identical' || (r.changed_attributes || '').split(';').some(k => k && k !== 'Dbxref'));
          assert.equal((await info()).recordsDisplay, review.length, 'default review row count');
          for (const [key,count] of Object.entries(gene.consequence_counts)) {
            assert.equal(await page.locator(`tr[data-consequence="${key}"] .gene-result-count strong`).innerText(), count.toLocaleString('en-US'));
          }
          await page.locator('#gene-view').selectOption('all');
          assert.equal((await info()).recordsDisplay, gene.relationships.length);
          for (const category of [...new Set(gene.relationships.map(r => r.protein_consequence))]) {
            await page.locator('#gene-consequence').selectOption(category);
            assert.equal((await info()).recordsDisplay, gene.relationships.filter(r => r.protein_consequence === category).length, 'consequence filter');
          }
          await page.locator('#gene-consequence').selectOption('');
          // Inspect one representative of every available relationship/evidence type.
          const selected = new Map();
          gene.relationships.forEach(r => {const k = r.status+'/'+r.relationship+'/'+r.protein_consequence; if (!selected.has(k)) selected.set(k,r);});
          for (const row of selected.values()) {
            await page.locator('input[type=search]').fill(row.id_a || row.id_b);
            const index = await page.evaluate(({a,b}) => {
              const table = new DataTable('#gene-table');
              return table.rows({search:'applied'}).indexes().toArray().find(i => table.row(i).data().id_a === a && table.row(i).data().id_b === b);
            }, {a:row.id_a,b:row.id_b});
            assert.notEqual(index,undefined);
            await page.evaluate(i => {const t = new DataTable('#gene-table');t.page(Math.floor(t.rows({search:'applied'}).indexes().toArray().indexOf(i)/t.page.len())).draw(false);},index);
            await page.locator(`#gene-table button[data-gene-index="${index}"]`).click();
            const detail = await page.locator('#gene-detail').innerText();
            for (const side of ['a','b']) if (row['id_'+side]) assert(detail.includes(row['id_'+side]));
            assert(detail.includes('1-based inclusive'), 'coordinate system not disclosed');
            if (row.status !== 'matched') assert.equal(await page.locator('.gene-detail-stats').count(),0);
            assert.deepEqual(await unclipped(page,'.gene-locus svg'),[], 'locus text clipping');
            const ticks = (await page.locator('.gene-locus svg text').allTextContents()).filter(t => /^[\d,]+$/.test(t));
            assert(ticks.length >= 2);
            for(const tick of ticks) assert.equal(tick,Number(tick.replaceAll(',','')).toLocaleString('en-US'),'coordinate grouping depends on browser locale');
            const context = page.locator('.gene-context-toggle input');
            if (!(await context.isDisabled())) {
              await context.check(); assert.deepEqual(await unclipped(page,'.gene-locus svg'),[],'context text clipping'); await context.uncheck();
            }
          }
          await page.locator('input[type=search]').fill('');
          await page.locator('#gene-view').selectOption('review');
          assert.equal((await info()).recordsDisplay,review.length);
          // Save one representative per assembly and viewport.
          if (name.includes('Bakta-vs-RefSeq') && name.includes('forward')) {
            const id = name.split(path.sep)[1];
            await page.locator('#gene-detail').screenshot({path:path.join(images,id+'-'+width+'.png'),style:'.gene-nav{visibility:hidden}'});
            const download = page.waitForEvent('download');
            await page.getByRole('button',{name:'Download locus SVG',exact:true}).click();
            await (await download).saveAs(path.join(images,id+'-'+width+'.svg'));
          }
        }
        assert.equal(await page.evaluate(() => document.documentElement.scrollWidth > innerWidth+1),false,'page horizontal overflow');
        assert.deepEqual(errors,[]);
        checks.push({file:name,width,passed:true});
      } catch(e) {failures.push({file:name,width,error:String(e),browserErrors:errors});}
      fs.writeFileSync(path.join(images,'progress.json'),JSON.stringify({checks:checks.length,failures:failures.length,lastFile:name,width}));
      await page.close();
    }
    // Inspect every standalone scientific SVG, including figures hidden in HTML details.
    for (const file of files(root,'.svg').filter(f => !f.includes(path.sep+'browser'+path.sep))) {
      const page = await browser.newPage({viewport:{width:1600,height:1300}});
      try {
        await page.goto(pathToFileURL(file).href);
        await page.evaluate(() => document.fonts.ready);
        assert.deepEqual(await unclipped(page,'svg'),[],'figure text outside exported canvas');
        checks.push({file:path.relative(root,file),passed:true});
      } catch(e) {failures.push({file:path.relative(root,file),error:String(e)});}
      await page.close();
    }
  } finally {await browser.close();}
  fs.writeFileSync(path.join(images,'verification.json'),JSON.stringify({node:process.version,browser:browser.version(),locale:'en-IN',assetHashes:Object.fromEntries(['gene_report.js','gene_report.css'].map(name=>[name,crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,'..','annostat','assets',name))).digest('hex')])),checks,failures},null,2));
  console.log(JSON.stringify({checks:checks.length,failures},null,2));
  if(failures.length) process.exitCode=1;
})().catch(e => {console.error(e);process.exitCode=1;});
