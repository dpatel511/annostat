/* Annostat gene-level report. Source attributes are rendered as text, never HTML. */
(() => {
    'use strict';
    const gene = JSON.parse(document.getElementById('gene-data').textContent);
    const rows = gene.relationships.map((row, index) => ({...row, index}));
    const labels = Object.keys(gene.cds_counts);
    const loci = JSON.parse(document.getElementById('gene-loci').textContent);
    const locusByRow = new Map();
    loci.forEach(group => group.rows.forEach(index => locusByRow.set(index, group)));
    const modelKey = (row, side) => JSON.stringify([side, row.seqid, row.strand, row[`id_${side}`]]);
    const inventory = new Map();
    rows.forEach(row => ['a', 'b'].forEach(side => {
        if (row[`id_${side}`] === null) return;
        const parts = row[`segments_${side}`].split(';').map(block => block.split('-').map(Number));
        inventory.set(modelKey(row, side), {key: modelKey(row, side), side, seqid: row.seqid, strand: row.strand,
            id: row[`id_${side}`], gene: row[`gene_${side}`], parts,
            low: Math.min(...parts.flat()), high: Math.max(...parts.flat()),
            wraps: parts.some((part, i) => i && (row.strand === '+' ? part[0] < parts[i-1][0] : part[0] > parts[i-1][0]))});
    }));
    const needsReview = row => row.status !== 'matched' || row.relationship !== 'exact' ||
        row.protein_consequence !== 'identical' || (row.changed_attributes || '').split(';').some(key => key && key !== 'Dbxref');
    const statusNames = {matched: 'Accepted pair', candidate: 'Candidate pair', ambiguous: 'Ambiguous correspondence', unmatched: 'Unmatched call'};
    const names = {identical: 'Identical protein', terminal_only: 'Terminal-only difference',
        initiator_recoding: 'Initiation-context difference', translation_difference: 'Other translated-sequence difference',
        limited: 'Limited protein evidence', unresolved: 'Unresolved correspondence'};
    const relationshipNames = {exact: 'Exact CDS structure', start_shift: 'Alternative start boundary', stop_shift: 'Alternative stop boundary',
        frame_disagreement: 'Frame disagreement', overlap: 'Other overlap', unmatched: 'Unmatched call'};
    const startContexts = {initiator: 'Initiator rule at annotated complete start',
        incomplete_start: 'Incomplete 5′ boundary; elongation translation',
        translation_exception: 'Explicit first-codon translation exception',
        unrecognized_codon: 'First codon is not a recognized initiator', unavailable: 'No complete first codon'};
    const shown = value => value === null || value === undefined ? 'Not available' : value === '' ? 'Not provided' : String(value);
    const signed = value => value === null || value === undefined ? 'Not available' : `${value >= 0 ? '+' : ''}${value}`;
    const measured = value => typeof value === 'number' && !Number.isInteger(value) ? value.toFixed(2) : shown(value);
    const el = (tag, value, cls) => {
        const node = document.createElement(tag);
        if (value !== undefined) node.textContent = value;
        if (cls) node.className = cls;
        return node;
    };
    const text = DataTable.render.text();
    const column = (title, data) => ({title, data, defaultContent: '—', render: text});
    const table = new DataTable('#gene-table', {
        data: rows, pageLength: 25, lengthMenu: [10, 25, 50, 100], order: [], deferRender: true,
        columns: [
            {title: 'Evidence', data: 'index', orderable: false, searchable: false,
                render: (index, type) => type === 'display' ? `<button type="button" data-gene-index="${index}">Inspect</button>` : index},
            ...['a', 'b'].map((side, i) => column(`${side.toUpperCase()}: ${text.display(labels[i])}`, row =>
                row[`id_${side}`] === null ? 'No counterpart' : `${row[`gene_${side}`] || 'No gene symbol'}\n${row[`id_${side}`]}\n${shown(row[`protein_length_${side}`])} aa`)),
            column('Correspondence', row => `${statusNames[row.status]}\n${relationshipNames[row.relationship]}`),
            column('Protein evidence', row => names[row.protein_consequence]),
            {title: 'Δ aa B − A', data: 'protein_length_delta_aa', defaultContent: '—', type: 'num',
                render: (value, type) => type === 'display' ? value === null ? '—' : signed(value) : value},
            {...column('Source context', row => ['seqid', 'product_a', 'product_b', 'gene_b', 'changed_attributes', 'segments_a', 'segments_b',
                'partial_a', 'partial_b', 'pseudo_a', 'pseudo_b', 'pseudogene_a', 'pseudogene_b', 'exception_a', 'exception_b',
                'dbxref_a', 'dbxref_b', 'status', 'relationship'].map(key => row[key] || '').join(' ')), visible: false}
        ], language: {emptyTable: 'No CDS relationships in these annotations.', search: 'Search current view:', zeroRecords: 'No relationships match this view. Try All relationships or clear the filters.'}
    });
    function updateViewSummary() {
        const info = table.page.info();
        document.getElementById('gene-view-summary').textContent = `Showing ${info.recordsDisplay.toLocaleString('en-US')} of ${info.recordsTotal.toLocaleString('en-US')} relationship rows. Overview counts and exports retain all rows. The review view includes structural, protein and source-attribute changes except database references alone.`;
    }
    table.on('draw', updateViewSummary);
    table.search.fixed('gene-view', (search, row) => needsReview(row)).draw();
    document.getElementById('gene-view').addEventListener('change', event => {
        table.search.fixed('gene-view', event.target.value === 'review' ? (search, row) => needsReview(row) : null).draw();
    });
    for (const [id, key] of [['gene-consequence', 'protein_consequence'], ['gene-structure', 'relationship']]) {
        document.getElementById(id).addEventListener('change', event => {
            const value = event.target.value;
            table.search.fixed(id, value ? (search, row) => row[key] === value : null).draw();
        });
    }
    function svgNode(tag, attrs, value) {
        const node = document.createElementNS('http://www.w3.org/2000/svg', tag);
        Object.entries(attrs).forEach(([key, val]) => node.setAttribute(key, val));
        if (value !== undefined) node.textContent = value;
        return node;
    }
    function locusModels(row) {
        const group = locusByRow.get(row.index);
        const keys = new Set();
        (group ? group.rows : [row.index]).forEach(index => ['a', 'b'].forEach(side => {
            if (rows[index][`id_${side}`] !== null) keys.add(modelKey(rows[index], side));
        }));
        return [...keys].map(key => inventory.get(key)).sort((a, b) => a.side.localeCompare(b.side) || a.low-b.low || a.id.localeCompare(b.id));
    }
    function locus(row, withContext = false) {
        const members = locusModels(row), selected = new Set(['a', 'b'].map(side => modelKey(row, side)));
        const membership = new Set(members.map(model => model.key));
        let models = [...members];
        const lowMember = Math.min(...members.map(model => model.low)), highMember = Math.max(...members.map(model => model.high));
        if (withContext && !members.some(model => model.wraps)) {
            for (const side of ['a', 'b']) {
                const others = [...inventory.values()].filter(model => model.side === side && model.seqid === row.seqid && !membership.has(model.key));
                models.push(...others.filter(model => model.high >= lowMember && model.low <= highMember));
                const left = others.filter(model => model.high < lowMember).sort((a, b) => b.high-a.high)[0];
                const right = others.filter(model => model.low > highMember).sort((a, b) => a.low-b.low)[0];
                if (left) models.push(left);
                if (right) models.push(right);
            }
        }
        const low = Math.min(...models.map(model => model.low)), high = Math.max(...models.map(model => model.high));
        const height = 150 + 51 * (models.length + ['a', 'b'].filter(side => !models.some(model => model.side === side)).length);
        const scale = pos => 75 + 750 * (pos-low)/Math.max(1, high-low+1);
        const svg = svgNode('svg', {viewBox: `0 0 900 ${height}`, role: 'img', 'aria-label': 'CDS locus comparison, 1-based inclusive'});
        svg.append(svgNode('title', {}, `CDS locus: ${labels[0]} / ${labels[1]}`));
        svg.append(svgNode('desc', {}, `${row.seqid}; genomic coordinates 1-based inclusive; one arrow lane per CDS. Source colours identify annotations, not confidence. Grey arrows are neighbouring context, not additional accepted matches.`));
        svg.append(svgNode('rect', {width: 900, height, fill: '#fff'}));
        svg.append(svgNode('text', {x: 20, y: 23, fill: '#17231e', 'font-size': 14, 'font-family': 'Arial'}, `${row.seqid} | genomic coordinates (1-based inclusive)`));
        let y = 50;
        for (const [side, i] of [['a', 0], ['b', 1]]) {
            const color = i ? '#2563eb' : '#087f5b';
            svg.append(svgNode('text', {x: 20, y, fill: color, 'font-size': 14, 'font-family': 'Arial', 'font-weight': 'bold'}, `${side.toUpperCase()}: ${labels[i]}`));
            y += 22;
            const source = models.filter(model => model.side === side).sort((a, b) => a.low-b.low || a.high-b.high || a.id.localeCompare(b.id));
            if (!source.length) {
                svg.append(svgNode('text', {x: 75, y: y+8, fill: '#64746c', 'font-size': 12, 'font-family': 'Arial'}, 'No CDS call in this group'));
                y += 51;
            }
            for (const model of source) {
                const isMember = membership.has(model.key), isSelected = selected.has(model.key);
                const fill = isMember ? color : '#9ca9b3';
                const caption = `${isSelected ? 'Selected link · ' : ''}${model.id}${model.gene ? ' | '+model.gene : ''} | ${model.parts.map(part => part.join('–')).join('; ')} (${model.strand})${isMember ? '' : ' · context only'}`;
                svg.append(svgNode('text', {x: 75, y, fill: '#253b4a', 'font-size': 12, 'font-family': 'Arial'}, caption.length > 103 ? caption.slice(0,100)+'…' : caption));
                model.parts.forEach(([start, end]) => {
                    const x = scale(start), width = Math.max(1, scale(end+1)-x), tip = Math.min(8, width/2), top = y+7;
                    const shape = model.strand === '+' ? `${x},${top} ${x+width-tip},${top} ${x+width},${top+8} ${x+width-tip},${top+16} ${x},${top+16}`
                        : `${x+tip},${top} ${x+width},${top} ${x+width},${top+16} ${x+tip},${top+16} ${x},${top+8}`;
                    const arrow = svgNode('polygon', {points: shape, fill, stroke: isSelected ? '#17231e' : fill, 'stroke-width': isSelected ? 1.5 : 0.5});
                    arrow.append(svgNode('title', {}, `${caption}; block ${start}–${end} (1-based inclusive)`));
                    svg.append(arrow);
                });
                y += 51;
            }
            y += 10;
        }
        const axisY = height-28;
        svg.append(svgNode('line', {x1:75, x2:825, y1:axisY-14, y2:axisY-14, stroke:'#a8b6c1'}));
        for (const [position, anchor] of [[low, 'start'], [Math.floor((low+high)/2), 'middle'], [high, 'end']])
            svg.append(svgNode('text', {x:scale(position), y:axisY+3, 'text-anchor':anchor, fill:'#566777', 'font-size':12, 'font-family':'Arial'}, position.toLocaleString('en-US')));
        return svg;
    }
    function groupPanel(row) {
        const group = locusByRow.get(row.index);
        if (!group || group.kind !== 'ambiguous') return null;
        const panel = el('section', undefined, 'gene-locus-group');
        panel.append(el('h4', 'Ambiguous overlap group'));
        panel.append(el('p', `${labels[0]}: ${group.cds_counts.a} distinct CDS models; ${labels[1]}: ${group.cds_counts.b} distinct CDS models; ${group.rows.length} relationship rows.`));
        panel.append(el('p', group.reason));
        panel.append(el('p', 'These links are unresolved correspondence, not confirmed split/merge events. No pairwise protein consequence is assigned.', 'note'));
        const wrapper = el('div', undefined, 'table-wrap'), table = el('table', undefined, 'gene-group-members');
        const head = el('tr');
        for (const text of ['Source', 'CDS ID / gene', 'Segments (1-based inclusive)', 'Strand']) head.append(el('th', text));
        const thead = el('thead'); thead.append(head); table.append(thead);
        const body = el('tbody');
        for (const model of locusModels(row)) {
            const tr = el('tr');
            for (const value of [labels[model.side === 'a' ? 0 : 1], `${model.id}${model.gene ? ' / '+model.gene : ''}`, model.parts.map(part => part.join('–')).join('; '), model.strand]) tr.append(el('td', value));
            body.append(tr);
        }
        table.append(body); wrapper.append(table); panel.append(wrapper);
        const links = el('div', undefined, 'gene-group-links');
        group.rows.forEach(index => {
            const candidate = rows[index], button = el('button', `${candidate.id_a} ↔ ${candidate.id_b}`);
            button.type = 'button'; button.dataset.geneIndex = index;
            button.setAttribute('aria-pressed', String(index === row.index)); links.append(button);
        });
        panel.append(el('p', 'Inspect a candidate link; this does not accept it as a match.', 'note'), links);
        return panel;
    }
    function alignment(row, full = false) {
        const a = row.alignment_a, b = row.alignment_b;
        const output = el('pre');
        let first = 0;
        while (first < a.length && a[first] === b[first]) first++;
        let start = full ? 0 : Math.max(0, first - 16);
        const end = full ? a.length : Math.min(a.length, start + 64);
        let posA = a.slice(0, start).replace(/-/g, '').length;
        let posB = b.slice(0, start).replace(/-/g, '').length;
        for (; start < end; start += 64) {
            const endBlock = Math.min(end, start + 64);
            const left = a.slice(start, endBlock), right = b.slice(start, endBlock);
            const countA = left.replace(/-/g, '').length, countB = right.replace(/-/g, '').length;
            for (const [side, seq, other, pos, count] of [['A', left, right, posA, countA], ['B', right, left, posB, countB]]) {
                output.append(document.createTextNode(`${side} ${String(count ? pos + 1 : '—').padStart(6)} `));
                for (let i = 0; i < seq.length; i++) {
                    const changed = seq[i] !== other[i];
                    output.append(el(changed ? 'mark' : 'span', seq[i], changed ? seq[i] === '-' || other[i] === '-' ? 'gene-gap' : 'gene-mismatch' : undefined));
                }
                output.append(document.createTextNode(` ${count ? pos + count : '—'}\n`));
            }
            output.append(document.createTextNode('\n'));
            posA += countA; posB += countB;
        }
        return output;
    }
    function showDetail(row) {
        const detail = document.getElementById('gene-detail');
        const header = el('div', undefined, 'gene-detail-header');
        header.append(el('span', 'Selected CDS relationship', 'gene-eyebrow'));
        header.append(el('h3', `${row.gene_a || row.id_a || 'No counterpart'} → ${row.gene_b || row.id_b || 'No counterpart'}`));
        header.append(el('p', `${labels[0]}: ${row.id_a || 'No counterpart'} → ${labels[1]}: ${row.id_b || 'No counterpart'}`));
        header.append(el('p', `${statusNames[row.status]} · ${relationshipNames[row.relationship]} · ${names[row.protein_consequence]}`));
        detail.replaceChildren(header);
        const group = groupPanel(row);
        if (group) detail.append(group);
        const length = side => row[`id_${side}`] === null ? 'no counterpart' : `${shown(row[`protein_length_${side}`])} aa`;
        const lengths = `Predicted protein lengths: A ${length('a')}; B ${length('b')}. `;
        detail.append(el('p', lengths + row.consequence_note, 'gene-interpretation'));
        if (row.start_codon_a && row.start_codon_a === row.start_codon_b &&
            row.first_codon_positions_a !== row.first_codon_positions_b) {
            detail.append(el('p', `Both models begin with ${row.start_codon_a}, but at different genomic positions. These are separate codon occurrences. ` + (row.status === 'matched' ? 'Compare the aligned residues at shared genomic codons, not the two first-codon labels.' : 'Equal first-codon strings do not establish CDS correspondence; no protein alignment is assigned to this link.'), 'gene-interpretation'));
        }
        const stats = el('div', undefined, 'gene-detail-stats');
        for (const [title, value] of [['Protein length Δ (B−A, aa)', signed(row.protein_length_delta_aa)],
            ['N-terminal Δ (B−A, aa)', signed(row.n_terminal_delta_aa)], ['C-terminal Δ (B−A, aa)', signed(row.c_terminal_delta_aa)],
            ['Identity (including gaps, %)', measured(row.alignment_identity_percent)]]) {
            const stat = el('div'); stat.append(el('span', title), el('strong', value)); stats.append(stat);
        }
        if (row.status === 'matched') detail.append(stats);
        const flags = [];
        for (const side of ['a', 'b']) for (const key of ['partial', 'pseudo', 'pseudogene', 'transl_except', 'exception']) {
            const value = row[`${key}_${side}`];
            if (value !== null && value !== '') flags.push(`${labels[side === 'a' ? 0 : 1]}: ${key}=${value}`);
        }
        if (flags.length) detail.append(el('p', 'Source annotation context · ' + flags.join(' · '), 'gene-flags'));
        const plot = el('div', undefined, 'gene-locus');
        let svg = locus(row);
        const contextLabel = el('label', undefined, 'gene-context-toggle'), contextToggle = el('input');
        contextToggle.type = 'checkbox';
        contextToggle.disabled = locusModels(row).some(model => model.wraps);
        if (contextToggle.disabled) contextToggle.title = 'Flanking context is not inferred for circular-origin models; inspect the separate coding blocks and segment coordinates.';
        contextLabel.append(contextToggle, document.createTextNode(' Show neighbouring CDS calls (grey; context only)'));
        contextToggle.addEventListener('change', () => {svg = locus(row, contextToggle.checked); plot.replaceChildren(svg);});
        detail.append(contextLabel);
        plot.append(svg); detail.append(plot);
        detail.append(el('p', 'Shared axis increases left to right; arrows indicate strand. Source colours identify annotations, not confidence; outlined arrows belong to the selected link. Grey neighbours are context only. Discontinuous segments are separate blocks, not continuous coding DNA. Use the segment table for circular-origin models.', 'note'));
        const download = el('button', 'Download locus SVG');
        download.type = 'button';
        download.addEventListener('click', () => {
            const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(svg)], {type: 'image/svg+xml'}));
            const link = el('a'); link.href = url; link.download = `gene-locus-${row.index+1}.svg`; link.click();
            setTimeout(() => URL.revokeObjectURL(url), 1000);
        });
        detail.append(download);
        for (const end of ['n', 'c']) {
            if (row[`${end}_terminal_delta_aa`]) {
                const side = row[`${end}_terminal_residues_a`] ? 'a' : 'b';
                detail.append(el('p', `${end.toUpperCase()}-terminal extra residues in ${side.toUpperCase()}:`));
                detail.append(el('pre', row[`${end}_terminal_residues_${side}`]));
            }
        }
        if (row.alignment_a && row.alignment_b) {
            detail.append(el('h4', 'Changed-sequence window'));
            detail.append(el('p', 'The window shows up to 64 columns around the first difference, including up to 16 preceding columns. Protein positions are 1-based; “—” means no residue in that window.', 'note'));
            const legend = el('p', undefined, 'gene-alignment-legend');
            for (const [label, cls] of [['Codon present in one model only', 'gene-gap'], ['Non-gap amino-acid mismatch', 'gene-mismatch']]) {
                const item = el('span'); item.append(el('b', label, cls)); legend.append(item);
            }
            detail.append(legend);
            detail.append(alignment(row));
            const full = el('details'); full.append(el('summary', 'Full protein alignment'));
            // Render long alignments only when requested, not for every row.
            full.addEventListener('toggle', () => {if (full.open && full.childElementCount === 1) full.append(alignment(row, true));});
            detail.append(full);
        }
        const measurements = el('details'); measurements.append(el('summary', 'Codon-mapped alignment measurements'));
        const metrics = el('table', undefined, 'gene-measurements'), metricBody = el('tbody');
        for (const [key, title] of [['shared_codon_count', 'Shared genomic codons'], ['identical_aa', 'Identical aa'],
            ['substitution_count', 'All non-gap mismatches'], ['initiator_mismatch_count', 'Initiation-context mismatches'],
            ['internal_mismatch_count', 'Other codon-mapped mismatches'], ['gap_aa', 'Gap aa']]) {
            const tr = el('tr'), label = el('th', title); label.scope = 'row';
            tr.append(label, el('td', measured(row[key]))); metricBody.append(tr);
        }
        metrics.append(metricBody); measurements.append(metrics); if (row.status === 'matched') detail.append(measurements);
        function sourceTable(fields) {
            const wrapper = el('div', undefined, 'table-wrap'), models = el('table', undefined, 'gene-models');
            const titles = el('tr');
            for (const title of ['Source evidence', `A: ${labels[0]}`, `B: ${labels[1]}`]) titles.append(el('th', title));
            const head = el('thead'); head.append(titles); models.append(head);
            const body = el('tbody');
            for (const [key, title] of fields) {
                const tr = el('tr'), label = el('th', title); label.scope = 'row'; tr.append(label);
                for (const side of ['a', 'b']) {
                    const value = row[`${key}_${side}`];
                    tr.append(el('td', row[`id_${side}`] === null ? 'No counterpart' :
                        key === 'start_context' ? startContexts[value] || shown(value) :
                        value === null && ['gene','product','dbxref','partial','pseudo','pseudogene','exception','transl_except','start_range','end_range'].includes(key) ? 'Not provided' : shown(value)));
                }
                body.append(tr);
            }
            models.append(body); wrapper.append(models); return wrapper;
        }
        detail.append(el('h4', 'Source CDS models'));
        detail.append(sourceTable([['id', 'CDS ID'], ['gene', 'Gene symbol'], ['product', 'Product label'],
            ['segments', 'Segments (1-based inclusive; biological order)'], ['phase', 'Initial GFF3 phase'],
            ['cds_length', 'CDS length (nt; all segments)'], ['protein_length', 'Protein length (aa)'],
            ['genetic_code', 'Translation table'], ['start_codon', 'First translated codon (DNA)'],
            ['first_codon_positions', 'First codon genomic bases (1-based; translation order)'],
            ['first_amino_acid', 'First predicted amino acid'], ['start_context', 'First-codon translation context']]));
        detail.append(el('p', 'First codons belong to each model’s own coding boundary. A shared codon can be an initiator in one model and internal in the other. Initiator rules describe the predicted translation, not experimental proof of initiation.', 'note'));
        const context = el('details'); context.append(el('summary', 'Source annotations and database references'));
        context.append(sourceTable([['partial', 'Partial attribute'], ['start_range', 'Start range'], ['end_range', 'End range'],
            ['pseudo', 'Pseudo attribute'], ['pseudogene', 'Pseudogene attribute'], ['transl_except', 'Translation exceptions'],
            ['unapplied_transl_except', 'Unapplied translation exceptions (translation provisional)'],
            ['exception', 'Other exceptions'], ['dbxref', 'Database cross-references']])); detail.append(context);
        detail.append(el('p', `${row.status === 'matched' ? 'Boundary deltas' : 'Selected-call boundary offsets (descriptive only)'} (B−A, bp, positive downstream): 5′ ${shown(row.start_delta_bp)}; 3′ ${shown(row.stop_delta_bp)}. Changed source attributes: ${row.changed_attributes || 'none / not applicable'}. Alignment: ${shown(row.alignment_method)}.`, 'note'));
        detail.append(el('p', 'These measurements do not identify the correct gene model or demonstrate altered function. Source annotations of partial/pseudogene status are context, not validation.', 'note'));
    }
    document.getElementById('gene-report').addEventListener('click', event => {
        const category = event.target.closest('button[data-consequence]');
        if (category) {
            const value = category.dataset.consequence;
            document.getElementById('gene-consequence').value = value;
            document.getElementById('gene-structure').value = '';
            document.getElementById('gene-view').value = 'all';
            table.search.fixed('gene-view', null).search('').search.fixed('gene-structure', null).search.fixed('gene-consequence', (search, row) => row.protein_consequence === value).draw();
            document.getElementById('gene-evidence').scrollIntoView({block: 'start'});
            document.querySelector('#gene-evidence input[type=search]').focus({preventScroll: true});
        }
        const button = event.target.closest('[data-gene-index]');
        if (button) {showDetail(rows[Number(button.dataset.geneIndex)]); document.getElementById('gene-detail').scrollIntoView({block: 'start'});}
    });
    const example = rows.find(row => row.status === 'ambiguous') || rows.find(row => row.status === 'matched' && row.protein_consequence === 'terminal_only') || rows.find(needsReview);
    if (example) showDetail(example);
    else if (rows.length) document.getElementById('gene-detail').append(el('p', 'No differences or review findings in the default view. Select All relationships to inspect unchanged pairs.'));
})();
