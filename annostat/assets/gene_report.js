/* Annostat gene-level report. Source attributes are rendered as text, never HTML. */
(() => {
    'use strict';
    const gene = JSON.parse(document.getElementById('gene-data').textContent);
    const rows = gene.relationships.map((row, index) => ({...row, index}));
    const labels = Object.keys(gene.cds_counts);
    const names = {identical: 'Identical protein', terminal_only: 'Terminal-only difference',
        initiator_recoding: 'Initiation-context difference', translation_difference: 'Other translated-sequence difference',
        limited: 'Limited protein evidence', unresolved: 'Unresolved correspondence'};
    const relationshipNames = {exact: 'Exact CDS structure', start_shift: 'Alternative start', stop_shift: 'Alternative stop',
        frame_disagreement: 'Frame disagreement', overlap: 'Other overlap', unmatched: 'Unmatched call'};
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
            column('Gene / CDS A → B', row => `${row.gene_a || row.gene_b || 'unnamed'}\n${row.id_a || '—'} → ${row.id_b || '—'}`),
            column('Correspondence', row => `${row.status === 'matched' ? 'Resolved pair' : row.status}\n${relationshipNames[row.relationship]}`),
            column('Protein evidence', row => names[row.protein_consequence]),
            column('Protein aa A → B', row => `${row.protein_length_a ?? '—'} → ${row.protein_length_b ?? '—'}`),
            {title: 'Δ aa B − A', data: 'protein_length_delta_aa', defaultContent: '—', type: 'num',
                render: (value, type) => type === 'display' ? value === null ? '—' : signed(value) : value},
            {...column('Source context', row => ['seqid', 'product_a', 'product_b', 'gene_b', 'changed_attributes', 'segments_a', 'segments_b',
                'partial_a', 'partial_b', 'pseudo_a', 'pseudo_b', 'pseudogene_a', 'pseudogene_b', 'exception_a', 'exception_b',
                'dbxref_a', 'dbxref_b', 'status', 'relationship'].map(key => row[key] || '').join(' ')), visible: false}
        ], language: {emptyTable: 'No CDS relationships in these annotations.', search: 'Search all relationships:'}
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
    function locus(row) {
        const blocks = ['a', 'b'].map(side => (row[`segments_${side}`] || '').split(';').filter(Boolean).map(block => block.split('-').map(Number)));
        const points = blocks.flat(2);
        const low = Math.min(...points), high = Math.max(...points);
        const scale = pos => 70 + 760 * (pos - low) / Math.max(1, high - low + 1);
        const svg = svgNode('svg', {viewBox: '0 0 900 225', role: 'img', 'aria-label': 'CDS genomic coordinates, 1-based inclusive'});
        svg.append(svgNode('title', {}, `${row.gene_a || row.gene_b || 'CDS'}: ${row.id_a || '—'} / ${row.id_b || '—'}`));
        svg.append(svgNode('desc', {}, `${row.seqid}; strand ${row.strand}; genomic coordinates 1-based inclusive. A ${labels[0]}: ${row.segments_a || 'no counterpart'}. B ${labels[1]}: ${row.segments_b || 'no counterpart'}. Positive boundary deltas point downstream; 5-prime delta ${shown(row.start_delta_bp)} bp; 3-prime delta ${shown(row.stop_delta_bp)} bp. Not evidence of functional impact.`));
        svg.append(svgNode('rect', {width: 900, height: 225, fill: '#fff'}));
        svg.append(svgNode('text', {x: 20, y: 22, fill: '#17231e', 'font-size': 14, 'font-family': 'Arial'}, `${row.seqid} | genomic coordinates (1-based inclusive); strand ${row.strand}`));
        blocks.forEach((segments, side) => {
            const y = 68 + side * 70, color = side ? '#2563eb' : '#087f5b';
            const suffix = side ? 'b' : 'a';
            const caption = `${labels[side]} | ${row[`id_${suffix}`] || 'No counterpart'} | ${row[`segments_${suffix}`] || '—'}`;
            svg.append(svgNode('text', {x: 70, y: y-9, fill: color, 'font-size': 12, 'font-family': 'Arial'}, caption.length > 105 ? caption.slice(0, 102) + '…' : caption));
            svg.append(svgNode('text', {x: 20, y: y + 17, fill: color, 'font-size': 16, 'font-family': 'Arial'}, side ? 'B' : 'A'));
            segments.forEach(([start, end]) => {
                const x = scale(start), width = Math.max(1, scale(end + 1) - x), tip = Math.min(9, width / 2);
                const shape = row.strand === '+' ? `${x},${y} ${x+width-tip},${y} ${x+width},${y+12} ${x+width-tip},${y+24} ${x},${y+24}`
                    : `${x+tip},${y} ${x+width},${y} ${x+width},${y+24} ${x+tip},${y+24} ${x},${y+12}`;
                const arrow = svgNode('polygon', {points: shape, fill: color});
                arrow.append(svgNode('title', {}, `${side ? 'B' : 'A'}: ${start}–${end} (inclusive)`));
                svg.append(arrow);
            });
        });
        svg.append(svgNode('text', {x: 70, y: 201, fill: '#64746c', 'font-size': 13, 'font-family': 'Arial'}, low.toLocaleString()));
        svg.append(svgNode('text', {x: 830, y: 201, 'text-anchor': 'end', fill: '#64746c', 'font-size': 13, 'font-family': 'Arial'}, high.toLocaleString()));
        return svg;
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
        header.append(el('h3', row.gene_a || row.gene_b || 'CDS evidence'));
        header.append(el('p', `${labels[0]}: ${row.id_a || 'No counterpart'} → ${labels[1]}: ${row.id_b || 'No counterpart'}`));
        header.append(el('p', `${row.status === 'matched' ? 'Resolved pair' : row.status} · ${relationshipNames[row.relationship]} · ${names[row.protein_consequence]}`));
        detail.replaceChildren(header);
        const length = side => row[`id_${side}`] === null ? 'no counterpart' : `${shown(row[`protein_length_${side}`])} aa`;
        const lengths = `Predicted protein lengths: A ${length('a')}; B ${length('b')}. `;
        detail.append(el('p', lengths + row.consequence_note, 'gene-interpretation'));
        const stats = el('div', undefined, 'gene-detail-stats');
        for (const [title, value] of [['Protein length Δ (B−A, aa)', signed(row.protein_length_delta_aa)],
            ['N-terminal Δ (B−A, aa)', signed(row.n_terminal_delta_aa)], ['C-terminal Δ (B−A, aa)', signed(row.c_terminal_delta_aa)],
            ['Identity (including gaps, %)', measured(row.alignment_identity_percent)]]) {
            const stat = el('div'); stat.append(el('span', title), el('strong', value)); stats.append(stat);
        }
        detail.append(stats);
        const flags = [];
        for (const side of ['a', 'b']) for (const key of ['partial', 'pseudo', 'pseudogene', 'transl_except', 'exception']) {
            const value = row[`${key}_${side}`];
            if (value !== null && value !== '') flags.push(`${labels[side === 'a' ? 0 : 1]}: ${key}=${value}`);
        }
        if (flags.length) detail.append(el('p', 'Source annotation context · ' + flags.join(' · '), 'gene-flags'));
        const plot = el('div', undefined, 'gene-locus');
        const svg = locus(row);
        plot.append(svg); detail.append(plot);
        detail.append(el('p', 'Axis increases left to right; arrows indicate strand. Discontinuous segments are separate blocks, not continuous coding DNA. Use the segment table for circular-origin models.', 'note'));
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
        metrics.append(metricBody); measurements.append(metrics); detail.append(measurements);
        function sourceTable(fields) {
            const wrapper = el('div', undefined, 'table-wrap'), models = el('table', undefined, 'gene-models');
            const titles = el('tr');
            for (const title of ['Source evidence', `A: ${labels[0]}`, `B: ${labels[1]}`]) titles.append(el('th', title));
            const head = el('thead'); head.append(titles); models.append(head);
            const body = el('tbody');
            for (const [key, title] of fields) {
                const tr = el('tr'), label = el('th', title); label.scope = 'row'; tr.append(label);
                for (const side of ['a', 'b']) tr.append(el('td', row[`id_${side}`] === null ? 'No counterpart' : row[`${key}_${side}`] === null && ['gene','product','dbxref','partial','pseudo','pseudogene','exception','transl_except','start_range','end_range'].includes(key) ? 'Not provided' : shown(row[`${key}_${side}`])));
                body.append(tr);
            }
            models.append(body); wrapper.append(models); return wrapper;
        }
        detail.append(el('h4', 'Source CDS models'));
        detail.append(sourceTable([['id', 'CDS ID'], ['gene', 'Gene symbol'], ['product', 'Product label'],
            ['segments', 'Segments (1-based inclusive; biological order)'], ['phase', 'Initial GFF3 phase'],
            ['cds_length', 'CDS length (nt; all segments)'], ['protein_length', 'Protein length (aa)'],
            ['genetic_code', 'Translation table'], ['start_codon', 'First translated codon (not proof of initiation)']]));
        const context = el('details'); context.append(el('summary', 'Source annotations and database references'));
        context.append(sourceTable([['partial', 'Partial attribute'], ['start_range', 'Start range'], ['end_range', 'End range'],
            ['pseudo', 'Pseudo attribute'], ['pseudogene', 'Pseudogene attribute'], ['transl_except', 'Translation exceptions'],
            ['exception', 'Other exceptions'], ['dbxref', 'Database cross-references']])); detail.append(context);
        detail.append(el('p', `Boundary deltas (B−A, bp, positive downstream): 5′ ${shown(row.start_delta_bp)}; 3′ ${shown(row.stop_delta_bp)}. Changed source attributes: ${row.changed_attributes || 'none / not applicable'}. Alignment: ${shown(row.alignment_method)}.`, 'note'));
        detail.append(el('p', 'These measurements do not identify the correct gene model or demonstrate altered function. Source annotations of partial/pseudogene status are context, not validation.', 'note'));
    }
    document.getElementById('gene-report').addEventListener('click', event => {
        const category = event.target.closest('button[data-consequence]');
        if (category) {
            const value = category.dataset.consequence;
            document.getElementById('gene-consequence').value = value;
            document.getElementById('gene-structure').value = '';
            table.search('').search.fixed('gene-structure', null).search.fixed('gene-consequence', (search, row) => row.protein_consequence === value).draw();
            document.getElementById('gene-evidence').scrollIntoView({block: 'start'});
            document.querySelector('#gene-evidence input[type=search]').focus({preventScroll: true});
        }
        const button = event.target.closest('[data-gene-index]');
        if (button) {showDetail(rows[Number(button.dataset.geneIndex)]); document.getElementById('gene-detail').scrollIntoView({block: 'start'});}
    });
    const example = rows.find(row => row.status === 'matched' && row.protein_consequence === 'terminal_only') || rows[0];
    if (example) showDetail(example);
})();
