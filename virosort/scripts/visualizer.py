import os
from typing import Dict, List, Any, Optional
from collections import defaultdict


def _generate_svg_probes(cluster_probes: List[tuple], consensus_len: int, probe_colors: List[str]) -> List[str]:
    """
    Generate SVG elements for probe visualization.
    
    Args:
        cluster_probes: list of (probe_id, probe_result) tuples
        consensus_len: length of consensus sequence
        probe_colors: list of colors for probes
    
    Returns:
        list of SVG string elements
    """
    svg_elements = []
    
    # Add consensus backbone
    svg_elements.append(f'''
        <line x1="10" y1="30" x2="{consensus_len + 10}" y2="30" 
              stroke="#333" stroke-width="3" stroke-linecap="round"/>
    ''')
    
    # Add each probe
    for i, (probe_id, probe_result) in enumerate(cluster_probes):
        color = probe_colors[i % len(probe_colors)]
        start = probe_result.get('start_pos', 0) + 10
        end = probe_result.get('end_pos', start + 50) + 10
        y = 35 + (i % 10) * 35  # Wrap after 10 probes
        
        # Check for adapters
        adapter_info = probe_result.get('adapter_info', {})
        has_left = adapter_info.get('has_adapter_left', False)
        has_right = adapter_info.get('has_adapter_right', False)
        
        # Draw probe body
        svg_elements.append(f'''
            <title>{probe_id}: {probe_result.get("core_sequence", "")[:50]}...</title>
            <line x1="{start}" y1="{y}" x2="{end}" y2="{y}" 
                  stroke="{color}" stroke-width="4" stroke-linecap="round"
                  opacity="0.7"
                  class="probe-line"
                  data-probe-id="{probe_id}"
                  data-start="{start}"
                  data-end="{end}"/>
        ''')
        
        # Mark adapters if present
        if has_left:
            svg_elements.append(f'<rect x="{start}" y="{y-5}" width="5" height="5" fill="#ff0000" opacity="0.5"/>')
        if has_right:
            svg_elements.append(f'<rect x="{end-5}" y="{y-5}" width="5" height="5" fill="#0000ff" opacity="0.5"/>')
    
    return svg_elements


def _generate_table_rows(cluster_probes: List[tuple]) -> List[str]:
    """
    Generate HTML table rows for probe mapping results.
    
    Args:
        cluster_probes: list of (probe_id, probe_result) tuples
    
    Returns:
        list of HTML table row strings
    """
    rows = []
    
    for probe_id, probe_result in cluster_probes:
        start = probe_result.get('start_pos', 'N/A')
        end = probe_result.get('end_pos', 'N/A')
        coverage = probe_result.get('coverage', 0)
        score = probe_result.get('score', 0)
        adapter_info = probe_result.get('adapter_info', {})
        has_left = adapter_info.get('has_adapter_left', False)
        has_right = adapter_info.get('has_adapter_right', False)
        
        rows.append(f'''
            <tr>
                <td title="{probe_result.get("core_sequence", "")[:80]}...">{probe_id}</td>
                <td>{start}</td>
                <td>{end}</td>
                <td>{coverage:.3f}</td>
                <td>{score:.3f}</td>
                <td>{'Yes' if has_left else 'No'}</td>
                <td>{'Yes' if has_right else 'No'}</td>
            </tr>
        ''')
    
    return rows


def _generate_cluster_html_content(cluster_num: str,
                                   consensus_seq: str,
                                   cluster_probes: List[tuple],
                                   taxonomy: Dict[str, str],
                                   analysis_stats: Dict[str, float],
                                   width: int,
                                   height: int,
                                   length_disparity: Dict = None,
                                   iupac_analysis: Dict = None,
                                   gap_analysis: Dict = None) -> str:
    """
    Generate HTML content for cluster visualization with extended analysis.
    
    Args:
        cluster_num: cluster number
        consensus_seq: consensus sequence
        cluster_probes: list of (probe_id, probe_result) tuples
        taxonomy: taxonomy info
        analysis_stats: analysis stats
        width: visualization width
        height: visualization height
        length_disparity: length disparity analysis for this cluster
        iupac_analysis: IUPAC analysis for this cluster
        gap_analysis: gap analysis for this cluster
    
    Returns:
        HTML string
    """
    consensus_len = len(consensus_seq)
    
    # Color scheme for probes
    probe_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    # Prepare analysis data
    ld_info = length_disparity if length_disparity else {}
    iupac_info = iupac_analysis if iupac_analysis else {}
    gap_info = gap_analysis if gap_analysis else {}
    
    # Build analysis cards
    length_card = ''
    if ld_info:
        disparity_status = 'Yes' if ld_info.get("has_disparity") else 'No'
        disparity_class = 'yes' if ld_info.get("has_disparity") else 'no'
        length_card = f'''<div class="analysis-card">
<h3>Length Disparity Analysis</h3>
<table class="analysis-table">
  <tr><th>Has Disparity</th><td><span class="disparity-badge disparity-{disparity_class}">{disparity_status}</span></td></tr>
  <tr><th>Num Sequences</th><td>{ld_info.get("num_sequences", "N/A")}</td></tr>
  <tr><th>Shortest Length</th><td>{ld_info.get("shortest_length", "N/A")} bp</td></tr>
  <tr><th>Longest Length</th><td>{ld_info.get("longest_length", "N/A")} bp</td></tr>
  <tr><th>Reference Length</th><td>{ld_info.get("reference_length", "N/A")} bp</td></tr>
  <tr><th>Reference ID</th><td>{ld_info.get("reference_id", "N/A")}</td></tr>
  <tr><th>Length Ratio</th><td>{format(ld_info.get('length_ratio', 0), '.4f')}</td></tr>
  <tr><th>Extended Leading</th><td>{ld_info.get("extended_leading_length", "N/A")} bp</td></tr>
  <tr><th>Extended Trailing</th><td>{ld_info.get("extended_trailing_length", "N/A")} bp</td></tr>
</table>
</div>'''
    
    iupac_card = ''
    if iupac_info:
        iupac_card = f'''<div class="analysis-card">
<h3>IUPAC Analysis</h3>
<table class="analysis-table">
  <tr><th>Status</th><td>{iupac_info.get("status", "N/A")}</td></tr>
  <tr><th>Num Sequences</th><td>{iupac_info.get("num_sequences", "N/A")}</td></tr>
  <tr><th>Consensus Length</th><td>{iupac_info.get("consensus_length", "N/A")} bp</td></tr>
  <tr><th>Total IUPAC</th><td>{iupac_info.get("total_iupac", "N/A")}</td></tr>
  <tr><th>IUPAC Fraction</th><td>{iupac_info.get("iupac_fraction", "N/A")}</td></tr>
</table>'''
        if iupac_info.get("iupac_count_by_type"):
            iupac_card += '''<h4>IUPAC Codes Breakdown</h4>
<div class="iupac-summary">
'''
            iupac_card += ''.join(f'  <div class="iupac-box"><span class="code">{code}</span>: <span class="count">{count}</span></div>\n' 
                                   for code, count in sorted(iupac_info.get("iupac_count_by_type", {}).items()) if count > 0)
            iupac_card += '</div>\n</div>'
        else:
            iupac_card += '''</div>
</div>'''
    
    gap_card = ''
    if gap_info:
        gap_card = f'''<div class="analysis-card">
<h3>Gap Analysis</h3>
<table class="analysis-table">
  <tr><th>Status</th><td>{gap_info.get("status", "N/A")}</td></tr>
  <tr><th>Num Sequences</th><td>{gap_info.get("num_sequences", "N/A")}</td></tr>
  <tr><th>Original Length</th><td>{gap_info.get("orig_length", "N/A")} bp</td></tr>
  <tr><th>Alignment Length</th><td>{gap_info.get("aln_length", "N/A")} bp</td></tr>
  <tr><th>Total Gaps</th><td>{gap_info.get("gaps", "N/A")}</td></tr>
  <tr><th>Gap Percentage</th><td>{format(gap_info.get('gap_percent', 0), '.2f')}%</td></tr>
  <tr><th>Insertions</th><td>{gap_info.get("insertions", "N/A")}</td></tr>
</table>
</div>'''
    
    # Build analysis grid
    analysis_grid = '<div class="analysis-grid">\n'
    if ld_info:
        analysis_grid += length_card + '\n'
    if iupac_info:
        analysis_grid += iupac_card + '\n'
    if gap_info:
        analysis_grid += gap_card + '\n'
    analysis_grid += '</div>'
    
    # Build HTML
    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width={width},height={height}">
    <title>Cluster {cluster_num} - {taxonomy.get('family', 'Unknown')} / {taxonomy.get('species', 'Unknown')}</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 20px;
            background-color: #f5f5f5;
        }}
        .header {{
            background-color: #ffffff;
            padding: 15px;
            border-radius: 8px;
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .header h1 {{
            margin: 0 0 10px 0;
            color: #333;
        }}
        .header .taxonomy {{
            color: #666;
            font-size: 14px;
        }}
        .stats {{
            display: flex;
            flex-wrap: wrap;
            gap: 15px;
            margin-top: 10px;
        }}
        .stat-box {{
            background-color: #e8f4f8;
            padding: 8px 12px;
            border-radius: 4px;
            font-size: 12px;
        }}
        .stat-label {{
            font-weight: bold;
            color: #006699;
        }}
        .visualization {{
            background-color: #ffffff;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 20px;
        }}
        .consensus-container {{
            overflow-x: auto;
            margin: 20px 0;
            padding: 10px;
            background-color: #fafafa;
            border-radius: 4px;
            font-family: 'Courier New', monospace;
            font-size: 11px;
            word-break: break-all;
        }}
        .consensus-sequence {{
            color: #333;
            line-height: 1.4;
        }}
        .consensus-position {{
            color: #999;
            font-size: 10px;
            margin-bottom: 5px;
        }}
        .analysis-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(350px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }}
        .analysis-card {{
            background-color: #f9f9f9;
            border: 1px solid #ddd;
            border-radius: 6px;
            padding: 15px;
        }}
        .analysis-card h3 {{
            margin: 0 0 15px 0;
            color: #006699;
            font-size: 14px;
            border-bottom: 2px solid #006699;
            padding-bottom: 8px;
        }}
        .analysis-table {{
            width: 100%;
            font-size: 11px;
            border-collapse: collapse;
        }}
        .analysis-table th {{
            text-align: left;
            padding: 6px;
            color: #666;
            font-weight: normal;
        }}
        .analysis-table td {{
            padding: 6px;
            color: #333;
        }}
        .analysis-table tr:nth-child(even) {{
            background-color: #f0f0f0;
        }}
        .disparity-badge {{
            display: inline-block;
            padding: 3px 8px;
            border-radius: 3px;
            font-size: 10px;
            font-weight: bold;
        }}
        .disparity-yes {{
            background-color: #ffe0b2;
            color: #e65100;
        }}
        .disparity-no {{
            background-color: #c8e6c9;
            color: #2e7d32;
        }}
        .iupac-summary {{
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
            margin-top: 10px;
        }}
        .iupac-box {{
            background-color: #fff3e0;
            border: 1px solid #ffb74d;
            border-radius: 4px;
            padding: 6px 10px;
            font-size: 11px;
        }}
        .iupac-box .code {{
            font-weight: bold;
            color: #e65100;
        }}
        .iupac-box .count {{
            color: #666;
        }}
        .probe-track {{
            position: relative;
            height: 40px;
            margin: 5px 0;
            background-color: #f0f0f0;
            border-radius: 4px;
        }}
        .probe-rect {{
            position: absolute;
            height: 30px;
            top: 5px;
            border-radius: 3px;
            cursor: pointer;
            transition: opacity 0.2s;
        }}
        .probe-rect:hover {{
            opacity: 0.7;
        }}
        .probe-rect.adapter-left {{
            border-left: 3px solid rgba(255,0,0,0.6);
        }}
        .probe-rect.adapter-right {{
            border-right: 3px solid rgba(0,0,255,0.6);
        }}
        .probe-tooltip {{
            position: absolute;
            background-color: rgba(0,0,0,0.85);
            color: white;
            padding: 8px 12px;
            border-radius: 4px;
            font-size: 11px;
            white-space: nowrap;
            z-index: 100;
            pointer-events: none;
            display: none;
        }}
        .probe-list {{
            margin-top: 20px;
        }}
        .probe-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 12px;
        }}
        .probe-table th {{
            background-color: #006699;
            color: white;
            padding: 10px;
            text-align: left;
        }}
        .probe-table td {{
            padding: 8px 10px;
            border-bottom: 1px solid #ddd;
        }}
        .probe-table tr:hover {{
            background-color: #f0f7fb;
        }}
        .adapter-badge {{
            display: inline-block;
            padding: 2px 6px;
            border-radius: 3px;
            font-size: 10px;
            margin: 0 2px;
        }}
        .adapter-left-badge {{
            background-color: rgba(255,0,0,0.2);
            border: 1px solid rgba(255,0,0,0.4);
        }}
        .adapter-right-badge {{
            background-color: rgba(0,0,255,0.2);
            border: 1px solid rgba(0,0,255,0.4);
        }}
        .legend {{
            display: flex;
            gap: 20px;
            margin-top: 15px;
            font-size: 12px;
        }}
        .legend-item {{
            display: flex;
            align-items: center;
            gap: 5px;
        }}
        .legend-color {{
            width: 15px;
            height: 15px;
            border-radius: 3px;
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Cluster {cluster_num} Visualization</h1>
        <div class="taxonomy">
            <strong>Family:</strong> {taxonomy.get('family', 'Unknown')} | 
            <strong>Genus:</strong> {taxonomy.get('genus', 'Unknown')} | 
            <strong>Species:</strong> {taxonomy.get('species', 'Unknown')}
        </div>
        <div class="stats">
            <div class="stat-box"><span class="stat-label">Consensus Length:</span> {consensus_len} bp</div>
            <div class="stat-box"><span class="stat-label">Mapped Probes:</span> {len(cluster_probes)}</div>
            <div class="stat-box"><span class="stat-label">Avg Coverage:</span> {analysis_stats.get('avg_coverage_over_unambig', 0):.3f}x</div>
            <div class="stat-box"><span class="stat-label">Frac Covered:</span> {analysis_stats.get('frac_bases_covered', 0):.4f}</div>
            {'<div class="stat-box"><span class="stat-label">Length Disparity:</span> ' + ('Yes' if ld_info.get("has_disparity") else 'No') + '</div>' if ld_info else ''}
            {'<div class="stat-box"><span class="stat-label">Total IUPAC:</span> ' + str(iupac_info.get("total_iupac", 0)) + '</div>' if iupac_info else ''}
            {'<div class="stat-box"><span class="stat-label">Gap %:</span> ' + format(gap_info.get('gap_percent', 0), '.2f') + '%' + '</div>' if gap_info else ''}
        </div>
    </div>
    
    <div class="visualization">
        <h2>Consensus Sequence</h2>
        <div class="consensus-position">
            Positions: 1 - {consensus_len}
        </div>
        <div class="consensus-container">
            <div class="consensus-sequence">{consensus_seq}</div>
        </div>
        
        <h2>Probe Mapping</h2>
        <svg id="probe-viz" width="{width}" height="{max(80 + len(cluster_probes) * 40, height - 200)}">
            <defs>
                <style>
                    .probe-tooltip {{
                        font-family: Arial, sans-serif;
                        font-size: 11px;
                    }}
                </style>
            </defs>
            {"".join(_generate_svg_probes(cluster_probes, consensus_len, probe_colors))}
        </svg>
        
        <div class="legend">
            <div class="legend-item">
                <div class="legend-color" style="background-color: #ff000055;"></div>
                <span>Adapter (5')</span>
            </div>
            <div class="legend-item">
                <div class="legend-color" style="background-color: #0000ff55;"></div>
                <span>Adapter (3')</span>
            </div>
            <div class="legend-item">
                <div class="legend-color" style="background: linear-gradient(to right, '#1f77b4', '#ff7f0e', '#2ca02c');"></div>
                <span>Probes (color-coded)</span>
            </div>
        </div>
        
        <div class="probe-list">
            <h3>Probe Details</h3>
            <table class="probe-table">
                <thead>
                    <tr>
                        <th>Probe ID</th>
                        <th>Start</th>
                        <th>End</th>
                        <th>Coverage</th>
                        <th>Score</th>
                        <th>5' Adapter</th>
                        <th>3' Adapter</th>
                    </tr>
                </thead>
                <tbody>
                    {"".join(_generate_table_rows(cluster_probes))}
                </tbody>
            </table>
        </div>
    </div>
    
    {analysis_grid}
    
    <script>
        // Tooltip functionality
        const svg = document.getElementById('probe-viz');
        const tooltip = document.createElement('div');
        tooltip.className = 'probe-tooltip';
        document.body.appendChild(tooltip);
        
        document.querySelectorAll('.probe-rect').forEach(rect => {{
            rect.addEventListener('mouseenter', (e) => {{
                tooltip.style.display = 'block';
                tooltip.textContent = e.target.title;
            }});
            
            rect.addEventListener('mousemove', (e) => {{
                tooltip.style.left = (e.pageX + 10) + 'px';
                tooltip.style.top = (e.pageY - 10) + 'px';
            }});
            
            rect.addEventListener('mouseleave', () => {{
                tooltip.style.display = 'none';
            }});
        }});
    </script>
</body>
</html>'''
    
    return html


def generate_all_cluster_visualizations(
    mapping_results: Dict[str, List[tuple]],
    consensus_clusters: Dict[str, List[str]],
    cluster_taxonomy_summary: Dict[str, Dict[str, str]],
    analysis_stats: Dict[str, Dict[str, float]],
    output_dir: str,
    width: int = 1200,
    height: int = 600,
    length_disparity: Dict = None,
    iupac_analysis: Dict = None,
    gap_analysis: Dict = None
) -> List[str]:
    """
    Generate HTML visualizations for all clusters.
    
    Args:
        mapping_results: probe mapping results per cluster
        consensus_clusters: consensus sequences per cluster
        cluster_taxonomy_summary: taxonomy info per cluster
        analysis_stats: analysis statistics per cluster
        output_dir: output directory for HTML files
        width: visualization width
        height: visualization height
        length_disparity: length disparity analysis
        iupac_analysis: IUPAC analysis
        gap_analysis: gap analysis
    
    Returns:
        List of generated HTML file paths
    """
    os.makedirs(output_dir, exist_ok=True)
    html_files = []
    
    # Group probes by cluster_num
    probes_by_cluster = {}
    for probe_id, result in mapping_results.items():
        cluster_num = result.get('cluster_num')
        if cluster_num is not None:
            if cluster_num not in probes_by_cluster:
                probes_by_cluster[cluster_num] = []
            probes_by_cluster[cluster_num].append((probe_id, result))
    
    # Generate visualization for each cluster that has consensus sequences
    for cluster_num in consensus_clusters.keys():
        consensus_dict = consensus_clusters.get(cluster_num, {})
        consensus_seq = consensus_dict.get('sequence') if isinstance(consensus_dict, dict) else (consensus_dict[0] if consensus_dict else None)
        cluster_probes = probes_by_cluster.get(cluster_num, [])
        taxonomy = cluster_taxonomy_summary.get(cluster_num, {})
        stats = analysis_stats.get(cluster_num, {})
        
        # Get cluster-specific analysis data
        ld_data = length_disparity.get(cluster_num, {}) if length_disparity else {}
        iupac_data = iupac_analysis.get(cluster_num, {}) if iupac_analysis else {}
        gap_data = gap_analysis.get(cluster_num, {}) if gap_analysis else {}
        
        html_content = _generate_cluster_html_content(
            cluster_num,
            consensus_seq,
            cluster_probes,
            taxonomy,
            stats,
            width,
            height,
            length_disparity=ld_data,
            iupac_analysis=iupac_data,
            gap_analysis=gap_data
        )
        
        # Write HTML file
        html_filename = f"cluster_{cluster_num}.html"
        html_path = os.path.join(output_dir, html_filename)
        with open(html_path, 'w') as f:
            f.write(html_content)
        
        html_files.append(html_path)
    
    return html_files


def generate_all_cluster_visualizations_revised(
    degapped_consensus: Dict[str, str],
    probe_details: Dict[str, Dict],
    output_dir: str,
    metadata_map: Dict[str, Dict],
    cluster_accessions: Dict[str, List[str]],
    adapter_seq: str = 'GTGGAGGTCGCTAGATGGTC',
    width: int = 1200,
    height: int = 600
) -> List[str]:
    """
    Generate HTML visualizations for all clusters with new probe mapping format.
    
    Args:
        degapped_consensus: dict mapping cluster_num to consensus sequence
        probe_details: dict mapping probe_id to probe details
        output_dir: output directory for HTML files
        metadata_map: dict mapping accession_id to taxonomy info
        cluster_accessions: dict mapping cluster_num to list of accession IDs
        adapter_seq: adapter sequence
        width: visualization width
        height: visualization height
    
    Returns:
        List of generated HTML file paths
    """
    os.makedirs(output_dir, exist_ok=True)
    html_files = []
    
    # Group probes by cluster_num
    probes_by_cluster = {}
    for probe_id, details in probe_details.items():
        cluster_num = details.get('cluster_num')
        if cluster_num is not None:
            if cluster_num not in probes_by_cluster:
                probes_by_cluster[cluster_num] = []
            probes_by_cluster[cluster_num].append((probe_id, details))
    
    # Generate visualization for each cluster
    def _get_cluster_num(key):
        if isinstance(key, int):
            return key
        try:
            return int(key.split('|')[0].replace('CONS', ''))
        except (ValueError, IndexError):
            return 0
    
    for cluster_num in sorted(degapped_consensus.keys(), key=_get_cluster_num):
        consensus_seq = degapped_consensus.get(cluster_num, '')
        cluster_probes = probes_by_cluster.get(cluster_num, [])
        
        # Get taxonomy from first accession
        taxonomy = {}
        for acc in cluster_accessions.get(cluster_num, []):
            if acc in metadata_map:
                taxonomy = metadata_map[acc]
                break
        
        if not taxonomy:
            taxonomy = {'family': 'Unknown', 'genus': 'Unknown', 'species': 'Unknown'}
        
        html_file = _generate_cluster_visualization(
            cluster_num,
            consensus_seq,
            cluster_probes,
            taxonomy,
            output_dir,
            width,
            height,
            adapter_seq
        )
        html_files.append(html_file)
    
    return html_files


def _generate_cluster_visualization(
    cluster_num: int,
    consensus_seq: str,
    cluster_probes: List[tuple],
    taxonomy: Dict[str, str],
    output_dir: str,
    width: int,
    height: int,
    adapter_seq: str
) -> str:
    """
    Generate single cluster visualization.
    
    Args:
        cluster_num: cluster number
        consensus_seq: consensus sequence
        cluster_probes: list of (probe_id, probe_details) tuples
        taxonomy: taxonomy info
        output_dir: output directory
        width: visualization width
        height: visualization height
        adapter_seq: adapter sequence
    
    Returns:
        Path to generated HTML file
    """
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>Cluster {cluster_num} - {taxonomy.get('family', 'Unknown')}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1 {{ color: #333; }}
        .info {{ background: #f0f0f0; padding: 15px; border-radius: 5px; margin-bottom: 20px; }}
        .visual {{ border: 1px solid #ccc; padding: 10px; margin: 10px 0; }}
        .probe {{ margin: 5px 0; padding: 5px; background: #e6f3ff; border-radius: 3px; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #4CAF50; color: white; }}
    </style>
</head>
<body>
    <h1>Cluster {cluster_num}</h1>
    
    <div class="info">
        <h2>Cluster Information</h2>
        <table>
            <tr><th>Property</th><th>Value</th></tr>
            <tr><td>Family</td><td>{taxonomy.get('family', 'Unknown')}</td></tr>
            <tr><td>Genus</td><td>{taxonomy.get('genus', 'Unknown')}</td></tr>
            <tr><td>Species</td><td>{taxonomy.get('species', 'Unknown')}</td></tr>
            <tr><td>Consensus Length</td><td>{len(consensus_seq)}</td></tr>
            <tr><td>Number of Probes</td><td>{len(cluster_probes)}</td></tr>
        </table>
    </div>
    
    <h2>Probe Mappings</h2>
    <table>
        <tr>
            <th>Probe ID</th>
            <th>Start</th>
            <th>End</th>
            <th>Length</th>
            <th>Sequence</th>
        </tr>
"""
    
    for probe_id, details in sorted(cluster_probes, key=lambda x: x[1].get('start', 0)):
        html_content += f"""        <tr>
            <td>{probe_id}</td>
            <td>{details.get('start', '-')}</td>
            <td>{details.get('end', '-')}</td>
            <td>{len(details.get('original_sequence', ''))}</td>
            <td style="font-family: monospace; font-size: 11px;">{details.get('original_sequence', '-')}</td>
        </tr>
"""
    
    html_content += """    </table>
    
    <div class="visual">
        <h2>Consensus Sequence Preview</h2>
        <p style="font-family: monospace; font-size: 10px; word-break: break-all;">
"""
    
    # Show first 500 bases of consensus
    preview_seq = consensus_seq[:500] if len(consensus_seq) > 500 else consensus_seq
    html_content += f"            {preview_seq}\n"
    
    if len(consensus_seq) > 500:
        html_content += f"            ... ({len(consensus_seq) - 500} more bases)\n"
    
    html_content += """        </p>
    </div>
</body>
</html>"""
    
    html_path = os.path.join(output_dir, f'cluster_{cluster_num}.html')
    with open(html_path, 'w') as f:
        f.write(html_content)
    
    return html_path
