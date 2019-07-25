<%def name="composition(comp, id_)">
    <div class="composition row" id="${id_}">
        <%
            from copy import copy
            comp = copy(comp)
            if ('total' in comp) or ('all' in comp):
                k_ = 'total' if 'total' in comp else 'all'
                total = float(comp.pop(k_))
            else:
                total = sum([float(v) for k, v in comp.items()])

            if 'long-range' in comp:
                intra_long = int(comp.pop('long-range'))
                intra = int(comp.pop('intra-chromosome'))
                comp['intra-chromosome(long-range)'] = intra_long
                comp['intra-chromosome(short-range)'] = intra - intra_long
        %>

        ${composition_table(comp, total, id_, flow=True)}
        ${composition_piechart(comp, total, id_)}

    </div>
    <div class="total_reads">
        <p class="total-reads"> Total reads: ${int(total)} </p>
    </div>
</%def>

<%def name="composition_table(qc_dict, total, id_, flow=False)">
    % if flow:
    <div class="comp_table column left" id="${id_}">
    % else:
    <div class="comp_table" id="${id_}">
    % endif
        <table>
            <tr>
                <th></th>
                <th>Reads number</th>
                <th>Ratio</th>
            </tr>

            % for item, val in qc_dict.items():
                <%
                    ratio = 0 if total == 0 else float(val)/float(total)
                    ratio_str = "{:.2%}".format(ratio)
                %>
                <tr>
                    <td>${item}</td>
                    <td>${val}</td>
                    <td>${ratio_str}</td>
                </tr>
            % endfor

        </table>
    </div>
</%def>

<%def name="composition_piechart(qc_dict, total, id_)">
    <%
        # Convert Python dict to JS objects Array literal.
        dataset = [{"item": item, "number": n} for item, n in qc_dict.items()]
        dataset = str(dataset)
    %>
    <div class="piechart column right" id="${id_}">
        <script> 
            var total = ${total}
            var dataset = ${dataset}
            var pieID = "${id_}"

            plotPieChart(total, dataset, pieID)

        </script>
    </div>
</%def>


<%def name="chr_interactions(interactions)">
    <%
        # Convert Python dict to JS objects Array literal.
        dataset = str(interactions)
    %>
    <div class="chr_interactions heatmap">
        <script>
            var dataset = ${dataset}

            heatmap(dataset)

        </script>
    </div>
</%def>


<%def name="barchart(qc_dict, id_, x_numeric, width=400, height=300, color='#000000', ymax=None, xLabel='', yLabel='')">
    <%
        if x_numeric:
            dataset = [{'x': int(k), 'y': v} for k, v in qc_dict.items()]
        else:
            dataset = [{'x': k, 'y': v} for k, v in qc_dict.items()]
        dataset = sorted(dataset, key=lambda i:i['x'] )
        dataset = str(dataset)
        if ymax is None:
            ymax = 'null'
    %>

    <div class="barchart" id="${id_}">

    <script>
        var dataset = ${dataset}
        var barID = "${id_}"

        barchart(dataset, barID, ${width}, ${height}, "${color}", ${ymax}, "${xLabel}", "${yLabel}")
    </script>
</%def>


<%def name="reads_counts_table(qc_contents)">
    <%
        raw_reads    = int(qc_contents['extract_PET']['main']['flag_stat']['total'])
        linker_intra = int(qc_contents['extract_PET']['main']['flag_stat']['intra-molecular linker'])
        linker_inter = int(qc_contents['extract_PET']['main']['flag_stat']['inter-molecular linker'])
        if linker_inter == 0:
            linker_reads = linker_intra
            linker_inter = 'NA'
            linker_intra = 'NA'
        else:
            linker_reads = linker_inter + linker_intra
        unique_mapped_reads = int(qc_contents['build_bedpe']['main']['paired']['total'])
        non_redundant_reads = int(qc_contents['bedpe2pairs']['comp']['total'])
        inter_chr_reads = int(qc_contents['bedpe2pairs']['comp']['inter-chromosome'])
        intra_chr_reads = int(qc_contents['bedpe2pairs']['comp']['intra-chromosome'])
        counts_table = [
            ('Raw reads', raw_reads),
            ('Linker reads', linker_reads),
            ('Linkers (A-A) and (B-B)', linker_intra),
            ('Linkers (A-B) and (B-A)', linker_inter),
            ('Uniquely mapped reads', unique_mapped_reads),
            ('Non-redundant mapped reads', non_redundant_reads),
            ('Interchromosomal contacts', inter_chr_reads),
            ('Intrachromosomal contacts', intra_chr_reads),
        ]
    %>
    <div class="counts_table">
        <table>
            <tr>
                <th></th>
                <th>Reads number</th>
                <th>Keep ratio</th>
            </tr>
            % for key, count in counts_table:
                <% 
                    if not isinstance(count, str):
                        ratio = 0 if raw_reads == 0 else count / raw_reads
                        ratio = "{0:.2%}".format(ratio)
                    else:
                        ratio = 'NA'
                %>
                <tr>
                    <td>${key}</td>
                    <td>${count}</td>
                    <td>${ratio}</td>
                </tr>
            % endfor
        </table>
    </div>
</%def>

<%def name="pet_span_stats_table(stats_list)">
    <table>
        <tr>
            % for item, _ in stats_list:
                <th>${item}</th>
            % endfor
        </tr>

        <tr>
            % for _, val in stats_list:
                <%
                    if val - int(val) != 0:
                        val_str = "{:.2}".format(val)
                    else:
                        val_str = "{}".format(int(val))
                %>
                <td>${val_str}</td>
            % endfor
        </tr>
    </table>
</%def>


<%def name="include_svg(svg)">
    ## include svg
    ${svg}
</%def>


<%def name="pet_span(qc_contents_res)">
    <div class="pet_span">
        <div class="pet_span_stats_table">
            ${pet_span_stats_table(qc_contents_res['pet_span_stats'])}
        </div>
        <div class="pet_span_svg" style="text-align:center">
            ${include_svg(qc_contents_res['pet_span_fig'])}
        </div>
    </div>
</%def>


##
## Page structure
##


<html>
    <head>
        <style>
            <%include file="/qc_report.css" />
        </style>
        <script>
            <%include file="/d3.v4.min.js" />
        </script>
        <script>
            <%include file="/piechart.js" />
        </script>
        <script>
            <%include file="/heatmap.js" />
        </script>
        <script>
            <%include file="/barchart.js" />
        </script>

    </head>

## Report Header

    <body>

        <div class="header">
            <div>
            <h1>DLO-HiC-Tools</h1>

            <p class="version">v${version}</p>
            <p class="sample-id datetime">
                sample: ${sample_id}, create: ${datetime}
            </p>
            </div>
        </div>

        <div class="content">

## Overview counts table:

            <div class="counts_table">
                <h2>Reads counts</h2>
                ${reads_counts_table(qc_contents)}
            </div>

## QC of each steps:

            <div class="steps">
                <h2>Steps:</h2>

                <div class="PET_extract">
                    <h3>1. PET extract</h3>
    ## Step 1
                    <div class="PET_len_stat">
                        <h4>(1). PET length distribution</h4>
                        <%
                            pet1_len_dist = qc_contents['extract_PET']['main']['PET1_len_dist']
                            pet2_len_dist = qc_contents['extract_PET']['main']['PET2_len_dist']
                            len_dist = [int(i) for i in pet1_len_dist.values()] + [int(i) for i in pet2_len_dist.values()]
                            barmax = max(len_dist)
                        %>
                        <table><tr>
                            <td>
                                <div class="PET1_len_dist">
                                    <h5> PET1 </h5>
                                    ${barchart(pet1_len_dist, "PET1_len_dist", True, xLabel="PET length", yLabel="Count", ymax=barmax)}
                                </div>
                            </td>
                            <td>
                                <div class="PET2_len_dist">
                                    <h5> PET2 </h5>
                                    ${barchart(pet2_len_dist, "PET2_len_dist", True, xLabel="PET length", yLabel="Count", ymax=barmax)}
                                </div>
                            </td>
                        </tr></table>
                    </div>

                    <div class="flag_stat">
                        <h4>(2). Flag counts</h4>
                        <%
                            raw = int(qc_contents['extract_PET']['main']['flag_stat']['total'])
                            from copy import copy
                            flag_stat  = qc_contents['extract_PET']['main']['flag_stat']
                            flag_stat = copy(flag_stat)
                            total_flag = flag_stat.pop('total')
                        %>
                        ${composition_table(flag_stat, total_flag, "flag_stat")}
                        <p> Total (raw reads): ${raw} </p>
                    </div>

                    <div class="match_score_dist">
                        <h4>(3). Match score distribution</h4>
                        <%
                            linker_m_score = qc_contents['extract_PET']['main']['linker_match_score_dist']
                            adapter_m_score = qc_contents['extract_PET']['main']['adapter_match_score_dist']
                        %>
                        <div class="linker_match_score_dist">
                            <h5> Linker <h5>
                            ${barchart(linker_m_score, "linker_match_score_dist", True, width=800, height=350, xLabel="Match Score", yLabel="Count", color="#ff9d9d")}
                        </div>
                        % if adapter_m_score:
                        <div class="adapter_match_score_dist">
                            <h5> Adapter <h5>
                            ${barchart(adapter_m_score, "adapter_match_score_dist", True, width=800, height=350, xLabel="Match Score", yLabel="Count", color="#928287")}
                        </div>
                        % endif
                    </div>

                    <%
                        has_adapter_infer = 'adapter' in qc_contents['extract_PET']
                    %>
                    % if has_adapter_infer:
                    <div class="infer_adapter">
                        <h4>(4). Adapter inference</h4>
                        
                        <h5> Adapter decision graph </h5>
                        <div class="adapter_svg">
                            ${include_svg(qc_contents['extract_PET']['adapter_svg'])}
                        </div>

                        <h5> Adapter sequence </h5>
                        <div class="adapter">
                            <p>${qc_contents['extract_PET']['adapter']['adapter_seq']}</p>
                        </div>
                        
                    </div>
                    % endif

                </div>

    ## Step 2
                <div class="build_bedpe">
                    <h3>2. build BEDPE</h3>
                    <%
                        from copy import copy
                        comp = copy(qc_contents['build_bedpe']['main'])
                        total = int(qc_contents['extract_PET']['main']['flag_stat']['valid reads'])
                        pet1 = comp['PET1']
                        pet2 = comp['PET2']
                        paired = int(comp['paired']['total'])
                        p_r = "{:.2%}".format(paired/total) if total != 0 else 0
                    %>
                    <h4>BWA alignment</h4>
                    <table>
                        <tr>
                        <td>
                            <p>PET1:</p>
                            ${composition_table(pet1, total, 'pet1_align')}
                        </td>
                        <td>
                            <p>PET2:</p>
                            ${composition_table(pet2, total, 'pet2_align')}
                        </td>
                        </tr>
                    </table>
                    <p>Total reads(valid PET extracted): ${total}</p>
                    <p>Paired unique mapped reads: ${paired}&nbsp;&nbsp;&nbsp; Ratio: ${p_r}</p>


                </div>

    ## Step 3
                <div class="noise_reduce">
                    <h3>3. noise reduce</h3>
                    <h4>Ratio of self-ligation and re-ligation</h4>
                    ${composition(qc_contents['noise_reduce']['main']['type'], "ratio_noise")}

                    <h4>Count of reads orientation and position relative to fragment<h4>
                    <div class="position_table">
                    <table>
                        <%
                            ocs = ["++", "+-", "-+", "--"]  # orientation combinations
                            pcs = ["ss", "st", "ts", "tt"]  # position compositions
                        %>
                        <tr>
                            <th></th>
                            % for oc in ocs:
                            <th>${oc[0]+'/'+oc[1]}</th>
                            % endfor
                        </tr>

                        % for pc in pcs:
                        <tr>
                            <td>${pc[0]+'/'+pc[1]}</td>
                            % for oc in ocs:
                            <td>${qc_contents['noise_reduce']['main']['position'][oc+pc]}</td>
                            % endfor
                        </tr>
                        % endfor
                    </table>
                    <div>

                </div>

    ## Step 4
                <div class="bedpe2pairs">
                    <h3>4. remove duplication</h3>
                    <h4>Duplicate ratio</h4>
                    <%
                        raw = int(qc_contents['extract_PET']['main']['flag_stat']['total'])
                        final = int(qc_contents['bedpe2pairs']['comp']['total'])
                        nr = int(qc_contents['noise_reduce']['main']['type']['normal'])

                        dup = nr - final
                        d_r_r = "{:.2%}".format(dup/raw) if raw != 0 else 0 
                        d_n_r = "{:.2%}".format(dup/nr) if nr != 0 else 0

                    %>
                    <table>
                    <tr><th>duplication</th><th>duplication / raw</th><th>duplication / noise-reduced</th></tr>
                    <tr><td>${dup}</td><td>${d_r_r}</td><td>${d_n_r}</td></tr>
                    </table>

                </div>
            </div>

## Final results:

            <div class="final_results">
                <h2>Results:</h2>

                <div class="inter_intra_chr">
                    <h3>Ratio of intra and inter chromosome interactions</h3>
                    ${composition(qc_contents['bedpe2pairs']['comp'], "inter_intra_chr")}
                </div>

                <div class="pet_span_dist">
                    <h3>PET span distribution</h3>
                    ${pet_span(qc_contents['bedpe2pairs'])}
                </div>

                <div class="chr_interactions">
                    <h3>Interaction between chromosomes</h3>
                    ${chr_interactions(qc_contents['bedpe2pairs']['chr_interactions'])}
                </div>
            </div>

            <div class="run_time">
                <h2>Run times:</h2>
                <div class="run_time_table">
                <table>
                    <tr>
                        % for item in ["", "start", "end", "duration"] :
                            <th>${item}</th>
                        % endfor
                    </tr>

                    % for pname in qc_contents['other']['time_points'].keys():
                        <tr>
                            <td>${pname}</td>
                            % for key in ["start", "end", "diff"]:
                            <td>${qc_contents['other']['time_points'][pname][key]}</td>
                            % endfor
                        </tr>
                    % endfor
                </table>
                </div>
            </div>

            <div class="pipe_config">
                <h2>Pipeline configuration:</h2>
                <div>
                    <textarea rows="20" cols="98" readonly>${qc_contents['other']['pipe_config']}</textarea>
                </div>
            </div>

        </div>

        <div class="footer">
            <div>
            <p class="footer">
                Created by DLO-HiC-Tools, project page: <a href="https://github.com/GangCaoLab/DLO-HiC-Tools">GitHub</a>
            </p>
            </div>
        </div>

    </body>

</html>