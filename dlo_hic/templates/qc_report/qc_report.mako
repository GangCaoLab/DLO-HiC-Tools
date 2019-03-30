<%def name="composition(comp, id_)">
    <div class="composition row" id="${id_}">
        <%
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
    <div class="comp_table left" id="${id_}">
    % endif
        <table>
            <tr>
                <th></th>
                <th>Reads number</th>
                <th>Ratio</th>
            </tr>

            % for item, val in qc_dict.items():
                <%
                    ratio = float(val)/float(total)
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
        raw_reads    = int(qc_contents['extract_PET']['main']['flag_stat']['all'])
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
                        ratio = count / raw_reads
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
        <div class="pet_span_svg">
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

            <p class="sample-id">Sample: ${sample_id}</p>
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
                            flag_stat  = qc_contents['extract_PET']['main']['flag_stat']
                            total_flag = flag_stat.pop('all')
                        %>
                        ${composition_table(flag_stat, total_flag, "flag_stat")}
                    </div>
                </div>

                <div class="build_bedpe">
                    <h3>2. build BEDPE</h3>

                </div>

                <div class="noise_reduce">
                    <h3>3. noise reduce</h3>

                </div>

                <div class="bedpe2pairs">
                    <h3>4. remove duplication</h3>

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

        </div>

    </body>

</html>