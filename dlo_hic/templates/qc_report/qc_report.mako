<%def name="composition(step, qc_dict)">
    <div class="composition row" id="${step}">
        <%
            total_key = 'total' if 'total' in qc_dict else 'all'
            total = qc_dict.pop(total_key)

            if step == "extract_PET":
                intra = int(qc_dict['intra-molecular linker'])
                inter = int(qc_dict['inter-molecular linker'])
                qc_dict = {
                    "intra-molecular linker": intra,
                    "inter-molecular linker": inter,
                }
                total = inter + intra
            else:
                intra = int(qc_dict.pop("intra-chromosome"))
                long_range = int(qc_dict.pop("long-range"))
                qc_dict["intra-chromosome (long-range)"] = long_range
                qc_dict['intra-chromosome (short-range)'] = (intra - long_range)
        %>

        ${composition_table(qc_dict, total, step)}
        ${composition_piechart(qc_dict, total, step)}

    </div>
    <div class="total_reads">
        <p class="total-reads"> Total reads: ${total} </p>
    </div>
</%def>

<%def name="composition_table(qc_dict, total, step)">
    <div class="qc_table column left" id="${step}">
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

<%def name="composition_piechart(qc_dict, total, step)">
    <%
        # Convert Python dict to JS objects Array literal.
        dataset = [{"item": item, "number": n} for item, n in qc_dict.items()]
        dataset = str(dataset)
    %>
    <div class="piechart column right" id="${step}">
        <script> 
            var total = ${total}
            var dataset = ${dataset}
            var stepName = "${step}"

            plotPieChart(total, dataset, stepName)

        </script>
    </div>
</%def>

<%def name="reads_compositions(qc_contents)">
    <div class="reads_compositions">
        % for step, qc_dict in qc_contents.items():
        <h3>Step: ${step}</h3>
        <% step = step.replace(".", "_") %>
            <div class="reads_compositions">
                ${composition(step, qc_dict)}
            </div>
        % endfor
    </div>
</%def>

<%def name="chr_interactions(qc_contents)">
    <%
        # Convert Python dict to JS objects Array literal.
        dataset = str(qc_contents)
    %>
    <div class="chr_interactions heatmap">
        <script>
            var dataset = ${dataset}
            <%include file="/heatmap.js" />
        </script>
    </div>
</%def>

<%def name="reads_counts_table(qc_contents)">
    <%
        raw_reads    = int(qc_contents['extract_PET']['all'])
        linker_intra = int(qc_contents['extract_PET']['intra-molecular linker'])
        linker_inter = int(qc_contents['extract_PET']['inter-molecular linker'])
        if linker_inter == 0:
            linker_reads = linker_intra
            linker_inter = 'NA'
            linker_intra = 'NA'
        else:
            linker_reads = linker_inter + linker_intra
        unique_mapped_reads = int(qc_contents['build_bedpe']['total'])
        non_redundant_reads = int(qc_contents['bedpe2pairs']['total'])
        inter_chr_reads = int(qc_contents['bedpe2pairs']['inter-chromosome'])
        intra_chr_reads = int(qc_contents['bedpe2pairs']['intra-chromosome'])
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
    ${svg}
</%def>

<%def name="pet_span(qc_contents)">
    <div class="pet_span">
        <div class="pet_span_stats_table">
            ${pet_span_stats_table(qc_contents['stats'])}
        </div>
        <div class="pet_span_svg">
            ${include_svg(qc_contents['svg'])}
        </div>
    </div>
</%def>

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

    </head>

    <body>
        <div class="header">
            <div>
            <h1>DLO-HiC-Tools</h1>

            <p class="sample-id">Sample: ${sample_id}</p>
            </div>
        </div>
        <div class="content">

        <div class="counts_table">
            <h2>Reads counts</h2>
            ${reads_counts_table(qc_contents['reads_counts'])}
        </div>

        <div class="pet_span">
            <h2>PET span distribution</h2>
            ${pet_span(qc_contents['pet_span'])}
        </div>

        <div class="compositions">
            <h2>Reads compositions of each step's result</h2>
            ${reads_compositions(qc_contents['reads_counts'])}
        </div>

        <div class="chr_interactions">
            <h2>Interaction between chromosomes</h2>
            ${chr_interactions(qc_contents['chr_interactions'])}
        </div>

        </div>
    </body>

</html>