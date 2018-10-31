<%def name="composition(step, qc_dict)">
    <div class="composition row" id="${step}">
        <%
            total_key = 'total' if 'total' in qc_dict else 'all'
            total = qc_dict.pop(total_key)

            if step not in ["extract_PET"]:
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

        <div class="compositions">
            <h2>Reads compositions of each step's result</h2>
            ${reads_compositions(qc_contents['reads_counts'])}

            <h3>Tendency</h3>
        </div>

        <div class="chr_interactions">
            <h2>Interaction between chromosomes</h2>
            ${chr_interactions(qc_contents['chr_interactions'])}
        </div>

        </div>
    </body>

</html>