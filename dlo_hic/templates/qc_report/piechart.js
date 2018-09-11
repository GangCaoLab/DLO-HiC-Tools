var plotPieChart = function(total, dataset, stepName) {

    var w = 160,
        h = w

    var formatPercentage = function(portion) {
        return ((portion / total) * 100).toFixed(2) + '%'
    }

    var svg = d3.select("div.piechart#" + stepName)
                .append("svg")
                .attr("id", stepName)
                .attr("width", w)
                .attr("height", h)

    var pie = d3.pie()

    var outerRadius = w / 2
    var innerRadius = w / 5

    var arc = d3.arc()
        .innerRadius(innerRadius)
        .outerRadius(outerRadius)

    var arcs = svg.selectAll("svg#" + stepName + " g.arc")
                  .data(pie(dataset.map(
                      function(d){return +d.number})))
                  .enter()
                  .append("g")
                  .attr("class", "arc")
                  .attr("transform", "translate(" + outerRadius + ", " + outerRadius + ")")

    var colorPattle = ["#fcd6cf", "#6ec9f7", "#ffcc99", "#9efadb", "#9efafa"]
    var color = d3.scaleOrdinal(colorPattle)

    arcs.append("path")
        .attr("fill", function(d, i) {
            return color(i)
        })
        .attr("d", arc)

    arcs.append("text")
        .attr("transform", function(d) {
            return "translate(" + arc.centroid(d) + ")"
        })
        .attr("text-anchor", "middle")
        .text(function(d) {
            return formatPercentage(d.data)
        })
        .attr("fill", "white")

    arcs.append("title")
        .data(dataset)
        .text(function(d) {
            return d.item
        })

}