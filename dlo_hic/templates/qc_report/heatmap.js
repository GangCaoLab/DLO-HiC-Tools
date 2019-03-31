function heatmap(dataset) {

    var 
        margin = { top: 0, right: 0, bottom: 50, left: 0 },
        web_width = parseInt(d3.select(".content").style("width")),
        w = web_width - margin.left - margin.right,
        h = web_width / 2
        dim = { w: web_width, h: h + margin.top + margin.bottom };

    var
        n = dataset.chromosomes.length,
        gridSpan = w / n,
        gridWidth = Math.sqrt(2) * (gridSpan / 2);

    var 
        data = dataset.data.filter(function(d){ if (d.pos[0] != d.pos[1]) {return d} }),
        diagonal_data = dataset.data.filter(function(d){ if (d.pos[0] == d.pos[1]) {return d} });

    var
        colorPattle = ["#ffffd9","#edf8b1","#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494","#081d58"],
        colorScale = d3.scaleQuantile()
            .domain([0, d3.max(data, function (d) { return d.value })])
            .range(colorPattle);


    var svg = d3.select("div.heatmap")
                .append("svg")
                .attr("width", dim.w)
                .attr("height", dim.h)
                .append("g")
                .attr("transform", "translate(" + margin.left + "," + margin.top + ")")


    var scale = function(posX, posY) {
        return (Math.abs(posY - posX) + 1) * (gridSpan / 2)
    }

    var xScale = function(posX, posY) {
        return gridSpan * posX + scale(posX, posY)
    }

    var yScale = function(posX, posY) {
        return h - scale(posX, posY)
    }


    var grids = svg.selectAll("rect")
        .data(data)
        .enter()
        .append("rect")
        .attr("x", function(d) {
            return xScale(d.pos[0], d.pos[1])
        })
        .attr("y", function(d) {
            return yScale(d.pos[0], d.pos[1])
        })
        .attr("width", gridWidth)
        .attr("height", gridWidth)
        .attr("transform", function(d) {
            var x = xScale(d.pos[0], d.pos[1])
            var y = yScale(d.pos[0], d.pos[1])
            return "rotate(45," + x + "," + y + ")"
        })
        .style("fill", function(d) {
            return colorScale(d.value)
        })
        .on('mouseenter', function (actual, i) {
            d3.select(this).attr('opacity', 0.5)
        })
        .on('mouseleave', function (actual, i) {
            d3.select(this).attr('opacity', 1)
        })
        .append("title")
        .text(function(d) {
            return d.chr[0] + " - " + d.chr[1] + ": " + d.value
        })



    var legendMargin = {top: 10, bottom: 10, left: 20, right: 20},
        legendElementWidth = ((w - legendMargin.left - legendMargin.right) / (colorScale.quantiles().length + 1)),
        legendHeight = 20,
        legendData = [0].concat(colorScale.quantiles())

    var legend = svg.selectAll(".legend")
        .data(legendData)

    legend.enter().append("g")
        .attr("class", "legend")
        .attr("transform", "translate(" + legendMargin.left + ",0)")
        .append("rect")
        .attr("x", function(d, i) {
            return legendElementWidth * i
        })
        .attr("y", h + legendMargin.top)
        .attr("width", legendElementWidth)
        .attr("height", legendHeight)
        .style("fill", function(d, i) {return colorPattle[i]; })

    legend.enter().append("text")
        .attr("class", "mono")
        .attr("transform", "translate(" + legendMargin.left + ",0)")
        .text(function(d) {return "\u2265 " + Math.round(d);})
        .attr("x", function(d, i) {return legendElementWidth * i})
        .attr("y", h + legendElementWidth/2)

}    