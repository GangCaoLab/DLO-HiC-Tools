function barchart(data, chartID, width=400, height=300,
                  color="#000000", ymax=null,
                  xLabel="x", yLabel="y") {

  // set the dimensions and margins of the graph
  var margin = {top: 20, right: 20, bottom: 40, left: 70}

  width = width - margin.left - margin.right
  height = height - margin.top - margin.bottom;

  // set the ranges
  var x = d3.scaleBand()
            .range([0, width])
            .padding(0.1);
  var y = d3.scaleLinear()
            .range([height, 0]);

  // append the svg object to the body of the page
  // append a 'group' element to 'svg'
  // moves the 'group' element to the top left margin
  var svg = d3.select("div.barchart#"+chartID).append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
      .append("g")
      .attr("transform", 
            "translate(" + margin.left + "," + margin.top + ")");


  // format the data
  data.forEach(function(d) {
    d.y = +d.y;
  });

  // Scale the range of the data in the domains
  x.domain(data.map(function(d) { return d.x; }));
  if (ymax) {
    y.domain([0, ymax]);
  } else {
    ymax = d3.max(data, function(d) { return d.y; });
    y.domain([0, ymax]);
  }
  var xmax = d3.max(data, function(d) { return d.x });

  // append the rectangles for the bar chart
  var rects = svg.selectAll(".bar")
      .data(data)
      .enter().append("rect")
      .attr("class", "bar")
      .attr("fill", color)
      .attr("x", function(d) { return x(d.x); })
      .attr("width", x.bandwidth())
      .attr("y", function(d) { return y(d.y); })
      .attr("height", function(d) { return height - y(d.y); });

  // add the x Axis
  svg.append("g")
      .attr("transform", "translate(0," + height + ")")
      .call(d3.axisBottom(x));

  // add the y Axis
  svg.append("g")
      .call(d3.axisLeft(y));

  svg.append('text')
    .attr('x', -y(ymax/2))
    .attr('y', -margin.left/1.2)
    .attr('transform', 'rotate(-90)')
    .attr('text-anchor', 'middle')
    .text(yLabel)

  svg.append('text')
    .attr('x', width / 2)
    .attr('y', height + margin.bottom/1.2)
    .attr('text-anchor', 'middle')
    .text(xLabel)

  rects.append("title")
       .data(data)
       .text(function(d) {
         return d.y
       });

  rects.on('mouseenter', function (actual, i) {
        d3.select(this).attr('opacity', 0.5)
    })
    .on('mouseleave', function (actual, i) {
        d3.select(this).attr('opacity', 1)
    })


}  