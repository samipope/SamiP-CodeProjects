d3.csv('COD_ISOA3.csv').then(data => {
    // console.log("Data loaded:", data);

    const diseases = [
        'Meningitis', 
        "Alzheimer's Disease and Other Dementias", 
        "Parkinson's Disease",
        "Nutritional Deficiencies",
        'Malaria',
        'Drowning',
        "Interpersonal Violence",
        "Maternal Disorders",
        "HIV/AIDS",
        "Drug Use Disorders",
        'Tuberculosis',
        "Cardiovascular Diseases",
        "Lower Respiratory Infections",
        "Neonatal Disorders",
        "Alcohol Use Disorders",
        "Self-harm",
        "Exposure to Forces of Nature",
        "Diarrheal Diseases",
        "Environmental Heat and Cold Exposure",
        'Neoplasms',
        "Conflict and Terrorism",
        "Diabetes Mellitus",
        "Chronic Kidney Disease",
        'Poisonings',
        "Protein-Energy Malnutrition",
        "Road Injuries",
        "Chronic Respiratory Diseases",
        "Cirrhosis and Other Chronic Liver Diseases",
        "Digestive Diseases",
        "Fire, Heat, and Hot Substances",
        "Acute Hepatitis"
    ];

    // ------------- GRAPHING: FIRST LINE GRAPH --------------------------
   // prep data --> make empty yearDiseaseMap
   const yearDiseaseMap = {};

   data.forEach(d => {
       const year = d.Year;
       if (!yearDiseaseMap[year]) {
           yearDiseaseMap[year] = {};
           diseases.forEach(disease => {
               yearDiseaseMap[year][disease] = 0;
           });
       }

       diseases.forEach(disease => {
           yearDiseaseMap[year][disease] += +d[disease];
       });
   });
//    console.log("Year-disease map:", yearDiseaseMap);

   const processedData = [];
   for (const year in yearDiseaseMap) {
       const entry = { Year: +year };
       diseases.forEach(disease => {
           entry[disease] = yearDiseaseMap[year][disease];
       });
       processedData.push(entry);
   }

    //console.log("Processed data:", processedData);
    // create chart with unscaled data
    // console.log("Creating unscaled chart...");
   createUnscaledChart(processedData, diseases);

// ------------ GRAPHING: SECOND (SCALED) LINE GRAPH ------------

    // normalize the data relative to the first year in the dataset
    // console.log("Normalizing data for scaled chart...");
   const firstYear = processedData[0].Year;
   const series = diseases.map(disease => {
       const values = processedData.map(d => ({ Year: d.Year, value: d[disease] }));
       const v0 = values[0].value;
       return { name: disease, values: values.map(d => ({ Year: d.Year, value: d.value / v0 })) };
   });

//    console.log("Series for scaled chart:", series);
//    console.log("Creating scaled chart...");
   createScaledChart(series);



// ------------- GRAPHING: ZOOMABLE BAR CHART -------------

//    console.log("Filtering data for 'Cardiovascular Diseases' for zoomable bar chart...");
   const heartDiseaseData = data.map(d => ({
       Year: +d.Year,
       value: +d["Cardiovascular Diseases"]
   }));

//    console.log("Filtered data for zoomable chart:", heartDiseaseData);
//    console.log("Creating zoomable bar chart...");
   createZoomableBarChart(heartDiseaseData);

   // make dropdown selection
   const dropdown = d3.select("#disease")
   .selectAll("option")
   .data(diseases)
   .enter()
   .append("option")
   .text(d => d)
   .attr("value", d => d);

   // function to UPDATE map based on disease chosen by user
  function updateChoropleth(selectedDisease) {
    const choroplethData = data
        .filter(d => d.Year === "2019")
        .reduce((acc, curr) => {
            acc[curr["alpha-3"]] = +curr[selectedDisease];
            return acc;
        }, {});
 
    // clear the previous map
    d3.select("#world-map").select("svg").remove();
 
    // create the updated choropleth map
    createWorldChoropleth(
        "#world-map",
        'https://raw.githubusercontent.com/holtzy/D3-graph-gallery/master/DATA/world.geojson',
        choroplethData,
        selectedDisease,
        d3.interpolateReds
    );
 }


// function to update the choropleth map based on selected disease
updateChoropleth(diseases[0]);

// update the map when the dropdown selection changes
d3.select("#disease").on("change", function() {
   const selectedDisease = d3.select(this).property("value");
   updateChoropleth(selectedDisease);
});

}).catch(error => {
console.error("Error loading or processing data:", error);
}); // end of running like end of main
  

// -------------------------- BEGIN FUNCTIONS FOR GRAPHS !! --------------
  

  // function to create a world choropleth map
function createWorldChoropleth(mapContainerId, geoDataUrl, dataMap, dataKey, colorScheme) {
    const width = 960;
    const height = 600;

    // create SVG
    const svg = d3.select(mapContainerId)
        .append("svg")
        .attr("width", width)
        .attr("height", height);

    // color scale
    const colorScale = d3.scaleSequential(colorScheme)
        .domain([0, d3.max(Object.values(dataMap))]);

    // projection and path
    const projection = d3.geoNaturalEarth1()
        .scale(160)
        .translate([width / 2, height / 2]);

    const path = d3.geoPath().projection(projection);

    // geoJSON data and render map
    d3.json(geoDataUrl).then(geoData => {
        // console.log("GeoJSON Data:", geoData);

        svg.append("g")
            .selectAll("path")
            .data(geoData.features)
            .enter()
            .append("path")
            .attr("d", path)
            .attr("fill", d => {
                const country = d.id;  // matching so map knows country
                const value = dataMap[country];
                // console.log(`Country: ${country}, Value: ${value}`);
                return value ? colorScale(value) : "#ccc";
            })
            .style("stroke", "black")
            .style("opacity", 0.8)
            .on("mouseover", function (event, d) {
                const country = d.properties.name;
                const deaths = dataMap[d.id] || "No data";
                tooltip.style("visibility", "visible")
                    .html(`Country: ${country}<br>Deaths: ${deaths}`);
            })
            .on("mousemove", function (event) {
                tooltip.style("top", (event.pageY - 10) + "px")
                    .style("left", (event.pageX + 10) + "px");
            })
            .on("mouseout", function () {
                tooltip.style("visibility", "hidden");
            });

        const tooltip = d3.select("body").append("div")
            .attr("class", "tooltip")
            .style("position", "absolute")
            .style("visibility", "hidden");
    }).catch(error => {
        console.error("Error loading GeoJSON data:", error);
    });
}

// function to create the unscaled chart
function createUnscaledChart(data, diseases) {
    // console.log("In createUnscaledChart function");
    const width = 1200;
    const height = 800;
    const marginTop = 20;
    const marginRight = 40;
    const marginBottom = 50;
    const marginLeft = 60;

    const x = d3.scaleLinear()
        .domain(d3.extent(data, d => d.Year))
        .range([marginLeft, width - marginRight])
        .clamp(true);

    const y = d3.scaleLinear()
        .domain([0, d3.max(data, d => d3.max(diseases, disease => d[disease]))])
        .range([height - marginBottom, marginTop]);

    const color = d3.scaleOrdinal(d3.schemeCategory10)
        .domain(diseases);

    const svg = d3.select("#unscaled-chart").append("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("viewBox", [0, 0, width, height])
        .attr("style", "max-width: 100%; height: auto; -webkit-tap-highlight-color: transparent;");

    svg.append("g")
        .attr("transform", `translate(0,${height - marginBottom})`)
        .call(d3.axisBottom(x).ticks(width / 80).tickFormat(d3.format("d")).tickSizeOuter(0))
        .call(g => g.select(".domain").remove());

    svg.append("g")
        .attr("transform", `translate(${marginLeft},0)`)
        .call(d3.axisLeft(y))
        .call(g => g.select(".domain").remove());

    svg.append("text")
        .attr("class", "x-axis-label")
        .attr("text-anchor", "end")
        .attr("x", width / 2)
        .attr("y", height - 10)
        .text("Year");

    svg.append("text")
        .attr("class", "y-axis-label")
        .attr("text-anchor", "end")
        .attr("x", -height / 2)
        .attr("y", 20)
        .attr("transform", "rotate(-90)")
        .text("Number of Deaths");

    const line = d3.line()
        .x(d => x(d.Year))
        .y(d => y(d.value));

    const tooltip = d3.select("body").append("div")
        .attr("class", "tooltip")
        .style("position", "absolute")
        .style("visibility", "hidden")
        .style("background", "white")
        .style("border", "solid 1px black")
        .style("padding", "5px")
        .style("border-radius", "5px");

    const serie = svg.append("g")
        .style("font", "bold 10px sans-serif")
        .selectAll("g")
        .data(diseases.map(disease => {
            return {
                name: disease,
                values: data.map(d => ({ Year: d.Year, value: d[disease] }))
            };
        }))
        .join("g");

    serie.append("path")
        .attr("fill", "none")
        .attr("stroke-width", 1.5)
        .attr("stroke-linejoin", "round")
        .attr("stroke-linecap", "round")
        .attr("stroke", d => color(d.name))
        .attr("d", d => line(d.values))
        .on("mouseover", function(event, d) {
            tooltip.style("visibility", "visible")
                .text(d.name);
            d3.select(this).attr("stroke-width", 3);
        })
        .on("mousemove", function(event, d) {
            tooltip.style("top", (event.pageY - 10) + "px")
                .style("left", (event.pageX + 10) + "px")
                .html(`Disease: ${d.name}`);
        })
        .on("mouseout", function() {
            tooltip.style("visibility", "hidden");
            d3.select(this).attr("stroke-width", 1.5);
        });

    serie.append("text")
        .datum(d => ({ name: d.name, value: d.values[d.values.length - 1].value }))
        .attr("fill", d => color(d.name))
        .attr("paint-order", "stroke")
        .attr("stroke", "white")
        .attr("stroke-width", 3)
        .attr("x", width - marginRight + 3)
        .attr("y", d => y(d.value))
        .attr("dy", "0.35em")
        .text(d => d.name);
}

// function to create the scaled chart
function createScaledChart(series) {
    console.log("In createScaledChart function");
    const width = 1200;
    const height = 800;
    const marginTop = 20;
    const marginRight = 40;
    const marginBottom = 50;
    const marginLeft = 60;

    const x = d3.scaleLinear()
        .domain(d3.extent(series[0].values, d => d.Year))
        .range([marginLeft, width - marginRight])
        .clamp(true);

    const k = d3.max(series, ({values}) => d3.max(values, d => d.value) / d3.min(values, d => d.value));
    const y = d3.scaleLog()
        .domain([1 / k, k])
        .rangeRound([height - marginBottom, marginTop]);

    const color = d3.scaleOrdinal(d3.schemeCategory10)
        .domain(series.map(d => d.name));

    const svg = d3.select("#scaled-chart").append("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("viewBox", [0, 0, width, height])
        .attr("style", "max-width: 100%; height: auto; -webkit-tap-highlight-color: transparent;");

    svg.append("g")
        .attr("transform", `translate(0,${height - marginBottom})`)
        .call(d3.axisBottom(x).ticks(width / 80).tickFormat(d3.format("d")).tickSizeOuter(0))
        .call(g => g.select(".domain").remove());

    svg.append("g")
        .attr("transform", `translate(${marginLeft},0)`)
        .call(d3.axisLeft(y).ticks(null, x => +x.toFixed(6) + "×"))
        .call(g => g.selectAll(".tick line").clone()
            .attr("stroke-opacity", d => d === 1 ? null : 0.2)
            .attr("x2", width - marginLeft - marginRight))
        .call(g => g.select(".domain").remove());

    svg.append("text")
        .attr("class", "x-axis-label")
        .attr("text-anchor", "end")
        .attr("x", width / 2)
        .attr("y", height - 10)
        .text("Year");

    svg.append("text")
        .attr("class", "y-axis-label")
        .attr("text-anchor", "end")
        .attr("x", -height / 2)
        .attr("y", 20)
        .attr("transform", "rotate(-90)")
        .text("Scaled Deaths Relative to 1990");

    const line = d3.line()
        .x(d => x(d.Year))
        .y(d => y(d.value));

    const tooltip = d3.select("body").append("div")
        .attr("class", "tooltip")
        .style("position", "absolute")
        .style("visibility", "hidden")
        .style("background", "white")
        .style("border", "solid 1px black")
        .style("padding", "5px")
        .style("border-radius", "5px");

    const serie = svg.append("g")
        .style("font", "bold 10px sans-serif")
        .selectAll("g")
        .data(series)
        .join("g");

    serie.append("path")
        .attr("fill", "none")
        .attr("stroke-width", 1.5)
        .attr("stroke-linejoin", "round")
        .attr("stroke-linecap", "round")
        .attr("stroke", d => color(d.name))
        .attr("d", d => line(d.values))
        .on("mouseover", function(event, d) {
            tooltip.style("visibility", "visible")
                .text(d.name);
            d3.select(this).attr("stroke-width", 3);
        })
        .on("mousemove", function(event, d) {
            const [year] = d3.pointer(event, this);
            const yearValue = Math.round(x.invert(year));
            tooltip.style("top", (event.pageY - 10) + "px")
                .style("left", (event.pageX + 10) + "px")
                .html(`Disease: ${d.name}<br>Year: ${yearValue}`);
        })
        .on("mouseout", function() {
            tooltip.style("visibility", "hidden");
            d3.select(this).attr("stroke-width", 1.5);
        });

    serie.append("text")
        .datum(d => ({ name: d.name, value: d.values[d.values.length - 1].value }))
        .attr("fill", d => color(d.name))
        .attr("paint-order", "stroke")
        .attr("stroke", "white")
        .attr("stroke-width", 3)
        .attr("x", width - marginRight + 3)
        .attr("y", d => y(d.value))
        .attr("dy", "0.35em")
        .text(d => d.name);
}

// function to create the zoomable bar chart
function createZoomableBarChart(data) {
    // console.log("Creating zoomable bar chart with data:", data);

    const width = 1200;
    const height = 800;
    const marginTop = 20;
    const marginRight = 40;
    const marginBottom = 50;
    const marginLeft = 60;

    const x = d3.scaleBand()
        .domain(data.map(d => d.Year))
        .range([marginLeft, width - marginRight])
        .padding(0.1);

    const y = d3.scaleLinear()
        .domain([0, d3.max(data, d => d.value)]).nice()
        .range([height - marginBottom, marginTop]);

    const xAxis = d3.axisBottom(x).tickSizeOuter(0);
    const yAxis = d3.axisLeft(y);

    const svg = d3.select("#zoomable-bar-chart").append("svg")
        .attr("width", width)
        .attr("height", height)
        .attr("viewBox", [0, 0, width, height])
        .attr("style", "max-width: 100%; height: auto;");

    const bars = svg.append("g")
        .attr("class", "bars")
        .attr("fill", "steelblue")
        .selectAll("rect")
        .data(data)
        .join("rect")
        .attr("x", d => x(d.Year))
        .attr("y", d => y(d.value))
        .attr("height", d => y(0) - y(d.value))
        .attr("width", x.bandwidth());

    svg.append("g")
        .attr("class", "x-axis")
        .attr("transform", `translate(0,${height - marginBottom})`)
        .call(xAxis);

    svg.append("g")
        .attr("class", "y-axis")
        .attr("transform", `translate(${marginLeft},0)`)
        .call(yAxis)
        .call(g => g.select(".domain").remove());

    // zoom behavior
    const zoom = d3.zoom()
        .scaleExtent([1, 8])
        .translateExtent([[marginLeft, marginTop], [width - marginRight, height - marginTop]])
        .extent([[marginLeft, marginTop], [width - marginRight, height - marginTop]])
        .on("zoom", zoomed);

    svg.call(zoom);

    function zoomed(event) {
        x.range([marginLeft, width - marginRight].map(d => event.transform.applyX(d)));
        svg.selectAll(".bars rect").attr("x", d => x(d.Year)).attr("width", x.bandwidth());
        svg.selectAll(".x-axis").call(xAxis);
    }

    // console.log("Zoomable bar chart function executed");
}
