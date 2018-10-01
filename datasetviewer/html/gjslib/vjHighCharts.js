
/////////////////////////////////////////////////////
function vjHighChartsSingleSeries(plotData1, title, divName, xTitle, yTitle, tickFreq, chartWidth) {
        var rows = [];
        var maxValue = 0;
	var cols = [];
        rows.push([null, 'Frequency']);
        for (var i in plotData1){
		if (plotData1[i]["y1"] <= 0 && divName == "chart_div2"){
			var a = [plotData1[i]["x"].trim(),0];
	        	rows.push(a);
        	        maxValue = Math.max(maxValue,plotData1[i]["y1"]);
		} else{
			var a = [plotData1[i]["x"].trim(),plotData1[i]["y1"]];
                        rows.push(a);
                        maxValue = Math.max(maxValue,plotData1[i]["y1"]);
		}
        }

        chart  = new Highcharts.Chart({

                colors: ["#7cb5ec", "#f7a35c", "#90ee7e", "#7798BF", "#aaeeee", "#ff0066", "#eeaaee",
                        "#55BF3B", "#DF5353", "#7798BF", "#aaeeee"],
                chart: {
                        renderTo: divName,
                        type: 'column'
                },
                title: {
                        text: title
                },
                yAxis: {
                        max: maxValue,
                        title: {
                                text: yTitle
                        }
                },
                tooltip: {
                        headerFormat: '<span>{point.key}</span><br>',
                        pointFormat: '<span style="color:{series.color};">{series.name}: ' +
                                '</span><span><b>{point.y}</b>',
                        footerFormat: '</span>',
                        shared: true,
                        useHTML: true
                },
                xAxis: {
                        tickPositioner: function() {
                                var result = [];
				if (divName == "chart_div1") {
                                	for(i = 0; i < rows.length; i++){
                                    		result.push(i);
					}
				} else {
					result = null;
				}
				return result;
                        },
                        title: {
                                text: xTitle
                        },labels: {
                                rotation:-50,
                                useHTML: true,
				style:{
					width:'70px',
					whiteSpace:'normal'//set to normal
				},
                                formatter: function () {
					var values = this.value;
					if (divName == "chart_div1") {
						//var cancer = values.split('(')[0];
						//var doid = values.split('(')[1].slice(0, -1);
						//var allHover = doid + '/ ' + cancer;
                                	        var allHover = values;
						return '<div class="hastip" title="' + allHover + '">' + values + '</div>';
					} else {
						return '<div class="hastip" title="' + values + '">' + values + '</div>';
					}
                                }
                        },
                },
                data: {
                        columns: rows,
                        switchRowsAndColumns: true

                }
});

}
