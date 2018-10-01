var libFuncs = {
  	
//////////////////////////////////////
rndrHeaderRowOne: function (configObj, containerId){

	var cn = '';
	for (var i=0; i < configObj.length; i++){
		cn += '<img src="' + configObj[i]["href"] + '" class="' + configObj[i]["class"] + '">';
	}
	$(containerId).html(cn);
	return;
},

//////////////////////////////////////
rndrHeaderRowTwo: function (configObj, containerId){

	var cn = '<img src="' + configObj[0]["href"] + '" width=60>';
	cn += '&nbsp; <u>H</u>igh performance <u>I</u>ntegrated <u>V</u>irtual <u>E</u>nvironment';
	$(containerId).html(cn);
	return;

},

//////////////////////////////////////
rndrHeaderRowThree: function (configObj, containerId){

	var cn = '<table width=100% height=100% class=nav cellspacing=0 border=0 cellpadding=0>';
	cn += '<tr><td style="padding:0 0 5 40;" valign=bottom>';
	cn += '<table cellpadding=0 cellspacing=0 border=0><tr>';

	for (var i=0; i < configObj.length-1; i++){
		cn += '<td class=gsection><a href="' + configObj[i]["href"];
		cn += '" class=' + configObj[i]["class"] + '> ' + configObj[i]["label"]  + '</a></td>';
	}
	cn += '</tr></table></td>';
	cn += '<td align=right class=glogin valign=bottom>';
//How to deal with this part -- for login
	cn += '<a href="' + configObj[configObj.length-1]["href"] + '" class=';
	cn += configObj[configObj.length-1]["class"];
	cn += '>' + configObj[configObj.length-1]["label"] + '</a></td></tr></table>';

	$(containerId).html(cn);
	return;
},


//////////////////////////////////////
rndrHeaderRowFour: function (configObj, containerId){

	var cn = '';
	for (var i=0; i < configObj.length; i++){
		cn += '<a href="' + configObj[i]["href"] + '" class=';
		cn += configObj[i]["class"] + '>' + configObj[i]["label"] + '</a> >> ';
	}
	$(containerId).html(cn);
	return;
},


//////////////////////////////////////
rndrHeaderRowFive: function (configObj, containerId){
//can I put all class right here or should they pass through configObj?
	var cn = '<DIV style="left:0px;top:0px;width:100%;height:100%;">';
	cn += '<DIV class=moduleicon1></DIV>';
	cn += '<DIV class=moduleicon2></DIV>';
	cn += '<DIV class=moduleicon3></DIV>';
	cn += '</DIV>';

	$(containerId).html(cn);
	return;

},

//////////////////////////////////////
rndrHeaderRowSix: function (configObj, containerId){

	var configText = configObj["title"];
	var cn = '<DIV class=modulename>';
	cn += '<DIV class=logotext1>' + configText + '</DIV></DIV>';

	$(containerId).html(cn);
	return;
},


//////////////////////////////////////
rndrFooter: function (configObj, containerId){

	var cn = '<table border=0 cellpadding=0 cellspacing=0>';
	cn += '<tr><td width=50></td>';

	for (var i=0; i < configObj.length; i++){
		cn += '<td class=rightbordered width=100>&nbsp; <a class="' + configObj[i]["class"];
		cn += '" href="' + configObj[i]["href"] + '">' + configObj[i]["label"] + '</a></td>';
	}
	cn += '</table>';

	$(containerId).html(cn);
	return;

},

//////////////////////////////////////
rndrModuleSections: function (configObj, containerId){

	var cn = '<table border=0>';
	for (var i=0; i < configObj.length; i++){
		var cls = (configObj[i]["id"] == SELECTEDMODULEID ? configObj[i]["classes"][1] : 
				configObj[i]["classes"][0]);
		var link = '<a href="' + configObj[i]["href"];
		link += '" class="' + cls + '" id="'+configObj[i]["id"] + '">';
		link += configObj[i]["label"] + '</a>';
		cn += '<tr height=20><td>' + link + '</td></tr>';
	}
	cn += '</table>';

	$(containerId).html(cn);
	return;
},


//////////////////////////////////////
fetchConfigJson: function (svcPath, svcName, jsonPath, configJson, funcInputList){

	var url = svcPath + svcName + '.py'
	alert(url);

	var xhr = new XMLHttpRequest();
	xhr.open("POST", url, true);
	xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
	xhr.onreadystatechange = function() {
		if (xhr.readyState === 4 && xhr.status == "200") {
			configJson = JSON.parse(xhr.responseText);
			for (var i=0; i < funcInputList.length; i++){
				libFuncs[funcInputList[i]["funcName"]](configJson[funcInputList[i]["configSection"]], funcInputList[i]["containerId"]);

			}
			libFuncs['getWaitMsg'](configJson["loadingimage"], waitMsg);
		}
	}

	var postData = 'injsonpath=' + jsonPath;
	xhr.send(postData);
	console.log('request='+postData);
	return;

},

//////////////////////////////////////
rndrErrorMsg: function (msg, containerId){

	$(containerId).html('<table width=90% height=200><tr><td valign=middle align=center>' +
			msg + '</td></tr></table>');
        return;
},

/////////////////////////////
jsonToHtmlTable: function (jsonObj, containerId, styleJsonObj){
	
	if (jsonObj["output"].length == 1){
		var table = '<table width=90% height=400><tr><td valign=middle align=center>' +
			'Sorry, your search did not return any result.' +
			'</td></tr></table>';
		$(containerId).html(table);
		return;
	}

	var table = '<table style="' + styleJsonObj["tableStyleWhole"] + '">';
	table += '<tr class=evenrow>';
	for (var j=0; j < jsonObj["output"][0].length; j++){
		table += '<td align=center nowrap style="' + styleJsonObj["headerRowStyle"] + '">' + jsonObj["output"][0][j] + '</td>';
	}
	table += '</tr>';
	for (var i=1; i < jsonObj["output"].length; i++){
		table += '<tr>';
		for (var j=0; j < jsonObj["output"][i].length; j++){
			table += '<td nowrap style="' + styleJsonObj["rowStyle"] + '">' + jsonObj["output"][i][j] + '</td>';
		}
		table += '</tr>';
	}
	table += '</table>';
	$(containerId).html(table);

	return;

},

/////////////////////////////
styleHash: function (tableStyleIndex, figureStyleIndex){

	var tableStyleList = [
		{"tableStyle":"position:absolute;left:0px;top:0px;width:100%;height:530px;overflow:auto;",
		"navStyle1":"position:absolute;left:0px;top:-20px;width:100%;height:15px;border-bottom:1px solid #ccc;",
		"navStyle2":"position:absolute;left:0px;bottom:-20px;width:100%;height:17px;border-top:1px solid #ccc;",
		"tableStyleWhole":"font-size:13px;width=90% border=0",
		"headerRowStyle":"padding:0 10 0 10;",
		"rowStyle":"padding:0 10 0 10;border:1px solid #eee;"},
		{"tableStyle":"", "headerRowStyle":"", "rowStyle":""},
		{"tableStyle":"", "headerRowStyle":"", "rowStyle":""},
		{"tableStyle":"", "headerRowStyle":"", "rowStyle":""}
	];

	var figureStyleList = [
		{"plotStyle":"position:absolute;left:0px;top:0px;width:100%;height:370px;border:1px solid #ccc;"},
		{"plotStyle":"position:absolute;left:0px;top:0px;width:100%;height:370px;border:1px solid #ccc;"}
	];

	var styleJsonObj = {"tableStyleList":tableStyleList[tableStyleIndex], "figureStyleList":figureStyleList[figureStyleIndex]};

	return styleJsonObj;

},


//////////////////////////////////////////
getJsonForSubmittedValues: function (className){

	var jqClass = "."+className;
	var emList = $(jqClass);
	var myHash = {};
	for(var i=0; i < emList.length; i++){
		myHash[emList[i].name] = emList[i].value;
	}
	return myHash;
},


//////////////////////////////////////////
drawingSeriesBarplot: function (jsonObj, containerId, styleJsonObj){

	if (jsonObj["output"].length == 1){
		return;
	}
	var rows = [];
	var maxValue = jsonObj["output"][1][1];
	var legend = jsonObj["description"]["legend"];
	if (legend[1] == ""){
		rows.push([null, legend[0]]);
	} else {
		rows.push([null, legend[0], legend[1]]);
	}

	for (var i = 1; i < jsonObj["output"].length; i++){
		var a = [];
		for (var j in jsonObj["output"][i]){
			a.push(jsonObj["output"][i][j]);
			if (j != 0) {
				maxValue = (maxValue < jsonObj["output"][i][j]) ? jsonObj["output"][i][j] : maxValue;
			}
		}
		rows.push(a);
	}

	var chart  = new Highcharts.Chart({

		colors: jsonObj["description"]["colors"],
		chart: {
			renderTo: containerId.substring(1),
			type: 'column'
		},
		title: {
			text: jsonObj["description"]["title"]
		},
		yAxis: {
			max: maxValue,
			title: {
				text: jsonObj["description"]["yAxis"]
			},labels: {
				style: {
					fontSize: '13px'
				}
			}
		},
		xAxis: {
			tickPositioner: function() {
				var result = [];
				for(i = 0; i < rows.length; i++){
					result.push(i);
				}
				return result;
			},
			title: {
				text: jsonObj["description"]["xAxis"]
			},labels: {
				rotation:-60,
				useHTML: true,
				formatter: function () {
					//console.log(this);
					var xLabel = jsonObj["xLabel"];
					if (xLabel != undefined){
						return '<div class="hastip" title="' + xLabel[this.value] + '">' + this.value + '</div>';
					} else {
						return '<div class="hastip" title="' + this.value + '">' + this.value + '</div>';
					}
				},
				style: {
					fontSize:'13px'
				}
			}
		},
		data: {
			columns: rows,
			switchRowsAndColumns: true
		},
		tooltip: {
			pointFormat: '<span style="color:{series.color}">{series.name}</span>: <b>{point.y}</b> ({point.percentage:.0f}%)<br/>',
			shared: true
		},
		plotOptions: {
			column: {
				stacking: jsonObj["description"]["stack"]
			}
		}
	});

},


//////////////////////////////////////////
drawingSeriesBoxplot: function (jsonObj, containerId, styleJsonObj){

	if (jsonObj["output"].length == 1){
        	return;
	}
	var rows = [];
	var quantiles = [];
	var maxValue = jsonObj["output"][1][1];
	var minValue = jsonObj["output"][1][1];

	for (var i=1; i < jsonObj["output"].length; i++){
		var a = [];
		for (var j in jsonObj["output"][i]){
			a.push(jsonObj["output"][i][j]);
			if (j != 0) {
				maxValue = (maxValue < jsonObj["output"][i][j]) ? jsonObj["output"][i][j] : maxValue;
				minValue = (minValue > jsonObj["output"][i][j]) ? jsonObj["output"][i][j] : minValue;
			}
		}
		quantiles.push(a);
		rows.push(jsonObj["output"][i][0]);
	}

	var chart = Highcharts.chart({

		chart: {
			renderTo: containerId.substring(1),
			type: 'boxplot'
		},
		title: {
			text: jsonObj["description"]["title"]
		},
		yAxis: {
			min: minValue,
			max: maxValue,
			title: {
				text: jsonObj["description"]["yAxis"]
			},labels: {
				style: {
					fontSize:'13px'
				}
			}
		},
		xAxis: {
			categories:rows,
			title: {
				text: jsonObj["description"]["xAxis"]
			},labels: {
				rotation:-60,
				useHTML: true,
				formatter: function () {
					var xLabel = jsonObj["xLabel"];
					if (xLabel != undefined){
						return '<div class="hastip" title="' + xLabel[this.value] + '">' + this.value + '</div>';
					} else {
						return '<div class="hastip" title="' + this.value + '">' + this.value + '</div>';
					}
				},
				style: {
					fontSize:'13px'
				}
			}
		},
		legend: {
			enabled: false
		},
		series: [{
			data: quantiles,
			tooltip: {
				headerFormat: '<em>Cancer Type {point.key}</em><br/>'
			}
			}]
		});
},


//////////////////////////////////////////
//For url link to search
urlParam: function(name){

	var parts1 = window.location.href.split("?");
	if(parts1.length > 1){
		var parts2 = parts1[1].split("&");
		for (var i in parts2){
			var parts3 = parts2[i].split("=");
			if (parts3[0] == name){
				return parts3[1];
			}
		}
	}
	return '';
},

////////////////////////////////////////
setContent: function (containerId, funcName, svcName, inJson, styleJsonObj, svcPath){

	var url = svcPath + svcName + '.py';
	var xhr = new XMLHttpRequest();
	xhr.containerId = containerId;
	xhr.funcName = funcName;
	xhr.open("POST", url, true);
	xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
	xhr.onreadystatechange = function() {
		if (xhr.readyState == 4 && xhr.status == 200) {
			console.log('response='+xhr.responseText);
			var resJson = JSON.parse(xhr.responseText);
			if (resJson["taskStatus"] == 1){
				libFuncs[this.funcName](resJson, this.containerId, styleJsonObj);
			}
			else {
				libFuncs["rndrErrorMsg"](resJson["errorMsg"], this.containerId);
			}
		}
	}
	var inJsonText = JSON.stringify(inJson);
	//var postData = 'mode=json&svc='+svcName+'&injson=' + inJsonText;
	var postData = 'injson=' + inJsonText;

	xhr.send(postData);
	console.log('request='+postData);
	return;
},

//////////////////////////////////////////
setHtmlContent: function (htmlPath, fileName, containerId){

	var url = htmlPath + '/' + fileName;
	var reqObj = new XMLHttpRequest();
	reqObj.containerId = containerId;
	reqObj.open("POST", url, true);
	reqObj.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
	reqObj.onreadystatechange = function() {
		if (reqObj.readyState == 4 && reqObj.status == 200) {
			$(reqObj.containerId).html(reqObj.responseText);
		}
	};
	var postData = '';
	reqObj.send(postData);
},


//////////////////////////////////////////
getInputElement: function (emName, emClass, emValue, emObj){

	emValue = (emValue == 'None' ? '' : emValue);
	var elm = '';
	if (emObj.type == 'textarea'){
		elm += '<textarea name="'+emName+' class="'+emClass+'" rows="'+ emObj.rows +
			'" cols="'+ emObj.cols +'" style="'
		elm +=  emObj.style +'">' + emValue+'</textarea>';
	}
	else if (emObj.type == 'select'){
		var parts1 = emObj.list.split(",");
		elm += '<select name="'+emName+'" class="'+emClass+'" style="'+emObj.style+'" onchange="'+emObj.onchange+'">';
		for (var i=0; i < parts1.length; i++){
			var parts2 = parts1[i].split("=");
			var s = (parts2[0] == emValue ? 'selected' : '');
			elm += '<option value="'+parts2[0]+'" '+s+'>'+parts2[1]+'</option>';
		}
		elm += '</select>';
	}
	else if (emObj.type == 'checkbox'){
		var s = (emValue == 1 ? 'checked' : '');
		elm += '<input type="'+emObj.type+'" name="'+emName+'" '+s+' style="' +
			emObj.style+'" class="'+emClass+'" >';
	}
	else{
		elm += '<input type="'+emObj.type+'" name="'+emName+'" value="'+emValue+'" style="' +
			emObj.style+'" class="'+emClass+'" >';
	}
	return elm;
},


//////////////////////////////////////////
getWaitMsg: function (configObj, waitMsg){

	var imgUrl = configObj["loadingimage"];
	var imgObj = '<img src="'+imgUrl+'">';
	waitMsg = '<table width=100% height=400><tr><td valign=middle align=center>'+imgObj+'</td></tr></table>';
	return;
},


//////////////////////////////////////////
getBatchingBtns: function (batchInfo){

	var style21 = 'position:absolute;height:15;bottom:1px;background:#ccc;color:#fff;text-align:center;';
	var style22 = 'position:absolute;height:13;bottom:1px;text-align:center;color:#777;font-size:11px;';
	var btns = '<a id=start class=batchnavlink href="">' +
		'<DIV style="'+style21+';width:30;right:185;"><<</DIV>' +
		'</a>' +
		'<a id=prev class=batchnavlink href="">' +
		'<DIV style="'+style21+';width:20;right:163;"><</DIV>' +
		'</a>' +
		'<DIV style="'+style22+';width:107;right:53;">page '+batchInfo+'</DIV>' +
		'<a id=next class=batchnavlink href="">' +
		'<DIV style="'+style21+';width:20;right:32;">></DIV>' +
		'</a>' +
		'<a id=end class=batchnavlink href="">' +
		'<DIV style="'+style21+';width:30;right:0;">>></DIV>' +
		'</a>';

	return btns;
}

};




