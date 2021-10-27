var bcoPrefix = '';
var dsPrefix = '';

var childList = {};
var isMain = {};
var nodeid2name = {};
var nodeid2index = {};
var treeJson = [];
var toggleCount = {};
var currentBatch = 1;
var batchSize = 100;
var batchCount = 0;
var datasetName = '';
var readMeFile = '';
var objId = '';
var objVer = '';
var resJson = {};
var filterState = {};
var pageId = 'home';
var categoryList = [];
var seen = {};


////////////////////////////////
$(document).ready(function() {




    var url = '/cgi-bin/init.py';
    var xhr = new XMLHttpRequest();
    xhr.open("POST", url, true);
    xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
    xhr.onreadystatechange = function() {
        if (xhr.readyState == 4 && xhr.status == 200) {
            resJson = JSON.parse(xhr.responseText);
            if(resJson["taskstatus"] != 1){
                $("#resultscn").html("<br><br>" + resJson["errormsg"]);
            }
            else{
                $("#moduleversioncn").html(resJson["moduleversion"]);
                setNavigation(resJson["domains"]);
                dsPrefix = resJson["dsprefix"];
                bcoPrefix = resJson["bcoprefix"];
                setPageLinks();
                fillFrameCn();
            }
        }
    };
    var postData = '';
    xhr.send(postData);
});






/////////////////////////
function modifyMenuLinks(){

    $('.glygenmenu').each(function() {

        var lastPart = this.href.split("/").pop();
        var url_one = "http://glygen.org/" + lastPart;
        var url_two = "http://" + server + ".glygen.org/" + lastPart;
        var menuUrl = (server == "prd" ? url_one : url_two);
        this.href = menuUrl;
    });
    return;
}



////////////////////////////
function setPageLinks(){
    

    var pageLinks = []
    var s = 'font-size:13px;';
    for (var i in resJson["pagelinks"]){
        var obj = resJson["pagelinks"][i];
        var l = '<a class=pagelink href="'+obj["url"]+'" style="'+s+'">'+obj["label"]+'</a>';
        pageLinks.push(l);
    }
    var links = pageLinks.join(" | ");
    $("#pagelinkscn").html(links);
    return;
}


////////////////////////////
function fillFrameCn(){


        var imgFile =  "/content/loading.gif";
        var gifImage = '<img src='+imgFile+' style="width:20%;margin-top:2%;">';        
        $("#pagecn").html(gifImage);

        var urlStr = location.href.substring(0, location.href.length);
        var parts = urlStr.split("/")
        var htmlRoot = parts[0] + '/' + parts[1] + '/' +parts[2];
        pageId = (parts[3] == '' ? 'home' : parts[3]);
        objId = (pageId.indexOf(dsPrefix) != -1 ? pageId : "");
        pageId = (pageId.indexOf(dsPrefix) != -1 ? "view" : pageId);
        


        if(pageId == 'home'){
		fillGridViewCn("1");
        }
        else if(pageId == 'view'){
            fillEntryViewCn();
        }
        else if(pageId == 'datamodel'){
		fillDataModelViewCn();
	}
        else if(pageId == 'history'){
            fillReleaseHistoryViewCn({});
        }
	else{
                fillStaticHtmlCn(htmlRoot + '/content/page.'+pageId+'.html', '#pagecn');
        }
	return;
}


////////////////////////////
function fillGridViewCn(){


        var imgFile =  "/content/loading.gif";
        var gifImage = '<img src='+imgFile+' style="width:20%;margin-top:2%;">';
        $("#pagecn").html(gifImage);

	var url = '/cgi-bin/search_objects.py';
	var xhr = new XMLHttpRequest();
	xhr.open("POST", url, true);
	xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
        xhr.onreadystatechange = function() {
                if (xhr.readyState == 4 && xhr.status == 200) {
                        //console.log('response='+xhr.responseText);
                        resJson = JSON.parse(xhr.responseText);
			rndrGridContent();
                    }
        };

        var inJson = {}
        var queryValue = $("#queryvalue").val().trim();
        inJson = {"queryvalue":queryValue};
        
        var postData = 'injson=' + JSON.stringify(inJson);
        xhr.send(postData);
	console.log(postData);

}

//////////////////////////////
function fillReleaseHistoryViewCn(catJson){

        var url = '/cgi-bin/get_releasehistory.py';
        var xhr = new XMLHttpRequest();
        xhr.open("POST", url, true);
        xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
        xhr.catJson = catJson;
        xhr.onreadystatechange = function() {
                if (xhr.readyState == 4 && xhr.status == 200) {
                        //console.log('response='+xhr.responseText);
                        resJson = JSON.parse(xhr.responseText);
                        if(resJson["taskstatus"] != 1){
                            $("#pagecn").html("<br><br>" + resJson["errormsg"]);
                        }
                        else{
                            var cn = '<table style="width:100%;">';
                            cn += '<tr><td id="filterscn"></td></tr>';
                            cn += '<tr><td id="releasehistorycn"></td></tr>';
                            cn += '</table>';
                            $("#pagecn").html(cn);
                            var catSelector = '<select id=catcombo>';
                            catSelector += '<option value="">--</option>';
                            for (var catName in resJson["categories"]){
                                var col = '<b>' + catName + '</b><br>';
                                for (var j in resJson["categories"][catName]){
                                    var catValue = resJson["categories"][catName][j];
                                    var catCombo = catName + ':' + catValue
                                    var s = (this.catJson["category_value"] == catValue ? "selected": "");

                                    catSelector += '<option value="'+catCombo+'" '+s+'>'+catCombo+'</option>';
                                }
                            }
                            catSelector += '</select>';
                            var filtersCn = '<br><table style="width:100%;">';
                            filtersCn += '<tr><td align=left><b>Filter by category</b><br>';
                            filtersCn += catSelector + '</td></tr>';
                            filtersCn += '</table><br>';
                            $("#filterscn").html(filtersCn);
                            drawTable(resJson["dataframe"], "releasehistorycn", {"pagesize":1000});
                        }
                }
        };
        objId = objId.replace("/", "");
        var inJson = catJson;
        var postData = 'injson=' + JSON.stringify(inJson);
        xhr.send(postData);
        console.log(postData);


}

////////////////////////////
function fillEntryViewCn(){

        var url = '/cgi-bin/get_dataset.py';
        var xhr = new XMLHttpRequest();
        xhr.open("POST", url, true);
        xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
        xhr.onreadystatechange = function() {
                if (xhr.readyState == 4 && xhr.status == 200) {
                        //console.log('response='+xhr.responseText);
                        resJson = JSON.parse(xhr.responseText);
                        if(resJson["status"] != 1){
                            $("#pagecn").html("<br><br>" + resJson["errormsg"]);
                        }
                        else{
                            rndrEntryContent();
                        }
                }
        };
        objId = objId.replace("/", "");
        var inJson = {"objid":objId, "objver":objVer};
        var postData = 'injson=' + JSON.stringify(inJson);
        xhr.send(postData);
        console.log(postData);
}


/////////////////////
function rndrEntryContent(){


        var latestObjVer = resJson["versions"][0].split(" ")[0];

        var versionSelector = '<select class=versionselector id=versioncn>';
        for (var i in resJson["versions"]){
            var parts = resJson["versions"][i].split(" ");
            var s = (parts[0] == resJson["selectedversion"] ? "selected" : "");
            versionSelector += '<option value="'+parts[0]+'" '+s+'>'+resJson["versions"][i]+'</option>';
        }
        versionSelector += '</select>';

        var dsId = resJson["info"]["objid"].replace(bcoPrefix, dsPrefix);

        var bcoId = resJson["info"]["objid"];
        var datasetName = resJson["info"]["filename"];
        var downloadUrl = '/ln2wwwdata/reviewed/'+datasetName
        var readmeUrl = '/' + dsId + '/'  + latestObjVer  + '/readme';


        var bcoUrl = '/' + bcoId + '/'  + latestObjVer;
        if (objVer != ''){
                downloadUrl =  '/ln2wwwdata/reviewed/' + objVer+ '/' +datasetName
                readmeUrl = '/' + dsId + '/'  + objVer  + '/readme';
                bcoUrl = '/' + bcoId + '/'  + objVer;
        }
       
        var links = '<a href="'+bcoUrl+'" style="font-size:12px;" target=_>BCO JSON</a>';
        links += ' | <a href="'+readmeUrl+'" style="font-size:12px;" target=_>README</a>';
        links += ' | <a href="'+downloadUrl+'" download="'+datasetName+'" style="font-size:12px;">DOWNLOAD</a>';
        //var linkId = "obj_" + parseInt(resJson["info"]["objid"].replace("ONCOMX", ""));
	//links += ' | <a id="'+linkId+'" class=commentlink style="font-size:12px;">COMMENT</a>';
        


        var s1 = 'font-size:16px;font-weight:bold;color:#004065;'; 
        var d1 = 'position:static;width:100%;background:#fff;';
        var divcn = '<div style="'+d1+'" id=datasetcn></div>';
        var cn = '<table width=100% style="font-size:13px;" border=0>';
        cn += '<tr height=30><td colspan=2 >&nbsp;</td></tr>';
        cn += '<tr><td colspan=2 align=left>'+dsId+' sample view</td></tr>';
        cn += '<tr><td style="'+s1+'" colspan=2 align=left>'+resJson["info"]["title"]+'</td></tr>';
        cn += '<tr><td align=left colspan=2>'+resJson["info"]["description"]+'</td></tr>';
        cn += '<tr><td style="border-bottom:1px solid #ccc;" align=left><br>Version<br>'+versionSelector+'</td><td align=right style="border-bottom:1px solid #ccc;" valign=bottom><br><br>'+links+'</td></tr>';
        cn += '</table><br>';
        cn += divcn;
        $("#pagecn").html(cn);
        
        if (["csv", "tsv"].indexOf(resJson["info"]["filetype"]) != -1){
            drawTable(resJson["dataframe"], "datasetcn", {"pagesize":1000});
        }
        else if (resJson["info"]["filetype"] == "fasta"){
            drawFasta(resJson["seqobjects"], "#datasetcn");
        }
        else if (resJson["info"]["rndrtype"] == "text"){
            drawPlainText(resJson["txtbuffer"], "#datasetcn");
        }
        else{
            drawHtmlText(resJson["txtbuffer"], "#datasetcn");
        }
        return
}


////////////////////////////////////////
function drawFasta(seqObjs, containerId){

    var cn = '<table width=100% style="font-size:12px;" border=0>';
    for (var i in seqObjs){
        var obj = seqObjs[i];
        cn += '<tr><td align=left>>'+ obj.seqid + ' ' + obj.seqdesc + '</td></tr>';
        var seq = "";
        var lineLen = 100;
        for (var j=0; j < parseInt(obj.seqbody.length/lineLen) + 1; j++){
            var startPos = j*lineLen;
            endPos = startPos + lineLen;
            endPos = (endPos > obj.seqbody.length - 1 ? obj.seqbody.length - 1: endPos);
            seq += obj.seqbody.substring(startPos, endPos) + '\n';
        }
        cn += '<tr><td align=left><pre>'+seq+'</pre></td></tr>';
        cn += '<tr height=30><td>&nbsp;</td></tr>';
    }
    cn += '</table>';
    $(containerId).html(cn);

    return;
}


////////////////////////////////////
function drawPlainText(txtBuffer, containerId){

    var cnBody = '<textarea style="width:100%;" rows=45>'+txtBuffer+'</textarea>';
    var cn = '<table width=100% style="font-size:12px;" border=0>';
    cn += '<tr><td><pre>'+cnBody+'<pre></td></tr>';
    cn += '</table>';
    $(containerId).html(cn);
    return;
}

////////////////////////////////////
function drawHtmlText(txtBuffer, containerId){

    var cn = '<table width=100% style="font-size:12px;" border=0>';
    cn += '<tr><td><pre>'+txtBuffer+'<pre></td></tr>';
    cn += '</table>';
    $(containerId).html(cn);
    return;
}




///////////////////////////
function fillJsonTextCn(){

	var link = '<a class=pagelink id=saveobject href="#" style="font-size:13px;">Save Object</a>'; 
	var cn = '<table width=100%>';
	cn += '<tr><td align=right>'+link+'</td></tr>';
	cn += '<tr><td><textarea id=jsontext rows=40 style="width:100%;padding:10;"></textarea></td></tr>';
	cn += '</table>';
	$("#pagecn").html(cn);

	var url = '/cgi-bin/get_single_object.py';
       	var xhr = new XMLHttpRequest();
       	xhr.open("POST", url, true);
       	xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
       	xhr.onreadystatechange = function() {
               	if (xhr.readyState == 4 && xhr.status == 200) {
                       	$("#jsontext").val(xhr.responseText);
               	}
       	};
        var inJson = {"objid":objId};
        var postData = 'injson=' + JSON.stringify(inJson);
       	xhr.send(postData);
       	console.log(postData);
}


///////////////////////////
function fillDataModelViewCn(){


	var url =  '/content/page.datamodel.html';
	var xhr = new XMLHttpRequest();
	xhr.open("GET", url, true);
	xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
	xhr.onreadystatechange = function() {
		if (xhr.readyState == 4 && xhr.status == 200) {
			//console.log('response='+xhr.responseText);
			$("#pagecn").html(xhr.responseText);
			var cn = '<br><table width=100% border=0>' + 
				'<tr><td id=predicatescn></td></tr></table>';
			$("#pagecn").append(cn);
			rndrDataModelTable();
		}
	};
	xhr.send();
			        

}


///////////////////////////
function rndrDataModelTable(){

	var url = '/cgi-bin/get_data_model_table.py';
	var xhr = new XMLHttpRequest();
	xhr.open("POST", url, true);
	xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
	xhr.onreadystatechange = function() {
		if (xhr.readyState == 4 && xhr.status == 200) {
                        var resObj = JSON.parse(xhr.responseText);
			drawTable(resObj["dataframe"], "predicatescn", {"pagesize":1000});
		}
	};
        var inJson = {};
        var postData = 'injson=' + JSON.stringify(inJson);
	xhr.send(postData);
}


//////////////////////////////////////////
function rndrGridContent(){

	var gridWidth = 250;
	var gridHeight = 140;
	var npass = 0;
	
        categoryList = resJson["categorylist"];
        seen = {"total":0, "datasetcount":{}};
        for (var i in categoryList){
            var catName = categoryList[i].toLowerCase();
            seen[catName] = {};
            if (!(catName in filterState)){
                filterState[catName] = {};
            }
            seen["datasetcount"][catName] = {};
        }

        var gridDivs = '';
        for (var i in resJson["datasets"]){
                        var obj = resJson["datasets"][i];
			var ii = parseInt(i) + 1;
                        var filterOutFlag = false;
                        for (var j in categoryList){
                            var catName = categoryList[j];
                            var catValue = "Other " + catName.substring(0, 2) + "...";
                            if (catName in obj["categories"]){
                                var catValue = obj["categories"][catName];
                            }
                            
                            if(!(catValue in filterState[catName])){
                                filterState[catName][catValue] = true;
                            }
                            if (filterState[catName][catValue] == false){
                                filterOutFlag = true;
                            }
                            seen[catName][catValue] = true;
                            seen["datasetcount"][catName][catValue] = (!(catValue in seen["datasetcount"][catName]) ? 1 : seen["datasetcount"][catName][catValue] + 1);
                        }
                        seen["total"] += 1;

                        if(filterOutFlag == true){
                            continue;
                        }
                        var xPos = i*(gridWidth + 80);
			var iconFileName = obj["iconfilename"];
                        var miniTable = '';
			if ("minitable" in obj){
				miniTable = rndrMiniTable(obj["minitable"]);
			}
                        npass += 1;
                        var moleculeType = ("molecule" in obj["categories"] ? 
                            obj["categories"]["molecule"] : "");
                        var speciesType = ("species" in obj["categories"] ?
                            obj["categories"]["species"] : "");
                        var fileType = ("file_type" in obj["categories"] ?
                            obj["categories"]["file_type"] : "");
                        var statusType = ("status" in obj["categories"] ?
                            obj["categories"]["status"] : "");

                        gridDivs += '<div class="gridcn">';
			var s = 'display:block;float:left;margin:20px 1% 0px 1%;';
                        s += 'color:#333;font-size:11px;text-align:center;width:98%;border:0px dashed orange;';
                        objId = bcoPrefix + "0000".substring(0, 10 - String(obj["_id"]).length) + String(obj["_id"]);
                        objId = obj["_id"].replace(bcoPrefix, dsPrefix);
                        var titleText = statusType + ' ' + moleculeType.toLowerCase();
                        titleText += ' dataset ' + objId + ' in ';
                        titleText += fileType.toUpperCase() + ' format.';
                        titleText += '<br>' + speciesType;
                        gridDivs += '<div style="'+s+'">';
                        gridDivs += titleText;
                        gridDivs += '</div>';
                        var s = 'display:block;float:left;width:90%;font-size:16px;margin:20px 5% 0px 5%;';
			s += 'font-weight:bold;color:#004065;vertical-align:bottom;text-align:center;';
                        s += 'border:0px dashed orange;';
			var detailsUrl = '/' + objId ;
			gridDivs += '<a href="'+detailsUrl+'" style="font-size:12px;">';
			gridDivs += '<div style="'+s+'">';
			gridDivs += obj["title"];
			gridDivs += '</div>';
			gridDivs += '</a>';
		
			var iconUrl = '/content/' + iconFileName;
			var s = 'display:block;float:left;width:80%;height:120px;margin:20px 10% 0px 10%;';
			s += 'font-size:12px;text-align:center;border:0px dashed orange;';
			gridDivs += '<div style="'+s+'">';
			gridDivs += miniTable;
			gridDivs += (miniTable == '' ? '<img src="'+iconUrl+'" height=90%>' : "");
			gridDivs += '</div>';
		
			var s = 'display:block;float:left;width:80%;margin:0px 10% 0px 10%;';
                        s += 'font-size:12px;text-align:center;border:0px dashed orange;';
			gridDivs += '<div style="'+s+'">';
			gridDivs += obj["description"] + ' ...<br>';
                        gridDivs += '<a href="'+detailsUrl+'" style="font-size:12px;">view details</a>';
			gridDivs += '</div>';	
			gridDivs += '</div>';
	}




        var filters = '<table width=100% border=0 style="font-size:13px;"><tr>';
        for (var i in categoryList){
	    var catName = categoryList[i];
	    var cn_one = '';
            var cn_two = '';
            for (var x in seen[catName]){
	    	var n = (x in seen["datasetcount"][catName] ? seen["datasetcount"][catName][x] : 0);
                var label = x + ' (' + n + ')';
		var chkd = "checked";
                if (x in filterState[catName]){
		    chkd = (filterState[catName][x] == true ? "checked" : "");
		}
                var chkbox = '<input class=filtercheckbox type=checkbox name='+catName+' '+chkd+' width=15 value="'+x+'">';
                if (x.substring(0, 5) != "Other"){
                    cn_one += '&nbsp;&nbsp;&nbsp;' + chkbox + " " + label + '<br>';
	        }
                else{
                    cn_two += '&nbsp;&nbsp;&nbsp;' + chkbox + " " + label + '<br>';
                }
            }
            var cn = '<b>Filter by '+catName+'</b><br>' + cn_one + cn_two;
            filters += '<td valign=top align=left>'+cn+'</td>';
        }
        filters += '</tr></table>';



	var s =  'width:80px;height:25px;';
        var applybtn = '<input type=submit class=filterbtn id=apply style="'+s+'" value=" Apply ">';
        var resetbtn = '<input type=submit class=filterbtn id=reset style="'+s+'" value=" Reset ">';
	var btns = '';
	var filterLink = '<a id=filterlink href="" style="font-size:13px;">Filters</a>';
	var filterTable = '<table width=100% border=0 style="font-size:13px;">' +
			'<tr height=20 style="border-bottom:1px solid #eee;">' + 
			'<td align=left>Total of '+ seen["total"]+ ' datasets ('+npass +' passed filter)</td>' + 
        		'<td align=right>'+filterLink+'</td>' + 
			'</tr>' + 
        		'<tr><td class=filtercontainer colspan=2 style="padding:20px;">'+filters+'</td></tr>' +
			'</table>';
        


        var style = 'display:block;float:left;width:100%;margin:20px 0px 20px 0px;border:1px solid #eee;';
	style += 'padding:10px;background:#fff;';
        var gridTable = '<div style="'+style+'">';
	gridTable +=  filterTable ;
	gridTable += '</div>';
        gridTable += gridDivs;
	$("#pagecn").html(gridTable);
	return
}



////////////////////////
function popMessage(msg){

	event.preventDefault(); 
	$("html, body").animate({ scrollTop: 0 }, "0");
	
	var s = "padding:2px 20px 2px 20px;";
        var closebtn = '<input name=btn2 id=closewindow type=submit style="'+s+'" value="&times;">';
        var table = '<table width=100% style="font-size:13px;" border=0>' +
                        '<tr height=25><td align=right>'+ closebtn+'</td></tr>' +
                        '</table>';

	var s = 'position:absolute;left:1%;top:5px;width:98%;height:25px;';
        s += 'filter: alpha(opacity=100);opacity: 1.00;z-index:1003;border:0px solid;';
        var div1 = '<DIV id=popdiv1 style="'+s+'">'+table+'</DIV>';

        var s = 'position:absolute;left:10%;top:50;width:80%;height:400;overflow:auto;';
        s += 'background:#fff;filter: alpha(opacity=100);opacity: 1.00;z-index:1002;padding:10px;';
        var div2 = '<DIV id=popdiv2 style="'+s+'">'+msg+'</DIV>';

        var s = 'position:absolute;left:10%;width:80%;height:500;top:5%;';
        s += 'background:#eee;filter: alpha(opacity=100);opacity: 1.00;z-index:1001;';
        var popdiv = '<DIV id=popdiv style="'+s+'">'+ div1 + div2+'</DIV>';

        var s = 'position:absolute;left:0px;top:0px;width:100%;height:3000px;background:#000;color:#fff;';
        s += 'filter: alpha(opacity=75);z-index:1000;';
        s += 'opacity: 0.75;';
        var bgdiv = '<DIV id=bgdiv style="'+s+'"></DIV>';

        $("body").append(bgdiv + popdiv);

	return;
}


//////////////////////////////////
function setReadmeContent(fileName, containerId){


        var url = '/' + fileName;
        var xhr = new XMLHttpRequest();
        xhr.containerId = containerId;
        xhr.open("GET", url, true);
        xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
        xhr.onreadystatechange = function() {
                if (xhr.readyState == 4 && xhr.status == 200) {
                        $(xhr.containerId).html('<pre>'+xhr.responseText + '</pre>');
                }
                else{
                        var msg = fileName + ' does not exist!';
                        var table = '<table width=100%>' +
                                '<tr height=400><td style="color:red;" align=center> ' + msg + '</td></tr>' +
                                '</table>';
                        $(xhr.containerId).html(table);
                }
        };
        xhr.send();
}






///////////////////////
function isValidJson(str) {
	    
	try {
		JSON.parse(str);
	} catch (e) {
		return e;
	}
	return true;
}


////////////////////////////////////
$(document).on('click', '.titlecn', function (event) {
    event.preventDefault();
    var jqId = "#" + this.id.replace("title", "")
    $(jqId).toggle();
    var jqId = "#" + this.id.replace("title", "download")
    $(jqId).toggle();

    if (this.id.indexOf("idmapping") != -1){
        $(".idmappingcn").toggle();
    }
    else if (this.id.indexOf("sitemapping") != -1){
        $(".sitemappingcn").toggle();
    }

});



$(document).on('click', '#savefile', function (event) {
    event.preventDefault();


    var inJson = {"action":"save_file", "filename":resJson["fileinfo"]["filename"]}; 
    $('.savefilefield').each(function() {
        inJson[$(this).attr("id")] = $(this).val();
    });
    
    for (var f in inJson){
        if (inJson[f].trim() == ""){
            alert("Please provide value for " + f );
            return false;
        }
    }
    
    var imgFile = "/content/loading.gif";
    var gifImage = '<img src='+imgFile+' style="width:20%;margin-top:2%;">';
    $("#progresscn").html(gifImage);

    let xhr = new XMLHttpRequest();
    xhr.onloadend = function() {
        if (xhr.status == 200) {
            console.log(xhr.responseText);
            resJson = JSON.parse(xhr.responseText);
            $("#allwrappercn").css("display", "none");
            $("#progresscn").css("display", "block");
            $("#progresscn").html("<br><br><br>File saved successfully!")
        } else {
            console.log("error " + this.status);
        }
    };
    var url = "/cgi-bin/upload.py?injson=" + JSON.stringify(inJson);
    xhr.open("GET", url);
    xhr.send();
    console.log(url);

});
    

$(document).on('click', '#submitsitemapping', function (event) {
    event.preventDefault();
    var imgFile = "/content/loading.gif";
    var gifImage = '<img src='+imgFile+' style="width:20%;margin-top:2%;">';
    $("#sitemappingresultcn").html(gifImage);

    var idField = $("#sitemapidfield option:selected").val();
    var posField = $("#sitemapposfield option:selected").val();
    var residueField = $("#sitemapresiduefield option:selected").val();
    var inJson = {
        "action":"perform_sitemapping",
        "filename":resJson["fileinfo"]["filename"],
        "idfield":idField,
        "posfield":posField,
        "residuefield":residueField
    };


    let xhr = new XMLHttpRequest();
    xhr.onloadend = function() {
        if (xhr.status == 200) {
            //console.log(xhr.responseText);
            resJson = JSON.parse(xhr.responseText);
            drawTable(resJson["mappingrows"], "sitemappingresultcn", {"pagesize":1000});
            drawTable(resJson["previewrows"], "previewcn", {"pagesize":1000});
            var fileUrl = '/ln2uploads/tmp/' + resJson["fileinfo"]["filename"]
            var downloadLink = '<a href="'+fileUrl+'" download>Download</a>';
            $("#previewdownloadcn").html(downloadLink);
        } else {
            console.log("error " + this.status);
        }
    };
    var url = "/cgi-bin/upload.py?injson=" + JSON.stringify(inJson);
    xhr.open("GET", url);
    xhr.send();
    console.log(url);

});



$(document).on('click', '#submitidmapping', function (event) {
    event.preventDefault();

    var imgFile = "/content/loading.gif";
    var gifImage = '<img src='+imgFile+' style="width:20%;margin-top:2%;">';
    $("#idmappingresultcn").html(gifImage);

    var selectedType = $("#idmapidtype option:selected").val();
    var selectedField = $("#idmapidfield option:selected").val();
    var inJson = {
        "action":"perform_idmapping", 
        "filename":resJson["fileinfo"]["filename"], 
        "fieldname":selectedField,
        "fieldtype":selectedType
    };

    let xhr = new XMLHttpRequest();
    xhr.onloadend = function() {
        if (xhr.status == 200) {
            //console.log(xhr.responseText);
            resJson = JSON.parse(xhr.responseText);
            drawTable(resJson["mappingrows"], "idmappingresultcn", {"pagesize":1000});
            drawTable(resJson["previewrows"], "previewcn", {"pagesize":1000});
            var fileUrl = '/ln2uploads/tmp/' + resJson["fileinfo"]["filename"]
            var downloadLink = '<a href="'+fileUrl+'" download>Download</a>';
            $("#previewdownloadcn").html(downloadLink);

        } else {
            console.log("error " + this.status);
        }
    };
    var url = "/cgi-bin/upload.py?injson=" + JSON.stringify(inJson);
    xhr.open("GET", url);
    xhr.send();
    console.log(url);

});



//////////////////////////////////
$(document).on('click', '#submitfile', function (event) {
    event.preventDefault();
 
    $("#allwrappercn").css("display", "block");
    $("#savecn").css('display', 'none');
    $("#summarywrappercn").css('display', 'none');
    $("#sanitywrappercn").css('display', 'none');
    $("#previewwrappercn").css('display', 'none');
    $("#savewrappercn").css('display', 'none');
    $("#idmappingwrappercn").css('display', 'none');
    $("#sitemappingwrappercn").css('display', 'none');

    $("#progresscn").css('display', 'block');
                                                    

    var imgFile = "/content/loading.gif";
    var gifImage = '<img src='+imgFile+' style="width:20%;margin-top:2%;">';
    $("#progresscn").html(gifImage);

    var file = $('#userfile')[0].files[0];
    var formData = new FormData();
    formData.append("userfile", file);

    var sizeLimit = 1000000000;
    if (file.size > sizeLimit){ 
        $("#summarywrappercn").css('display', 'block');
        $("#summarywrappercn").css('height', '200px');
        $("#summarywrappercn").css('padding', '80px 0px 0px 0px');
        var msg = 'Your submitted file is ' + file.size + ' Bytes big. ';
        msg += 'This exceeds maximum allowed file size of ' + sizeLimit + ' Bytes.';
        $("#summarywrappercn").html('<font color=red>' + msg + '</font>');
        $("#progresscn").css('display', 'none');
        return;
    }


    let xhr = new XMLHttpRequest();
    // track upload progress
    //xhr.upload.onprogress = function(event) {
        //console.log(`Uploaded ${event.loaded} of ${event.total}`);
    //};

    // track completion: both successful or not
    xhr.onloadend = function() {
        if (xhr.status == 200) {
            //console.log(xhr.responseText);
            resJson = JSON.parse(xhr.responseText);
            $("#progresscn").css('display', 'none');
            if (resJson["sanityrows"].length > 2){
                $("#sanitywrappercn").css('display', 'block');
                $("#qcreportwrappercn").css('display', 'block');
                drawTable(resJson["qcreportrows"], "qcreportcn", {"pagesize":1000});
                drawTable(resJson["sanityrows"], "sanitycn", {"pagesize":1000});
            }
            else{
                $("#qcreportwrappercn").css('display', 'block');
                $("#savewrappercn").css('display', 'block');
                $("#idmappingwrappercn").css('display', 'block');
                $("#sitemappingwrappercn").css('display', 'block');

                $("#previewwrappercn").css('display', 'block');
               
                drawTable(resJson["qcreportrows"], "qcreportcn", {"pagesize":1000});
                drawTable(resJson["previewrows"], "previewcn", {"pagesize":1000});
                var fileUrl = '/ln2uploads/tmp/' + resJson["fileinfo"]["filename"]
                var downloadLink = '<a href="'+fileUrl+'" download>Download</a>';
                $("#previewdownloadcn").html(downloadLink);


                var cn_one = 'GlyGen ID Type<br>';
                cn_one += '<select id=idmapidtype style="width:250px;">';
                cn_one += '<option value="protein">Protein</option>';
                cn_one += '<option value="glycan">Glycan</option>';
                cn_one += '</select>';

                var fieldList = resJson["fileinfo"]["fieldlist"];
                var options = ''
                for (var j in fieldList){
                    options += '<option value="'+fieldList[j]+'">'+fieldList[j]+'</option>';
                }
                var cn_two = 'User ID Field<br>';
                cn_two += '<select id=idmapidfield style="width:250px;">';
                cn_two += options + '</select>';                                
               
                var cn_three = 'User Protein ID Field<br>';
                cn_three += '<select id=sitemapidfield style="width:250px;">';
                cn_three += options + '</select>';
        
                var cn_four = 'User Protein Position Field<br>';
                cn_four += '<select id=sitemapposfield style="width:250px;">';
                cn_four += options + '</select>';

                var cn_five = 'User Protein Residue Field<br>';
                cn_five += '<select id=sitemapresiduefield style="width:250px;">';
                cn_five += options + '</select>';

                cn_two += '<input id="submitidmapping" type="submit" value="Submit">';
                cn_five += '<input id="submitsitemapping" type="submit" value="Submit">';

                $("#idmapselectorone").html(cn_one);
                $("#idmapselectortwo").html(cn_two);
               
                $("#sitemapselectorone").html(cn_three);
                $("#sitemapselectortwo").html(cn_four);
                $("#sitemapselectorthree").html(cn_five);

            }
        } else {
            console.log("error " + this.status);
        }
    };
    xhr.open("POST", "/cgi-bin/upload.py");
    xhr.send(formData);

});



$(document).on('change', '#submitfile', function (event) {
    event.preventDefault();


    var fileObj = document.getElementById("userfile");
    var fileList = fileObj.files;
    var formData = new FormData();
    formData.append("userfile", fileList[0]);

    var url = '/cgi-bin/upload.py';
    var xhr = new XMLHttpRequest();
    xhr.open("POST", url, true);
    xhr.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
    xhr.onreadystatechange = function() {
        if (xhr.readyState == 4 && xhr.status == 200) {
            //console.log(xhr.responseText);
            $("#pagecn").html(xhr.responseText);
        }
        else{
            $("#pagecn").html("error");
        }
    };
    //var inJson = {"objid":objId};
    //var postData = 'injson=' + JSON.stringify(inJson);
    xhr.send(formData);
    console.log(formData);

});



$(document).on('click', '#searchbtn', function (event) {
    event.preventDefault();
    $('html').animate({scrollTop:0}, 'fast');
    $('body').animate({scrollTop:0}, 'fast');
    fillGridViewCn("1");
});





$(document).on('click', '#filterlink', function (event) {
    event.preventDefault();
    $(".filtercontainer").toggle();
});






$(document).on('click', '.xxxx', function (event) {
        event.preventDefault();

	pageId = this.id;
	setPageLink();
        fillFrameCn();

});





///////////////////////////////////////////////////
$(document).on('click', '.filtercheckbox', function (event) {


	event.preventDefault();
        var clickedCatValue = $(this).val();
        for (var j in categoryList){
            var catName = categoryList[j];
            $("input[type=checkbox][name="+catName+"]").each(function () {
                var catValue = $(this).val();
                if (catValue == clickedCatValue){
                    var currentState = filterState[catName][catValue];
                    var newState = (currentState == false ? true : false)
                    filterState[catName][$(this).val()] = newState;
                }
            });
        }

	rndrGridContent();

});



///////////////////////////////////////////////////
$(document).on('change', '#catcombo', function (event) {

    event.preventDefault();
    var catCombo = $("#catcombo option:selected").val();
    var catJson = {}
    if (catCombo != ''){
        var parts = catCombo.split(":");
        catJson = {"category_name":parts[0], "category_value":parts[1]};
    }
    fillReleaseHistoryViewCn(catJson);

});




///////////////////////////////////////////////////
$(document).on('change', '.versionselector', function (event) {
        
    event.preventDefault();
    objVer = $("#versioncn option:selected").val();
    var imgFile = "/content/loading.gif";
    var gifImage = '<img src='+imgFile+' style="width:20%;margin-left:40%;margin-top:2%;">';
    $("#pagecn").html(gifImage);

    fillEntryViewCn();

});




///////////////////////////////////////////////////
$(document).on('click', '.readmelink', function (event) {
 
	event.preventDefault(); 
	var i = parseInt(this.id.split("_")[1]);

	$("html, body").animate({ scrollTop: 0 }, "0");
	
 
	var s = "padding:2px 20px 2px 20px;";
        var closebtn = '<input name=btn2 id=closewindow type=submit style="'+s+'" value="&times;">';
        var table = '<table width=100% style="font-size:13px;" border=0>' +
                        '<tr height=25><td align=right>'+ closebtn+'</td></tr>' +
                        '</table>';

	var s = 'position:absolute;left:1%;top:5px;width:98%;height:25px;';
        s += 'filter: alpha(opacity=100);opacity: 1.00;z-index:1003;border:0px solid;';
        var div1 = '<DIV id=popdiv1 style="'+s+'">'+table+'</DIV>';

        var s = 'position:absolute;left:5%;top:5%;width:90%;height:90%;overflow:auto;';
        s += 'background:#fff;filter: alpha(opacity=100);opacity: 1.00;z-index:1002;padding:10px;';
        var div2 = '<DIV id=popdiv2 style="'+s+'"></DIV>';

        var s = 'position:absolute;left:10%;width:80%;height:90%;top:5%;';
        s += 'background:#eee;filter: alpha(opacity=100);opacity: 1.00;z-index:1001;';
        var popdiv = '<DIV id=popdiv style="'+s+'">'+ div1 + div2+'</DIV>';

        var s = 'position:absolute;left:0px;top:0px;width:100%;height:3000px;background:#000;color:#fff;';
        s += 'filter: alpha(opacity=75);z-index:1000;';
        s += 'opacity: 0.75;';
        var bgdiv = '<DIV id=bgdiv style="'+s+'"></DIV>';

        $("body").append(bgdiv + popdiv);



	var readmeFile = '/ln2wwwdata/reviewed/' + resJson["datasets"][i]["readmefilename"];
	setReadmeContent(readmeFile, '#popdiv2');

});

///////////////////////////////////////
$(document).on('click', '#closewindow', function (event) {
        event.preventDefault();
        $("#popdiv1").remove();
        $("#popdiv2").remove();
        $("#popdiv").remove();
        $("#bgdiv").remove();
});


///////////////////////////////////////
$(document).on('click', '#btn', function (event) {
    event.preventDefault();
    file_size = document.getElementById("my_file").files[0].size;
    alert(file_size);
    
});



function setNavItemAsCurrent(itemText) {
     $('.nav > li > a').each(function () {
        if ($(this).text().indexOf(itemText) >= 0) {
                $(this).parent().addClass('current');
        }
    });
}


function setNavigation(domainUrls){

    var url = window.location.href;
    var fullFilename = url.substring(url.lastIndexOf('/') + 1);
    var filename = fullFilename.substring(0, fullFilename.lastIndexOf('.'));
    var navItemText = filename.replace(/_/g, ' ').toUpperCase();
    var glygen_url = window.location.origin;
    if (glygen_url.indexOf('beta-') >= 0) {
        glygen_url = glygen_url.replace("beta-", "beta.");
    }
    if (glygen_url.indexOf('data.') >= 0) {
       glygen_url = glygen_url.replace("data.", "");
    } else if (glygen_url.indexOf('sparql.') >= 0) {
       glygen_url = glygen_url.replace("sparql.", "");
    }
   
   var domain = glygen_url + "/";
    
    if (navItemText == '') {
        navItemText = 'HOME';
    } else if (navItemText == 'INDEX') {
        navItemText = 'HOME';
    } else if (navItemText == 'CONTACT') {
        navItemText = 'HELP';
    } else if (navItemText == 'HOW TO CITE') {
        navItemText = 'HELP';
    } else if (navItemText == 'ABOUT') {
        navItemText = 'HELP';
    } else if (navItemText == 'RESOURCES') {
        navItemText = 'MORE';
    } else if (navItemText == 'MEDIA') {
        navItemText = 'MORE';
    } else if (navItemText == 'FRAMEWORKS') {
        navItemText = 'MORE';
    } else if (navItemText == 'GLYGEN SETTINGS') {
        navItemText = 'MY GLYGEN';
    }

    if (url.indexOf('data.') >= 0) {
        navItemText = 'DATA';
    } else if (url.indexOf('sparql.') >= 0) {
        navItemText = 'SPARQL';
    }
    else if (url.indexOf('api.') >= 0) {
        navItemText = 'API';
    }

    $('.nav > li').removeClass('current');
    setNavItemAsCurrent(navItemText);
    
    $("#a_portal").attr('href', domainUrls["portal"]);
    $("#a_data").attr('href', domainUrls["data"]);
    $("#a_api").attr('href', domainUrls["api"]);
    $("#a_sparql").attr('href', domainUrls["sparql"]);

    
    $.each($(".a_header"), function(i, v) {
        var nav_url = $(v).attr('href');
        $(v).attr('href', domain + nav_url);
    });


}



