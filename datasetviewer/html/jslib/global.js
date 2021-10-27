var boxListIndex = 0;




/////////////////////////
function getSections(){

    	var style1 = 'background:'+moduleMenuBg+';' + 'color:'+moduleMenuFg+';';
    	var style2 = 'background:'+moduleMenuFg+';' + 'color:'+moduleMenuBg+';';


    	var rows = '';
    	for(var i=0;i< sectionObjs.length;i++){
        	var cls = (sectionObjs[i].id ==  pageId  ? "menucellselected" : "menucell");
		var style = (sectionObjs[i].id ==  pageId  ? style2 : style1);
        	var urlBase = (moduleVersion == 'tst' ? '/tst/' + moduleBase : '/'+moduleBase);
		var args = urlBase + '/' + sectionObjs[i].id;
        	args +=  (sectionObjs[i].action ? '&action='+ sectionObjs[i].action : '');
		//var link = '<a href="' + args +' " class="'+cls+'"  style="'+style+'" d="' + 
		//			sectionObjs[i].id +'">' ;
        	//link +=  sectionObjs[i].label + '</a>';
        	var link = sectionObjs[i].label ;
		if(sectionObjs[i].access != '000'){
                	rows += '<tr height=20><td id="'+sectionObjs[i].id+'" class="'+cls+'" style="'+style+'">'+link+'</td></tr>';
        	}
    	}

	var s = 'width:15px;height:15px;font-size:16px;border:1px solid #fff;color:'+moduleMenuFg+';';
	s += 'padding:0px 0px 0px 0px;text-align:center;vertical-align:middle;cursor:hand;';

	var closeBtn = '<DIV id=closemodulemenu style="'+s+'">&times;</DIV>'	


    	return '<table width=100% height=100% style="border:1px solid #fff;">' + 
		'<tr height=50><td align=right valign=top>'+closeBtn+'</td></tr>' + 
		rows + 
		'<tr height=100%></td></tr>' +
		'</table>';

}


/////////////////////////
function getSearchBoxCn(){

        var searchType = $("input[name=searchtype]").val();

	var formTable = "";
	if (queryFormJson["nosearch"] == true){
		formTable = getSubtitleTable();	
	}
	else if (searchType == "advanced"){
		formTable = getAdvancedFormTable(NFIELDS);
	}
	else{
	 	formTable = getGenericFormTable();
        }
	return formTable;
}


///////////////////////////
function getSubtitleTable(){

	var table = '<table width=70% cellspacing=2 cellpadding=0 style="border:1px solid #fff;">';
        table += '<tr>';
        table += '<td style="padding:0px 0px 0px 20px;">'+moduleRelease+'</td>';
        table += '</tr>';
        table += '</table>';
	return table;

}


///////////////////
function fillStaticHtmlCn(fileName, containerId){

        
        var url =  fileName;
        var reqObj = new XMLHttpRequest();
        reqObj.containerId = containerId;
        reqObj.open("GET", url, true);
        reqObj.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
        reqObj.onreadystatechange = function() {
                if (reqObj.readyState == 4 && reqObj.status == 200) {
                        var htmlText = reqObj.responseText;
                        $(reqObj.containerId).html(htmlText);
                }
                else{
                        var msg = fileName + ' does not exist!';
                        var table = '<table width=100%>' +
                                '<tr height=400><td style="color:red;" align=center> ' + msg + '</td></tr>' +
                                '</table>';
                        $(reqObj.containerId).html(table);
                }
        };
        reqObj.send();
}


////////////////////////
function getGenericFormTable(){
        

        var fieldValueList = $('.valuetxtboxcls');
        if (fieldValueList.length > 0 && fieldValueList[0].value.trim() != ''){
                queryFormJson["boxlist"][boxListIndex][0]["value"] = fieldValueList[0].value;
        }

        var s1 = 'padding:10 5 10 10;font-size:16px;';
        var s2 = "font-size:16px;padding:0 15 5 10;";
        var s3 = 'padding:0 5 0 5;font-size:16px;';
        var s4 = 'padding:0 0 0 5;font-size:16px;cursor:hand;';
        var s5 = 'padding:0 5 0 3;font-size:16px;cursor:hand;';
        var s = 'width:5px;height:5px;border-left:1px solid #777;border-bottom:1px solid #777;';
        s += '-ms-transform: rotate(-45deg);-webkit-transform: rotate(-45deg);transform: rotate(-45deg);';
        var arrowDown = '<div style="'+s+'"></div>';

	var advancedCols = '<td id=advancedform style="'+s4+'" class=advancedform>Advanced</td>';
	advancedCols += '<td style="'+s5+'" class=advancedform>'+arrowDown+'</td>';
	advancedCols = (queryFormJson["advanced"] == true ? advancedCols : "");

	var table = '<table width=70% cellspacing=2 cellpadding=0 style="border:1px solid #fff;">';
	table += '<tr>';
	for (var i in queryFormJson["boxlist"][boxListIndex]){
		var s = 'width:'+queryFormJson["boxlist"][boxListIndex][i]["containerwidth"]+';';
		table += '<td style="'+s+'">' + getElement(queryFormJson["boxlist"][boxListIndex][i]) + '</td>' ;
	}
	table += advancedCols;
	table += '</tr>';
	table += '</table>';
	table += '<table width=70% cellspacing=2 cellpadding=0>';
	table += '<tr>';
        for (var i in queryFormJson["boxlist"][boxListIndex]){
                var s = 'padding:5px 0px 0px 0px;font-size:16px;color:#777;' + 
			'width:'+queryFormJson["boxlist"][boxListIndex][i]["containerwidth"]+';';
 		var boxLabel = (queryFormJson["boxlist"][boxListIndex][i]["boxlabel"] != undefined ? 
					queryFormJson["boxlist"][boxListIndex][i]["boxlabel"] : "");
		table += '<td style="'+s+'" valign=top>' + boxLabel + '&nbsp;</td>' ;
        }
	table += '</tr>';
        table += '</table>';
	table += '<input type=hidden class=fieldselectorcls value="genericfield">';
	table += '<input type=hidden name=searchtype value="generic"><br><br>' ;

	return table;
}


//////////////////
function getAdvancedFormTable(nfields){

        var style = 'width:15px;height:15px;font-size:16px;border:1px solid #ccc;color:#ccc;text-align:center;';
        style += 'cursor:hand';
        var closeIcon = '<div id=genericform style="'+style+'">&times;</div>';

        var style = 'width:15px;height:15px;font-size:16px;border:1px solid #ccc;color:#ccc;text-align:center;';
        style += 'cursor:hand';
        var addFieldIcon = '<div id=addfieldicon style="'+style+'">+</div>';
        var delFieldIcon = '<div id=delfieldicon style="'+style+'">-</div>';

        var fList = queryFieldList.split(",");
        var optHash = {};
	for (var i in fList){
		var val = fList[i].split("=")[0];
		var lbl = fList[i].split("=")[1];
		optHash[val] = lbl;
	}

	var jList = junctionList.split(",");
        var junHash = {};
        for (var i in jList){
                var val = jList[i].split("=")[0];
                var lbl = jList[i].split("=")[1];
                junHash[val] = lbl;
        }
	


	var s = (deviceType == "pc" ? 'width:80px;height:30px;' : 'width:60px;');
        var searchbtn = '<input type=submit class=searchbtn name=searchbtn style="'+s+'" value=" Search ">';
        var rows = '';
        for (var i=0; i < nfields; i++){
                var emStyle = (deviceType == "pc" ? "width:100%;height:30px;" : "width:100px;");

                var j = i%2;
                var selectedValue = $('select[name=field_'+i+']').val();
                selectedValue = (selectedValue == undefined ? fList[j].split("=")[0] : selectedValue);
                
		var emInfoObj = {"name":"field_"+i, "class":"fieldselectorcls", "value":selectedValue,
			  type:"select", style:emStyle, readonly:false, options:optHash
		};
		fieldInput = getElement(emInfoObj);

                var emStyle = (deviceType == "pc" ? 'width:100%;height:30px;padding:3;border:1px solid #fff;' :
                                                        'width:100px;padding:3;');
                var selectedValue = $('input[name=value_'+i+']').val();
                selectedValue = (selectedValue == undefined ? '' : selectedValue);
                var emInfoObj = {type:"text", style:emStyle, readonly:false};
                
		var emInfoObj = {"name":"value_"+i, "class":"valuetxtboxcls", "value":selectedValue,
                          type:"text", style:emStyle, readonly:false
                };
		valueInput = getElement(emInfoObj);

                var selectedValue = $('select[name=junction_'+i+']').val();
                selectedValue = (selectedValue == undefined ? '' : selectedValue);
                var emStyle = (deviceType == "pc" ? "width:60px;height:30px;" : "width:100px;");
                var emInfoObj = {"name":"junction_"+i, "class":"junctionselectorcls", "value":selectedValue,
			type:"select", style:emStyle, readonly:false, options:junHash
		};
                junctionInput = getElement(emInfoObj);

                junctionInput = (i == nfields -1 ? searchbtn : junctionInput);
                closeIcon = (i == 0 ? closeIcon : '');
                rows += '<tr>' +
                        '<td style="border:0px solid;" width=30%>&nbsp;<br>Field<br>'+fieldInput+'</td>' +
                        '<td style="border:0px solid;" width=50%>&nbsp;<br>Term<br>'+valueInput +'</td>' +
                        '<td style="border:0px solid;" nowrap>&nbsp;<br>&nbsp;<br>'+junctionInput+'</td>' +
                        '<td style="border:0px solid;" width=20% align=right valign=top>'+closeIcon+'</td>' +
                        '</tr>';
        }
        rows += '<tr><td colspan=2 style="border:0px solid;" align=right>' + addFieldIcon + '</td>' +
                '<td>'+delFieldIcon+'</td>' +
                '</tr>';
        var style = 'font-size:16px;border:1px solid #fff;padding:10px;';
        return '<table width=70% cellspacing=2 cellpadding=0 style="'+style+'" border=0>' + rows +  '</table>' +
                '<input type=hidden name=searchtype value="advanced">';

}




//////////////////////////////////////////////////////
function getBioxpressFormTableOne(){

        var transcriptType = $('select[name=transcripttype]').val();
        transcriptType = (transcriptType == undefined ? "featureType_mrna" : transcriptType);
        var emStyle = (deviceType == "pc" ? "width:200px;height:30px;" : "width:140px;");

        var slist = 'featureType_mrna=mRNA transcript, featureType_mirna=microRNA transcript';
        var emInfoObj = {"name":"transcripttype", "class":"searchElement1A", "value":transcriptType,
		type:"select", style:emStyle, readonly:false, 
		list:slist, onchange:"setSearchTypePage();"
	};
        var selector1 = getElement(emInfoObj);

        var slist = 'featureName=HGNC gene symbol,xrefSrc_uniProtAc=UniProt AC,xrefSrc_refSeqId=RefSeq AC';
        var emInfoObj = {"name":"searchfield1", "class":"searchElement1A","value":"",
		type:"select", style:emStyle, readonly:false, list:slist};
        var selector2A = getElement(emInfoObj);

        var slist = 'featureName=miRNA alias symbol,xrefSrc_refSeqId=RefSeq AC,xrefSrc_HGNCSymbol=HGNC miRNA symbol,xrefSrc_miRBaseId=miRBase AC,xrefSrc_Ensembl=Ensembl miRNA AC';
        var emInfoObj = {"name":"searchfield1", "class":"searchElement1A","value":"", 
		type:"select", style:emStyle, readonly:false, list:slist};
        var selector2B = getElement(emInfoObj);

        var s = (deviceType == "pc" ? 'width:300px;padding:3;height:30px;' : 'width:100px;padding:3;');
        var emInfoObj = {"name":"searchvalue1","class":"searchElement1A","value":"BRCA1",
		type:"text", style:s, readonly:false};
        var inBox1A = getElement(emInfoObj);
	var emInfoObj = {"name":"searchvalue1","class":"searchElement1A","value":"hsa-mir-21",
                type:"text", style:s, readonly:false};
        var inBox1B = getElement(emInfoObj);

        var bioxpressInfo = '  &nbsp; The BioXpress basic search is a gene or miRNA-centric search. ';
        var exampleA = 'Example: BRCA1, P38398, NP_009225';
        var exampleB = "Example: hsa-mir-21, NR_029493, MIR21, MI0000077, ENSG00000199004";

        var selector2 = (transcriptType == "featureType_mrna" ? selector2A : selector2B);
        var inBox1 = (transcriptType == "featureType_mrna" ? inBox1A : inBox1B);
        var example = (transcriptType == "featureType_mrna" ? bioxpressInfo + exampleA : bioxpressInfo + exampleB);
        var searchHGNC = '  &nbsp; Click ' + '<a href="http://www.genenames.org/cgi-bin/symbol_checker" target="_blank">here</a> ' +
                        'to find HGNC approved symbol for your gene of interest.';


	var s = (deviceType == "pc" ? 'width:80px;height:30px;' : 'width:60px;');
        var searchbtn = '<input type=submit class=searchbtncls name=searchbtn style="'+s+'" value=" Search ">';



        var s1 = 'padding:10 5 5 10;font-size:16px;';
        var s2 = "font-size:16px;padding:5 10 5 10;";
        var s3 = 'padding:10 5 5 10;font-size:16px;';

        var style = 'font-size:16px;border:1px solid #fff;color:'+moduleMenuFg+';';
        var formTable = '<table width=100% cellspacing=0 cellpadding=0 style="'+style+'">'+
                        '<tr height=30>'+
                        '<td nowrap style="'+s1 + '">&nbsp;Search for<br>'+ selector1+ '</td>' +
                        '<td style="'+s1+'">&nbsp;By<br>' + selector2+ '</td>' +
                        '<td style="'+ s1 + '" colspan=2>&nbsp;Query<br>'+ inBox1+'</td>' +
			'<td width=100%><br>'+searchbtn+'</td>' +
			'</tr>' +
                        '</table>' +
                        '<table width=100% cellspacing=0 cellpadding=0 style="color:'+moduleMenuFg+';">'+
                        '<tr>'+
                        '<td id=searchstatus style="'+s2+'" colspan=2>' + example+'<br><br>'+ searchHGNC +
                        '</td></tr>' +
                        '</table>';

        return formTable;
}

///////////////////////////////////////
function getBioxpressFormTableTwo(){

        var slist = 'sampleName=TCGA Cancer Type,sampleSynonym=Cancer Ontology (DOID)';
        var emStyle = (deviceType == "pc" ? "width:150px;height:30px;" : "width:120px;");
	var emInfoObj = {"name":"searchfield1", "class":"searchElement1A", "value":"",
		type:"select", style:emStyle, readonly:false, list:slist};
        var inBox1A = getElement(emInfoObj);

        var slist = 'transcriptType_mrna=mRNA,transcriptType_mirna=miRNA';
        var emStyle = (deviceType == "pc" ? "width:100px;height:30px;" : "width:120px;");
	var emInfoObj = {"name":"transcripttype", "class":"searchElement1A", "value":"",
		type:"select", style:emStyle, readonly:false, list:slist};
        var inBox1B = getElement(emInfoObj);

        var emStyle = (deviceType == "pc" ? "width:250px;height:30px;padding:0px 0px 0px 5px;" : "width:120px;");
	var emInfoObj = {"name":"searchvalue1", "class":"searchElement1A", "value":"BRCA",
		type:"text", style:emStyle, readonly:false};
        var inBox1C = getElement(emInfoObj);

        var slist = 'Expression_Up=Up,Expression_Down=Down,Expression_Both=Both';
        var emStyle = (deviceType == "pc" ? "width:70px;height:30px;" : "width:120px;");
	var emInfoObj = {"name":"searchfield3", "class":"searchElement1A", "value":"",
		type:"select", style:emStyle, readonly:false, list:slist};
        var inBox1E = getElement(emInfoObj);

        var slist = 'adjPValue_0.1=adjusted_p_value<0.1,adjPValue_0.05=adjusted_p_value<0.05,';
        slist += 'adjPValue_0.01=adjusted_p_value<0.01,adjPValue_0.001=adjusted_p_value<0.001,';
        slist += 'adjPValue_Manual=Manual Curation';
        
	var emStyle = (deviceType == "pc" ? "width:250px;height:30px;" : "width:120px;");
	var emInfoObj = {"name":"searchfield4", "class":"searchElement1A", "value":"",
		type:"select", style:emStyle, readonly:false, list:slist};
        var inBox1F = getElement(emInfoObj);


	var s = (deviceType == "pc" ? 'width:80px;height:30px;' : 'width:60px;');
        var searchbtn = '<input type=submit class=searchbtncls name=searchbtn style="'+s+'" value=" Search ">';




        var s1 = 'padding:10 5 5 10;font-size:11px;';
        var s2 = "font-size:10px;padding:5 10 5 10;";
        var s3 = 'padding:10 5 5 10;font-size:11px;';

        var example = 'Example: BRCA or DOID:1612 / Breast cancer [BRCA]';
        var style = 'font-size:13px;border:1px solid #fff;';
        var formTable = '<table width=100% cellspacing=0 cellpadding=0 style="'+style+'">'+
                        '<tr height=30 class=evenrow>'+
                        '<td nowrap style="'+s1 + '">&nbsp;Search for<br>'+ inBox1A+ '</td>' +
                        '<td style="'+s1+'">&nbsp;Transcript type<br>' + inBox1B+ '</td>' +
                        '<td style="'+s1+'">&nbsp;Query<br>' + inBox1C+ '</td>' +
                        '<td style="'+ s1 + '" >&nbsp;Trend<br>'+ inBox1E+'</td>' +
                        '<td style="'+ s1 + '" >&nbsp;Significance cutoff<br>'+ inBox1F+'</td>' +
			'<td width=100%><br>'+ searchbtn+'</td>' +
			'</tr>' +
                        '</table>' +
                        '<table width=100% cellspacing=0 cellpadding=0>'+
                        '<tr>'+
                        '<td id=searchstatus style="'+s2+'" colspan=2>' + example+'</td>' +
                        '</tr>'+
                        '</table>';

        return formTable;
}



///////////////////////
function getElement(emObj){

	var emName = ("name" in emObj ? emObj["name"] : "");
	var emValue = ("value" in emObj ? emObj["value"] : "");
	var emClass = ("class" in emObj ? emObj["class"] : "")

        emValue = (emValue == 'None' ? '' : emValue);
        var elm = '';
        if (emObj.type == 'textarea'){
                elm += '<textarea id="'+emObj.id+'" name="'+emName+'" class="'+emClass+'" rows="'+ emObj.rows +
                                '" cols="'+ emObj.cols +'" style="'
                elm +=  emObj.style +'">' + emValue+'</textarea>';
        }
        else if (emObj.type == 'select'){
                emObj.selected =  $('select[name='+emName+']').val();
		elm += '<select name="'+emName+'" class="'+emClass+'" style="'+emObj.style + '" ' ;
		elm +=  ' id="'+emObj.id+'" onchange="'+emObj.onchange+'">';
                for (var k in emObj["options"]){
                        var s = (k == emObj["selected"]  ? 'selected' : '');
                        elm += '<option value="'+k+'" '+s+'>'+emObj["options"][k]+'</option>';
                }
                elm += '</select>';
	}
        else if (emObj.type == 'checkbox'){
                var s = (emValue == 1 ? 'checked' : '');
                elm += '<input type="'+emObj.type+'" name="'+emName+'" '+s+' style="' +
                                        emObj.style+'" class="'+emClass+'" >';
        }
        else{
                elm += '<input ';
		for (var k in emObj){
			if(["containerwidth", "readonly"].indexOf(k) == -1){
				elm += k + '="' + emObj[k] + '" ';
			}
		}
		elm += '>';
	}
        return elm;

}


//////////////////////////////////
function rndrMiniTable(inObj){

	var style = "width:100%;font-size:16px;border:1px solid #ccc;"
	var table = '<table style="'+style+'" align=center cellspacing=1>';
	table += '<tr height=30>';
	for (var j in inObj["headers"]){
		var style = "text-align:center;font-weight:bold;padding:5px 0px 0px 0px;border:1px solid #ccc;"
		style += (j in inObj["colwidth"] ? "width:" + inObj["colwidth"][j] : "");
		table += '<td style="'+style+'" >'+inObj["headers"][j]+'</td>';
	}
	table += '</tr>';

	for (var i =0; i < inObj["content"].length; i ++){
		table += '<tr height=30>';
		for (var j in inObj["content"][i]){
                	var style = "text-align:center;padding:5px 0px 0px 0px;border:1px solid #ccc;"
			table += '<td style="'+style+'">'+inObj["content"][i][j]+'</td>';
       	 	}
        	table += '</tr>';
	}
	table += '</table>';
	return table;
}


////////////////////////////
function getWaitMsg(){

        var imgUrl =  '/imglib/loadingicon.gif';
        var imgObj = '<img src="'+imgUrl+'">';
        var cn = '<table width=100% height=400><tr><td valign=middle align=center>'+imgObj+'</td></tr></table>';
        return cn;
}


////////////////////////////
function getErrorMsg(msg){
        var cn = '<table width=100% height=200><tr><td valign=middle align=center style="color:red;">' +
                        msg + '</td></tr></table>';
        return cn;
}









/////// Event Handlers ///////////////////
//$(document).on('click', ':input', function (event) {
//        event.preventDefault();
//	handleEvents({"emclass":this.className, "emid":this.id});
//});


///////////////////////////////////////////
$(document).on('click', '.faqquestion', function (event) {
      	event.preventDefault();
      	var jqId = "#answer_" + this.id.split("_")[1];
	$(jqId).toggle();

});


//////////////////////////////////////
$(document).on('click', '.moduleiconcn', function (event) {
        event.preventDefault();
        handleEvents({"emclass":this.className, "emid":this.id});
});


//////////////////////////////////////
$(document).on('click', '#closemodulemenu', function (event) {
        event.preventDefault();
        handleEvents({"emclass":this.className, "emid":this.id});
});




/////////////////////////////////
function handleEvents(inObj){


	//$('html').animate({scrollTop:0}, 'fast');
        //$('body').animate({scrollTop:0}, 'fast');
        
	if (inObj["emclass"] == "nav-link"){
        	gpageId = inObj["emid"];
		setColumnOneCn();
        	fillStaticHtmlCn('/content/page.'+gpageId+'.html', '#clmn2cn');
	}
	else if (inObj["emclass"] == "searchbtncls"){
		var cn = '<table width=100%><tr height=200><td align=center>';
		cn += 'This functionality is still under development. Click ';
		cn += '<a href="/data">here</a> to return to the main page.'
		cn += '</td></tr></table>';
		$("#pagecn").html(cn);	
	}
	else if (inObj["emclass"] == "advancedform"){
		$("input[name=searchtype]").val("advanced");
	        NFIELDS = 1;
        	$("#modulesearchboxcn").html(getSearchBoxCn())
	}
	else if (inObj["emclass"] == "genericform"){
        	$("input[name=searchtype]").val("generic");
        	$("#modulesearchboxcn").html(getSearchBoxCn())
	}
	else if (inObj["emclass"] == "faqquestion"){
        	var jqId = "#answer_" + this.id.split("_")[1];
        	$(jqId).toggle();
        }
	else if (inObj["emclass"] == "moduleiconcn"){
		$("#modulesectionscn").toggle();
        }
	
	if (inObj["emid"] == "addfieldicon"){
		NFIELDS += 1;
	        NFIELDS = (NFIELDS > 4 ? 4 : NFIELDS);
        	$("#modulesearchboxcn").html(getSearchBoxCn());
	}
	else if (inObj["emid"] == "delfieldicon"){
		NFIELDS -= 1;
	        NFIELDS = (NFIELDS < 1 ? 1 : NFIELDS);
        	$("#modulesearchboxcn").html(getSearchBoxCn());
	}
	else if (inObj["emid"] == "closemodulemenu"){
                $("#modulesectionscn").css("display", "none");
        }
}



///////////////////////////////////////////////////
function getUrlParam(name){
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
}



