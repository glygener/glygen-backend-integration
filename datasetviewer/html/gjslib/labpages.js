$(document).on('click', '.navlinkmobile', function (event) {
        event.preventDefault();
        var emId = this.id;
        var emList = $(".sliderbox");
        var offset = 30;
        var boxWidth = 50;
        var leftJump = boxWidth + 5;
        var nBoxes = emList.length;
        var nVisibleBoxes = 5;
        var maxLeft = (nBoxes - 1)*leftJump + offset;
        var maxVisibleLeft = (nVisibleBoxes - 1)*leftJump + offset;
        leftJump = (emId == 'navleft' ? "+=" + leftJump : "-=" + leftJump)
        for (var i=0; i < emList.length; i++){
                var jqId = '#' + emList[i].id;
                $(jqId).animate({left: leftJump}, function() {
                        var left = parseInt($(this).css("left").replace("px", ""));
                        if(left < 0){
                                $(this).css("left", maxLeft);
                        }
                        if(left > maxLeft){
                                $(this).css("left", offset);
                                left = parseInt($(this).css("left").replace("px", ""));
                        }
                        $(this).css("display", "block");
                        if(left < 0 || left > maxVisibleLeft){
                                $(this).css("display", "none");
                        }
                });
        }
});

$(document).on('click', '.navlink', function (event) {
      	event.preventDefault();
      	var emId = this.id;
	var emList = $(".sliderbox");
	var offset = 25;
	var boxWidth = 100;
      	var leftJump = boxWidth + 10;
        var nBoxes = emList.length;
	var nVisibleBoxes = 7;
	var maxLeft = (nBoxes - 1)*leftJump + offset;
	var maxVisibleLeft = (nVisibleBoxes - 1)*leftJump + offset;
	leftJump = (emId == 'navleft' ? "+=" + leftJump : "-=" + leftJump)
	for (var i=0; i < emList.length; i++){
		var jqId = '#' + emList[i].id;
		$(jqId).animate({left: leftJump}, function() {
			var left = parseInt($(this).css("left").replace("px", ""));
			if(left < 0){
				$(this).css("left", maxLeft);
			}
			if(left > maxLeft){
				$(this).css("left", offset);
                                left = parseInt($(this).css("left").replace("px", "")); 
			}
			$(this).css("display", "block");
			if(left < 0 || left > maxVisibleLeft){
				$(this).css("display", "none");
			}
		});
	}
});

$(document).on('click', '.imglinkmobile', function (event) {
        event.preventDefault();
        var emId = this.id;



        var imgUrl = htmlRoot + '/content-mobile/' + emId;
        var parts = emId.split(".");
        var titledivId = "page_" + parts[1] + "_slidertitle_" + parts[3];
        var jqId = '#' + titledivId ;
        $("#selectedslidertitle").html($(jqId).html());
        var textdivId = "page_" + parts[1] + "_slidertext_" + parts[3];
        var jqId = '#' + textdivId ;
        $("#selectedslidertext").html($(jqId).html());
        $(".sliderbox").css("background", "#fff");
        var jqId = '#sliderbox_' + parts[3];
        $(jqId).css("background", "#000");
        $("#selectedsliderimage").html('<img src="'+imgUrl+'" width=100% height=150>');
        if(emId == 'grid11.left'){
                var emList = $("#grid11").children();
                for (var i=0; i < emList.length; i++){
                        if(emList[i].className == 'slidermovingbox'){
                                alert(i + ' ' + emId + ' ' + emList[i].className);
                        }
                }
        }
});

$(document).on('click', '.imglink', function (event) {
        event.preventDefault();
        var emId = this.id;
	var imgUrl = htmlRoot + '/content/' + emId;
	var parts = emId.split(".");
	var titledivId = "page_" + parts[1] + "_slidertitle_" + parts[3];
        var jqId = '#' + titledivId ;
        $("#selectedslidertitle").html($(jqId).html());
	var textdivId = "page_" + parts[1] + "_slidertext_" + parts[3];
	var jqId = '#' + textdivId ;
	$("#selectedslidertext").html($(jqId).html());
	$(".sliderbox").css("background", "#fff");
	var jqId = '#sliderbox_' + parts[3];
	$(jqId).css("background", "#000");
	$("#selectedsliderimage").html('<img src="'+imgUrl+'" width=100% height=200>');
        if(emId == 'grid11.left'){
                var emList = $("#grid11").children();
                for (var i=0; i < emList.length; i++){
                        if(emList[i].className == 'slidermovingbox'){
                                alert(i + ' ' + emId + ' ' + emList[i].className);
                        }
                }
        }
});


$(document).on('click', '#gmenuicon', function (event) {
        event.preventDefault();
        $("#gsectionscn").toggle();

});

  
