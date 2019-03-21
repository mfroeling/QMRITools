function setSearchTextField(paramname, field) {
	var passed = location.search.substring(1);
	var query = getParm(passed,paramname);
	if(typeof query != 'undefined') {
		query = query.replace(/\+/g," ");
		var loc = document.location;
		if(/.*documentation-search.html/.test(loc)) {
			document.title = decodeURIComponent(query) + ' - Wolfram Search';
		}
		$(field).attr('value', decodeURIComponent(query));
		$('#query').attr('class', 'searchboxsub-on');
	} 	
}

function getParm(string,parm) {
	// returns value of parm from string
	var startPos = string.indexOf(parm + "=");
	if (startPos > -1) {
	 	startPos = startPos + parm.length + 1;
	 	var endPos = string.indexOf("&",startPos);
	 	if (endPos == -1)
	 	 	endPos = string.length;
	 	return string.substring(startPos,endPos);
	}
}
/* new */
function searcher(query) {
	$.ajax({ 
		url: '/lucene/DirectHit?query=' + encodeURIComponent(query),
		type: 'GET',
		complete: searchComplete
	});
    return false;
}

function searchertwo(query, collection, undefined) {
	if(collection == undefined) {
		collection = 'system_modeler';
	}
	$.ajax({ 
		url: '/lucene/DirectHit?query=' + encodeURIComponent(query)  + '&collection=' + collection,
		type: 'GET',
		complete: searchComplete
	});
    return false;
}
function langsearcher(query, lang){
	$.ajax({ 
		url: '/lucene/DirectHit?query=' + encodeURIComponent(query) + '&lang=' + lang,
		type:'GET', 
		complete: searchComplete
	});
	return false;
}
function searchComplete(req) {
	var url = req.responseText;
	if (url.indexOf('/search.html') >= 0) {
		url = url.replace('/search.html', '/documentation-search.html');
	}
	window.location = url;
}
$(document).ready(function(){
	$('.search-wolfram-results-container a').live("click",function(e){
		e.preventDefault();
		var url = $(this).attr('href');
		if(url.indexOf('documentation-search.html') < 0) {
			url = url.replace('search.html', 'documentation-search.html');
		}
		window.location = url;
	});
	// placeholder for all browsers
	$('#query').keypress(function() {
			$('.placeholder').removeClass('show');
		}).blur(function() {
		  var input = $(this);
		  if(input.val() == '') {
			$('.placeholder').addClass('show');
			//$("#query").focus();
		  }
		});
	$('#query').bind('paste', function(){
			$('.placeholder').removeClass('show');
		}).blur(function() {
		  var input = $(this);
		  if(input.val() == '') {
			$('.placeholder').addClass('show');
			//$("#query").focus();
		  }
	});

	// have searchbar be on focus and cursor at the beginning
	var queryVal = $("#query").val();
	var bodyID = $('body').attr('id');
	//$('#query').change(function() {
	if(queryVal == "") {
		$('.placeholder').addClass('show');
	}
	if(queryVal == "" && (bodyID == "rootGuide" || bodyID == "languageRootGuide")) {
		//$("#query").focus();
	}
	//});
	// search page text and magnifying glass
	$('.search-all img').attr('src','../images/mathematica/small-magnifying-glass.png');
	$('#nextdemo').before('<strong>Next</strong> &raquo;');
	$('#nextdemo').css('display','none');
	if($('#prevdemo').length > 0 ){
		$('#prevdemo').before('&laquo; <strong>Prev</strong>');
		$('#prevdemo').css('display','none');
		$('#prevnextsep').before('|');
		$('#prevnextsep').css('display','none');
	}
});
