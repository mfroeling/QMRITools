function setSearchTextField(paramname, field) {
	var passed = location.search.substring(1);
	var query = getParm(passed,paramname);
	var query = getParm(passed,paramname);
	query = query.replace(/\+/g," ");
	var loc = document.location;
	if(/.*search.html/.test(loc)) {
	 	document.title = decodeURIComponent(query) + ' - Wolfram Search';
	}
	field.value = decodeURIComponent(query);
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
	return '';
}
