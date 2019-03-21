function getParameterByName(name) {
	var list = getParameterListByName(name);

	if(list.length > 0) {
		return list[0];
	} 
	
	return "";
}

function getParameterListByName(name) {
	name = name.replace(/[\[]/, "\\\[").replace(/[\]]/, "\\\]");
	var regexS = "[\\?&]" + name + "=([^&#]*)";
	var regex = new RegExp(regexS);
	var matches = window.location.search.match(new RegExp(regexS, 'g'));
	var results = [];
	
	if(matches == null) {
		return [];
	}
	
	for(var i = 0; i < matches.length; i++) {
		var match = regex.exec(matches[i]);
		
		if(match != null) {
			if(match[1]) {
				results.push(decodeURIComponent(match[1].replace(/\+/g, " ")));
			}
		}
	}
	
	return results;
}
