$(document).ready( function() { 
	var filelocation = 'Files/'+baselang+'/'+baselang.toLowerCase(); // example: Files/Plot.en/plot where plot is start of text file name
	var urlAnchor = window.location.hash; // get hash tag from url
	$("#Examples .expanded .expanded").each(function(i){ // append ids to inner divs, open and populate contents on click
		var expId = "#exp" + i;
		var subOuterId = "#subOuter" + i;
		var expLink = subOuterId + " h3 .expand";
		var file = filelocation + (i)+'.html'; // example: Files/Plot.en/plot0.html
		var files = [];
		var divContent = "";
		var numOfSubOuters = $('#Examples .expanded .expanded').length - 1;
		var last = "#exp"+numOfSubOuters;
		
		$(this).attr('id', 'exp' + i);
		$(this).parents('.expandSubOuter').attr('id', 'subOuter' + i);
		files[i] = filelocation + (i) + '.html';
		if(i != 0) {
			$(expLink).click(function(){
				if(!$(expId).html()) {$(expId).load(file);}
			});
		}
		function anchorload(anchor, link) {
			if (divContent.indexOf(anchor) >= 0) { // if data contains id with the anchor value minus # symbol
				var expHeading = "#" + $(expId).parents('.expandSubOuter').attr('id');
				$(expId).parents('.expandOuter').children('.expanded').removeClass('hide');
				$(expId).parents('.expandSubOuter').children('.expanded').removeClass('hide');
				$(expId).removeClass('hide');
				$(expHeading + ' h3 .expand').addClass('hide');
				$(expHeading + ' h3 .collapse').removeClass('hide');
				if($('.expandSubsubOuter').length){
					$('.expandSubsubOuter .expanded').each(function(a){
						var ssData = $(this).html();
						if(ssData.indexOf(anchor) >= 0) {
							var subExpHeading = "#" + $(this).parents('.expandSubsubOuter').attr('id');
							var subExp = "#" + $(this).attr('id');
							$(subExpHeading + ' h4 .expand').addClass('hide');
							$(subExpHeading + ' h4 .collapse').removeClass('hide');
							$(subExp).removeClass('hide');
						}
					});
				}
				if($('.expandSubsubsubOuter').length) {
					$('.expandSubsubsubOuter .expanded').each(function(b){
						var sssData = $(this).html();
						if(sssData.indexOf(anchor) >= 0) {
							var ssExpHeading = "#" + $(this).parents('.expandSubsubsubOuter').attr('id');
							var ssExp = "#" + $(this).attr('id');
							$(ssExpHeading + ' h5 .expand').addClass('hide');
							$(ssExpHeading + ' h5 .collapse').removeClass('hide');
							$(ssExp).removeClass('hide');
						}
					});
				}
				window.location = link;
			} 
		};
		if(urlAnchor !="") { // if anchor shows up in the url, then run this
			var link = urlAnchor;		
			var anchor = '"'+ link.replace('#','')+'"';
			if (i != 0 && !$(expId).html()) {
				$(expId).load(files[i], function(data){ // get data from each file that hasn't been opened
					divContent = data;
					anchorload(anchor, link);
				});
			} 
		}
		$('a.ExampleButtonLink').click(function() {
			var link = $(this).attr('href');
			var anchor = '"'+ link.replace('#','')+'"';
			if (!$(expId).html() && i != 0) {
				$(expId).load(files[i], function(data){ // get data from each file that hasn't been opened
					divContent = data;
					anchorload(anchor, link);
				});
			} else {
				divContent = $(expId).html();
				anchorload(anchor, link);
			}
		});
		$('a.opener').click(function() { 
			if(!$(expId).html()) { $(expId).load(files[i]); }
		});
		$(last).css({'margin-bottom':'-20px'});
	});	
});