$(document).ready( function() {

	// if javascript is disabled everything is expanded because USABILITY
	// but some things are closed by default
	$('.closed .expanded').addClass('hide');
	$('.closed .expand').removeClass('hide');
	$('.closed .collapse').addClass('hide');

	// do the expanding/collapsing
	$(".expandInner .expand").click(function(e) {
		if($('.expand a').attr('href')) {
			if($('.expand a').attr('href') !="") {
				$('.expand a').click( function() { window.location = $(this).attr('href') });
			}
		}
		else e.preventDefault();
		$(this).parents(".expandOuter").children('.expanded').removeClass('hide');
        $(this).parent().children('.collapse').removeClass('hide');
        $(this).addClass('hide');
		return false;
	});
	$(".expandInner .collapse").click(function(e) {
		if($('.collapse a').attr('href')) {
			if($('.collapse a').attr('href') !="") {
				$('.collapse a').click( function() { window.location = $(this).attr('href') });
			}
		}
		else e.preventDefault();
		$(this).parents(".expandOuter").children('.expanded').addClass('hide');
        $(this).parent().children('.expand').removeClass('hide');
        $(this).addClass('hide');
		return false;
	}); 

	// subs
	$(".expandSubInner .expand").click(function(e) {
		if($('.expand a').attr('href')) {
			if($('.expand a').attr('href') !="") {
				$('.expand a').click( function() { window.location = $(this).attr('href') });
			}
		}
		else e.preventDefault();
		$(this).parents(".expandSubOuter").children('.expanded').removeClass('hide');
        $(this).parent().children('.collapse').removeClass('hide');
        $(this).addClass('hide');
		return false;
	});
	$(".expandSubInner .collapse").click(function(e) {
		if($('.collapse a').attr('href')) {
			if($('.collapse a').attr('href') !="") {
				$('.collapse a').click( function() { window.location = $(this).attr('href') });
			}
		}
		else e.preventDefault();
		$(this).parents(".expandSubOuter").children('.expanded').addClass('hide');
        $(this).parent().children('.expand').removeClass('hide');
        $(this).addClass('hide');
		return false;
	});

	// subsubs
	$(".expandSubsubInner .expand").click(function(e) {
		if($('.expand a').attr('href')) {
			if($('.expand a').attr('href') !="") {
				$('.expand a').click( function() { window.location = $(this).attr('href') });
			}
		}
		else e.preventDefault();
		$(this).parents(".expandSubsubOuter").children('.expanded').removeClass('hide');
        $(this).parent().children('.collapse').removeClass('hide');
        $(this).addClass('hide');
		return false;
	});
	$(".expandSubsubInner .collapse").click(function(e) {
		if($('.collapse a').attr('href')) {
			if($('.collapse a').attr('href') !="") {
				$('.collapse a').click( function() { window.location = $(this).attr('href') });
			}
		}
		else e.preventDefault();
		$(this).parents(".expandSubsubOuter").children('.expanded').addClass('hide');
        $(this).parent().children('.expand').removeClass('hide');
        $(this).addClass('hide');
		return false;
	});

	// subsubsubs
	$(".expandSubsubsubInner .expand").click(function(e) {
		if($('.expand a').attr('href')) {
			if($('.expand a').attr('href') !="") {
				$('.expand a').click( function() { window.location = $(this).attr('href') });
			}
		}
		else e.preventDefault();
		$(this).parents(".expandSubsubsubOuter").children('.expanded').removeClass('hide');
        $(this).parent().children('.collapse').removeClass('hide');
        $(this).addClass('hide');
		return false;
	});
	$(".expandSubsubsubInner .collapse").click(function(e) {
		if($('.collapse a').attr('href')) {
			if($('.collapse a').attr('href') !="") {
				$('.collapse a').click( function() { window.location = $(this).attr('href') });
			}
		}
		else e.preventDefault();
		$(this).parents(".expandSubsubsubOuter").children('.expanded').addClass('hide');
        $(this).parent().children('.expand').removeClass('hide');
        $(this).addClass('hide');
		return false;
	});
	$(".expanded h6:first").addClass('first');

	// show/hide all examples
	$(".expandInner .opener").click(function(e) {
		e.preventDefault();
		$('#Examples .expanded .collapse').removeClass('hide');
		$('#Examples .expanded .expand').addClass('hide');
		$('#Examples .expanded .expanded').removeClass('hide');
        $(this).parent().children('.closer').removeClass('hide');
        $(this).addClass('hide');
		return false;
	});
	$(".expandInner .closer").click(function(e) {
		e.preventDefault();
		$('#Examples .expanded .collapse').addClass('hide');
		$('#Examples .expanded .expand').removeClass('hide');
		$('#Examples .expanded .expanded').addClass('hide');
        $(this).parent().children('.opener').removeClass('hide');
        $(this).addClass('hide');
		return false;
	});
	
	$('.overviews').last().css('margin-bottom','35px');
	
	// dropdown menus
    $("ul.dropdown li").hover(function(){
		$(this).addClass("hover");
		var mainWrapperHeight = $('.main-wrapper').height(); // height of main wrapper
		var dropdownHeight = $('li.hover ul').height()+291; // height of selected dropdown
        $('ul:first',this).css({'display':'block','visibility':'visible'});
		if (dropdownHeight > mainWrapperHeight) {  
			$('.main-wrapper').height(dropdownHeight);
		}
    }, function(){
        $(this).removeClass("hover");
        $('ul:first',this).css({'display':'none','visibility':'hidden'});
		$('.main-wrapper').css('height', 'auto');
    });
    $("ul.dropdown li ul li:has(ul)").find("a:first").append(" &raquo; ");
	
	
	// feedback
	var feedbackformhtml = $('#feedbackForm').html();
	var url = document.URL;
	$('#urlLabel').text("URL: " + url);
	$('#urlLabel').attr('title', url);
	$('#referer').attr('value', url);
	
	// feedback toggle
	$('.footer-give-feedback a').click(function(e) {
		e.preventDefault();
		$('#feedbackForm').html(feedbackformhtml);
		var url = document.URL;
	    $('#urlLabel').text("URL: " + url);
	    $('#urlLabel').attr('title', url);
	    $('#referer').attr('value', url);
		$('#feedbackForm').slideToggle().toggleClass('hide');
		$('div.footer-give-feedback').toggleClass('feedback-open');
		$('.footer-give-feedback a').toggleClass('close');
		$('.mainContent').toggleClass('feedback');
		$('div.largebutton').hover( function(){$(this).addClass('hover')}, function(){$(this).removeClass('hover')});
		return false;
	});
	$('.footer-give-feedback a').toggle(
		function() { $(this).text('') },
		function() { $(this).text('意见反馈') }
	);	
	
	// feedback submit
	$('div#submit').live("click", function() {
	 
	var feedbackval = $('#feedbackMessage').val();
	
	if(feedbackval === null || feedbackval === '')
	{
	   $('#feedbackMessageTable').addClass('errorHighlight');
	   $('#feedbackForm').errorBox();
	}
	else {
			$('div#thank_you').removeClass('hide');
			$('table#formTable').addClass('hide');
			$('#feedbackMessageTable').removeClass('errorHighlight');
			$.post("../includes/reference-feedback.cgi", { 
			       feedback: $('#feedbackMessage').val(), 
			       name: $('input#name').val(),
			       email: $('input#email').val(),
			       url: document.URL
			} );
		}
		});
	
	// Mathematica 9 is now available buttons hover
	$('div.button').hover( function(){$(this).addClass('hover')}, function(){$(this).removeClass('hover')});
	$(".podheader").click(function() {
	    $("+ .links-list", this).toggle();
	    $(".links-list").not($("+ .links-list", this)).hide();
	});
});

//FLASH MANIPULATES

function swap(obj, theWidth, theHeight, fileName, divId){
var flash1='<object classid="clsid:d27cdb6e-ae6d-11cf-96b8-444553540000" codebase="http://fpdownload.macromedia.com/pub/shockwave/cabs/flash/swflash.cab#version=7,0,0,0" width="'+ theWidth + '" height="'+ theHeight + '" id="benefits" align="middle"><param name="allowScriptAccess" value="sameDomain" /><param name="movie" value="'+ fileName + '" /><param name="loop" value="false" /><param name="menu" value="false" /><param name="quality" value="high" /><param name="bgcolor" value="#ffffff" /><embed src="'+ fileName + '" loop="false" menu="false" quality="high" bgcolor="#ffffff" width="'+ theWidth + '" height="'+ theHeight + '" name="benefits" align="middle" allowScriptAccess="sameDomain" type="application/x-shockwave-flash" pluginspage="http://www.macromedia.com/go/getflashplayer" /></object>';

	document.getElementById(divId).innerHTML = flash1;
}



