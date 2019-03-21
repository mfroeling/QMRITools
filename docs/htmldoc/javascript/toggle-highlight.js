$(document).ready(function() {
	$('a.highlight-link').click(function() {
	  $('highlighting').not(this).removeClass('highlighting');
	  $(this).toggleClass('highlighting');
	});
	$('a.highlight-link').toggle(
			function() {
				$('.modified-text').css({'background':'#fffcb8'});
				$('h2.modified-text').css({'border-width':'1px 0 1px 0','border-color':'#e3e3e3','border-style':'solid','height':'30px','line-height':'30px'});
				$('h2.modified-text a.expand').css({'background-position':'28px 8px'});
				$('h2.modified-text a.collapse').css({'background-position':'25px 10px'});
			},
			function() {
				$('.modified-text').css({'background':''});
				$('h2.modified-text').css({'border-width':'','border-color':'','border-style':'','line-height':'','height':''});
				$('h2.modified-text a').css({'background-position':''});
			}
	);
});