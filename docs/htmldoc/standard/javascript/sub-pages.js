$(document).ready( function() {
	// subsubs
	$(".expandSubsubInner a.expand").click(function(e) {
		e.preventDefault();
		$(this).parents(".expandSubsubOuter").children('.expanded').removeClass('hide');
        $(this).parent().children('.collapse').removeClass('hide');
        $(this).addClass('hide');
		return false;
	});
	$(".expandSubsubInner a.collapse").click(function(e) {
		e.preventDefault();
		$(this).parents(".expandSubsubOuter").children('.expanded').addClass('hide');
        $(this).parent().children('.expand').removeClass('hide');
        $(this).addClass('hide');
		return false;
	});

	// subsubsubs
	$(".expandSubsubsubInner a.expand").click(function(e) {
		e.preventDefault();
		$(this).parents(".expandSubsubsubOuter").children('.expanded').removeClass('hide');
        $(this).parent().children('.collapse').removeClass('hide');
        $(this).addClass('hide');
		return false;
	});
	$(".expandSubsubsubInner a.collapse").click(function(e) {
		e.preventDefault();
		$(this).parents(".expandSubsubsubOuter").children('.expanded').addClass('hide');
        $(this).parent().children('.expand').removeClass('hide');
        $(this).addClass('hide');
		return false;
	});
	$(".expanded h6:first").addClass('first');
	
	// append ids to subsub and subsubsub sections
	$('.expandSubsubOuter .expanded').each(function(i){
		$(this).attr('id', 'subExp' + i);
		$(this).parents('.expandSubsubOuter').attr('id', 'ssOuter' + i);
		if($('a.opener').hasClass('hide')) {
			$(this).removeClass('hide');
			$('.expandSubsubOuter .expand').addClass('hide');
			$('.expandSubsubOuter .collapse').removeClass('hide');
		}
	});
	$('.expandSubsubsubOuter .expanded').each(function(i){
		$(this).attr('id', 'ssExp' + i);
		$(this).parents('.expandSubsubsubOuter').attr('id', 'sssOuter' + i);
		if($('a.opener').hasClass('hide')) {
			$(this).removeClass('hide');
			$('.expandSubsubsubOuter .expand').addClass('hide');
			$('.expandSubsubsubOuter .collapse').removeClass('hide');
		}
	});
});