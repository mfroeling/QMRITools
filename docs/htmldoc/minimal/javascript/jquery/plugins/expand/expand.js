$(document).ready(function() {
	// firefox hack
	if ($.browser.mozilla && self.document.location.hash != "") {
		$.scrollTo(self.document.location.hash);
	}
	
	$(".featureDesc .expand").click(function(e) {
		e.preventDefault();
		$(this).parents(".featureCont").children('.expanded').removeClass('hide');
        $(this).parent().children('.collapse').removeClass('hide');
        $(this).addClass('hide');
		return false;
	});
	
	$(".featureDesc .collapse").click(function(e) {
		e.preventDefault();
		$(this).parents(".featureCont").children('.expanded').addClass('hide');
        $(this).parent().children('.expand').removeClass('hide');
        $(this).addClass('hide');
		return false;
	});
	
});
