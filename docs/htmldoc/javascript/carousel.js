$(document).ready(function() {
// carousel on load
	$('#mycarousel').tinycarousel({ start: 1, display: 4, animation: true });
	$('#mycarousel-minimize').css({'display':'none'});
// show all view
	$('#mycarousel-show-all').click(function(){
		$('#mycarousel').tinycarousel({ start: 1, display: 4, animation: true });
		$('#mycarousel').stop().tinycarousel();
		$('#mycarousel .overview').css({'width':'630px', 'height':'auto', 'left':'0px'});
		$('#mycarousel .overview').stop();
		$('.page-out-of, .divider, .total-pages ').css({'display':'none'});
		$('.viewport').css({'height':(($('.overview').height()))+10+'px'});
		$('#mycarousel').css({'overflow':'auto', 'height':'auto'});
		$('#mycarousel-show-all').css({'display':'none'});
		$('#mycarousel-minimize').css({'display':'inline-block'});
		$('#mycarousel .prev, #mycarousel .next').css({'display':'none'}); 
		$('.disable-next, .disable-prev').remove();
		return false;
	});
// minimized view
	$('#mycarousel-minimize').click(function(){	
		$('#mycarousel').tinycarousel({ start: 1, display: 4, animation: true });
		$('.page-out-of, .divider, .total-pages').css({'display':'inline-block'});
		clicks = 1;
		$('.page-out-of').html(clicks);
		$('.overview, .viewport').css({'height':'170px'});
		$('#mycarousel').css({'overflow':'hidden', 'height':'215px'});
		$('#mycarousel-minimize').css({'display':'none'});
		$('#mycarousel-show-all').css({'display':'inline-block'});
		$('#mycarousel .next, #mycarousel .prev').css({'display':'inline-block'});
		return false;
	});
// minimize entire carousel section
	$('#FeaturedExamples .collapse').click(function(){ 
		$('.viewport, .overview, .buttons, .jcarousel-scroll').css({'display':'none'});
		$('#mycarousel').css({'height':'32px'});
	});
// expand back to carousel view
	$('#FeaturedExamples .expand').click(function(){
		$('.viewport, .overview').css({'display':'block'});
		$('.buttons, .divider, .page-out-of, .total-pages, .jcarousel-scroll, #mycarousel-show-all').css({'display':'inline-block'});
		clicks = 1;
		$('.page-out-of').html(clicks);
		$('#mycarousel').tinycarousel({ start: 1, display: 4, animation: true });
		$('#mycarousel .overview').css('left','0px');
		$('#mycarousel .overview').stop();
		$('.overview, .viewport').css({'height':'170px'});
		$('#mycarousel').css({'overflow':'hidden'});
		$('#mycarousel').css({'height':'215px'});
		$('#mycarousel-minimize').css({'display':'none'});		
		$('.disable-next, .disable-prev').remove();
		$('.page-out-of').before('<span class="disable-prev"></span>');
	});
// find total number of pages
	function numPages() {
		var numPages = Math.ceil($('.overview li').length/4);
		return numPages;
	};
	var totalPages = numPages();
	$('.total-pages').text(totalPages);
	
// click counter for current page
	var clicks = 1; 
	$('.page-out-of').html(clicks);
	$(".next").click(function(){ 
		clicks++; 
		$('.page-out-of').html(clicks);
	});
	$(".prev").click(function(){
		clicks--; 
		$('.page-out-of').html(clicks);
	});
	
// append disable-prev and disable-next
	if(clicks == 1) {
		$('.page-out-of').before('<span class="disable-prev"></span>');
	} 
	$(".buttons").click(function(){
		if(clicks == 1) {
			$('.page-out-of').before('<span class="disable-prev"></span>');
			$('.disable-next').remove();
		} 
		if(clicks == totalPages) { 
			$('.total-pages').after('<span class="disable-next"></span>');
			$('.disable-prev').remove();
		} 
		if (clicks < totalPages && clicks > 1) { 
			$('.disable-next').remove();
			$('.disable-prev').remove();
		}
	});
	$('.overview a').each(function(){
		$(this).attr('title',$('.carousel-caption',this).text());
	});
// turncate text below images after two lines	
	$('.carousel-caption').dotdotdot({
		height: 28,
		tolerance: 0,
		wrap: 'letter'
	});
});