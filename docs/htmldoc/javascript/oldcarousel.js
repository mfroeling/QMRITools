$(document).ready(function() {
// carousel on load
	$('#mycarousel').tinycarousel({ display: 4 });
	$('#mycarousel-minimize').css({'display':'none'});
// show all view
	$('#mycarousel-show-all').click(function(){ 
		$('.page-out-of').css({'display':'none'});
		$('.divider').css({'display':'none'});
		$('.total-pages').css({'display':'none'});
		$('#mycarousel').tinycarousel({ start: 1 });
		$('.overview').css({'width':'630px'});
		$('.overview').css({'height':'auto'});
		$('.viewport').css({'height':(($('.overview').height()))+10+'px'});
		$('#mycarousel').css({'overflow':'auto'});
		$('#mycarousel').css({'height':'auto'});
		$('#mycarousel-show-all').css({'display':'none'});
		$('#mycarousel-minimize').css({'display':'inline-block'});
		$('#mycarousel .prev').css({'display':'none'}); 
		$('#mycarousel .next').css({'display':'none'});
		$('.disable-next').remove();
		$('.disable-prev').remove();
	});
// minimized view
	$('#mycarousel-minimize').click(function(){ 
		$('#mycarousel').tinycarousel({ start: 1, display: 4 });
		$('.page-out-of').css({'display':'inline'});
		$('.divider').css({'display':'inline'});
		$('.total-pages').css({'display':'inline'});
		clicks = 1;
		$('.page-out-of').html(clicks);
		$('.overview').css({'height':'170px'});
		$('.viewport').css({'height':'170px'});
		$('#mycarousel').css({'overflow':'hidden'});
		$('#mycarousel').css({'height':'215px'});
		$('#mycarousel-minimize').css({'display':'none'});
		$('#mycarousel-show-all').css({'display':'inline-block'});
		$('#mycarousel .next').css({'display':'inline-block'});
		$('#mycarousel .prev').css({'display':'inline-block'});
	});
// minimize entire carousel section
	$('.collapse').click(function(){ 
		$('.viewport').css({'display':'none'});
		$('.overview').css({'display':'none'});
		$('.buttons').css({'display':'none'});
		$('.jcarousel-scroll ').css({'display':'none'});
		$('#mycarousel').css({'height':'32px'});
	});
// expand back to carousel view
	$('.expand').click(function(){
		$('.viewport').css({'display':'block'});
		$('.overview').css({'display':'block'});
		$('.buttons').css({'display':'inline-block'});
		$('.page-out-of').css({'display':'inline'});
		$('.divider').css({'display':'inline'});
		$('.total-pages').css({'display':'inline'});
		$('.jcarousel-scroll').css({'display':'inline-block'});
		clicks = 1;
		$('.page-out-of').html(clicks);
		$('#mycarousel').tinycarousel({ start: 1, display: 4 });
		$('.overview').css({'height':'170px'});
		$('.viewport').css({'height':'170px'});
		$('#mycarousel').css({'overflow':'hidden'});
		$('#mycarousel').css({'height':'215px'});
		$('#mycarousel-minimize').css({'display':'none'});		
		$('.disable-next').remove();
		$('.disable-prev').remove();
		$('#mycarousel-show-all').css({'display':'inline-block'});
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
// turncate text below images after two lines	
	$('.carousel-caption').dotdotdot({
		height: 28,
		tolerance: 0,
		wrap: 'letter'
	});
/*
	$(".carousel-caption").each(function(i){
    len=$(this).text().length;
    if(len>50)
    {
      $(this).text($(this).text().substr(0,50)+'...');
    }
  });   */
});