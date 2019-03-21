$(document).ready(function() {
    $('.InCell img').each(function() {
        if($(this).is('[data-src]')) {
            var src = $(this).attr('data-src');
            var w = $(this).attr('data-big').split(' ')[0];
            var h = $(this).attr('data-big').split(' ')[1];
            $(this).prop('src', src);
            $(this).prop('width', w);
            $(this).prop('height', h);
        }
    });
});