(function($) {
    
    $.fn.wflipbook = function(_imageItems, _options) {
	var _opts = $.extend({}, $.fn.wflipbook.settings, _options);
	this.data('imageItems', _imageItems);
	this.data('opts', _opts);
        $.fn.wflipbook.initialize(this);
    };
    
    $.fn.wflipbook.settings = {
        initialStart: 0,
        transSpeed: 5000,
        fadeSpeed: 'fast',
        maxLoop: -1
    };

    $.fn.wflipbook.initialize = function(element) {
        element.initT = setTimeout(function() {
          $.fn.wflipbook.startFlip(element);
        }, element.data('opts').initialStart);
    };

    $.fn.wflipbook.startFlip = function(element) {
        element.wflipbook_offSet = 0;
        element.wflipbook_loopCount = 0;
        $.fn.wflipbook.flip(element);
    };
    
    $.fn.wflipbook.flip = function(element) {
        var current = element.data('imageItems')[element.wflipbook_offSet];
        $(element).fadeOut(element.data('opts').fadeSpeed, function() {
            if(isArray(current)) {
                if(current.length == 2) {
                    $(element).html('<img src="' + current[0] + '" alt="' + current[1] + '" />').fadeIn(element.data('opts').fadeSpeed);
                } else {
                    $(element).html('<a href="' + current[2] + '"><img src="' + current[0] + '" alt="' + current[1] + '" /></a>').fadeIn(element.data('opts').fadeSpeed);
                }
            } else {
                $(element).html('<img src="' + current + '" alt="" />').fadeIn(element.data('opts').fadeSpeed);
            }
        });
        
        if((element.wflipbook_offSet + 1) >= element.data('imageItems').length) {
            element.wflipbook_offSet = 0;
            element.wflipbook_loopCount++;
        } else {
            element.wflipbook_offSet++;
        }
        
        if(element.wflipbook_loopCount != element.data('opts').maxLoop) {
            element.wflipbook_T = setTimeout(function() {
              $.fn.wflipbook.flip(element);
            }, element.data('opts').transSpeed);
        }
    };

    function isArray(test) {
        if ( test instanceof Array ) {
            return true;
        } else { 
            return false;
        }
    };
})(jQuery);
