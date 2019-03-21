// modified to use ordered lists, to escape conflict

(function($){ 
     $.fn.extend({  
         accordion: function() {       
            return this.each(function() {
				if($(this).data('accordiated'))
					return false;									
				$.each($(this).find('ol, li>div'), function(){
					$(this).data('accordiated', true);
					$(this).hide();
				});
				$.each($(this).find('a'), function(){
					$(this).click(function(e){
						activate(e.target);
						return void(0);
					});
				});
				
				var active = false;
				if(location.hash)
					active = $(this).find('a[href=' + location.hash + ']')[0];
				
				if(active){
					activate(active, 'toggle','parents');
					$(active).parents().show();
				}
				
				function activate(el,effect,parents){
					$(el)[(parents || 'parent')]('li').toggleClass('active').siblings().removeClass('active').children('ol, div').slideUp('fast');
					$(el).siblings('ol, div')[(effect || 'slideToggle')]((!effect)?'fast':null);
				}
				
            });
        } 
    }); 
})(jQuery);