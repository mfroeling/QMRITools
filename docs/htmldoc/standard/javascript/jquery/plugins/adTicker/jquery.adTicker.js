(function($) {

    $.fn.adTicker = function(ad_settings) {

        var ad_timeout;

        return $(this).each(function() {

            var configuration = [],
                ad_list = [],
                the_ad,ad_show,
                i=0,j=0;

            for(i in ad_settings) {
                the_ad = $('#' + i);
                if( the_ad.length ) {
                    configuration.push(ad_settings[i]);
                    ad_list.push( the_ad[0] );
                }
            }

            ad_show = function(index) {

                //Hide all of the ads
                $(ad_list).hide();
    
                //Show specifically the ad we're after
                $(ad_list[index]).show();
    
                var duration = configuration[index]['duration'] || '5000';
                var next_index = (1+index) % ad_list.length;
                
                ad_timeout = setTimeout(function() {
                    ad_show(next_index); 
                }, duration)

            };

            //Select an ad at random to display
            ad_show(Math.floor( Math.random()*ad_list.length ));
        });    

    };

}(jQuery));
