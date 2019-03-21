(function($) {

    $.fn.newsTicker = function(news_settings) {

        var news_timeout;

        return $(this).each(function() {

            var configuration = [],
                news_list = [],
                the_news,news_show,
                i=0,j=0;

            for(i in news_settings) {
                the_news = $('#' + i);
                if( the_news.length ) {
                    configuration.push(news_settings[i]);
                    news_list.push( the_news[0] );
                }
            }

            news_show = function(index) {

                //Hide all of the news items
                $(news_list).hide();
    
                //Show specifically the news item we're after
                $(news_list[index]).show();
    
                var duration = configuration[index]['duration'] || '5000';
                var next_index = (1+index) % news_list.length;
                
                news_timeout = setTimeout(function() {
                    news_show(next_index); 
                }, duration)

            };

            //Select an news item at random to display
            news_show(Math.floor( Math.random()*news_list.length ));
        });    

    };

}(jQuery));
