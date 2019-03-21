(function($) {
    $.fn.errorBox = function() {

        var parent_object = this[0];
        var field_report = [];
 
        $('.errorbox', this).each( function(index, element) {

            var the_name = $(this).attr('for');
            var the_text = $(this).text();
            $(this).remove();
            field_report.push(the_name);
            //only create the div element if we can find the item described by 'the_name'

            if( parent_object[the_name] ) {
        
                //make new div
                var the_bubble = document.createElement('div');
                    the_bubble.className = 'errorBoxWrapper';

                var the_html = '<div class="errorTop"><!-- --></div>' +
                               '  <div class="errorMiddle boldtitle">' + the_text + '</div>' +
                               '<div class="errorBottom"><!-- --></div>';

                the_bubble.innerHTML = the_html;

                //bind hover event to the_name and make it toggle the_div..
                //and then there is no need for unique id tracking!
                
                var eh = dotclosest(parent_object[the_name], 'errorHighlight');
                $(eh).prepend(the_bubble);

                var object = ($(parent_object[the_name]).size() > 1) ? eh : parent_object[the_name];   

                //and hide the bubble
                $(the_bubble).hide();

                if( !navigator.userAgent.match(/Mobile|iPhone|iPad|Android/i) ) {
                    $(eh).bind('mouseover', {'field': parent_object[the_name],
                                             'object': object, 'error': the_bubble}, turn_on_errorbox);

                    $(eh).bind('mouseout', {'field': parent_object[the_name],
                                            'error': the_bubble}, turn_off_errorbox);
                }

                $(parent_object[the_name]).bind('focus', {'field': parent_object[the_name],
                                                'object': object, 'error': the_bubble}, turn_on_errorbox);

                $(parent_object[the_name]).bind('blur', {'field': parent_object[the_name],
                                                'error': the_bubble}, turn_off_errorbox);


                setTimeout(function() {
                    if(!window.visible_errorbox) {
                        $(parent_object[the_name]).trigger('focus');
                    }
                }, 500);

            }
        });

        if(field_report.length) {
            var report_str = '/common/includes/form_errors.txt?url=' + escape(document.URL) + '&errors=' + escape('{' + field_report.join('}_{') + '}');
            $.get(report_str);
        }
    };
})(jQuery);

function turn_on_errorbox(event) {
    if(event.data && event.data.error) {

        // if the error box that we're being asked to turn on is not the one already
        // being shown, then we hide the current one.
        // Instead of returning here because the box we want is alread showing, 
        // we run through the entire logic again, in case the box needs to be
        // relocated for some reason
        if(window.visible_errorbox && window.visible_errorBox != event.data.error) {
            $(window.visible_errorbox).hide();
        }

        window.visible_errorbox = event.data.error;

        var object = event.data.object;
        var the_bubble = event.data.error;

        var os = $(object).offset()

        //if the object is a TD, we want to move the bubble up by half of its height,
        //otherwise, we move it up by the full amount of its height.
        var is_td = (object.nodeName === 'TD') ? 2 : 1;

        //we can now calculate where the bubble will need to sit
        var oh = $(the_bubble).outerHeight() || 61;
        var ow = $(the_bubble).outerWidth() || 230;
        var top_offset = os.top - 61 / is_td;
        var left_offset = os.left + ($(object).outerWidth() / 2) - (ow / 2);

        $(the_bubble).css({
            'position':'absolute',
            'z-index' :'1000',
            'top'     : top_offset,
            'left'    : left_offset
        });
        $(event.data.error).show();
    }
}

function turn_off_errorbox(event) {

    if(window.visible_errorbox) {
        $(window.visible_errorbox).hide();
    }
    window.visible_errorbox = undefined;

}

function dotclosest(the_element, the_class) {
    if( the_element === null || the_element === undefined || the_element.length === 0) {
        return null;
    }

    if(the_element.length > 0 && the_element.nodeName !== 'SELECT') {
        the_element = the_element[0];
    }

    if( $(the_element).hasClass(the_class) ) {
        return the_element;
    } else {
        return dotclosest($(the_element).parent(), the_class);
    }
}
