(function($) {

    $.fn.placeholder = function(default_text, state_class) {
        
        // a custom class can be passed in, if placeholder is conflicting
        if(!state_class) { state_class = 'placeholder'; }

        var remove_text = function() {
            if( $(this).val() == default_text ) {
                $(this).attr('value', '');
                $(this).removeClass(state_class);
            }   
        };
        
        var add_text = function() {
            if( !$(this).val() ) {
                $(this).attr('value', '');
                $(this).removeClass('placeholder'); 
            } else if ( $(this).val() == default_text ) {
                $(this).addClass(state_class) // just to make sure
            }
        };

        var bind_events = function() {
            // on focus we clear the field if it contains the default text
            $(this).focus( remove_text );

            // on blur we put the text in the field if it is empty
            $(this).blur( add_text );
    
            //trigger the blur event, to add the text, if needed
            $(this).blur();
        };

        // bind the appropriate events for each of the matching elements
        $(this).each(bind_events);

        // when a form is submitted, go through and clean out
        // all of the fields that have their placeholder text
        $('form').submit( function() {
            $('.' + state_class,this).each(remove_text);
            return true;
        } );

    };

})(jQuery);
