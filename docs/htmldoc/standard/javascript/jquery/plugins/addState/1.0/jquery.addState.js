(function($) {
    $.fn.addState = function(state, on_event, on_callback, off_event, off_callback) {

        if(off_event && !off_callback) {
            //rec'd 4 parameters

            /**
             * two possibilities:
             *     str      str        fn          str
             * fn(state, on_event, on_callback, off_event) (which off_callback (the 5th param) was just omitted)
             *
             *                     (on_callback) (off_event)    (off_callback)
             *     str      str        str         fn
             * fn(state, on_event, off_event,    off_callback) (which on_callback (the 3rd param) was omitted)
             */
    
            if( typeof(on_callback) === 'string' ) {
                off_callback = off_event;
                off_event = on_callback;
                on_callback = undefined;
            }
        } else if (on_callback && !off_event && !off_callback) {
            //rec'd 3 parameters
            off_event = on_callback;
            on_callback = undefined;
        } else if(!on_callback) {
            //didn't receive enough parameters
            throw 'addState called with too few parameters!';
        }

        var toggle_class = function(obj, cls) {
            if( $(obj).hasClass(cls) ) {
                $(obj).removeClass(cls);
            } else {
                $(obj).addClass(cls);
            }
        };

        var do_bind = function(obj, event_a, callback_a, event_b, callback_b) {
            $(obj).bind(event_a, function() {
                toggle_class(obj,state);
                $(this).unbind(event_a);
                do_bind(this, event_b, callback_b, event_a, callback_a);
                if( $.isFunction(callback_a) ) {
                    callback_a.apply(this);
                }
            });
        };

        $(this).each(function() {
            if( $(this).hasClass(state) ) {
                do_bind(this, off_event, off_callback, on_event, on_callback);
            } else {
                do_bind(this, on_event, on_callback, off_event, off_callback);
            }
        });
    };
})(jQuery);
