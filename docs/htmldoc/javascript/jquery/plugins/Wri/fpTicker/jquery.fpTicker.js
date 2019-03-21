/**
 * The animation was removed from this function.
 * The parameters were left, in case they need to be used again at some point
 */

(function($) {
    $.fn.fpTicker = function(xmlfile,rate,fadein,fadeout,holdover) {
        $(this).each(function() {
            var tick = new Wri_FPTicker(this,xmlfile,rate,fadein,fadeout,holdover);
            tick.tock();
        });
    };
})(jQuery);

/**
 * Flag that keeps track if the mose is hovering over the news ticker.
 */
var hovering = false;

/**
 * Detects mouse hover ing on the ticker and switches flag on mouse in and 
 * mouse out
 */    
$(document).ready(function() {
    $('.fpticker').hover(function () {
        hovering = true;
    },
    function() {
        hovering = false;
    });
});

/**
 * Front page ticker for the WRI homepage.  The only required arguments are the
 * elem and xmlfile.  The ticker we be built inside the element with the id name
 * represented by elem, and the headlines will be pulled in from the xml file, 
 * located at where xmlfile is set to.
 */
var Wri_FPTicker = function(elem,xmlfile,rate,fadein,fadeout,holdover) {
    this.elem = elem;
    $(this.elem).html('<div class="inner"></div>');
    this.rate = rate || 5000; 
    this.fadein = fadein || 500; 
    this.fadeout = fadeout || 500; 
    this.holdover = holdover || 10;
    this.headlines = this.parse_xml(xmlfile);
    this.current_message = 0;
};

Wri_FPTicker.prototype = {
        


   /**
    * Calling this function starts an infinite loop of throwing messages on the
    * ticker.  
    */
    tock: function() {
        if(!hovering) {
            this.display_message(this.advance());
        }
        (function(fn,obj) {
            setTimeout(function() {
                fn.call(obj);
            },obj.rate);
        })(this.tock,this)   
    },

   /**
    * Has the new headline message built, and switches out the old message for
    * new message
    */
    display_message: function(i) {
        this.switch_message( 
            this.build_message( $('text', this.items[i-1]).text(), $('url', this.items[i-1]).text() ), this.fadein, this.fadeout, this.holdover);
    },

   /**
    * Builds the message as a pre-built url or plaintext string if no url is
    * provided
    */
    build_message: function(txt,url) {
        return (url) ? '<a href="' + url + '">' + txt + '</a>' : txt;
    },

   /**
    * Fades out the current message, replaces the text and fades back in
    */
    switch_message: function(txt,fi,fo,ho) {
        $(this.elem).children('.inner').html(txt);
    },

   /**
    * Pushes the current message forward, returns the new message index.
    */
    advance: function() {
        return this.current_message = this.next_message(this.current_message);
    },

   /**
    * Calculates the next message in the queue.
    * This rotates -- at the end of the message list it goes back to the front.
    */
    next_message: function(i) {
        return  ((i + 1) % this.items.length) || this.items.length;
    },

   /**
    * Requests the XML file, parses it and builds the list of headline items.
    */
    parse_xml: function(xmlfile) {
        var obj = this;
        $.ajax({
            url: xmlfile,
            async: false,
            success: 
                function(data) {
                    obj.items = $(data).find("hlist").find("headline");
                },
            error:
                function(e) {
                    throw 'Error loading XML: ' + e.status + ' ' + e.statusText;
                }
        });
        return obj.items;
    }

};
