var BoxWithShadow = function(configuration) {
    var BWS_Wrapper, 
        BWS_select, 
        BWS_parts, 
        id,
        image_types, 
        isImage,
        isNode,
        isMSIE,
        complianceMode,
        positionType,
        displayType,
        visibleHeight,
        visibleWidth,
        bodyOverflow,
        events,
        screen_components, 
        i;


    if(configuration === undefined) {
        configuration = {};
    }

    image_types = configuration.image_types || ['jpg','jpeg','png','gif'];

    BWS_parts = ['shadow','stage','header','exit','left','content', 'caption', 'right','footer'];

   /**
    * The screen is made of 9 parts:
    * 
    *   1) The wrapper - the outermost layer, everything is a child of this
    *   2) The shadow - a translucent layer that covers the entire page.
    *   3) The header - the top area of the stage
    *   4) The footer - the bottom area of the stage
    *   5) The exit - the area of the stage that closes the screen when clicked
    *   6) The stage - the center piece located on top of the shadow
    *   7) The left - the area of the stage that navigates left (smaller index)
    *   8) the right - the area of the stage that navigates right (larger index)
    *   9) The content - the section of the stage that holds content
    */

    id = (configuration.id) ? 'id="' + configuration.id + '"' : '';

    BWS_Wrapper = $('<div ' + id + ' class="BWS_Wrapper"></div>').appendTo('body');
    $(BWS_Wrapper).hide().html(
       '<div class="shadow"></div>' +
       '<div class="stage">' +
           '<div class="header"><!-- --></div>' +
           '<div class="exit"><!-- --></div>' +
           '<div class="content"><!-- --></div>' +
           '<div class="caption"><!-- --></div>' +
           '<div class="left"><!-- --></div>' +
           '<div class="right"><!-- --></div>' +
           '<div class="footer"><!-- --></div>' +
        '</div>'
    );

    // Building a map of the components, starting with the wrapper
    screen_components = { wrapper: BWS_Wrapper }; 
    
    // Helper function to find a component
    BWS_select = function(component) {
        return $('.' + component, BWS_Wrapper)[0];
    };

    // Loop over the list of parts, find each component and save a reference to it
    for(i = 0; i < BWS_parts.length; i++) {
        var obj = BWS_select( BWS_parts[i] );
        var default_text = configuration[ BWS_parts[i] ];
        if( default_text ) {
            $(obj).html(default_text);
        }
        screen_components[ BWS_parts[i] ] = obj;
    }

    isImage = function(element) {
        if( isNode(element) ) {

            switch( element.nodeName.toLowerCase() ) {
            case 'img':
                return true;
            case 'a':
                var img_re = new RegExp('.' + image_types.join('|') + '$', 'i');
                if( $(element).attr('href').match(img_re) ) {
                    return true;
                }         
                break;

            }
        }
        return false;
    };

    isNode = function(element) {
        return (element.nodeType && element.nodeType === 1);
    };

    events = {
        shadow: {click: ['close']},
        exit: {click: ['close']},
        left: {click: ['prev','activate']},
        right: {click: ['next','activate']}
    };

    isMSIE = (navigator.userAgent.toLowerCase().indexOf('msie') >= 0 && navigator.userAgent.toLowerCase().indexOf('opera') === -1);
    complianceMode = document.compatMode;
    if(isMSIE) {
    	var ver = getInternetExplorerVersion();
    	
    	if(ver > 8 || ver == -1) {
    		isMSIE = false;
    	} else {
    		isMSIE = true;
    	}
    }
    positionType = (isMSIE) ? 'absolute' : 'fixed';
    displayType = 'inline';
    bodyOverflow = $('body').css('overflow');

    visibleHeight = function(obj) {

        var padding = parseInt($(obj).css('padding-top'),10) +
                      parseInt($(obj).css('padding-bottom'),10) || 0;

        var margin = parseInt($(obj).css('margin-top'),10) +
                     parseInt($(obj).css('margin-bottom'),10) || 0;

        var border = parseInt($(obj).css('border-top-width'),10) +
                     parseInt($(obj).css('border-bottom-width'),10) || 0;

        var height = obj.clientHeight || $(obj).height();
        return padding + margin + border + height;

    };

    visibleWidth = function(obj) {

        var padding = parseInt($(obj).css('padding-left'), 10) +
                      parseInt($(obj).css('padding-right'), 10) || 0;

        var margin = parseInt($(obj).css('margin-left'), 10) +
                     parseInt($(obj).css('margin-right'), 10) || 0;

        var border = parseInt($(obj).css('border-left-width'), 10) +
                     parseInt($(obj).css('border-right-width'), 10) || 0;

        return padding + margin + border + $(obj).width(); 

    };

    return {

        components: screen_components,
        options: configuration,
        master_list: [], 
        cursor: 0,
        active: false,
        events_bound: false,

       /**
        * Add an item to the gallery/list
        */
        enqueue: function(element) {
            this.master_list.push(element);
        },
    
       /**
        * Remove an item from the gallery/list, return the removed item
        */
        pop: function(element) {
            return this.master_list.pop(element);
        },

       /**
        * Move the cursor to a specified index, if it is a valid index
        */
        seek: function(index) {
            index = Number(index);
            if( this.withinBounds(index) ) {
                this.cursor = index;
                return true;
            }
            return false;            
        },

       /**
        * Return the item at a specified index, if the index is valid.
        * If no index is provided, returns the item located under the cursor
        */
        peek: function(index) {
            if( index === undefined ) {
                index = this.cursor;
            }
            return (this.withinBounds(index)) ? this.master_list[ index ] : undefined;
        },

       /**
        * Opens up the box with shadow
        */
        open: function() {
            var BWS = this;

            if( this.open_enabled === false ) {
                return;
            }
            
            if( isMSIE ) {
                $('body').css('overflow','hidden');
            }

            $(this.components.wrapper).show();
            this.styleComponents(function() {
                BWS.resizeShadow();
            });
    
            if( this.master_list.length > 1 ) {
                $(this.components.left).show();
                $(this.components.right).show();
            } else {
                $(this.components.left).hide();
                $(this.components.right).hide();
            }
        
            if(this.events_bound === false) {
                this.bindEvents();
            }

            return this;
        },

       /**
        * Closes the box with shadow
        */
        close: function() {
            $(this.components.wrapper).hide();
            if( isMSIE ) {
                $('body').css('overflow',bodyOverflow);
            }
        },

       /**
        * Activate an item at an index, if the index is valid. Finishes by
        * calling the open function
        */
        activate: function(index) {
            var BWS = this;
            if(index !== undefined) {
                this.seek(index);
            }
            var element = this.peek();

            if( $(element).attr('title') ) {
                if( this.options.header ) {
                    this.setCaption( $(element).attr('title') );
                } else {
                    this.setHeaderText( $(element).attr('title') );
                }
            } else {
                if( this.options.header ) {
                    $(this.components.header).html( this.options.header );
                } else {
                    this.setHeaderText('');
                }
                this.setCaption('');
            } 
            $(this.components.stage).css('opacity', 0);
            this.open();
            this.activateFunction(element).call(this, element);

        },

       /**
        * Advance the cursor. Wraps around to 0 when it gets too high
        */
        next: function() {
            this.cursor = (this.cursor + 1) % this.master_list.length;
            return this;
        },
    
       /**
        * Retreats the cursor back, wraps around to the 
        * highest index when it gets below 0
        */
        prev: function() {
            this.cursor = (this.cursor - 1);
            if( this.cursor < 0 ) {
                this.cursor = this.master_list.length - 1;
            }
            return this;
        },

       /**
        * Determines if an index is within the bounds of the list
        */
        withinBounds: function(index) {
            index = Number(index);
            return ( this.master_list.length > 0 && 
                         index < this.master_list.length && index >= 0 );
        },

       /**
        * Determines which function to use to activate the element, based on
        * the settings, and the type of element being activated. Returns the
        * function.
        */
        activateFunction: function(element) {
            if (this.options.activate && this.options.activate instanceof Function) {
                return this.options.activate;
            }    

            if( isImage(element) ) {
                return (this.options.activate && this.options.activate.img) ? this.options.activate.img : this.stageImage;   
            } else if ( isNode(element) && element.nodeName.toLowerCase() === 'a') {
                return (this.options.activate && this.options.activate.a) ? this.options.activate.a : this.stageURL;
            } else if( isNode(element) && this.options.activate ) {
                if( this.options.activate[ element.nodeName.toLowerCase() ] ) {
                    return this.options.activate[ element.nodeName.toLowerCase() ];
                } else if (this.options.activate.inline) {
                    return this.options.activate.inline;
                }
            }
    
            // As a default, we just treat it like inline content.
            return this.stageInlineContent;
        },

       /**
        * Bind the necessary DOM events to the appropriate functions
        */
        bindEvents: function() {
            var BWS = this;
            var i,j,k;

            var resize_window = function() {
                BWS.resizeShadow();
                BWS.centerStage();
            };


            var attach_event = function(obj,evt,fn) {
                if( fn instanceof Function ) {
                    $(obj).bind(evt, function() {
                        fn.call(BWS);
                    });
                }
            };
        
            for(i in events) {
                for(j in events[i]) {
                    for(k = 0; k < events[i][j].length; k++) {
                        attach_event(this.components[i], j, this[ events[i][j][k] ]);
                    }
                }
            }  
    
            $(window).resize(resize_window);

          
            this.events_bound = true;
        },

       /**
        * Resize the shadow so that it takes up 100% of the viewport, and
        * matches the specified style (bgcolor, opacity)
        */
        resizeShadow: function() {

            var topAdj = 0, leftAdj = 0;
            // If the user is on IE, we have to adjust the location
            // of the shadow to where the user is scrolled to
            if(isMSIE) {
                topAdj = $(document).scrollTop();
                leftAdj = $(document).scrollLeft();
            }

            $(this.components.shadow).css({
                position: positionType,
                top: topAdj, left: leftAdj,
                height: $(window).height(),
                width: $(window).width(),
                zIndex: '1000'
            });
        },
    
       /**
        * Centers the stage in the viewport, makes sure it is the
        * appropriate size based on defaults/settings
        */
        centerStage: function(callback) {
    
            var topAdj = 0, leftAdj = 0;
            if( isMSIE ) {
                topAdj = $(document).scrollTop();
                leftAdj = $(document).scrollLeft();
            }
            var sHeight,sWidth;

            sHeight = visibleHeight( this.components.stage );

            if( !this.options.width ) {
                $(this.components.stage).css('width', $(this.components.content).children().outerWidth());
            } else if (isMSIE) {
                $(this.components.stage).css('width', this.options.width)
            }
            sWidth = visibleWidth( this.components.stage );

            if( isMSIE && complianceMode === 'BackCompat') {
                $(this.components.stage).css('width', sWidth);
            }

            topAdj += ($(window).height() - sHeight) / 2;
            leftAdj += ($(window).width() - sWidth) / 2;
    
            $(this.components.stage).css({
                zIndex: '2000',
                position: positionType,
                display: displayType,
                top: topAdj, left: leftAdj
            });

            if(callback && callback instanceof Function) {
                callback();
            }
        },

        styleComponents: function(callback) {

            var i,k;
            var styles = {
                shadow: {backgroundColor: '#222', opacity: '0.8'},
                stage: {backgroundColor: '#fff'},
                content: {height: '100%', width: '100%'}

            };
    
            if( this.options.height ) {
                styles.content.height = this.options.height; 
            }

            if( this.options.width ) {
                styles.content.width = this.options.width;
            }

            // Map for shortcut options to css parameters,
            // get overridden by actual css parameters
            var optionMap = {
                shadow: {shadowColor: 'backgroundColor', shadowOpacity: 'opacity'},
                stage: {
                    stageColor: 'backgroundColor',
                    height: 'height',
                    width: 'width'
                }
            };
            
            for(i in optionMap) {
                for(k in optionMap[i]) {
                    if( this.options[k] ) {
                        styles.shadow[ optionMap[i][k] ] = this.options[k];
                    }
                }
            }

            if( this.options.css ) {
                for(i in this.options.css) {
                    if( styles[i] ) {
                        for(k in this.options.css[i]) {
                            styles[i][k] = this.options.css[i][k];
                        }
                    } else {
                        styles[i] = this.options.css[i];
                    }
                }
            }

            for(i in styles) {
                if( this.components[i] ) {
                    $( this.components[i] ).css( styles[i] );
                }
            }

            if(callback && callback instanceof Function) {
                callback();
            }
        }, 

       /**
        * Send a URL to the stage, by creating an iframe. The default behavior for
        * an <a> element, which doesn't link to a picture
        */
        stageURL: function(element) {
            var BWS = this;
            this.loadStarted();
            var iframe = $('<iframe frameborder="0" width="100%" height="100%" allowtransparency="true"></iframe>');
            $(iframe).attr('src', $(element).attr('href'));
            this.loadComplete(iframe);
        },

       /**
        * Puts an image up on the stage
        */
        stageImage: function(element) {
            var BWS = this;
            this.loadStarted();
            var src = $(element).attr('src') || $(element).attr('href');
            var img = document.createElement('img');
            $(img).load(function() {
                BWS.loadComplete(img);
            }).attr('src', src);
        },

       /**
        * Take inline content, and put a copy of it up on the stage
        */
        stageInlineContent: function(element) {
            if( isNode(element) ) {
                element = element.cloneNode(true);
            }
            this.loadComplete(element);
        },
        
       /**
        * Generic function for putting content up on the stage
        */
        stageContent: function(content) {
            var BWS = this;    
            $(this.components.content).hide().html(content).show(1,function() {
                BWS.centerStage(function() {
                    $(window).resize();
                });
                var opacity = (isMSIE) ? 'none' : 1;
                $(BWS.components.stage).css('opacity', opacity);
            });
        },

        setHeaderText: function(str) {
            if( str ) {
                $(this.components.header).html('<p>' + str + '</p>').show();
            } else {
                $(this.components.header).html('').hide();
            }
        },

        setCaption: function(str) {
            if( str ) {
                $(this.components.caption).html('<p class="CaptionText">' + str + '</p>').show();
            } else {
                $(this.components.caption).html('').hide();
            }
        },    

       /**
        * Call when the element starts loading
        */
        loadStarted: function() {
            $(this.components.wrapper).addClass('loading');
        },

       /**
        * Call when the element on the stage has finished loading
        */
        loadComplete: function(content) {
            this.stageContent(content);
            $(this.components.wrapper).removeClass('loading');
        }

    };

};

(function($) {
    $.fn.boxWithShadow = function(config_params) {

        // If no elements were matched, return here
        // so nothing gets written to the DOM unnecessarily
        if(this.length === 0) {
            return;
        }
        
        var BWS = new BoxWithShadow(config_params);
        var non_displayed = 0;
        $(this).each( function(index,value) {
            if( ! $(this).hasClass('dont-display') ) {
                BWS.enqueue( this );
            } else {
                non_displayed++;
            }
            $(this).bind('click', {offset: non_displayed}, function(e) {
                var real_index = index - e.data.offset;
                BWS.activate(real_index);
                e.preventDefault();
            });
        } );
        return BWS;
    };
}(jQuery));

function getInternetExplorerVersion()
//Returns the version of Internet Explorer or a -1
//(indicating the use of another browser).
{
var rv = -1; // Return value assumes failure.
if (navigator.appName == 'Microsoft Internet Explorer')
{
 var ua = navigator.userAgent;
 var re  = new RegExp("MSIE ([0-9]{1,}[\.0-9]{0,})");
 if (re.exec(ua) != null)
   rv = parseFloat( RegExp.$1 );
}
return rv;
}

