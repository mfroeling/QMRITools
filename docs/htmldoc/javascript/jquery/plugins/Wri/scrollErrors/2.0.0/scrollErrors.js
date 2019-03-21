jQuery.WriValidator = function WriValidatorClass(conf){
    // Private attrinbutes ************************************************************************
    
    // the tests objects
    var _conf = conf;
    
    // the keys for the tests
    var _keys = [];
    
    // log of the objects with currently show errors grouped by test key
    var _errors =[];

    // initialize *********************************************************************************

    // get the keys once, so you don't have to loop over them several times
    for (var test in _conf){
        _keys.push(test);
    }

    //private methods *****************************************************************************

    // adds objects that have currently shown errors to the _errors log
    function _addShownErrors(key, objs){
        if (_errors.hasOwnProperty(key)){
            _errors[key] = $.merge(_errors[key], objs);
        }else{
            _errors[key] = objs;
        }
    }

     // removes objects from the _errors log for the given key and returns them
    function _removeShownErrors(key){
        var ret = $([]);
        
        if (_errors.hasOwnProperty(key)){
            ret = _errors[key];
            delete _errors[key];
        }
        
        return ret;
    }

    // returns the scroll top for an error associated with a specific object & test
    function _getscrollTop(key, obj){
            if($.isFunction(_conf[key].scrollTop)){
                return _conf[key].scrollTop(obj);    
            }else{
                return _conf[key].scrollTop;    
            }          
    }
    
    // the default callback called by scrollRequiredTestErrors
    function _defaultScrollRequiredTestErrorsCallBack(){;
        this.showRequiredTestErrors();
    }

    // the default callback called by scrollAllTestErrors
    function _defaultScrollAllTestErrorsCallBack(){
        this.showAllTestErrors();
    }
        
    //public methods ******************************************************************************
    return {
        // make sure each test object has a refrence to the current instance of the 
        // WriScrollErrorsClass class
        initialize: function(){
            for (var i = 0; i < _keys.length; i++){
                _conf[_keys[i]].v = this;
            }
        },
        
        // returns all the keys in the conf
        getKeys: function(){
            return _keys;
        },
        
        // returns all the selected jquery objects for a given key
        getObjects: function(key){
            if($.isFunction(_conf[key].selector)){
                return _conf[key].selector();    
            }else{
                return _conf[key].selector;    
            } 
        },

        // is a given test required 
        isRequired: function(key){
            if($.isFunction(_conf[key].required)){
                return _conf[key].required();    
            }else{
                return _conf[key].required;    
            }
        },

        // determines if a single jquery selection result object passes a given tests validation
        isObjectValid: function(key, obj){
            return _conf[key].valid(obj);
        },
        
        // determines if a test is internally valid. that is to say did every single field object 
        // returned by the selector passed the valid function
        isTestValid: function(key){
            var that = this;
            var objs = this.getObjects(key);
            var result = true;
            
            objs.each(function(){
                if(!that.isObjectValid(key,$(this))){
                    result = false;
                } 
            });
            
            return result;
        },

        // did every single required test pass Internal validation
        requiredTestsValid: function(){
            for (var i = 0; i < _keys.length; i++){
                if(this.isRequired(_keys[i])){
                    if(!this.isTestValid(_keys[i])){
                        return false;
                    }   
                }   
            }
            
            return true;
        },
        
        // did every single test pass internal validation
        allTestsValid: function(){ 
            for (var i = 0; i < _keys.length; i++){
                if(!this.isTestValid(_keys[i])){
                    return false;
                }  
            }
            
            return true;
        },

        // shows the error(s) associated with a given test
        showTestErrors: function(key){
            var that = this;
            var selected = this.getObjects(key);
            
            selected.each(function(){
                if(!that.isObjectValid(key,$(this))){
                    // the object was invalid so show an error for it, and log it
                    _addShownErrors.call(that, key, $(this));
                    _conf[key].show($(this));
                }  
            });            
        },
        
        // hides the error(s) associated with a given test that are currently shown
        hideTestErrors: function(key){
            var that = this;
            var selected = _removeShownErrors(key);
            
            selected.each(function(){
                    _conf[key].hide($(this));
            });            
        },

        // shows the error(s) associated with all test that are required
        showRequiredTestErrors: function(){
            for (var i = 0; i < _keys.length; i++){
                if(this.isRequired(_keys[i])){
                    this.showTestErrors(_keys[i]);
                }   
            }
        },

        // hides the error(s) associated with all test that are required
        hideRequiredTestErrors: function(){
            for (var i = 0; i < _keys.length; i++){
                if(this.isRequired(_keys[i])){
                    this.hideTestErrors(_keys[i]);
                }   
            }
        },

        // shows the error(s) associated with all tests
        showAllTestErrors: function(){
            for (var i = 0; i < _keys.length; i++){
                this.showTestErrors(_keys[i]);
            }
        },

        // hides the error(s) associated with all tests
        hideAllTestErrors: function(){
            for (var i = 0; i < _keys.length; i++){
                this.hideTestErrors(_keys[i]);
            }
        },
                
        // runs the required tests, and returns the scroll top for the highest field on the page that 
        // failed validation, or false if no errors are found. You can then use this value to just set scroll top, 
        // or you could use animate to go to it.
        getRequiredScrollTop: function(options){
            var hasError = false;
            var theScrollTop = 1000000;
            var that = this;
            
            for (var i = 0; i < _keys.length; i++){
                if(this.isRequired(_keys[i])){
                    // the key is required 
                                       
                    // get all the selectors for the key
                    var objects = this.getObjects(_keys[i]);
            
                    // loop through the field objects and see what ones are invalid
                    objects.each(function(){
                        if(!that.isObjectValid(_keys[i],$(this))){
                            // the field object was invalid so set the flag
                            hasError = true;
                            
                            // get the scropp top for the invalid field
                            var theTop = _getscrollTop.call(that, _keys[i], $(this));
                                if(theScrollTop > theTop){
                                    // the current scroll top is higher than the 
                                    // last so reset the place holder
                                    theScrollTop = theTop; 
                                }    
                        } 
                    });   
                }   
            }
            if(hasError){
                return theScrollTop;
            }else{
                return false;
            }              
        },        

        // runs all the tests, and returns the scroll top for the highest field on the page that 
        // failed validation, or false if no errors are found. You can then use this value to just set scroll top, 
        // or you could use animate to go to it.
        getAllScrollTop: function(options){
            var hasError = false;
            var theScrollTop = 1000000;
            var that = this;
            
            for (var i = 0; i < _keys.length; i++){
                // get all the selectors for the key
                var objects = this.getObjects(_keys[i]);
                
                // loop through the field objects and see what ones are invalid
                objects.each(function(){
                    
                    if(!that.isObjectValid(_keys[i],$(this))){
                        // the field object was invalid so set the flag
                        hasError = true;
                            
                        // get the scropp top for the invalid field
                        var theTop = _getscrollTop.call(that, _keys[i], $(this));
                            if(theScrollTop > theTop){
                                // the current scroll top is higher than the 
                                // last so reset the place holder
                                theScrollTop = theTop; 
                            }    
                    } 
                });   
            }
            if(hasError){
                return theScrollTop;
            }else{
                return false;
            }
        }
    };	
};