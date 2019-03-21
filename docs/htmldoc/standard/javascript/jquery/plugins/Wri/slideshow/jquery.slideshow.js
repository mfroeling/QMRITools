/*
* jquery.slideshow.js
* A simple plugin to create 'slideshows' of divs or a combination of divs
*	jschoenick@wolfram.com 3/26/09
*		- Created
*	nicholash@wolfram.com 06/23/09
*		- Use jquery's each function instead of a for loop on arrays,
*		some cleanup.
*		(also fixed usage comment)
*/


// elements - an array of selector-strings describing each 'slide', ex:
//	[ '#box1', '#box2,#box2caption', '.thirdslide box,#box3' ]
// options (optional) 
//	- a hash of options, see first bit of function for definable options.
jQuery.simpleSlide = function(elements, options)
{
	var obj = jQuery.extend(
	{
		autoTransition : 3000, // milliseconds, 0 for off
		loopForever : false, // Default behavior is to stop autoTransition-ing after one loop
		fadeTime : 1000, // milliseconds for a complete fade, 0 for instant
		switchBoxClass : null,
		switchBoxClassActive : null,
		switchBoxContent : '',
		switchBoxWrapper : null, // an object that the list of switchboxes will be prepended into.
		// The following are used internally and shouldn't be set as
		// options unless you know what you are doing
		
		hasLooped : false,
		switching : false,
		eleList : elements,
		switchBoxes : new Array(),
		eleIndex : 0,
		pendingThink : null,
		resetThink : function()
		{
			var obj = this;
			if (obj.pendingThink)
			{
				window.clearTimeout(obj.pendingThink);
				obj.pendingThink = null;
			}
			if (obj.autoTransition)
				obj.pendingThink = window.setTimeout(function() { obj.think() }, obj.autoTransition);
		},
		finishSwitch : function()
		{
			var obj = this;
			obj.switchBoxes[obj.eleIndex].addClass(obj.switchBoxClassActive);
			var ele = $(obj.eleList[obj.eleIndex]);
			ele.fadeIn(obj.fadeTime/2);
			obj.switching = false;
		},
		switchTo : function(index)
		{
			var obj = this;
			if (!obj.switching)
			{
				obj.switching = true;
				var ele = $(obj.eleList[obj.eleIndex]);
				obj.switchBoxes[obj.eleIndex].removeClass(obj.switchBoxClassActive);
				obj.eleIndex = index;
				ele.fadeOut(obj.fadeTime/2, function () { obj.finishSwitch() });
			}
		},
		think : function()
		{
			var obj = this;
			if (obj.autoTransition && ((obj.loopForever == true) || (obj.hasLooped == false)))
			{
				if (obj.eleList.length <= obj.eleIndex + 1)
				{
					obj.switchTo(0);
					obj.hasLooped = true;
				}
				else
					obj.switchTo(obj.eleIndex + 1);
			}
			obj.resetThink();
		}
	}, options);
	
	if (obj.switchBoxWrapper)
	{
		// create switchBoxes
		$.each(obj.eleList,function(i,img)
		{
			var box = $('<div>' + obj.switchBoxContent + '<!-- -->' + '</div>');
			box.click(function ()
			{
				if (i != obj.eleIndex)
					obj.switchTo(i);
				obj.resetThink();
			});
			box.addClass(obj.switchBoxClass);
			if (i == 0)
				box.addClass(obj.switchBoxClassActive);
			box.appendTo(obj.switchBoxWrapper);
			obj.switchBoxes.push(box);
		});
	}
	
	obj.resetThink();
	return obj;
};

