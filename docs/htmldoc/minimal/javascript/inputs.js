//THIS IS ALL SOMEHOW FOR THE INPUT POPUPS

var w3 = document.getElementById ? 1 : 0;
var ns = document.layers ? 0 : 0;
var ie = document.all ? 1 : 0;
var t_out = 0;
var t_out2 = 0;
var current_layer = 0;

function showLayer(thelayer, rol_imagename, rol_imagesrc) {
	if (t_out)		{ clearTimeout(t_out); }
	if (t_out2)		{ clearTimeout(t_out2); }
	if (current_layer)	{ hideLayer(current_layer); }
	if (w3) {
		document.getElementById(thelayer).style.display = 'block';
		if (document.getElementById('searchselect')) {
			document.getElementById('searchselect').style.display = 'none';
		}
		if ((rol_imagename) && (rol_imagesrc)) {
			document.getElementById(rol_imagename).src = rol_imagesrc;
		}
	}
	else if (ns) {
		eval("document.layers." + thelayer + ".display = 'block';");
		if ((rol_imagename) && (rol_imagesrc) && (document.images)) {
			eval("document.images." + rol_imagename + ".src = " + rol_imagesrc);
		}
	}
	else if (ie) {
		eval("document.all." + thelayer + ".style.display = 'block';");
		if (document.all.searchselect) {
			document.all.searchselect.style.display = 'none';
		}
		if ((rol_imagename) && (rol_imagesrc)) {
			eval("document.all." + rol_imagename + ".src = " + rol_imagesrc + ";");
		}
	}
	current_layer = thelayer;
}

function hideLayer(thelayer, rol_imagename, rol_imagesrc) {
	if (w3) {
		document.getElementById(thelayer).style.display = 'none';
		if ((rol_imagename) && (rol_imagesrc)) {
			document.getElementById(rol_imagename).src = rol_imagesrc;
		}
		if (document.getElementById('searchselect')) {
			document.getElementById('searchselect').style.display = 'block';
		}
	}
	else if (ns) {
		if ((rol_imagename) && (rol_imagesrc) && (document.images)) {
			eval("document.images." + rol_imagename + ".src = " + rol_imagesrc + ";");
		}
		eval("document.layers."+thelayer+".display = 'none';");
	}
	else if (ie) {
		eval("document.all."+thelayer+".style.display = 'none';");
		if ((rol_imagename) && (rol_imagesrc)) {
			eval("document.all." + rol_imagename + ".src = " + rol_imagesrc + ";");
		}
		if (document.all.searchselect) {
			document.all.searchselect.style.display = 'block';
		}
	}
}


function cons_input (i) {
  id_out = i + '_out';
  id_in = i + '_in';
  el_out = document.getElementById (id_out);
  el_in = document.getElementById (id_in);
    if (el_in.firstChild.nodeType == 1) {
      el_in = el_in.firstChild;
    }
    var innertxt = el_in.innerHTML.replace(/(\r|\n)+/g, '\n');
    rows = innertxt.split('\n').length + 1;
    var url = '';
    for (ch = el_out.firstChild; ch; ch = ch.nextSibling) {
      if (ch.className == 'IFU') {
        url = strValue (ch);
      }
    }
    /* We use .innerHTML when we want &lt; and .firstChild.nodeValue when we want < */
    div = '\n\t\t\t<div id="' +id_in+ '" class="IFF">\n'
        + '\t\t\t\t<form><textarea rows="' +rows+ '">' + innertxt + '</textarea></form>\n'
        + '\t\t\t\t<div class="closePopup"><a href="javascript:hideLayer(' +"'" +id_in+ "'"+ ')">x</a></div>\n'
        + '\t\t\t\t<div class="clear"></div>\n'
        + '\t\t\t\t<div class="bottomStuff">\n'
        + '\t\t\t\t\t<div class="mailLink"><a href="' + 'mailto:?Subject=Wolfram%20Language%20Example&amp;Body=' + url + '%0d%0a%0d%0a%0d%0a' + escape(el_in.firstChild.nodeValue)+ '">email</a></div>\n'
    if (url) { div += '\t\t\t\t\t<div class="IFU">' + url + '</div>\n' }
    div += '\t\t\t\t\t<div class="clear"></div>\n\t\t\t\t</div>\n\t\t\t</div>\n';
    el_out.innerHTML = div;
    showLayer(id_in);
}
/*  showInputForm */
function input(i) {
  id_out = i + '_out';
  id_in = i + '_in';
  el_out = document.getElementById (id_out);
  el_in = document.getElementById (id_in);
  if (!el_in) {
    el_url = "Files/" + baselang + "/" + i + ".txt";
    req_func = function (str) {
      el_out.innerHTML = str;
      cons_input (i);
    }
    makeRequest (el_url, req_func);
  }
  else if (el_in.className == 'InputFormText' || el_in.className == 'IFT') {
    cons_input (i);
  }
  else {
    return el_in.style.display == 'block' ? hideLayer(id_in) : showLayer(id_in);
  }
}

function strValue (el) {
  s = '';
  if (el.nodeType == 3) {
    s = s + el.nodeValue;;
  } else {
    for (var i = el.firstChild; i; i = i.nextSibling) {
      s = s + strValue(i);
    }
  }
  return s;
}


function makeRequest(url, func) {

 var http_request = false;

 if (window.XMLHttpRequest) // Mozilla, Safari,...
 {
   http_request = new XMLHttpRequest();
   if (http_request.overrideMimeType)
   {
     http_request.overrideMimeType('text/plain');
   }
 }
 else if (window.ActiveXObject) // IE
 {
   try {
     http_request = new ActiveXObject("Msxml2.XMLHTTP");
   } catch (e) {
     try {
       http_request = new ActiveXObject("Microsoft.XMLHTTP");
     } catch (e) {}
   }
 }

 if (!http_request)
 {
   alert('This feature requires an AJAX capable browser,\r such as Microsoft Internet Explorer version 5.0 and above,\rany version of Mozilla or Firefox, Netscape\rversion 7.1 and above, Apple Safari version 1.2\rand above, or Opera version 8.0 and above.');
   return false;
 }

 http_request.onreadystatechange = function()
  {   
    if (http_request.readyState == 4)
    {
      if (http_request.status == 200)
      {
        func (http_request.responseText);
      }
    } 
  };
 http_request.open('GET', url, true);
 http_request.send(null);

}