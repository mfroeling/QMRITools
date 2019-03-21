#!/usr/local/bin/perl

print qq(Content-type: text/html\n\n);

use strict;
use Mail::RFC822::Address qw(valid);
use lib '/www/proj/libs';
use WRI::WebForms;
use CGI;

my $vars = new CGI->Vars();


if($vars->{'status'} eq "send") {
    my $errors = validate();

    if($errors->{'errors'} == 1) {
        printForm($errors);
    } else {
        submitEmail();
        printThanks();
    }
} else {
    printForm(); 
}



sub printForm() {

    my $errors = shift @_;
    print <<EOF;
    <!DOCTYPE html>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html;">
<link rel="stylesheet" type="text/css" href="/common/includes/CSS.cgi?b=2009nonHPheader;p=css2009/;v=en">
<link rel="stylesheet" type="text/css" href="/common/includes/CSS.cgi?b=Main;p=css2003/;v=en">
<link rel="stylesheet" type="text/css" href="/common/includes/CSS.cgi?b=M6;p=css2003/www.wolfram.com/;v=en">
<link rel="stylesheet" type="text/css" href="/common/includes/CSS.cgi?b=M6;p=css2003/www.wolfram.com/;v=en">
<link rel="shortcut icon" href="/common/images2003/favicon.ico" type="image/x-icon">
<script type="text/javascript" src="/common/javascript/jquery/core/1.2.6/jquery.js"></script>
<script type="text/javascript" src="/common/javascript/jquery/plugins/addState/1.0/jquery.addState.js"></script>
<script language="JavaScript" type="text/javascript" src="/mathematica/javascript/button.js"></script>
<script language="JavaScript" type="text/javascript" src="/common/javascript/placeholder.js"></script>
<link rel="stylesheet" type="text/css" href="/common/css/css2010/m8.css">


<script type="text/javascript">

    \$(document).ready( function() {

        var the_defaults = {'sign-up-email':'leave your email address'};
        var the_form = document.getElementById('the-form');

        var make_it_say = function(text, is_default) {
            this.value = text;
            if(is_default) {
                \$(this).addClass('placeholder');
            } else {
                //remove the class..
                \$(this).removeClass('placeholder');
            }
        };

        for(var i in the_defaults) {
            if( the_form[i].value == '' ) {
                make_it_say.call(the_form[i], the_defaults[i], true);
            } else if (the_form[i].value == the_defaults[i]) {
                \$(the_form[i]).addClass('placeholder');
            }

            the_form[i].onfocus = function() {
                if( this.value == the_defaults[this.name] ) {
                    make_it_say.call(this, '', false);
                }
            };
        
            the_form[i].onblur = function() {
                if( this.value == '' ) {
                    make_it_say.call(this, the_defaults[this.name], true);
                }
            };
        }

    });
</script>
<style>
/*sign up form*/
input.text { border:1px solid #c1c1c1; height: 33px; line-height:33px; width:300px; color:#000; box-shadow: inset 2px 2px #eee; margin: 1px 8px 0 0; padding: 0 12px; float: left; }
input.text.placeholder { font-type:italic; color:#666666 }
#share { font:12px/17px Arial,Verdana,Geneva,sans-serif; margin-top:5px; }
#thanks { font:14px/17px Arial,Verdana,Geneva,sans-serif; margin-bottom:30px; }
p.formtext { margin:0 0 15px 0; font-size: 14px; }
.errors { color:#e00400; display:block; margin:0; }
.errors img { vertical-align: middle; margin: 0 5px 0 0; }
.unfinished-image { float: left; }
.unfinished-form { float: left; margin: 15px 0 0 20px; }
.text { float: left; }
#thanks { width: 300px; float: left; }
</style>

</head>
<body style="background-image:none; background-color:#fff; padding: 20px 30px;">

<img src="../images/unfinished.png" alt="unfinished" class="unfinished-image">
<form action="" method="POST" id="the-form" class="unfinished-form">
    <p class="formtext"><strong>This page isn't quite ready yet.</strong><br>Sign up to be notified about progress:</p>
    <input type="text" class="text placeholder" id="sign-up-email" name="sign-up-email" value="$vars->{'sign-up-email'}">
    <div class="largebutton" style="float:left; ">
        <span class="largebuttonLeft"><!-- --></span>
        <span class="largebuttonRight"><input id="submit" type="submit" value="Send"></span>
        <div class="clearingFloats"></div>
    </div>
    <div class="clearingFloats"></div>
    $errors->{'sign-up-email'}
    <input type="hidden" name="status" value="send">
    <input type="hidden" name="url" value="$vars->{'url'}">
</form>


</body>
</html>
EOF
}

sub printThanks() {
    
    print <<EOF;
    <!DOCTYPE html>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html;">
<link rel="stylesheet" type="text/css" href="/common/includes/CSS.cgi?b=2009nonHPheader;p=css2009/;v=en">
<link rel="stylesheet" type="text/css" href="/common/includes/CSS.cgi?b=Main;p=css2003/;v=en">
<link rel="stylesheet" type="text/css" href="/common/includes/CSS.cgi?b=M6;p=css2003/www.wolfram.com/;v=en">
<link rel="stylesheet" type="text/css" href="/common/includes/CSS.cgi?b=M6;p=css2003/www.wolfram.com/;v=en">
<link rel="shortcut icon" href="/common/images2003/favicon.ico" type="image/x-icon">
<script type="text/javascript" src="/common/javascript/jquery/core/1.2.6/jquery.js"></script>
<script type="text/javascript" src="/common/javascript/jquery/plugins/addState/1.0/jquery.addState.js"></script>
<script language="JavaScript" type="text/javascript" src="/mathematica/javascript/button.js"></script>

<link rel="stylesheet" type="text/css" href="/common/css/css2010/m8.css">
<style>
/*sign up form*/
input.text { border:1px solid #c1c1c1; height:25px; width:255px; color:#000;}
input.text.placeholder { font-type:italic; color:#666666 }
#share { font:12px/17px Arial,Verdana,Geneva,sans-serif; margin-top:5px; }
#thanks { font:14px/17px Arial,Verdana,Geneva,sans-serif; margin: 30px 20px; width: 300px; float: left; }
p.formtext { margin:0 0 10px 0; }
.errors { color:#e00400; display:block; margin:0;}
.unfinished-image { float: left; }
</style>

</head>
<body style="background-image: none; background-color: #fff; padding: 20px 30px;">
    <img src="../images/unfinished.png" alt="unfinished" class="unfinished-image">
    <p id="thanks"><img src="../images/thanks.png" alt="Thanks!"><br><br>We'll keep you posted.</p>
</body>
</html>  
EOF
}



sub submitEmail() {
    my $formname = 'wolfram-language-hyperlinks';
    my $wf   = WRI::WebForms->new('prd');
    my $dbh  = $wf->connect_sqlDB($wf->{'connection'});
    $wf->insertSubmissions($formname, $vars);
}




sub validate() {

    my $errors;
    $errors->{'errors'} = 0;

    
    if($vars->{'sign-up-email'} eq "" || $vars->{'sign-up-email'} eq "leave your email address" 
       || !valid($vars->{'sign-up-email'})) {
        $errors->{'errors'} = 1;

        $errors->{'sign-up-email'} = "<span class=\"errors\" /><img src=\"../images/warning.png\" alt=\"warning\">Please enter a valid email address.</span>"
    }

    return $errors;
}
