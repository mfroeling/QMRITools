#!/usr/local/bin/perl

use strict;
use CGI;
use lib '/www/proj/libs';
use WRI::WebForms;

my $EMAIL_TO = $ENV{'WRI_IS_TESTING'}
		? 't-testmessages@wolfram.com' # Send here on dev/test
		: 'documentationfeedback@wolfram.com'; # Send here on production. Separate multiple emails with commas.

my $CONF = {
    webform_name => 'reference-feedback'
};

print qq(Content-Type: text/html\n\n);
my $form = CGI->new->Vars;
my $webform = WRI::WebForms->new('prd');
$webform->insertSubmissions($CONF->{'webform_name'}, $form);
send_email();


sub send_email {
    
    my $SENDMAIL = '/usr/sbin/sendmail';
    open(MAIL, "| $SENDMAIL -ti") || die "Error opening sendmail: $!";
    my $sel = select(MAIL);


    print <<EOF;
To: $EMAIL_TO
From: Reference Website <info\@wolfram.com>
Subject: Mathematica Documentation Feedback ($form->{'url'}) 

Someone filled out the feedback form on this page:
$form->{'url'}

Feedback: 
$form->{'feedback'}

Name (optional): 
$form->{'name'}

Email (optional): 
$form->{'email'}


All submissions can be found at:
https://bizi.wolfram.com/wbi/webforms/formData.cgi?formname=reference-feedback

EOF


    select ($sel);
    close(MAIL);

}
