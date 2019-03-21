#! /usr/local/bin/perl
use CGI; 
BEGIN {
        print qq(Content-type: text/html\n\n);
        open(STDERR,'>&STDOUT'); $| = 1;
        $SIG{'__DIE__'} = sub {
                print qq(<p class="brightred"><b>Error: @_</b></p>);
                exit;};}
		
my $query = CGI->new();
my $q = CGI::escape($query->param('q'));
my $unescapedQ = CGI::unescape($query->param('q'));
my $lang = $query->param('lang');
my $request_uri = $ENV{'REQUEST_URI'};

if ($q && !($request_uri =~ '/search/')) {
print <<EOF;
	<div id="searchlink">
		<p>
		Search for all pages containing <a href="/search/?q=$q&forward=0">$unescapedQ</a>
		</p>
	</div>
EOF
} elsif($q && ($request_uri =~ '/search/')) {
print <<EOF;
<div id="searchlink">
   <p>
	<span id="spelling_container"></span>
  </p>
</div>	
EOF
}
