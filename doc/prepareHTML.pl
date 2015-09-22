#!/usr/bin/perl
use strict;

# This is a total hack, but ran into problems with MathML using
# HTML::DOM and problems with HTML5 using XML::DOM. Since it's just a
# few simple things that needs to be done on one particular file, this
# should do

while (<>) {
    if (/<\/head>/) {
        print '<link href="stylesheets/styles.css" rel="stylesheet" type="text/css">'."\n";
        print '<link href="stylesheets/github-light.css" rel="stylesheet" type="text/css">'."\n";
        print '<link href="stylesheets/article.css" rel="stylesheet" type="text/css">'."\n";
    }
    $_ =~ s/<nav .*>/<header>/;
    $_ =~ s/<\/nav>/<\/header>/;
    $_ =~ s/style=".*"/style=""/;
    $_ =~ s/ltx_page_main/ltx_page_main wrapper/;
    print;
}
exit;
