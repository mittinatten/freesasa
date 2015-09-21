#!/usr/bin/perl
use strict;
use XML::DOM;

while (<>) {
    if (/<\/head>/) {
        print '<link href="stylesheets/styles.css" rel="stylesheet" type="text/css">'."\n";
        print '<link href="stylesheets/github-light.css" rel="stylesheet" type="text/css">'."\n";
        print '<link href="stylesheets/article.css" rel="stylesheet" type="text/css">'."\n";
    }
    $_ =~ s/style=".*"/style=""/;
    $_ =~ s/ltx_page_main/ltx_page_main wrapper/;
    $_ =~ s/--/&en;/;
    print;
}
exit;
