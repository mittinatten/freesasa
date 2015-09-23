#!/usr/bin/perl
use strict;

# This is a total hack, but ran into problems with MathML using
# HTML::DOM and problems with HTML5 using XML::DOM. Since it's just a
# few simple things that needs to be done on one particular file, this
# should do

my $tracking = <<HERE;
<script type="text/javascript"
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  MathML: {
    useMathMLspacing: false
  }
});
MathJax.Hub.Config({
  "HTML-CSS": {
    availableFonts: ["Latin-Modern"],
    webFont: "Latin-Modern",
    scale: 90
  }
});
</script>
<script>
    (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
        (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
        })(window,document,'script','//www.google-analytics.com/analytics.js','ga');

     ga('create', 'UA-66397346-1', 'auto');
     ga('send', 'pageview');
</script>
HERE

while (<>) {
    if (/<\/head>/) {
        print '<link href="stylesheets/styles.css" rel="stylesheet" type="text/css">'."\n";
        print '<link href="stylesheets/github-light.css" rel="stylesheet" type="text/css">'."\n";
        print '<link href="stylesheets/article.css" rel="stylesheet" type="text/css">'."\n";
    }
    print $tracking if (/<\/body>/);
    $_ =~ s/<nav .*>/<header>/;
    $_ =~ s/<\/nav>/<\/header>/;
    $_ =~ s/style=".*"/style=""/;
    $_ =~ s/ltx_page_main/ltx_page_main wrapper/;
    print;
}
exit;
