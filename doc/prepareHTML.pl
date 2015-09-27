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
     var manualLink = document.getElementById('manual');
     addListener(manualLink, 'click', function() {ga('send', 'event', 'manual', 'click');});
     function addListener(element, type, callback) {
         if (element.addEventListener) element.addEventListener(type, callback);
         else if (element.attachEvent) element.attachEvent('on' + type, callback);
     }
</script>
HERE

my $style = <<HERE;
  <link href="stylesheets/styles.css" rel="stylesheet" type="text/css">
  <link href="stylesheets/github-light.css" rel="stylesheet" type="text/css">
  <link href="stylesheets/article.css" rel="stylesheet" type="text/css">
HERE

my $pdf = '<p>[ <a href="manual.pdf" id="manual">Download PDF</a> ]</p>';

while (<>) {
    print $style if (/<\/head>/);
    print $tracking if (/<\/body>/);
    $_ =~ s/<nav .*>/<header>/;
    $_ =~ s/<\/nav>/$pdf\n<\/header>/;
    $_ =~ s/style=".*"/style=""/;
    $_ =~ s/ltx_page_main/ltx_page_main wrapper/;
    print;
}
exit;
