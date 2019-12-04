# FDS Summary Web Page

This directory contains a webpage template for summarizing the results  of a firebot run. 
The summary consists of png images created by firebot for the FDS user and verification guides
and links to all of the guides created by firebot.
The summary may viewed (after firebot is run) by opening the web page `index.html`
found in this directory.
This summary pages may also be viewd from other computers if the -w option is invoked when running firebot.  
At NIST, firebot is run using `-w /var/www/html/firebot` and then accessed using the URL http://blaze.el.nist.gov/firebot

firebot copies all png image files generated for the FDS User and Verification guides
to the images/user and images/verification subdirectories of FDS_Summary.
Many but not all of these images are displayed in the `index.html`
web page. To display other images, add desired links to `index_template.html` and rerun firebot.
In addition to copying image files and guides, 
firebot copies `index_template.html` to `index.html` while adding the current date and FDS and SMV repo revisions.
Note, firebot does not write to index.html directly to avoid modifying a file that is in the fds repo.
