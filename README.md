Xeva

Integration of molecular and pharmacological profiles of patient-derived models

How to install:

At present there are 2 ways to use/test Xeva:

1. Download zip file from Xeva github page. Unzip it. and install 
    install("Xeva-master")

2. Using install_github : To do this you will required devtools. Install devtools in R as:
    install.packages("devtools") 
    
    Once devtools installed load library:
    library(devtools)
    
    As Xeva is yet in privet reposeatry, you will required a a personal access tokens on github. 
    Go to the url https://github.com/settings/tokens and login to your account. Type a Token description
    and click on Generate Token. Copy the token text.
    Now in R:
    install_github("bhklab/Xeva", auth_token = "your access token")
    
    It will install Xeva to your system. 


3. Other way to use Xeva is to copy the Xeva as R project on your local machine. 
    On the Xeva github page click on the 


To install Xeva 
Download the zip
https://github.com/bhklab/Xeva/archive/master.zip
