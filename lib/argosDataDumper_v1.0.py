import getpass
import json
import os
import requests

basename = "https://hive3.biochemistry.gwu.edu/"

def apiLogin():
   global cookies

   user = input("\nEnter your username on HIVE3 (email address): ")
   password  = getpass.getpass()
   loginParams = {'api': '0', 'cmdr': 'login', 'login': user, 'pswd': password}
   response = requests.get(basename, params=loginParams)
   if response.text != '':
      print(response.text)
      apiLogin()  # Try again
   else:
      cookies = response.cookies

apiLogin()




