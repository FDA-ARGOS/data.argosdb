#!/bin/bash

echo "Finished downlaoding data from ARGOS API. Attempting to push data sets to HIVE 3...\n"
read -p "Enter your login for HIVE 3: " user
read -p "Enter your password" pass
curl -k -c ~/hive3.cookie "https://hive3.biochemistry.gwu.edu/dna.cgi?cmdr=login&login=${user}&pswd=${pass}"

# Log into HIVE 3:
# curl -k -c ~/hive3.cookie 'https://hive3.biochemistry.gwu.edu/dna.cgi?cmdr=login&login=argos@argos.org&pswd=qQ123123'
# Remember to get input from user for login/pass!!!

for DATASET in ~/scripts/*.json;
#for DATASET in ~/api/data/*.json;
do
        curl -v -b  ~/hive3.cookie  -X POST -F "content=@${DATASET}" -F 'cmd=objSetFile' -F 'raw=1'  -F 'bin=0' -F 'type=u-file' -F 'filename=_.json' -F 'name=apiTestFile'  -F 'ext=json'  'https://hive3.biochemistry.gwu.edu/dna.cgi';
done
#        curl -v -b  ~/hive3.cookie  -X POST -F 'content=@/home/vboxuser/scripts/test.json' -F 'cmd=objSetFile' -F 'raw=1'  -F 'bin=0' -F 'type=u-file' -F 'filename=_.json' -F 'name=apiTestFile'  -F 'ext=json'  'https://hive3.biochemistry.gwu.edu/dna.cgi';

# The @content tag is where you specify the file.

