# Redirect the output of Robel's dataset-maker/check-dataset-sanity.py to a file and run this script on it to get a count of how many errors there are by table.

import json

with open ("dataset-sanity-log-10.24.2024", "r") as inFile:
    myObj = json.load(inFile)

myList = []

for key in myObj:
    values = myObj[key]
    for value in values:
        myList.append(value)
    print (key + ", Number of errors: " + str(len(myList)))

