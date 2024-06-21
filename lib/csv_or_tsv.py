# Christie Woodside
# June 17, 2024

'''Code test to see if the output file is a csv or tsv. If there is one column then it is a csv and 
if there are more than one then it was successfully changed to a tsv'''

''' Taken fro github: https://stackoverflow.com/questions/77331752/check-if-a-file-is-a-tsv-file'''

import csv
# paste in the filename url in each spot. Leave the "" blank at the end
for filename in ["in.csv", "in.tsv", ""]:
    try:
        with open(filename, "r", newline="") as file_in:
            reader = csv.reader(file_in, delimiter="\t")
            headers = next(reader)

            column_count = len(headers)
            if column_count == 1:
                raise ValueError("Column Count is 1")
            
            ## -------------------
            ## inspect some of the rows?
            ## -------------------
            ## for row in reader:
            ##     ....
            ## -------------------

    except Exception as e:
        print(f"file \"{filename}\" is probably not a tsv")
        # print(e)
        continue

    print(f"file \"{filename}\" appears to have {column_count} columns")


