# FDA-ARGOS Data Dictionary

List of controlled vocabulary terms for ARGOSdb QC and biosample metadata table data properties.

This sheet was created to aid in the integration of ARGOS data from many disparate sources. Each of the column headers in each of the respective data sheets displayed on [data.argosdb.org](https://data.argosdb.org) was recommended by project members and collaborators from the FDA. The resulting list was manually curated to combine similar terms, and provide a consistent representation across all datasets in ARGOSdb.

The primary use case for the Data Dictionary is to ensure all data submitted to data.argosdb.org is following a consistent representation of the data properties.

## To generate the Schema files

The draft dictionary is maintained as a Excell file for project members to contribute to. After a final review that sheet is downloaded and each of the respective tabs are converted to a `tsv` and are saved in the respective versiond directory here. Each version should contain the following files:

- `release_notes.tsv` => An itemized list of the changes implemented in the current version.
- `README.tsv` => A summary of each sheet and the column headers
- `property_definition.tsv` => List of all controlled vocabulary terms for ARGOSdb
- `annotation_property_list.tsv` => List of all non-core table data properties
- `core_property_list.tsv` => List of all core table data properties

### Running the `write_schema` function

Ensure that the virtual enviornment is activated and updated. From the repository root run:
```shell
sh env/bin/activate
pip install -r requirements.txt
```

The `dictionary_utils.py` script contains a few different functions. They are listed by simply invoking the script with out any additional parameters:
```shell
(env) data.argosdb$: python lib/dictionary_utils.py 
usage: argosdb_dict_utils [options]

positional arguments:
  {functions,write_schema,validate_columns}
    functions           List of all available functions
    write_schema        Used to convert a TSV into a JSNO schema. If no mapping file
                        is provided, performs default conversions.
    validate_columns    Validates columns in a list of files in a directory using
                        provided column headers

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```

To prepare the new schema directories run the following commands from the project root:

```shell
mkdir schema/[CURRENT NEW VERSION]
mkdir schema/[CURRENT NEW VERSION]/core
mkdir schema/[CURRENT NEW VERSION]/annotation
```

To generate the schema files you need to run the `dictionary_utils.py` twice. Once for the CORE and once for the ANNOTATION sheet. 
```shell
'-m', '--multi',
        'Flag to indicate if multiple items are being processed'
'-i', '--input',
        required=True,
        "Data shet to validate. This can be a 'tsv' or 'csv'. "
            "Other file types are not permitted."
            "If the 'multi' flag is supplied then this should be a directory"
            "path for a group of files"
'-d', '--directory',
        "Directory"
'-f', '--definitions',
        "Definitions"
'-s', '--schema',
        "Root json schema to validate against."
```

For the Core tables run: 

`python lib/dictionary_utils.py write_schema -i data_dictionary/[current version]/core_property_list.tsv -d schema/[current version]/core/ -f data_dictionary/[current version]/property_definition.tsv`

For the Annotation tables run: 

`python lib/dictionary_utils.py write_schema -i data_dictionary/[current version]/annotation_property_list.tsv -d schema/[current version]/annotation/ -f data_dictionary/[current version]/property_definition.tsv`

## Schema Tests (using `data_sheet_validator.py`)

The schema tests take a `CSV` or `TSV` and convert them to a `.json` file. Each line in the data sheet becomes a JSON object that is validated by the supplied schema. 

**The tests are run in the ARGOS server because that is where the data files live, but the tests can be used on any individual file**

Ensure that the virtual enviornment is activated and updated. From the repository root run:
```shell
sh env/bin/activate
pip install -r requirements.txt
```

### Setting up the Test Directories

We need a versioned directory for the test results and test files. run the following commands to create them:

	mkdir tests/[current version]
        mkdir tests/[current version]/test_results
        mkdir tests/[current version]/test_files

### Gathering Test files

For the testing a new dictionary version it is best to use an abriveated table for each data type. To do this we run [test_file_generation.sh](/lib/shell/test_file_generation.sh) This is a simple script to take the HEAD of each data file to make the test files. You need to provide 2 directories as the arguments: 
1. first directory is the location of files to copy to test.
2. second directory is the location you want the test files to be generated in.

Example command:

        sh lib/shell/test_file_generation.sh /data/shared/argosdb/generated/datasets/reviewed/ tests/v1.4/test_files/

### Running the `data_sheet_validator.py` function

once the testing files are inplace run the following: 

        python3 lib/data_sheet_validator.py -m -i tests/v1.4/test_files/ -s schema/v1.4/ -o tests/v1.4/test_results/ > tests/v1.4/test_results/summary.txt
