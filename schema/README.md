# ARGOS DB Schema versions

1) copy excel tabs as TSV files to data_dictionary/<version>/:
    a) release_notes.tsv
    b) README.tsv
    c) property_definition.tsv
    d) non-core_property_list.tsv
    e) core_property_list.tsv
2) from root run `python lib/dictionary_utils.py write_schema -i data_dictionary/<version>/non-core_property_list.tsv -d schema/<version>/non-core/ -f data_dictionary/<version>/property_definition.tsv`

3) from root run `python lib/dictionary_utils.py write_schema -i data_dictionary/<version>/core_property_list.tsv -d schema/<version>/core/ -f data_dictionary/<version>/property_definition.tsv`

4) from root run `python lib/dict_release_notes.py -o data_dictionary/<old_version> -n data_dictionary/<new_version>` and add output to the release notes. 

5) 