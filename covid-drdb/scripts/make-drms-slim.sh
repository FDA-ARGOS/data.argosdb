#!/bin/bash

DB_FILE=$1

sqlite3 $DB_FILE <<EOF
DELETE FROM articles AS ref WHERE NOT EXISTS (SELECT 1 FROM resistance_mutation_articles rma WHERE ref.ref_name = rma.ref_name);
PRAGMA writable_schema = 1;
DELETE FROM sqlite_master WHERE name NOT LIKE '%resistance_mutation%' AND name NOT LIKE '%ref_amino_acid%' AND name NOT IN ('antibodies', 'compounds', 'articles');
PRAGMA writable_schema = 0;
VACUUM;
EOF
