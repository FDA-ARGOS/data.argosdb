#! /bin/bash

set -e

VERSION=$1

pipenv run python scripts/db_to_sqlite.py "postgresql://postgres@covid-drdb-devdb:5432/postgres" /dev/shm/covid-drdb-$VERSION.db --all
echo "Written build/covid-drdb-$VERSION.db"
ln -s covid-drdb-$VERSION.db /dev/shm/covid-drdb-latest.db
ls ./views/*.sql | sort -h | while read filepath; do
    sqlite3 /dev/shm/covid-drdb-latest.db < $filepath
done
echo "build/covid-drdb-latest.db -> covid-drdb-$VERSION.db"

cp /dev/shm/covid-drdb-$VERSION.db /dev/shm/covid-drdb-$VERSION-slim.db
./scripts/make-slim-version.sh /dev/shm/covid-drdb-$VERSION-slim.db
sqlite3 /dev/shm/covid-drdb-$VERSION-slim.db < drop_views.sql
echo "Written build/covid-drdb-$VERSION-slim.db"

cp /dev/shm/covid-drdb-$VERSION.db /dev/shm/covid-drdb-$VERSION-variants.db
./scripts/make-variants-slim.sh /dev/shm/covid-drdb-$VERSION-variants.db
echo "Written build/covid-drdb-$VERSION-variants.db"

cp /dev/shm/covid-drdb-$VERSION.db /dev/shm/covid-drdb-$VERSION-drms.db
./scripts/make-drms-slim.sh /dev/shm/covid-drdb-$VERSION-drms.db
echo "Written build/covid-drdb-$VERSION-drms.db"

mkdir -p build/
mv /dev/shm/*.db build/
