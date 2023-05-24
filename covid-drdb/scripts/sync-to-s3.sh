#! /bin/bash

set -e

VERSION=$1
S3_BUCKET=s3://cms.hivdb.org
S3_PREFIX=covid-drdb

mkdir /github-assets
pushd /github-assets >/dev/null
gh release download --repo hivdb/covid-drdb-payload $VERSION
ls -1 | while read name; do
  echo "compress: $name"
  pigz -9 $name
  mv $name.gz $name
done
popd >/dev/null

aws s3 sync /github-assets/ "${S3_BUCKET}/${S3_PREFIX}/" \
  --content-encoding gzip --cache-control max-age=2592000
