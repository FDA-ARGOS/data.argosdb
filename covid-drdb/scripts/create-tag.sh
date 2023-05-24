#! /bin/bash

set -e
export TZ=America/Los_Angeles

PRE_RELEASE=
DRYRUN=

while (($#)); do
  if [[ "$1" == "--prerelease" ]]; then
    PRE_RELEASE=$1
  elif [[ "$1" == "--dryrun" ]]; then
    DRYRUN=$1
  fi
  shift
done

if [[ "$PRE_RELEASE" == "--prerelease" ]]; then
  VERSION=$(date +"%Y%m%d-%H%M%S")
else
  VERSION=$(date +"%Y%m%d")
fi

GIT="git -C payload/"

if gh release view --repo hivdb/covid-drdb-payload $VERSION > /dev/null 2>&1; then
  echo "Release abort: version $VERSION is already released." 1>&2
  if [[ "$PRE_RELEASE" != "--prerelease" ]]; then
    echo "You may want to 'make pre-release' for deploy a testing version." 1>&2
  fi
  exit 3
fi
builder_local_commit=$(git rev-parse HEAD)

if [ -z "$PRE_RELEASE" ]; then
  builder_remote_commit=$(git rev-parse HEAD --branches=origin/master)
  remote_commit=$($GIT rev-parse HEAD --branches=origin/master)
  local_commit=$($GIT rev-parse HEAD)
  if [[ "$builder_remote_commit" != "$builder_local_commit" ]]; then
    echo "Release abort: the local repository covid-drdb seems not up-to-date.  Forgot running 'git pull --rebase' and 'git push'?" 1>&2
    exit 1
  fi
  if [[ "$remote_commit" != "$local_commit" ]]; then
    echo "Release abort: the local repository covid-drdb-payload seems not up-to-date. Forgot running 'git pull --rebase' and 'git push'?" 1>&2
    exit 1
  fi

  if [ -n "$(git status -s .)" ]; then
    git status
    echo "Release abort: uncommitted changes are found. Please submit them & run 'git push' first." 1>&2
    exit 1
  fi

  if [ -n "$($GIT status -s .)" ]; then
    $GIT status
    echo "Release abort: uncommitted changes are found under payload/ directory. Please submit them & run 'git push' first." 1>&2
    exit 1
  fi

fi

if [[ "$DRYRUN" == "--dryrun" ]]; then
  echo "Dry run: skip creating tag $VERSION" 1>&2
  exit 1
fi

echo -n $VERSION
