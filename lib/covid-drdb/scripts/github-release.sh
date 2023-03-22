#! /bin/bash

set -e
export TZ=America/Los_Angeles

git config --global --add safe.directory /covid-drdb
git config --global --add safe.directory /covid-drdb-payload

PRE_RELEASE=
VERSION=
DRYRUN=
REPO=hivdb/covid-drdb-payload

while (($#)); do
  if [[ "$1" == "--prerelease" ]]; then
    PRE_RELEASE=$1
  elif [[ "$1" == "--dryrun" ]]; then
    DRYRUN=$1
  else
    VERSION=$1
  fi
  shift
done

if [ -z "$VERSION" ]; then
  if [[ "$PRE_RELEASE" == "--prerelease" ]]; then
    VERSION=$(date +"%Y%m%d-%H%M%S")
  else
    VERSION=$(date +"%Y%m%d")
  fi
fi

if [[ "$VERSION" =~ ^[0-9]{8}-[0-9]{6}$ ]]; then
  PRE_RELEASE="--prerelease"
  TODAY=$(echo "$VERSION $(date +%Z)" | sed -E 's/^([0-9]{4})([0-9]{2})([0-9]{2})-([0-9]{2})([0-9]{2})([0-9]{2})/\1-\2-\3 \4:\5:\6/g')
elif [[ "$VERSION" =~ ^[0-9]{8}$ ]]; then
  PRE_RELEASE=
  TODAY=$(echo "$VERSION" | sed -E 's/^([0-9]{4})([0-9]{2})([0-9]{2})/\1-\2-\3/g')
else
  echo "Release abort: invalid version format. Please use YYYYMMDD or YYYYMMDD-HHMMSS." 1>&2
  exit 1
fi

GIT="git -C payload/"

if gh release view --repo $REPO $VERSION > /dev/null 2>&1; then
  echo "Release abort: version $VERSION is already released." 1>&2
  if [[ "$PRE_RELEASE" != "--prerelease" ]]; then
    echo "You may want to 'make pre-release' for deploy a testing version." 1>&2
  fi
  exit 3
fi
builder_local_commit=$(git rev-parse HEAD)

if [[ "$PRE_RELEASE" == "--prerelease" ]]; then
  title="Pre-release $VERSION"
  description="Pre-release date: $TODAY\n\n
Usage of the \`.db\` files:\n\n
The \`.db\` files are SQLite3 databases. You can simply use this command to open them (if SQLite3 is installed):\n\n
\`\`\`bash\n
sqlite3 covid-drdb-$VERSION.db\n
\`\`\`\n\n
Once the SQLite shell prompts, type \`.tables\` to list all tables.\n\n
You can also use any SQLite viewer to open the database, e.g.:\n\n
- https://inloop.github.io/sqlite-viewer/\n
- https://sqlitebrowser.org/\n\n
Built with hivdb/covid-drdb@${builder_local_commit}\n"
else
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

  prev_tag=$(gh release list --exclude-pre-releases --repo $REPO -L 1 | cut -d$'\t' -f3)
  prev_commit=$($GIT rev-list -n 1 $prev_tag)

  title="COVID-DRDB $VERSION"
  description="Release date: $TODAY"

  if [[ "$prev_commit" != "null" ]]; then
    description="Release date $TODAY\n\nChanges since previous release ($prev_tag):\n
$($GIT log --pretty=format:'- %s (%H, by %an)\n' --abbrev-commit $prev_commit..$local_commit)\n\n
Usage of the \`.db\` files:\n\n
The \`.db\` files are SQLite3 databases. You can simply use this command to open them (if SQLite3 is installed):\n\n
\`\`\`bash\n
sqlite3 covid-drdb-$VERSION.db\n
\`\`\`\n\n
Once the SQLite shell prompts, type \`.tables\` to list all tables.\n\n
You can also use any SQLite viewer to open the database, e.g.:\n\n
- https://inloop.github.io/sqlite-viewer/\n
- https://sqlitebrowser.org/\n\n
Built with hivdb/covid-drdb@${builder_local_commit}\n"
  fi
fi

if [[ "$DRYRUN" == "--dryrun" ]]; then
  echo "Dry run: will not create a release."
  echo
  echo "Title: $title"
  echo "Tag: $VERSION"
  echo -e "Description:\n$description"
  exit 0
fi

scripts/export-sqlite.sh $VERSION

if [ ! -f "build/covid-drdb-$VERSION.db" ]; then
  echo "Release abort: file 'build/covid-drdb-$VERSION.db' is not found. Something is wrong, please contact Philip." 1>&2
  exit 2
fi

if [ ! -f "build/covid-drdb-$VERSION-slim.db" ]; then
  echo "Release abort: file 'build/covid-drdb-$VERSION-slim.db' is not found. Something is wrong, please contact Philip." 1>&2
  exit 2
fi

if [ ! -f "build/covid-drdb-$VERSION-variants.db" ]; then
  echo "Release abort: file 'build/covid-drdb-$VERSION-variants.db' is not found. Something is wrong, please contact Philip." 1>&2
  exit 2
fi

if [ ! -f "build/covid-drdb-$VERSION-drms.db" ]; then
  echo "Release abort: file 'build/covid-drdb-$VERSION-drms.db' is not found. Something is wrong, please contact Philip." 1>&2
  exit 2
fi

echo -e $description | gh release create --repo $REPO --title "$title" $PRE_RELEASE --notes-file - $VERSION
ls -1 build/covid-drdb-${VERSION}*.db | while read name; do
  gh release upload --repo $REPO $VERSION $name
done

if [[ "PRE_RELEASE" == "--prerelease" ]]; then
  echo "Pre-release $VERSION created: https://github.com/$REPO/releases/tag/$VERSION"
else
  echo "Release $VERSION created: https://github.com/$REPO/releases/tag/$VERSION"
fi
