#!/bin/sh
# Check out the git submodules listed in tools/submodule-pins.
#
# GitHub source tarballs (remotes, renv, pak, and codeload archives) contain
# empty submodule directories. remotes::install_github() then shallow-clones
# each submodule's default branch, which is not the pinned commit for the
# optimotu forks, and R CMD build drops the clone's .git before configure
# runs. This script replaces a checkout that is missing or has the wrong
# content, and leaves a correct vendored tree or a matching git checkout
# alone.
#
# Usage:
#   sh tools/fetch-submodules.sh
#   sh tools/fetch-submodules.sh --check
#   sh tools/fetch-submodules.sh --print-digests

set -u

root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd) || exit 1
cd "$root" || exit 1
pins=$root/tools/submodule-pins
# /tmp on this host is a small volume and may be full. Prefer TMPDIR,
# then ~/tmp, which is where R builds are pointed.
blobfile=
for scratch in ${TMPDIR:+"$TMPDIR"} "$HOME/tmp" /tmp; do
  [ -n "$scratch" ] || continue
  mkdir -p "$scratch" 2>/dev/null || continue
  if blobfile=$(mktemp "$scratch/optimotu-blob.XXXXXX" 2>/dev/null); then
    break
  fi
done
if [ -z "$blobfile" ]; then
  echo "fetch-submodules: need a writable temp directory (set TMPDIR)" >&2
  exit 1
fi
trap 'rm -f "$blobfile"' EXIT
mode=fetch

case "${1:-}" in
  "") mode=fetch ;;
  --check) mode=check ;;
  --print-digests) mode=print ;;
  *)
    echo "usage: sh tools/fetch-submodules.sh [--check|--print-digests]" >&2
    exit 2
    ;;
esac

if [ ! -f "$pins" ]; then
  echo "missing $pins" >&2
  exit 1
fi

# Paths omitted from an R source tarball (.Rbuildignore, plus the .gitignore
# files R CMD build always drops). The digest covers every other tracked file.
skip_path() {
  case "$1" in
    .git | .git/* | .github | .github/* | .vscode | .vscode/* | .travis.yml | \
      .gitignore | */.gitignore)
      return 0
      ;;
  esac
  return 1
}

# R CMD build drops CR from text and appends a final newline. Hash the
# same normalized bytes so a built tarball matches the git commit.
# Binary files (embedded NUL) are hashed unchanged.
hash_normalized() {
  file=$1
  raw=$(wc -c < "$file" | tr -d ' ')
  stripped=$(tr -d '\0' < "$file" | wc -c | tr -d ' ')
  if [ "$raw" != "$stripped" ]; then
    git hash-object "$file"
    return
  fi
  # Drop a trailing CR on each line and always end with LF, matching
  # R CMD build's text cleanup.
  awk '{ sub(/\r$/, ""); print }' "$file" | git hash-object --stdin
}

# Digest of the committed tree (not a dirty worktree).
digest_commit() {
  dir=$1
  git -C "$dir" ls-tree -r --name-only HEAD |
    while IFS= read -r path; do
      if skip_path "$path"; then
        continue
      fi
      git -C "$dir" cat-file blob "HEAD:$path" > "$blobfile" || return 1
      blob=$(hash_normalized "$blobfile") || return 1
      printf '%s %s\n' "$blob" "$path"
    done |
    LC_ALL=C sort |
    git hash-object --stdin
}

# Digest of files on disk, using the same exclusions as a source tarball.
digest_worktree() {
  dir=$1
  find "$dir" \
    \( -name .git -o -name .github -o -name .vscode \) -prune -o \
    -type f -print |
    while IFS= read -r file; do
      path=${file#"$dir"/}
      if skip_path "$path"; then
        continue
      fi
      blob=$(hash_normalized "$file") || return 1
      printf '%s %s\n' "$blob" "$path"
    done |
    LC_ALL=C sort |
    git hash-object --stdin
}

fetch_sha() {
  dest=$1
  url=$2
  sha=$3
  rm -rf "$dest"
  mkdir -p "$dest"
  git -c init.defaultBranch=master init -q "$dest" || return 1
  git -C "$dest" remote add origin "$url" || return 1
  if git -C "$dest" fetch -q --depth 1 origin "$sha"; then
    git -C "$dest" checkout -q --detach FETCH_HEAD || return 1
  else
    echo "shallow fetch of $sha failed; cloning $url" >&2
    rm -rf "$dest"
    git clone -q "$url" "$dest" || return 1
    git -C "$dest" checkout -q --detach "$sha" || return 1
  fi
  got=$(git -C "$dest" rev-parse HEAD) || return 1
  if [ "$got" != "$sha" ]; then
    echo "checkout of $dest is $got, expected $sha" >&2
    return 1
  fi
}

fail() {
  echo "fetch-submodules: $*" >&2
  exit 1
}

# Windows checkouts may deliver this pin file with CRLF. A trailing
# CR would make an otherwise matching digest compare unequal.
cr=$(printf '\r')
while IFS= read -r line || [ -n "$line" ]; do
  line=${line%"$cr"}
  case "$line" in
    "" | \#*) continue ;;
  esac
  path= url= sha= sentinel= digest=
  IFS=' ' read -r path url sha sentinel digest <<EOF
$line
EOF
  if [ -z "$path" ] || [ -z "$url" ] || [ -z "$sha" ] || [ -z "$sentinel" ]; then
    fail "bad pin line: $line"
  fi
  case "$path" in
    src/*) ;;
    *) fail "refusing submodule path $path" ;;
  esac
  case "$path" in
    *..*) fail "refusing submodule path $path" ;;
  esac

  if [ "$mode" = "print" ]; then
    if [ ! -e "$root/$path/.git" ]; then
      fail "$path is not checked out; cannot print its digest"
    fi
    head=$(git -C "$root/$path" rev-parse HEAD) || fail "cannot read $path"
    if [ "$head" != "$sha" ]; then
      fail "$path is at $head, pins have $sha"
    fi
    got=$(digest_commit "$root/$path") || fail "digest failed for $path"
    printf '%s %s %s %s %s\n' "$path" "$url" "$sha" "$sentinel" "$got"
    continue
  fi

  if [ -z "$digest" ]; then
    fail "$path is missing its content digest"
  fi

  if [ -e "$root/.git" ]; then
    recorded=$(git -C "$root" ls-tree HEAD -- "$path" | awk '{print $3}')
    if [ "$recorded" != "$sha" ]; then
      fail "$path gitlink is ${recorded:-missing}, pins have $sha"
    fi
  elif [ "$mode" = "check" ]; then
    fail "--check requires a git checkout of optimotu"
  fi

  dest=$root/$path
  if [ -e "$dest/.git" ]; then
    head=$(git -C "$dest" rev-parse HEAD 2>/dev/null || true)
    if [ "$head" = "$sha" ]; then
      got=$(digest_commit "$dest") || fail "could not digest $path"
      if [ "$got" != "$digest" ]; then
        fail "$path digest is $got, pins have $digest"
      fi
      continue
    fi
    if [ "$mode" = "check" ]; then
      fail "$path is at ${head:-missing}, expected $sha"
    fi
  else
    if [ -f "$dest/$sentinel" ]; then
      got=$(digest_worktree "$dest") || fail "could not digest $path"
      if [ "$got" = "$digest" ]; then
        continue
      fi
    fi
    if [ "$mode" = "check" ]; then
      fail "$path content does not match $sha"
    fi
  fi

  echo "fetching $path at $sha" >&2
  fetch_sha "$dest" "$url" "$sha" || fail "could not fetch $url at $sha"
  if [ ! -f "$dest/$sentinel" ]; then
    fail "$path is missing $sentinel after fetch"
  fi
  got=$(digest_commit "$dest") || fail "could not digest fetched $path"
  if [ "$got" != "$digest" ]; then
    fail "fetched $path digest is $got, pins have $digest"
  fi
done < "$pins"
