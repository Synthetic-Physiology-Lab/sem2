#!/bin/bash

# FIX IT
# -> ignore LAMMPS updates -> different .sh
# -> overwrite src with src_changed -> + WARNING


# this is called by Install.sh when trying to include changes in src
# -> it will check that we don't screw up things

# init (all paths relative to curdir)
curdir=`pwd`
srcdir=..
origdir="src_orig"
changeddir="src_changed"

# check files in src
# already_patched == src is same as src_changed
# orig_diff == src is different from src_orig
# -> triggers LAMMPS-update handling iff new_sem=false
# new_sem == it seems that USER-SEM is there (triggered by "//@@@SEM" comment)
# -> only set if files differ from src_changed (i.e. already_patched = false)
already_patched=true
orig_diff=false
new_sem=false
OLDIFS=$IFS
IFS=$'\n'
for i in $( ls "$changeddir"/* ); 
do
  filename=$(basename $i)
  
  filegrep_sem=`grep "//@@@SEM" "$srcdir/$filename"`
  if [[ "$filegrep_sem" != "" ]]; then
    new_sem=true
  fi
  
  filediff_orig=`diff -Na "$origdir/$filename" "$srcdir/$filename"`
  filediff_changed=`diff -Na "$changeddir/$filename" "$srcdir/$filename"`  
  if [[ "$filediff_changed" != "" ]]; then
    already_patched=false
    if [[ "$filediff_orig" != "" ]]; then
      if [[ "$new_sem" != "true" ]]; then
        echo "!!! WARNING: $filename differs from expected one !!!"
      fi
      orig_diff=true
    fi
  fi
done
IFS=$OLDIFS

# only continue if we are not already patched
if [ "$already_patched" == "true" ]; then
  echo "src files are already patched => nothing to be done"
else
  # careful when files differ

  # guess1: user changed src_changed, src_orig is correct but USER-SEM is already installed
  if [ "$new_sem" == "true" ]; then
    echo "I assume that you had USER-SEM installed and just changed src_changed."
    echo "I will overwrite your src-content with src_changed."
    echo "A backup of the files in src will be created in old/src."
    echo "Proceed? (y/n)"
    read proceed
    if [[ "$proceed" != "y" ]]; then
      echo "ABORTING"
      exit
    fi
    # backup
    rm -rf "old/src"
    mkdir -p "old/src"
    # get new files
    OLDIFS=$IFS
    IFS=$'\n'
    for i in $( ls "$changeddir"/* ); 
    do
      filename=$(basename $i)
      cp "$srcdir/$filename" "old/src"
      if [ ! -e "$origdir/$filename" ]; then
        echo "WARNING: $filename was missing in src_orig (I added it)!!"
        cp "$srcdir/$filename" "$origdir"
      fi
    done
    IFS=$OLDIFS

  # guess2: there was a LAMMPS update -> try to patch it
  elif [ "$orig_diff" == "true" ]; then
    echo "!!! Original src-files not as expected => probably new LAMMPS version used !!!"
    echo "I can try to apply patch to new files (backup kept in old-folder)."
    echo "Proceed? (y/n)"
    read proceed
    if [[ "$proceed" != "y" ]]; then
      echo "ABORTING"
      exit
    fi
    # keep backup
    mkdir old
	cd old
    rm -rf "$origdir"
    rm -rf "$changeddir"
    rm src_patch
    mkdir "$origdir"
    mkdir "$changeddir"
	cd ..
	cp "$origdir"/* "old/$origdir" 2>/dev/null
	cp "$changeddir"/* "old/$changeddir" 2>/dev/null
    # get patch
    LC_ALL=C
    TZ=UTC0
    cd old
    diff -Naur "$origdir" "$changeddir" > src_patch
    cd ..
    # get new files
    OLDIFS=$IFS
    IFS=$'\n'
    for i in $( ls "old/$changeddir"/* ); 
    do
      filename=$(basename $i)
      cp "$srcdir/$filename" "$origdir"
      cp "$srcdir/$filename" "$changeddir"
    done
    IFS=$OLDIFS
    # apply patch
    cd "$changeddir"
    patch -p1 --no-backup-if-mismatch < "$curdir/old/src_patch"
    cd "$curdir"
  fi
  
  # now the easy part -> overwrite src-files with changed ones
  cp "$changeddir"/* "$srcdir"
fi
