#!/bin/bash
#
# This script is used to update the version number of MPSolve,
# and/or comments regarding the copyright years.
#

function update_year {
  # Get the year currently in the header files.
  source_year=$(grep "2001-" include/mps/mps.h | cut -d '-' -f2 | cut -d ',' -f1)
  current_year=$(date +'%Y')

  echo "YEAR UPDATE"

  echo "> Year in source files: ${source_year}"
  echo "> Current year: ${current_year}"

  if [ "$source_year" != "$current_year" ]; then
    echo "> Updating the Copyright statement with years 2001-${current_year}"
  fi
}

function update_version {
  source_version=$(grep "This file is part of MPSolve" include/mps/mps.h | cut -d ' ' -f9)
  current_version=$1

  echo "VERSION UPDATE"

  echo "> Version in source files: ${source_version}"
  echo "> Current version: ${current_version}"

  if [ "$1" = "" ]; then
    echo "> Please specify a valid version to migrate to"
  else
    echo "> Updating MPSolve version from ${source_version} to ${current_version}"

    # Update version in source and header files
    for file in $(find src include -name \*.c -or -name \*.cpp -or -name \*.h -or -name \*.h.in); do
      sed -i "s/MPSolve ${source_version}/MPSolve ${current_version}/g" $file
    done
    sed -i "s/MPSOLVE ${source_version}/MPSOLVE ${current_version}/g" README

    echo "> Updated source file, please take care of configure.ac manually"
  fi
}

update_year

if [ "$1" != "" ]; then
  update_version $1
fi
