#!/bin/bash

set -euo pipefail
IFS=$'\n\t'

rm -f sources.tar.gz
tar -czf sources.tar.gz \
  LICENCE \
  models/ \
  python/main.py python/eval.py \
  java/src/ java/lib/ \
  java/gradle/ java/gradlew java/gradlew.bat java/build.gradle java/settings.gradle