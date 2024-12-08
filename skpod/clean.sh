#!/bin/sh

pushd "$(dirname -- "$0")"

rm -f build.sh
[ -d build_prog ] && cd build_prog && ( rm -f * || true ) && cd ..
[ -d gen_prog ] && cd gen_prog && ( rm -f * || true ) && cd ..
[ -d configs ] && cd configs && ( rm -f * || true ) && cd ..
rm -f run.sh

popd
