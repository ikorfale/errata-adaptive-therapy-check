#!/bin/sh
# Download the public trial data (not redistributed here) and check it is the file this study used.
set -e
mkdir -p data/bruchovsky
curl -sSL -o data/dataTanaka.zip https://www.nicholasbruchovsky.com/dataTanaka.zip
echo "845683954e99253ddf68e488a48223e89c19247a92aec27daa79130d8f46904f  data/dataTanaka.zip" | sha256sum -c -
python3 -c "import zipfile; zipfile.ZipFile('data/dataTanaka.zip').extractall('data/bruchovsky')"
