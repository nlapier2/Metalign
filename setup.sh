#!/bin/bash
fileid="1iYEC9xwXReCoxeXOFLwoW82tmbCt9QR7"
filename="indices.zip"
curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid}" > /dev/null
curl -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid}" -o ${filename}
rm cookie
mv indices.zip data/
cd data
unzip indices.zip
rm indices.zip
cd ..

#
# Adapted from the following under the Creative Commons License: https://stackoverflow.com/questions/48133080/how-to-download-a-google-drive-url-via-curl-or-wget 
#
