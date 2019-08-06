#!/bin/bash
fileid="1f_zbd9JqVsIgOe-wwUvytfIBoiTOa9T2"
filename="db_info.txt.zip"
curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid}" > /dev/null
curl -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid}" -o ${filename}
rm cookie

fileid="1NXxC7uQ2xc7mIQ1G8f-2ZVW7RXqdM32p"
filename="cmash_db.zip"
curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid}" > /dev/null
curl -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid}" -o ${filename}
rm cookie

mv db_info.txt.zip cmash_db.zip data/
cd data
unzip db_info.txt.zip
unzip cmash_db.zip
rm db_info.txt.zip cmash_db.zip
cd ..

git clone https://github.com/lh3/minimap2.git && cd minimap2 && make && cd ..
git clone https://github.com/dkoslicki/CMash.git && cd CMash && pip install --user -r requirements.txt && cd ..

#
# Adapted from the following under the Creative Commons License: https://stackoverflow.com/questions/48133080/how-to-download-a-google-drive-url-via-curl-or-wget
#
