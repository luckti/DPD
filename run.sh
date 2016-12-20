#! /bin/bash
echo 'Start time:' |tee log.txt
date  |tee -a log.txt

chmod +x solver preproces postproces

./preproces |tee -a log.txt

echo '==========================start solver===================================='  |tee -a log.txt
./solver |tee -a log.txt

echo '======================start postprocess==================================='  |tee -a log.txt
./postproces |tee -a log.txt

echo 'End time:' |tee -a log.txt
date  |tee -a log.txt