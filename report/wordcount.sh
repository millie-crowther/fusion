DUMMY=`./build.sh`
WORDS=`pdftotext report.pdf - | wc -w`
echo "wordcount: ${WORDS}"
echo "target is 30000"
