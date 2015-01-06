# Extract theory SD code and turn it into a CMSSW library.
# Only needed to do this when getting a new version of SD.

echo "Extracting Code"
tar xfvz SD-tosvn-v0.6.tar.gz 


echo "Moving Code to right place"
cd SD-tosvn
mv *.h ../SDAlgorithm/interface/
mv *.cxx ../SDAlgorithm/src/
cd ..
rm -rf SD-tosvn


echo "Renaming cxx -> cc"
cd SDAlgorithm/src
for i in *.cxx; do mv "$i" "`basename $i .cxx`.cc"; done


echo "Removing TestAnalysis (we can't use it at the moment)"
rm -v TestAnalysis.cc


echo "Fix all include statements (this may take some time)"
source fix_includes.sh
cd ../interface
source fix_includes.sh


echo "We're done"
cd ../..