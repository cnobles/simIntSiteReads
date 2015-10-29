for i in $(ls -d GTSP*/); do
    echo ${i%%//};
    Rscript ~/makeUCSCHubFromBedFiles/GTSP2BED_fromFolder.R ${i%%//};
done

mkdir bed
mv *.bed bed

cd bed
cat *.bed > caller.bed
cd ..

