for i in $(ls -d GTSP*/); do
    echo ${i%%//};
    Rscript ~/makeUCSCHubFromBedFiles/GTSP2BED_fromFolder.R ${i%%//};
done

