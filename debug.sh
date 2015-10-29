for i in $(ls -d GTSP*/); do
    echo ${i%%//};
    Rscript ~/makeUCSCHubFromBedFiles/GTSP2BED_fromFolder.R -t all ${i%%//};
    Rscript ~/makeUCSCHubFromBedFiles/GTSP2BED_fromFolder.R -t uniq ${i%%//};
    Rscript ~/makeUCSCHubFromBedFiles/GTSP2BED_fromFolder.R -t multi ${i%%//};
done

cat GTSP*.all.bed > caller.all.bed
cat GTSP*.uniq.bed > caller.uniq.bed
cat GTSP*.multi.bed > caller.multi.bed

rm GTSP*.bed

cat <<EOT > bed.csv
Sample,bedfile,Notes,freeze,hub
truth,truth.bed,truth,hg18,simint
callerall,callerall.bed,all,hg18,simint
calleruniq,calleruniq.bed,uniq,hg18,simint
callermulti,callermulti.bed,multi,hg18,simint
EOT

## may need to move to microb244 to make UCSC tracks
## bed.csv truth.bed caller.bed

Rscript ~/intSiteCaller/check_stats.R | cut -f1-30 | head -1 > tmp.txt
Rscript ~/intSiteCaller/check_stats.R | cut -f1-30 | awk '{for(i=1;i<=NF;i++)a[i]+=$i} END{for(i=1;i<=NF;i++)printf "%d%s", a[i], (i==NF?"\n":",")}' >> tmp.txt

cat tmp.txt | tr "\t" "," | tr "," "\t" | transpose.pl

