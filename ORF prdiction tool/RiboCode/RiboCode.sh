for i in 0 4 8
do
for j in a b c
do
RiboCode -a RiboCode_annot -c /Data_3/Suna/RiboCode/output/metaplot/D${i}${j}_pre_config.txt -l no -g -o /Data_3/Suna/RiboCode/output/D${i}${j}/RiboCode_ORFs_result
done
done
