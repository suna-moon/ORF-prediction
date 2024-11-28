for type in uORF dORF lncRNA
do
    path=`find /Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/finalinput/${type} -type f`
    for filename in ${path}
        do
            name=${filename##*/}
            clustalo -i ${filename} -o /Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/MSA_clustalOmega/${type}/${name} --outfmt=fa --force
        done
done
