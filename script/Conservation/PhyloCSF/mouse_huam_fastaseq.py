m_write = []
h_write = []
for type in ['uORF', 'lncRNA', 'dORF']:
    mouse = open(f'/Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/mouse/{type}.fasta', 'r')
    human = open(f'/Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/human/{type}.fasta', 'r')
    for i, val in enumerate(zip(mouse, human)):
        m = val[0].rstrip('\n').lstrip('>')
        h = val[1].rstrip('\n').lstrip('>')
        if i % 2 != 0:
            m_write.append(m)
            h_write.append(h)
            out = open(f'/Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/finalinput/{type}/{name}.fasta', 'w')
            out.write('\n'.join(m_write) + '\n')
            out.write('\n'.join(h_write) + '\n')
            out.close()
            m_write = []
            h_write = []
        else:
            m_write.append('>Mouse') # +m
            h_write.append('>Human') # +h
            name = m.lstrip('>').split(':')[0] # +f'_{i}'

    

mouse.close()
human.close()
