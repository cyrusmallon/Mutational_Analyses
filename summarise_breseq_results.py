import sys
import glob

# sample             x
# timepoint          x
# evidence           x 
# seq_id             x
# position           x
# mutation           x
# frequency          x
# coverage           x    only for RA, in output.gd, e.g. total_cov=66/101
# score              x    only for (all) RA, in output.gd, e.g. polymorphism_score=25.1 
# annotation         -    for all, in annotated.gd  
# mutation_type      x    first col of output.gd
# Mutation_category  x    looks like for all, in annotated_1.gd
# gene               x    looks like for all, in annotated_1.gd, e.g. gene_name=rocD-C-Prokka_02755
# description        x    looks like gene_product Ribosome binding ATPase YchF/Vitamin B12 transporter BtuB

# http://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/gd_format.html

# breseq_outputs = glob.glob('output_*/*_*')

breseq_outputs = glob.glob('output_*')


print("hello!")


for b in breseq_outputs:

    sname = b.split('/')[-1]
    print(sname + "\n")
    mut_d = {}
    
    ogd = b + '/output/output.gd'
    
    with open(ogd) as fin:
        outputgd = fin.readlines()
    
    for line in outputgd:
        
        if line.startswith('#'):
            continue
        
        line = line.rstrip()
        
        llist = line.split('\t')
        parent_id = llist[2]
        
        if parent_id != '.':
            
            # entry type
            e_type = llist[0]
            e_nr = llist[1]
            p_id = parent_id
            
            # if the mutation has more than one index
            p_id_list = p_id.split(',')
            
            # for each index
            for id in p_id_list:
                
                seq_id = llist[3]
                pos = llist[4]
                new_seq = llist[5]
                insert_pos = ''       
                freq = 'NA'
              
                for c in llist:
                    if c.startswith('frequency='):
                        freq = float( c.split('frequency=')[1])
                        # print(freq, "\n")
                    if c.startswith('insert_position='):
                        insert_pos = c.split('insert_position=')[1]
                
                if e_type == 'SUB':
                    new_seq = llist[6]
                    
                if id not in mut_d:
                    mut_d[id] = {}
                    
                mut_d[id]['seq_id'] = seq_id
                mut_d[id]['position'] = pos
                
                if insert_pos != '':
                    mut_d[id]['position'] = pos + ':' + insert_pos
                    
                mut_d[id]['mutation'] = new_seq
                mut_d[id]['frequency'] = freq        
                mut_d[id]['coverage'] = 'NA'
                mut_d[id]['score'] = 'NA'
                mut_d[id]['annotation'] = 'NA'
                mut_d[id]['mutation_type'] = e_type
                mut_d[id]['mutation_category'] = 'NA'
                mut_d[id]['gene_name'] = 'NA'
                mut_d[id]['desciption'] = 'NA'
        else:
            evid = llist[0]
            evnt_idx = llist[1]
            
            if evnt_idx in mut_d:
                
                mut_d[evnt_idx]['evidence'] = evid
                cov = 'NA'
                lastcol = llist[-1]
                if lastcol.startswith('total_cov'):
                    t = lastcol.split('total_cov=')[-1]
                    covlist = t.split('/')
                    cov = str(int(covlist[0]) + int(covlist[1]))
                    
                poly_score = 'NA'
                cons_score = 'NA'
                for c in llist:
                    if c.startswith('consensus_score='):
                        cons_score = c.split('consensus_score=')[1]
                    if c.startswith('polymorphism_score='):
                        poly_score = c.split('polymorphism_score=')[1]
                    if c.startswith('pos_hash_score'):
                        cons_score = c.split('pos_hash_score=')[1]
                    if c.startswith('max_pos_hash_score'):
                        poly_score = c.split('max_pos_hash_score=')[1]  
                    if c.startswith('coverage_minus'):
                        cov = int(c.split('coverage_minus=')[1])
                    if c.startswith('coverage_plus'):
                        cov += int(c.split('coverage_plus=')[1])
                score = cons_score + '/' + poly_score
                mut_d[evnt_idx]['score'] = score
                mut_d[evnt_idx]['coverage'] = str(cov)
    
    agd = b + '/output/evidence/annotated.gd'
    with open(agd) as fin_annot:
        annotatedgd = fin_annot.readlines()
     
    for line in annotatedgd:
    
        if line.startswith('#'):
            continue
        
        line = line.rstrip()
    
        llist_annot = line.split('\t')
        col3 = llist_annot[2]
        if col3 != '.':
            e_nr_annot = llist_annot[1]
            evt_idx_annt = llist_annot[2]
            print(e_nr_annot, evt_idx_annt)
            
            evt_idx_annt_list = evt_idx_annt.split(',')
            
            for i in evt_idx_annt_list:
                evt_idx_annt = i
                intergen = True
                for col in llist_annot:
                    # print(col)
                    if col.startswith('mutation_category='):
                        mut_category = col.split('mutation_category=')[1]
                        mut_d[evt_idx_annt]['mutation_category'] = mut_category
                    if col.startswith('gene_name='):
                        # print("line 157")
                        gene_name = col.split('gene_name=')[1]
                        mut_d[evt_idx_annt]['gene_name'] = gene_name
                    if col.startswith('gene_product='):
                        description = col.split('gene_product=')[1]
                        mut_d[evt_idx_annt]['desciption'] = description
                    if col.startswith('aa_ref_seq='):
                        aa_ref_seq = col.split('aa_ref_seq=')[1]
                        intergen = False
                    if col.startswith('aa_new_seq='):
                        aa_new_seq = col.split('aa_new_seq=')[1]
                        intergen = False
                    if col.startswith('aa_position='):
                        aa_position = col.split('aa_position=')[1]
                        intergen = False
                    if col.startswith('codon_new_seq='):
                        codon_new_seq = col.split('codon_new_seq=')[1]
                        intergen = False
                    if col.startswith('codon_ref_seq='):
                        codon_ref_seq = col.split('codon_ref_seq=')[1]
                        intergen = False
                    if col.startswith('gene_position='):
                        gene_position = col.split('gene_position=')[1]
                        
                
                if intergen == True:
                    print(gene_position)
                    mut_d[evt_idx_annt]['annotation'] = gene_position
                else:
                    mut_d[evt_idx_annt]['annotation'] = ''.join([aa_ref_seq,
                                                                   aa_position,
                                                                   aa_new_seq,
                                                                   ' (',
                                                                   codon_ref_seq,
                                                                   '->',
                                                                   codon_new_seq,
                                                                   ')'])
                
    
    with open('summary/' + sname + '_breseq_table_v3.txt', 'w') as fout:
        header = ['evidence',
                  'seq_id',
                  'position',
                  'mutation',
                  'frequency',
                  'coverage',
                  'score',
                  'annotation',
                  'mutation_type',
                  'mutation_category',
                  'gene',
                  'description'
                  ]
        fout.write('\t'.join(header) + '\n')
        
        for m in mut_d:
            
            fr = mut_d[m]['frequency']
            if fr != 'NA':
                fr = str(float(mut_d[m]['frequency'])*100)
            
            l = [mut_d[m]['evidence'],
                 mut_d[m]['seq_id'],
                 mut_d[m]['position'],
                 mut_d[m]['mutation'],
                 fr,
                 mut_d[m]['coverage'],
                 mut_d[m]['score'],
                 mut_d[m]['annotation'],
                 mut_d[m]['mutation_type'],
                 mut_d[m]['mutation_category'],
                 mut_d[m]['gene_name'],
                 mut_d[m]['desciption']
                 ]
            fout.write('\t'.join(l) + '\n')
