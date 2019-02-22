#!/usr/bin/env python3
# coding: utf-8

from Bio import SeqIO
from Bio.SeqUtils import GC,GC_skew
from subprocess import PIPE,Popen
import os
import re
import colour
from collections import OrderedDict
import svgwrite
import math
import argparse

class SumGene:
    def __init__(self, name=None, tlen=0, featurelen=0, tseq=''):
        self.name = name
        self.len = tlen
        self.featurelen = featurelen
        self.seq = tseq

class Gene:
    def __init__(self, name=None, locus_tag=None, seq=None, genetype=None, position=None, strand=None,
                 genlen=None, group=None):
        self.name = name
        self.locus_tag = locus_tag
        self.seq = seq
        self.type = genetype
        self.position = position
        self.strand = strand
        self.genlen = genlen
        self.group = group
        
    def __str__(self):
        return "{0} {1} {2} {3} {4} {5} {6}".format(
            self.type,self.name,self.locus_tag,self.position[0],self.position[1],self.strand,self.genlen)
       
def parse_gb2class(gbfile,genedict,sumgene):
    with open("tmp.fasta","w+") as f:
        contigcount=0
        for record in SeqIO.parse(gbfile, "genbank"):
            sumgene.name = record.description
            sumgene.len += len(record)
            sumgene.featurelen += len(record.features)
            sumgene.seq += record.seq
            #extract features nucleotide data to fasta file
            for feature in record.features:
                if feature.type not in ("source", "gene", ""):
                    gene = Gene()
                    gene.name = feature.qualifiers.get('gene','N')[0]
                    gene.locus_tag = feature.qualifiers.get('locus_tag','N')[0]
                    gene.type = feature.type
                    location_str = re.sub('[<>]','',str(feature.location))
                    position_str = re.split(r'[\[\]\:\{\}]', location_str)
                    if position_str[0] == 'join':
                        jpos = [int(i) for i in position_str[2:4]+position_str[5:7]]
                        if math.fabs(jpos[1]-jpos[2]) < 50:
                            gene.position = [min(jpos)+contigcount,max(jpos)+contigcount]
                        else:
                            gene.position = [int(i)+contigcount for i in position_str[2:4]]
                    else:
                        gene.position = [int(i)+contigcount for i in position_str[1:3]]
                    gene.strand = feature.location.strand
                    gene.genlen = gene.position[1] - gene.position[0] + 1
                    gene.seq = sumgene.seq[gene.position[0]:gene.position[1]]
                    print(">"+gene.locus_tag, file=f)
                    print(gene.seq, file=f)
                    genedict[gene.locus_tag] = gene
            contigcount += len(record)

def run_blastx():
    db_name = "/bip5_disk/jiamin_105/COG/egg/NOGDB/COG"
    out_str = "6 qseqid sallseqid qcovs pident length slen mismatch gapopen qstart qend sstart send evalue bitscore"
    gene_fna = "tmp.fasta"
    header={}

    with Popen(['blastx','-query',gene_fna,'-db',db_name,'-num_threads','30','-outfmt',out_str,'-max_target_seqs','1'],
               stdout=PIPE,bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            headerlist = line.split('\t')
            header[headerlist[0]] = headerlist[1].split('_')[0]
    # with open('tmp_blastx.tsv','r') as bxf:
        # for line in bxf:
            # headerlist = line.split('\t')
            # header[headerlist[0]] = headerlist[1].split('_')[0]       
    p.wait()
    os.remove("tmp.fasta")
    return header    
    
def run_diamond():
    db_name = "/bip5_disk/fangren106/Biopy/diamond_COGdb/COG"
    gene_fna = "tmp.fasta"
    header={}
    with Popen(['diamond','blastx','-f','6','-d',db_name,'-q',gene_fna,'-k','1'],
               stdout=PIPE,bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            headerlist = line.split('\t')
            header[headerlist[0]] = headerlist[1].split('_')[0]
    p.wait()
    os.remove("tmp.fasta")
    return header

def run_debug():
    header={}
    with open('test_CCU063.tsv','r') as bxf:
        for line in bxf:
            headerlist = line.split('\t')
            header[headerlist[0]] = headerlist[1].split('_')[0]  
    return header

def extract_group(DB_annotate_file):
    dbtag2group_dict={}
    with open(DB_annotate_file, 'r') as dbf:
        for dbline in dbf:
            tmp = dbline.split('\t')
            dbtag2group_dict[tmp[1]] = tmp[-1].rstrip('\n')
    return dbtag2group_dict

def generate_colormap(genedict,gene_group_color):
    g2c = {}  
    for tag,g in genedict.items():
        if g.group not in g2c and g.group != None:
            g2c[g.group] = 0

    # pick color
    glen = len(g2c)
    if glen == 0:
        glen = 1
    color1 = colour.Color(gene_group_color[0].lower())
    color2 = colour.Color(gene_group_color[1].lower())
    c_list = list(color1.range_to(color2,glen))
    # print(c_list)
    cnum = 0
    sort_g2c = OrderedDict(sorted(g2c.items(),key=lambda x:x[0]))
    for g,c in sort_g2c.items():
        sort_g2c[g] = c_list[cnum]
        cnum+=1
    return sort_g2c

def position_mapping(sumgene,genepos,radius):
    # compute Polar coordinates and transform to Right angle coordinates
    s_theda = genepos[0]/sumgene.len*2*math.pi
    e_theda = genepos[1]/sumgene.len*2*math.pi
    start_x = round(radius*math.sin(s_theda),5)
    start_y = round(radius*math.cos(s_theda),5)
    end_x = round(radius*math.sin(e_theda),5)
    end_y = round(radius*math.cos(e_theda),5)
    map_pos=(start_x,start_y,end_x,end_y)
    return map_pos

def draw_number_mark(sumgene,outer_r,dwg):
    mark_count = 0
    while mark_count <= sumgene.len:
        mark_pos = position_mapping(sumgene,[mark_count,mark_count+5000],outer_r+5)
        dwg.add(dwg.line((1500+mark_pos[0],1500-mark_pos[1]),(1500+mark_pos[2],1500-mark_pos[3]),
                         stroke='black',stroke_width=10))
        m_legend_pos = position_mapping(sumgene,[mark_count+2500,0],outer_r+20)
        # adjust number mark text position
        if sumgene.len*9/16 < mark_count < sumgene.len*15/16:
            dwg.add(dwg.text(str(mark_count),insert=(1500+m_legend_pos[0],1500-m_legend_pos[1]),
                         stroke_width=10,style="text-anchor: end;"))
        elif sumgene.len*1/16 < mark_count < sumgene.len*7/16:
            dwg.add(dwg.text(str(mark_count),insert=(1500+m_legend_pos[0],1500-m_legend_pos[1]),
                         stroke_width=10,style="text-anchor: start;"))
        else:
            dwg.add(dwg.text(str(mark_count),insert=(1500+m_legend_pos[0],1500-m_legend_pos[1]),
                         stroke_width=10,style="text-anchor: middle;"))
        mark_count += 500000
    return dwg

def draw_cds_rna(sumgene,genedict,sort_g2c,dwg):
    for tag,gene in genedict.items():
        if gene.group != None:
            if gene.type == 'CDS':
                if gene.strand == 1: # positive strand
                    map_p = position_mapping(sumgene, gene.position, 760)
                    dwg.add(dwg.line((1500+map_p[0],1500-map_p[1]),(1500+map_p[2],1500-map_p[3]),
                                     stroke=str(sort_g2c[gene.group]).lower(),stroke_width=50))
                if gene.strand == -1: # negative strand
                    map_n = position_mapping(sumgene, gene.position, 690)
                    dwg.add(dwg.line((1500+map_n[0],1500-map_n[1]),(1500+map_n[2],1500-map_n[3]),
                                     stroke=str(sort_g2c[gene.group]).lower(),stroke_width=50))
        if gene.type == 'tRNA':
            map_t = position_mapping(sumgene, gene.position, 620)
            dwg.add(dwg.line((1500+map_t[0],1500-map_t[1]),(1500+map_t[2],1500-map_t[3]),
                             stroke='blue',stroke_width=50))
        if gene.type == 'rRNA':
            map_r = position_mapping(sumgene, gene.position, 550)
            dwg.add(dwg.line((1500+map_r[0],1500-map_r[1]),(1500+map_r[2],1500-map_r[3]),
                             stroke='blue',stroke_width=50))
    return dwg

def draw_gc(sumgene,dwg,gc_color):
    #gc kmer number
    unit_num = 300
    unit_len = int(sumgene.len/unit_num)
    #first coordinates
    start_unit = 0
    end_unit = unit_len-1
    mid_unit = (end_unit-start_unit)/2
    #radius
    gcc_r = 430
    gcs_r = 230
    gcc_range = (330,530)
    gcs_range = (130,330)
    gcc_mean = GC(sumgene.seq)
    gc_contents = []
    gc_skews = []
    #compute gc
    while(unit_num > 0):
        
        gc_content = GC(sumgene.seq[start_unit : end_unit])
        gcc_variance = gc_content-gcc_mean
        gc_contents.append(gcc_variance)
        gc_skew = GC_skew(sumgene.seq[start_unit : end_unit],window=unit_len)[0]
        gc_skews.append(gc_skew)
        start_unit += unit_len
        if end_unit+unit_len <= sumgene.len:
            end_unit += unit_len
        else:
            end_unit = sumgene.len
        unit_num -= 1 
    gcc_max = max(gc_contents)
    gcc_min = min(gc_contents)
    gcs_max = max(gc_skews)
    gcs_min = min(gc_skews)
    gcc_scal = (gcc_range[1]-gcc_range[0])/(gcc_max-gcc_min)
    gcs_scal = (gcs_range[1]-gcs_range[0])/(gcs_max-gcs_min)
    gc_count = 0
    
    unit_num = 300
    start_unit = 0
    end_unit = unit_len-1
    mid_unit = (end_unit-start_unit)/2
    while(unit_num > 0):
        # draw gc_content 
        c_pos = position_mapping(sumgene, [start_unit,end_unit], gcc_r)
        gcc_nr = (gc_contents[gc_count]-gcc_min)*gcc_scal+gcc_range[0]
        cm_pos = position_mapping(sumgene, [mid_unit,0], gcc_nr)
        points = [(1500+c_pos[0],1500-c_pos[1]),(1500+c_pos[2],1500-c_pos[3]),(1500+cm_pos[0],1500-cm_pos[1])]
        if gcc_nr >= gcc_r: 
            dwg.add(dwg.polygon(points, fill=gc_color[0], stroke_width=0))   
        else:
            dwg.add(dwg.polygon(points, fill=gc_color[1], stroke_width=0))

        # draw gc_skew
        s_pos = position_mapping(sumgene, [start_unit,end_unit], gcs_r)
        gcs_nr = (gc_skews[gc_count]-gcs_min)*gcs_scal+gcs_range[0]
        sm_pos = position_mapping(sumgene, [mid_unit,0], gcs_nr)
        points=[(1500+s_pos[0],1500-s_pos[1]),(1500+s_pos[2],1500-s_pos[3]),(1500+sm_pos[0],1500-sm_pos[1])]
        if gcs_nr >= gcs_r:
            dwg.add(dwg.polygon(points,fill=gc_color[2], stroke_width=0))
        else:
            dwg.add(dwg.polygon(points,fill=gc_color[3], stroke_width=0))

        start_unit += unit_len

        if end_unit+unit_len <= sumgene.len:
            end_unit += unit_len
        else:
            end_unit = sumgene.len

        mid_unit += unit_len
        unit_num -= 1 
        gc_count += 1
    return dwg
    
def draw_svg(sumgene,genedict,sort_g2c,dwg,gc_color):
    print("Draw svg")
    # outer circle radius
    outer_r = 800
    # draw outer circle
    dwg.add(dwg.circle(center=(1500,1500),r=outer_r,fill='none',stroke='black',stroke_width=5))
    # draw legend 
    addr = 0
    for group,color in sort_g2c.items():
        dwg.add(dwg.line((2400,1800+addr),(2420,1800+addr),stroke=str(color).lower(),stroke_width=20))
        dwg.add(dwg.text(group,insert=(2430,1805+addr),stroke_width=20))
        addr+=25
    dwg = draw_number_mark(sumgene,outer_r,dwg)
    dwg = draw_cds_rna(sumgene,genedict,sort_g2c,dwg)
    dwg = draw_gc(sumgene,dwg,gc_color)
    return dwg    


def main():
    # get argv variables
    gbfile = args.gbfile
    DB_annotate_file = args.feature_file 
    svg_file_name = args.output
    gene_group_color = (args.gg_c[0], args.gg_c[1])
    gc_color = (args.gccp_c, args.gccn_c, args.gcsp_c, args.gcsn_c) 
    
    print("\n")
    print("Data Extract")
    print('open genbank: "' + gbfile + '"')
    genedict={}
    sumgene = SumGene()
    parse_gb2class(gbfile,genedict,sumgene)
    print("read genbank done")
    
    # aliment
    if args.blastx:
        header = run_blastx()
        print("blastx done")
    elif args.debug:
        header = run_debug()
        print("debug mode: open tsv done")
    else:
        header = run_diamond()
        print("diamond done")

    # extract group data from annotation file 
    dbtag2group_dict = extract_group(DB_annotate_file)

    # map (feature tag & db tag) Dict and (db tag & tag group) Dict to genedict
    # add gene group to genedict
    for key,value in header.items():
        if key in genedict:
            if value in dbtag2group_dict:
                genedict[key].group = dbtag2group_dict[value]
    # generate gene_group colormap
    sort_g2c = generate_colormap(genedict, gene_group_color)        
    print("Data Extract done\n")
    
    # draw svg
    dwg = svgwrite.Drawing(svg_file_name,size=(u'5000', u'5000'))
    dwg = draw_svg(sumgene,genedict,sort_g2c,dwg,gc_color)
    dwg.save()
    print("Draw svg done")
    print("Program \"COGplot.py\" complete")
    print("COGplot has been generated to \""+ svg_file_name + "\"\n")
    

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument("gbfile", help="genbank file name")
    p.add_argument('-o','--output', help="output svg file name", default='COGplot.svg')
    p.add_argument('-ff','--feature_file', help="db annotation file",
                   default='/bip5_disk/fangren106/Biopy/diamond_COGdb/bactNOG.annotations.tsv')
    p.add_argument('--blastx', help="use blasx run aliment(diamond)", action='store_true')
    p.add_argument('-g','--debug', help="open test tsv file instead alignment(for debug)", action='store_true')
    p.add_argument('--gg_c', help="gene group color range(give (start,end) color)", nargs=2, default=['#800000','blue'])
    p.add_argument('--gccp_c', help="gc content positive triangle color(limegreen)", default='limegreen')
    p.add_argument('--gccn_c', help="gc content negative triangle color(mediumpurple)", default='mediumpurple')
    p.add_argument('--gcsp_c', help="gc skew positive triangle color(deepskyblue)", default='deepskyblue')
    p.add_argument('--gcsn_c', help="gc skew negative triangle color(gold)", default='gold')
    args = p.parse_args()
    
    main()
