import re
import csv
from tqdm import tqdm
import argparse


def parse_cigar(cigar):
    matches = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    length = 0
    for count, op in matches:
        if op in "M=XD":
            length += int(count)
    return length

def load_read_groups(rg_file):
    with open(rg_file, 'r') as f:
        rg_list = [line.strip() for line in f]
    return rg_list

def process_sam_file(sam_file, rg_list):
    paired_reads = {}
    coordinates = []
    wo_pair_count = 0

    chrs = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    
    with open(sam_file, 'r') as f:
        for line in tqdm(f):
            if line.startswith('@'):
                continue
            
            fields = line.split('\t')
            qname = fields[0]
            flag = int(fields[1]) 
            rname = fields[2] 
            pos = int(fields[3]) 
            cigar = fields[5]  
            tlen = int(fields[8]) 
            rg = fields[-1].strip()
            if rg not in rg_list:
                continue
            if rname not in chrs:
                continue
                
            strand = '+' if not flag & 0x10 else '-'

            length = parse_cigar(cigar)
            end_pos = pos + length - 1
            
            if qname not in paired_reads:
                paired_reads[qname] = {
                    'rname': rname,
                    'start1': pos,
                    'end1': end_pos,
                    'start2': None,
                    'end2': None,
                    'tlen': tlen,
                    'rg': rg,
                    'strand1': strand
                }
            else:
                paired_reads[qname]['start2'] = pos
                paired_reads[qname]['end2'] = end_pos
                paired_reads[qname]['strand2'] = strand

    headers = ['read_name', 'chromosome', 'start1', 'end1', 'start2', 'end2', 'tlen', 'rg', 'strand1', 'strand2']
    with open('reads.csv', mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=headers)
    
        writer.writeheader()
        
        for read_name, data in paired_reads.items():
            writer.writerow({
                'read_name': read_name,
                'chromosome': data['rname'],
                'start1': data['start1'],
                'end1': data['end1'],
                'start2': data.get('start2', 'None'),
                'end2': data.get('end2', 'None'),
                'tlen': data['tlen'],
                'rg': data['rg'],
                'strand1': data['strand1'],
                'strand2': data.get('strand2', 'None')
            })

    for qname, data in tqdm(paired_reads.items()):
        if data['start2'] is None:
            wo_pair_count += 1
        else:
            start = min(data['start1'], data['start2'])
            end = max(data['end1'], data['end2'])
            coordinates.append({
                'read_name': qname,
                'chromosome': data['rname'],
                'start': start,
                'end': end,
                'rg': data['rg'],
                'strand1': data['strand1'],
                'strand2': data['strand2']  
            })

    return coordinates, wo_pair_count

def save_coordinates_to_csv(coordinates, output_file):
    headers = ['read_name', 'chromosome', 'start', 'end', 'rg', 'strand1', 'strand2']
    
    with open(output_file, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=headers)
        writer.writeheader()
        for coord in tqdm(coordinates):
            writer.writerow(coord)

parser = argparse.ArgumentParser()
parser.add_argument("--name", type=str, required=True)
args = parser.parse_args()

GSM_id = args.name
name = f'{GSM_id}.sam'
path = 'data/'

GSM_id = name.split('.')[0]

rg_file = f'{path}/top_100_{GSM_id}.txt'
sam_file = f'{path}/{name}'

rg_list = load_read_groups(rg_file)

coordinates, wo_pair_count = process_sam_file(sam_file, rg_list)

#print(f"Количество одиночных прочтений (без пары): {wo_pair_count}")

save_coordinates_to_csv(coordinates, f'{path}/{GSM_id}_top_100_rg_coordinates_of_reads.csv')
