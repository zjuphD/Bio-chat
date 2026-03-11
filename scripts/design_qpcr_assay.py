import argparse
import requests
import json
import time
import math

def get_transcript(accession):
    if not accession: return None
    print(f"Fetching sequence for {accession}...")
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    for _ in range(3):
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                lines = response.text.strip().split('\n')
                return "".join(lines[1:])
        except:
            time.sleep(1)
    return None

def get_exon_junctions(accession):
    print(f"Fetching exon annotation for {accession}...")
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=gb&retmode=text"
    junctions = []
    for _ in range(3):
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                lines = response.text.split('\n')
                current_end = 0
                for line in lines:
                    if line.strip().startswith('exon'):
                        parts = line.split()
                        if len(parts) > 1:
                            coords = parts[1].replace('<', '').replace('>', '')
                            if '..' in coords:
                                try:
                                    start, end = coords.split('..')
                                    if current_end > 0:
                                        junctions.append(int(current_end))
                                    current_end = int(end)
                                except: pass
                return junctions
        except:
            time.sleep(1)
    return junctions

def map_homolog_junctions(target_seq, homolog_acc):
    print(f"Mapping homolog junctions from {homolog_acc} to target...")
    homolog_seq = get_transcript(homolog_acc)
    if not homolog_seq: return []
    
    h_juncs = get_exon_junctions(homolog_acc)
    target_juncs = []
    
    for j in h_juncs:
        if 20 < j < len(homolog_seq) - 20:
            anchor = homolog_seq[j-20 : j+20]
            idx = target_seq.find(anchor)
            if idx != -1:
                target_juncs.append(idx + 20)
            else:
                # Sliding window alignment (Simplified Smith-Waterman logic)
                for offset in range(-5, 6):
                    sub_anchor = homolog_seq[j-15+offset : j+15+offset]
                    idx_sub = target_seq.find(sub_anchor)
                    if idx_sub != -1:
                        target_juncs.append(idx_sub + 15 - offset)
                        break
    
    return sorted(list(set(target_juncs)))

def calculate_tm_nn(seq):
    if not seq: return 0
    # Values from SantaLucia (1998)
    nn_data = {
        'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2), 'AT': (-7.2, -20.4), 'TA': (-7.2, -21.3),
        'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7), 'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
        'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0), 'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
        'CG': (-10.6, -27.2), 'GC': (-9.8, -24.4), 'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9)
    }
    h, s = 0.2, -5.7 # Initial
    seq = seq.upper()
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        if dinuc in nn_data:
            dh, ds = nn_data[dinuc]
            h += dh
            s += ds
    if seq[0] in 'AT': h += 2.2; s += 6.9
    if seq[-1] in 'AT': h += 2.2; s += 6.9
    # Na+ = 50mM, Primer = 0.2uM
    s += 0.368 * (len(seq) - 1) * math.log(0.05)
    tm = (1000 * h) / (s + 1.987 * math.log(0.0000002 / 4.0)) - 273.15
    return round(float(tm), 1)

def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join([complement.get(b, 'N') for b in seq[::-1]])

def check_self_dimer(seq):
    rc = rev_comp(seq)
    for length in range(4, 8):
        for i in range(len(seq) - length + 1):
            if seq[i:i+length] in rc: return True
    return False

def check_3prime_quality(seq):
    last_5 = seq[-5:]
    return (last_5.count('G') + last_5.count('C')) <= 4

def spans_exon_junction(start, end, junctions):
    for j in junctions:
        if start + 4 <= j <= end - 4: return True
    return False

def design_assay(target_seq, offtarget_seq=None, junctions=None, amplicon_range=(70, 150)):
    candidates = []
    seq_len = len(target_seq)
    if not junctions: junctions = []
    
    for f in range(20, seq_len - amplicon_range[1] - 30):
        for length in range(17, 24):
            f_seq = target_seq[f : f+length]
            f_tm = calculate_tm_nn(f_seq)
            if 58 <= f_tm <= 65:
                if check_3prime_quality(f_seq) and not check_self_dimer(f_seq):
                    for r_len in range(17, 24):
                        for r in range(f + amplicon_range[0], f + amplicon_range[1]):
                            if r + r_len > seq_len: break
                            r_seq = rev_comp(target_seq[r : r+r_len])
                            r_tm = calculate_tm_nn(r_seq)
                            if 58 <= r_tm <= 65 and abs(f_tm - r_tm) <= 2.5:
                                if check_3prime_quality(r_seq) and not check_self_dimer(r_seq):
                                    is_gdna_safe = True if not junctions else (spans_exon_junction(f, f+length, junctions) or spans_exon_junction(r, r+r_len, junctions))
                                    if is_gdna_safe:
                                        candidates.append({
                                            "Amplicon": f"{r + r_len - f}bp",
                                            "F_Primer": {"seq": f_seq, "tm": f_tm},
                                            "R_Primer": {"seq": r_seq, "tm": r_tm}
                                        })
                                        if len(candidates) >= 5: return candidates
    return candidates

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--target")
    parser.add_argument("--homolog")
    args = parser.parse_args()
    target_seq = get_transcript(args.target)
    junctions = []
    if target_seq:
        if args.homolog: junctions = map_homolog_junctions(target_seq, args.homolog)
        else: junctions = get_exon_junctions(args.target)
        results = design_assay(target_seq, junctions=junctions)
        print(json.dumps(results, indent=2))
