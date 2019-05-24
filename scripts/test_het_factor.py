import sys
import numpy as np
from collections import Counter


GT_SNPS_FN = '../validation_results/output.vcf'


def conv_call_str(call_str):
    vals = sorted(call_str.replace('|','').replace('/',''))
    if vals == ['0','0']:
        return 'hom_ref'
    elif vals == ['0','1']:
        return 'het'
    return 'hom_alt'

gt_calls = {}
with open(GT_SNPS_FN) as fp:
    for line in fp:
        if line.startswith('#'): continue
        chrm, pos, _, ref, alt, _, _, _, _, gt, mega_call = line.split()
        if gt == '.' or mega_call == '.': continue
        # uncomment for SNPs
        #if len(ref) > 1 or len(alt) > 1: continue
        # uncomment for deletions
        if (len(ref) == 1 and len(alt) == 1) or len(ref) < len(alt): continue
        # uncomment for insertions
        #if (len(ref) == 1 and len(alt) == 1) or len(ref) > len(alt): continue
        # uncomment for indels
        #if (len(ref) == 1 and len(alt) == 1): continue
        gt_calls[(chrm, pos)] = conv_call_str(gt)

llrs_fp = open('llrs.txt', 'w')
calls, calls_full = [], []
with open(sys.argv[1]) as fp:
    for line in fp:
        if line.startswith('#'): continue
        chrm, pos, _, ref, alt, _, _, _, m_fmt, mega_call = line.split()
        try:
            gt = gt_calls[(chrm, pos)]
        except KeyError:
            continue
        mega_call_txt = conv_call_str(mega_call.split(':')[0])
        calls.append((gt, mega_call_txt))
        calls_full.append((gt, mega_call_txt, ref, alt))

        llrs_idx = next(i for i, fmt in enumerate(m_fmt.split(':')) if fmt == 'LLRS')
        pos_llrs = list(map(float, mega_call.split(':')[llrs_idx].split(',')))
        pos_mean_llhr = np.mean(pos_llrs)
        for llhr in pos_llrs:
            llrs_fp.write('\t'.join(map(
                str,
                (llhr, chrm + '_' + pos, pos_mean_llhr, gt, mega_call_txt))) + '\n')
llrs_fp.close()


count_calls = Counter(calls)
match, tot = 0, 0
print('Truth\tMegaCall\tCount')
for (gt, mega_call), count in count_calls.items():
    print('{}\t{}\t{}'.format(gt, mega_call, count))
    tot += count
    if gt == mega_call: match += count
print('{}\t{}\t{}'.format(match / tot, match, tot))


with open('calls_full.txt', 'w') as fp:
    for gt, mega_call, ref, alt in sorted(calls_full):
        fp.write('{}\t{}\t"{}"\t"{}"\n'.format(gt, mega_call, ref, alt))
