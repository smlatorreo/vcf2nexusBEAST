from sys import argv
import gzip

if len(argv) < 2:
    print('\nGiven a VCF file for diploid genotypes, this program converts it into a BEAST-suitable nexus format to the STDOUT\n')
    print('USAGE: python3 vcf2nexusBEAST.py [OPTIONS] input.vcf\n')
    print('OPTIONS:')
    print('--hets : Set this flag if your input file contains diploid heterozygous calls')
    print('--bin : REF = 0 ; ALT = 1\n')
    print('sergio.latorre _at_ tuebingen.mpg.de')
    print('https://github.com/smlatorreo')
    print('Sergio Latorre. Max Planck Institute for Developmental Biology.\n')
    exit()

def op(vcffile):
    if vcffile.endswith('gz'):
        vcf = gzip.open(vcffile, 'rt')
    else:
        vcf = open(vcffile, 'r')
    return vcf

vcffile = argv[-1]

linesheader = 0
with op(vcffile) as f:
    for line in f.readlines():
        if line.startswith('#CHROM'):
            linesheader += 1
            samples = line.strip().split('\t')[9:]
        elif line.startswith('#'):
            linesheader += 1
        else:
            break

# BINARY
def binary(vcffile, linesheader, samples):
    output = {}
    for s in samples:
        output[s] = []
    with op(vcffile) as f:
        for _ in range(linesheader):
            next(f)
        for line in f.readlines():
            for s in range(len(output)):
                A1 = line.split('\t')[s + 9][0]
                A2 = line.split('\t')[s + 9][2]
                output[samples[s]].append((A1, A2))
    return output

# SNPS
def snps(vcffile, linesheader, samples):
    output = {}
    for s in samples:
        output[s] = []
    with op(vcffile) as f:
        for _ in range(linesheader):
            next(f)
        for line in f.readlines():
            REF = line.split('\t')[3]
            ALT = line.split('\t')[4] 
            for s in range(len(output)): 
                if line.split('\t')[s + 9][0] == line.split('\t')[s + 9][2]: 
                    if line.split('\t')[s + 9][0] == '0': 
                        A1 = REF 
                        A2 = REF 
                    elif line.split('\t')[s + 9][0] == '1': 
                        A1 = ALT 
                        A2 = ALT 
                    else: 
                        A1 = '.' 
                        A2 = '.' 
                    output[samples[s]].append((A1, A2)) 
                else: 
                    if line.split('\t')[s + 9][0] == '0': 
                        A1 = REF 
                    elif line.split('\t')[s + 9][0] == '1': 
                        A1 = ALT 
                    else: 
                        A1 = '.' 
                    if line.split('\t')[s + 9][2] == '0': 
                        A2 = REF 
                    elif line.split('\t')[s + 9][2] == '1': 
                        A2 = ALT 
                    else: 
                        A2 = '.' 
                    output[samples[s]].append((A1, A2))
    return output

# Main Control
if '--hets' in argv:
    hets = True
else:
    hets = False

if '--bin' in argv:
    bina = True
    output = binary(vcffile, linesheader, samples)
else:
    bina = False
    output = snps(vcffile, linesheader, samples)

# Print
print('#NEXUS')
print('BEGIN DATA;')
if hets == False:
    print('\tDIMENSIONS NTAX={} NCHAR={};'.format((len(samples)), len(output[samples[0]])))
else:
    print('\tDIMENSIONS NTAX={} NCHAR={};'.format((len(samples) * 2), len(output[samples[0]])))
if bina == False:
    print('\tFORMAT DATATYPE=DNA GAP=.;')
else:
    print('\tFORMAT DATATYPE=BINARY SYMBOLS=\"01\" GAP=.;')
print('\tMATRIX')

if hets == False:
    for sample in output:
        out = ''
        for pos in range(len(output[sample])):
            out = out + (output[sample][pos][0])
        print('{0} {1}'.format(sample, out))
else:
    for sample in output:
        out1 = ''
        out2 = ''
        for pos in range(len(output[sample])):
            out1 = out1 + (output[sample][pos][0])
            out2 = out2 + (output[sample][pos][1])
        print('{0}_1 {1}'.format(sample, out1))
        print('{0}_2 {1}'.format(sample, out2))

print(';')
print('END;')
