import sys #reads through command line 
import gzip #reads unzipped file 
#filename = sys.argv[1]
f = gzip.open('small_testfile.chr22.vcf.gz', 'rb')
g = gzip.open('chr22.1000Genomes.SNP_info.gz', 'rb')
h = gzip.open('chromosome22.txt.gz', 'w')
i = open('samples.txt', 'w') 
rsID = {}
temp = []


for line in g:
        linelist = line.strip().split() #splits and strips list
        (chr, pos, id, ref, alt, qual, filter, info, format)= linelist[0:9]#define variables
        rsID[pos] = id #finds corresponding id to position
        
for line in f:
        if line.startswith ("#CHROM"):# checks if line starts with "#CHROM"
                temp = line.split() #splits line 
                for word in temp:
                        if word.startswith('LD'):#checks if it starts with 'LD'
                                i.write(word + '\n') #writes list of samples to file
        if line.startswith ('#'): #reads next lines 
                f.next()

        else:
                linelist = line.strip().split() #splits and strips list
                (chr, pos, id, ref, alt, qual, filter, info, format)= linelist[0:9]   #define variables
                r2 = info[15:22] #collects r^2 values
                r2 = float(r2)# convert to a float so you can filter
                if r2 >= 0.8:#if r^2 value greater than 0.8
                        gt_dosagelist = linelist[9:]
                        dosagelist = map(lambda x : float(x.split(":")[1]), gt_dosagelist)#lambda function to split each entry in gt_dosagelist and collect the dosage
                        freqalt = sum(dosagelist)/(len(dosagelist)*2)#calc ALT allele freq (I found that the ALT freq does not always = MAF), len(dosagelist)*2 because 2 chromosomes
                        dosagestring = ' '.join(str(e) for e in dosagelist) #converst to string
                        h.write('chr22 ' + rsID[pos] + ' ' + pos + ' ' + ref + ' ' + alt + ' ' + str(freqalt) + dosagestring + '\n') #writes to file 
    
            
           
h.close()        


        
    
