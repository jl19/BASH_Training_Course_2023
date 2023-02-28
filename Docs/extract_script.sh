gzcat Mus_musculus.GRCm38.cdna.all.fa.gz| \
grep '>' |
awk '{FS= " "}BEGIN{ print "TXNAME,GENEID"};{print substr($1,2) "," substr($4,6)};' > tx2gene.mm.GRCm38.cdna.csv 

