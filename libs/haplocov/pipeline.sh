perl addToTableNCBI.pl --metadata metadata.tsv --ref ref.fna --seq sequences.fasta --outfile out.HaploCoV --nproc 16
perl computeDefining.pl out.HaploCoV
perl assign.pl --dfile out.HaploCoV.definingVariants.txt --infile out.HaploCoV --outfile out.HaploCoV.assigned
perl augmentClusters.pl --metafile out.HaploCoV.assigned --deffile out.HaploCoV.definingVariants.txt --posFile out.HaploCoV.frequentVariants.txt --outfile out.HaploCoV.definingVariantsNew.txt --dist 5 --size 50
perl assign.pl --dfile out.HaploCoV.definingVariantsNew.txt --infile out.HaploCoV --outfile out.HaploCoV.assignedNew
cut -f 10 out.HaploCoV.assignedNew | sort | uniq -c