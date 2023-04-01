# Vic_probes
una colección de códigos para la identificación de regiones génicas únicas en especies a partir de secuencias genómicas

# Vic probes code 1 : Dedemos aegurarnos de que las especies que estamos trabajando presenten monofilia reciproca, esto sera inferido a traves de un analisis filogenómico
```r
setwd("/home/hp/Documentos/sondas_lbm/borreliella")
getwd()
dir()
#######################
### cambiar nombres 1 (no usar)###
#######################

library(seqinr) ;
r <- dir() ;
head <- gsub("*.fasta","",r) ;
a <- 0 ;
for (i in 1:length(head)){ ;
  a <- read.fasta(r[i]) ;
  names(a) <- paste0(rep(head[i],length(a)),"__",names(a)) ;
  write.fasta(a,names(a), file.out=paste0(head[i],".fasta")) ;
} ;
q("no") ;

#######################
### cambiar nombres 2 (usar)###
#######################

library(seqinr) ;
r <- dir() ;
a1 <- gsub(".1_.*",".1",r) ;
a2 <- gsub(".2_.*",".2",a1) ;
head <- gsub(".3_.*",".3",a2) ;
a <- 0 ;
for (i in 1:length(head)){ ;
  a <- read.fasta(r[i]) ;
  names(a) <- paste0(rep(head[i],length(a)),"_",names(a)) ;
  write.fasta(a,names(a), file.out=paste0(head[i],".fasta")) ;
} ;
q("no") ;

#### remove ###
rm *.fna

#######################
### cambiar nombres 3 (usar para ffn file)###
#######################

library(seqinr) ;
r <- dir() ;
a1 <- gsub(".1.ffn",".1",r) ;
a2 <- gsub(".2.ffn",".2",a1) ;
head <- gsub(".3.ffn",".3",a2) ;
a <- 0 ;
for (i in 1:length(head)){ ;
  a <- read.fasta(r[i]) ;
  names(a) <- paste0(rep(head[i],length(a)),"_",names(a)) ;
  write.fasta(a,names(a), file.out=paste0(head[i],".fasta")) ;
} ;
q("no") ;

#### remove ###
rm *.ffn ; 
cat *.fasta > input.fasta ;
ls -lh ; 

#anotar los contigs con PROKKA#
conda activate prokka_env
mkdir -p annotation
for r1 in *fasta
do
prefix=$(basename $r1 .fasta)
prokka --cpus 30 $r1 -o ${prefix} --prefix ${prefix} --kingdom Bacteria ; 
mv ${prefix}/*.gff annotation/${prefix}.gff
done
conda deactivate ;
cp */*.ffn ffn/

#extraer las regiones comunes: PANGENOME#
roary -p 4 -f roary -e -n -v -cd 95 -i 95 annotation/*gff ;
ls -lh ; 

## generar alineamiento de snps final (2 opciones) ##
snp-sites -m -o snp1.phy core_gene_alignment.aln ; 
snp-sites -m -c -o snp2.phy core_gene_alignment.aln ; 
ls -lh; 

## phylogeny ##
raxmlHPC-PTHREADS -p 1254512 -m GTRCAT -s ehrlichia_core.phy -n nwk -# 20 -T 30

## metadata ##
library(tidyr)

c <- read.csv("borreliella.tsv", header=TRUE, sep="\t")
head(c)
d <- data.frame(accession=c$Assembly.Accession, assembly=c$Assembly.Name, species=c$Organism.Name, 
                cepa=c$Organism.Infraspecific.Names.Strain, length=c$Assembly.Stats.Total.Sequence.Length,  
                level=c$Assembly.Level, date=c$Assembly.Submission.Date, genes=c$Annotation.Count.Gene.Total, prots=c$Annotation.Count.Gene.Protein.coding,method=c$Assembly.Sequencing.Tech)
head(d)

e <- separate(d,"species",c("a","b","c","d","e","f","g"), sep=" ")
head(e)
e <- data.frame(e[,1:2],species=paste0(e$a," ",e$b),e[,10:16])
head(e)
f <- as.data.frame(table(e$species))
f1 <- f[order(f$Var1),]
g <- e[e$species %in% f1$Var1[c(2:12,14:31,33:35)],]
dim(g)
dim(e)

write.table(e,"borreliella_metadata.tsv", row.names=FALSE, sep="\t")
write.table(e$accession,"borreliella_genomes.txt", row.names=FALSE, sep="\t",quot=F, col.names = F)
```

# Vic probes code 2 : Dedemos adjuntar la informacion de presencia u ausencia de genes a la metadata generada y plotearla en el árbol generado 
```r
##########################
### FINAL METADATA ##
##########################

m <- read.csv("borreliella_metadata.tsv", header=TRUE, sep="\t")
dim(m)
head(m)

library(ape)
tree <- read.tree("RAxML_bestTree_rooted.nwk") 
names(tree)
tl <- gsub("'","",tree$tip.label)

length(setdiff(tl,m$accession))
length(setdiff(m$accession,tl))

m2 <- m[m$accession %in% tl,]
dim(m2)

length(setdiff(tl,m2$accession))
length(setdiff(m2$accession,tl))

dat <- read.csv("gene_presence_absence.csv", header=TRUE)
dim(dat)
length(names(dat))
length(unique(names(dat)))
names(dat)[1:30]

head(dat[1:10,1:30])
dbt <- dat[,c(1,15:ncol(dat))]
dim(dbt)
x <- gsub("[[:alnum:]].*","1",as.matrix(dbt[,2:ncol(dbt)]))
y <- gsub("$.*","0",x)
z <- gsub("10","1",y)
w <- data.frame(Gene=dbt$Gene,z)
v <- as.data.frame(t(w))
p <- v[c(2:nrow(v)),]
names(p) <- v[1,]
class(p)
dim(p)
row(p)

q <- data.frame(sample=names(dbt[,2:ncol(dbt)]),p)
dim(q)
q[1:30,1:10]

dim(m)
dim(m2)
head(m2)
names(m2) <- c("sample",names(m2)[2:10])
q4 <- merge(m2,q,by="sample",all.x=FALSE)
dim(q4)
names(q4)[1:20]
q4[1:10,1:25]

##### write_metadata ######
dir()
write.csv(q4,"metadata_presabs.csv", row.names=FALSE)
```

# Vic probes code 3 : Ahora, necesitamos identificar los genes propios de una especie determinada 
```r
##################
##### assay ######
##################

genes <- names(q4)[11:length(names(q4))]
length(genes)
names(q4)[1:20]
unique(q4$species)

sp <- "Borreliella burgdorferi"

## loop para identificar genes unicos de B. burgdorferi ##
a <- 0
b <- 0
c <- 0
d <- 0
e <- 0
f <- 0 
for (i in 11:length(names(q4))) {
  a <- sum(as.numeric(q4[q4$species %in% sp ,as.numeric(i)]))
  b <- sum(as.numeric(q4[!q4$species %in% sp ,as.numeric(i)]))
  c <- append(c,a)
  d <- append(d,b)
  e <- append(e,i)
  f <- append(f,genes[i-10])
}

g <- data.frame(n_strains=c,otros=d,order=e, gene=f)
head(g)
dim(g)
as.data.frame(table(q4$species))
dim(q4)

h <- g[g$n_strains >  140 & g$otros < 1,]
hj <- data.frame(gene=dat[dat$Gene %in% h$gene,1],annotation=dat[dat$Gene %in% h$gene,3])
hi <- merge(h,hj,by="gene",all.x=TRUE)
head(hi)
dim(hi)

write.table(hi, "genes_tabla.tsv", row.names=F, quot=F, sep="\t")
write.table(hi$gene, "genes_extraer.txt", row.names=F, col.names=F, quot=F)

q6 <- q4[,c(1:9,hi$order)]
dim(q6)
head(q6[1:10,])

write.table(q6, "metadata_selecta.tsv", row.names=F, quot=F, sep="\t")

###############################################################
##### seleccionar solos las entradas acordes con el arbol #####
###############################################################

library(ape)

t <- read.tree("RAxML_bestTree_rooted.nwk")
plot(t)
tip <- gsub("'","",t$tip.label)
k <- q6[q6$sample %in% tip,]
k
dim(k)
l <- k[!duplicated(k$sample),]
dim(l)
write.table(l,"metadata_borreliella.tsv", row.names=F, sep="\t")
```

# Vic probes code 4 : Ahora, necesitamos extraer las secuencias unicas identificadas en el pangenoma de roary, blastearlas con los genomas anotados y tener multifastas con las secuencias totales por región génica 
```r
#############################
##### modificar headers (BASH) #####
#############################

for r1 in $(cat genes_extraer.txt)
do
prefix=$(basename $r1)
grep "$prefix" pan_genome_reference.fa > ${prefix}.vic
done
ls -lh *.vic
mkdir selected_genes
mv *.vic pan_genome_reference.fa genes_extraer.txt selected_genes
cd selected_genes/ ; 
cat *.vic > list.txt ; 
rm *.vic 
sed 's/>*//g' list.txt > seqs.txt
rm list.txt 
seqtk subseq pan_genome_reference.fa seqs.txt > genes_1.fasta ; 
aliview genes.fasta ;
sed 's/.*group/>group/g' genes_1.fasta > genes.fasta
sed 's/.* />/g' genes_1.fasta > genes2.fasta
seqtk subseq genes2.fasta genes_extraer.txt > genes3.fasta ; 
aliview genes3.fasta ;
ls -lh ; 

## blast BASH ##
makeblastdb -in input.fasta -dbtype nucl ;
blastn -db input.fasta -query genes3.fasta -perc_identity 80 -qcov_hsp_perc 90 -max_target_seqs 90 -outfmt 6 -num_threads 4 > blast.csv ;
ls -lh ;

## back to R ## 
blast <- read.csv("blast3/blast.csv", header=FALSE, sep="\t") ; 
names(blast) <- c("query.acc.ver","subject.acc.ver","perc.identity","alignment.length","mismatches","gap.opens","q.start","q.end","s.start","s.end","evalue","bit.score") ; 

library(tidyr)
bb <- separate(blast,"subject.acc.ver",c("a","b","c","d"),sep="_")[,2:3]
bc <- as.data.frame(paste0(bb[,1],"_",bb[,2]))
head(bc)

blast2 <- data.frame(query.acc.ver=blast[,1],sample=bc,blast[,2:12])
names(bc) <- "sample"
m3 <- m2[,c(1,3,4)]
names(m3) <- c("sample","species","cepa")
bd <- merge(bc,m3,by="sample",all.x=TRUE)
head(blast2)

names(blast2) <- c(names(blast2)[1],"sample",names(blast2)[3:13])
head(blast2)

write.table(blast2,"blast3/blast_renamed.tsv", sep="\t", row.names=F)

blast3 <- merge(blast2,bd,by="sample",all.x=TRUE)
blast4 <- blast3[!duplicated(blast3),]
blast5 <- blast4[order(blast4$query.acc.ver),]
head(blast5)

## loop ## 
getwd()
dir()
head(blast)
aa <- 0 
for (i in unique(blast$query.acc.ver)){
aa <- blast[blast$query.acc.ver == i, 2 ]
write.table(aa, paste0("blast3/",paste0(i,".txt")),row.names = F, col.names = F, quot=F)
}

write.table(blast5,"blast3/blast_renamed_2.tsv", sep="\t", row.names=F)
quit("no") ; 



## extract sequences from input.fasta ##
for r1 in $(ls *.txt)
do
prefix=$(basename $r1 .txt)
seqtk subseq input.fasta $r1 > ${prefix}.fasta ; 
done 
ls *.fasta ; 
```

# Vic probes code 5 : Ahora, generaremos plots y summary tables para una mejor identificacion de dichas regiones 
```r
## plots and statistics ## 
## boxplot ##
library(ggplot2)
a <- ggplot(blast2, aes(x=reorder(query.acc.ver,-perc.identity),y=perc.identity)) + 
geom_boxplot() + scale_y_discrete(limits=seq(from=90, to= 100, by = 0.5)) + theme_minimal() + 
theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + ylim(70,100)
a

## the loop for the summary data ##
head(blast2)
a0 <- 0
a1 <- 0
a2 <- 0
a3 <- 0 
a01 <- 0
a11 <- 0 
a21 <- 0 
a31 <- 0 
a41 <- 0 
a51 <- 0 

for (i in unique(blast$query.acc.ver)) {
  a0 <- mean(blast[blast$query.acc.ver == i , 3 ])
  a1 <- sd(blast[blast$query.acc.ver == i , 3 ])
  a2 <- unique(blast[blast$query.acc.ver == i , 1 ])
  a3 <- length(blast[blast$query.acc.ver == i , 1 ])
  a4 <- mean(blast[blast$query.acc.ver == i , 4 ])
  a5 <- sd(blast[blast$query.acc.ver == i , 4 ])
  a01 <- append(a01,a0)
  a11 <- append(a11,a1) 
  a21 <- append(a21,a2) 
  a31 <- append(a31,a3) 
  a41 <- append(a41,a4) 
  a51 <- append(a51,a5) 
}


b <- data.frame(gene=a21,mean=a01,sd=a11,strains=a31,mean_alin_length=a41,sd_alin_length=a51)
b1 <- b[2:length(b$gene),]
head(b1)
## in case ##
b2 <- data.frame(gene=rep(unique(blast$query.acc.ver),3),b1[,2:6])

write.table(b1,"blast3/blast_summary.tsv",row.names=F, sep="\t", quot=F)
write.table(b2,"blast3/blast_summary.tsv",row.names=F, sep="\t", quot=F)

head(blast5)
head(b1)

## check the results ##
head(b1)
b2 <- b1[b1$strains < 1, ]

library(plotly)
install.packages("plotly")
p1 <- ggplot(b2, aes(x = mean, y = sd, color=factor(strains),size=strains)) + 
geom_point(alpha = 0.5) + theme_minimal() + geom_text(aes(label=gene),hjust=0.5, vjust=0.5, size = 1) + 
labs(title="Ehrlichia specific CDSs") +
  scale_size(range = c(1,20), name="number of strains")
p1
ggplotly(p1)

head(b1)
dat[1:10,1:25]

############################## 
### filter more results #####

a <- read.csv("blast3/blast_summary.tsv", header=T, sep="\t")
dim(a)
head(a)

getwd()

sort(unique(a$strains))

st <- 90
a1 <- a[a$strains == st & a$mean > 98 & a$sd < 5 , ]
a1

b <- read.csv("blast3/blast_renamed_2.tsv", header=T, sep="\t")
dim(b)
head(b)

b0 <- 0
c0 <- 0
d0 <- 0
e0 <- 0
f0 <- 0

for(i in unique(a1$gene)){
b0 <- b[b$query.acc.ver == i , 14]
c0 <- unique(b0)
d0 <- paste(c0,collapse=' ')
e0 <- append(e0,i)
f0 <- append(f0,d0)
}

g0 <- data.frame(query.acc.ver=e0,species=f0)
names(g0) <- c("gene","species")

a2 <- a[a$gene %in% g0$gene , ]
a2

head(g0)
head(a2)

a3 <- merge(a2,g0,by="gene",all.x=T)
a3
head(a3)
write.table(a3,"blast3/blast_summary_species.tsv", sep="\t",row.names=F)

getwd()

library(plotly)
#install.packages("plotly")#
p1 <- ggplot(a3, aes(x = mean, y = sd, color=factor(strains),size=mean_alin_length)) + 
  geom_point(alpha = 0.5) + theme_minimal() + geom_text(aes(label=gene),hjust=0.5, vjust=0.5, size = 5) + 
  labs(title="Borreliella burgdorferi CDSs") + xlim(99,100) + 
  scale_size(range = c(1,20), name="number of strains")
p1
ggplotly(p1)

```
