########March 31, 2021 - https://www.qubioinfo.com/advanced-r-workshop led by Dr. Amber Paulson (Amber.Rose.Paulson@gmail.com / @Dragonflywingz / https://github.com/damselflywingz)

#DADA2: High-resolution sample inference from Illumina amplicon data - https://www.nature.com/articles/nmeth.3869

#The best place to start is the DADA2 Pipeline Tutorial (1.16) -https://benjjneb.github.io/dada2/tutorial.html

#the input is demultiplexed Illumina-sequenced 2x250 V4 region paired-end fastq files and the output is an amplicon sequence variant (ASV) table and taxonomic assignments

#https://benjjneb.github.io/dada2/dada-installation.html
library(dada2); packageVersion("dada2") #1.19.1

#install.packages('here')
library(here)
here::i_am("2021_03_31_webinar_APaulson.Rproj")

#I downloaded the mothur MiSeq SOP paired-end dataset into the /data directory of this R project from https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip
#Note: I only included four sequence libraries to cut down on required computing time for this webinar

path <- here("data/miseqsopdata/MiSeq_SOP")
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
# This does some string manupulations to get matched lists fo the forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#**Considerations for your own data: The string manipulations may have to be modified if your filename format is different.


#Inspect read quality profiles of forward reads
plotQualityProfile(fnFs[1:2])

#quality plots of the forward reads look good and the last 10 nucleotides to be trimmed off

#Inspect read quality profiles of reverse reads
plotQualityProfile(fnRs[1:2])

#quality of the reverse is worse - to be trimmed after 160.

#Filter & trim

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Filter and trim with standard parameters
#maxEE parameter sets the maximum number of 'expected errors' allowed in a read

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

#Considerations for your own data: The standard filtering parameters are starting points, not set in stone (see more information on DADA2 tutorial)
#there is a new commandline tool to determine optial triming parameters - https://github.com/Zymo-Research/figaro#figaro I have not tried

#Learn the Error Rates
#DADA2 algorithm makes use of a parametirc error model (err) and every amplicon dataset has a different set of error rates.

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE) #this plots the error rates for each possible transition (A->C, A->G,..) are shown.

#here the estimated error rates (black line) are a good fit to the points (observed rates), and the error rate (red line) drops with increased quality, which is expected.


#Sample Inference

##the core sample inference algorithm can now be applied to the filtered and trimmed sequence data

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#inspect the dada-class object
dadaFs[[1]]
help("dada-class")

#pool=TRUE, dada2 algorithim pools samples prior to sample inference ('pseudo-pooling')
#pseudo-pooling can increase sensitivity to sequence variants that might be present at low frequencies in multipel samples
dadaFs_pool <- dada(filtFs, err=errF, pool=TRUE)
dadaFs_pool[[1]]

#what changed?


#Merge paired reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

#mergers objects are a list of dataframes. Each datafram contains the merged $sequence, and information its $abundance, and indices of $forard and $reverse variants that were merged.

#construct sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#this will remove sequences that are much longer or much shorter than expected (could be result of non-specific priming)
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256] #expected length is ~250
table(nchar(getSequences(seqtab2)))

#Amber BONUS step (not in the dada2 tutorial) - collapseNoMismatch collapses identical sequences up to shifts or length variation that have no mismatches or indels when alighed
seqtab2.col <-collapseNoMismatch(seqtab2, minOverlap = 20, orderBy = "abundance",
                                 identicalOnly = FALSE, vec = TRUE, band = -1, verbose = TRUE)
#no sequences collapsed in this case

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#how many bimeras identified?
dim(seqtab)
dim(seqtab.nochim)


sum(seqtab.nochim)/sum(seqtab)
#bimeras only represent <1%. Note if large proportion of data is bimeras, consider revisiting upstream processing. This could be due to residual primer sequences with ambiguous nucleotides that were not previously removed.

#Track through pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), rowSums(seqtab2.col), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "collapsed","nonchim")
rownames(track) <- sample.names
track

#check sanity - does the filtering, denoising, merging and chimera removal make sense? Are there any big drops?

#assign taxonomy by naive Bayesian classifier method (uses a training set of references with known taxonomy)
#put the reference files in the 'reference/' sub-directory of this Rproj.
#see https://benjjneb.github.io/dada2/training.html, we used RDP trainset 16 (https://doi.org/10.5281/zenodo.801827)
#reference files formatted for assignTaxonomy - https://zenodo.org/record/2541239
#reference files formatted for assgnSpecies - https://zenodo.org/record/2658728#.XM7qtBNKhYh
#note many other reference databases exist, Silva, Greengenes, GenBank, etc. I have used Greengenes in the past also.

#this a slow part of the tutprial (~ 5 minutes) - genus level is considered at 97 % but this is somewhat arbitrary!
taxa <- assignTaxonomy(seqtab.nochim, here("reference/RefSeq-RDPv2_16S_species.fa"))

#add the species level assignment - this is for matches 100 % 
taxa <- addSpecies(taxa, here("reference/RefSeq-RDP_dada2_assignment_species.fa"))
head(taxa)

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#see Dada2 tutorial for altertnative taxonomy assignment by IDTAXA algorithim


#hand-off to phyloseq #https://joey711.github.io/phyloseq/

#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')

library(phyloseq)
packageVersion("phyloseq") #1.34.0

library(Biostrings); packageVersion("Biostrings") #2.58.0

library(ggplot2); packageVersion("ggplot2") #3.3.2

theme_set(theme_bw())

#now construct a data.frame from the information in the file names (this can be read in as a table or CSV too*)
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)

#this part pulls values from subject
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"

#renames the rows
rownames(samdf) <- samples.out
samples.out

#now create the phyloseq object - otu table (or in this case ASV table) is made from the count table and the 'sample_data' is made from samdf (metadata)

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

#check out the ps object
ps

head(otu_table(ps))
#this is ugly the sequences are the names!

#rename with short names "ASV1", "ASV2", etc.

#this pulls the sequences used as the names into an object called dna
dna <- Biostrings::DNAStringSet(taxa_names(ps))
head(dna)

names(dna) <- taxa_names(ps)
head(dna)

ps.1 <- merge_phyloseq(ps, dna)
ps.1
#now there is a refseq() included in the ps object

#create the short names
taxa_names(ps.1) <- paste0("ASV", seq(ntaxa(ps.1)))
ps.1

otu_table(ps.1)
#now it is easier to see the count data!

#meta-data
sample_data(ps.1)

#taxonomies
tax_table(ps.1)

##############end of dada2 tutorial - https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf




####Start of Phyloseq tutorial
##a small dataset is provided as a 'test set' that is made of four RDS objects found in the /data subdirectory

ASV_small <- readRDS(here("data/ASV_small.RDS"))
taxa2 <- readRDS(here("data/taxa.RDS"))
sample <- readRDS(here("data/sample.RDS"))
seq <- readRDS(here("data/seq.RDS"))


#import the phyloseq object by merging in each of the RDS objects
ps2 <- phyloseq(otu_table(ASV_small, taxa_are_rows =FALSE),
               sample_data(sample),
               tax_table(taxa2),
               refseq(seq))

##Now we will loosely follow the Phyloseq tutorial available on F1000 Reserach - https://f1000research.com/articles/5-1492
##data import tutorial is also a good overview to run through on your own - https://joey711.github.io/phyloseq/import-data.html

ps2

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 6295 taxa and 4 samples ]
#sample_data() Sample Data:       [ 4 samples by 11 sample variables ]
#tax_table()   Taxonomy Table:    [ 6295 taxa by 8 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 6295 reference sequences ]

#typical reproducible microbiome workflows will removal suspected non-bacterial taxa not assigned at phylum-level
ps3 <- subset_taxa(ps2, !is.na(Kingdom) & !Phylum %in% c("", "uncharacterized"))
ps3 <- subset_taxa(ps3, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

#similarly, remove contaminating chloroplast sequences (come from plants)
ps3 <- subset_taxa(ps3, (Family!="Chloroplast") | is.na(Family)) ##remove chloroplast-associated interferring sequences
#How many ASVs where removed?

plot_richness(ps3, x="Sample_AC", measures=c("Shannon", "Simpson"), color="Sample_AC")

#access the sample variables and the total library counts
sample_data(ps3)
sample_sums(ps3)

####Create a diagnostic biplot:

#ordination - cannot cover in webinar today please see - "Ten quick tips for effective dimensionality reduction" - https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006907
ps3.ord <- ordinate(ps3, "NMDS", "bray")

#bi-plot is 'double ordination' because we ordinate for taxa and sample at the same time to look for outliers
bi_plot <- 
  plot_ordination(
    physeq = ps3, 
    ordination = ps3.ord, 
    type="biplot",
    color="Phylum", 
    shape="Sample_full",
    label="Species")
bi_plot

#inspect the bi-plot for outliers - usually there are more libraries this is just an example test set

##inspect the rarefacation curves to ensure sampling depth was sufficient to capture most of the diversity in the community

###rarefaction wrapper ggrare of the ranacapa package, based on vegan rarecurve 
library(vegan)

#if (!requireNamespace("devtools", quietly = TRUE))
#  install.packages('devtools')

#library(devtools)
devtools::install_github("gauravsk/ranacapa")
library(ranacapa)

p <- ggrare(ps3, step = 50, color = "Sample_AC", label = "Sample_AC", se = FALSE)
p + theme_bw() + guides(fill=FALSE, color=FALSE)

# to remove low libraries
ps3 <- prune_samples(sample_sums(ps3) >= 5000, ps3)

#to sort for most abundant
top50 <- names(sort(taxa_sums(ps3), decreasing=TRUE))[1:50]
ps.top50 <- prune_taxa(top50, ps3)
ps.top50

#this transforms the counts to a relative value
ps.top50_rel = transform_sample_counts(ps.top50, function(otu) otu/sum(otu))
sample_sums(ps.top50)
sample_sums(ps.top50_rel)

library(scales)
plot <- plot_bar(ps.top50_rel, x="Sample_AC",fill="Genus") + 
  geom_bar(stat="identity", position="stack",color=NA) + 
  facet_grid(~Location1, scales="free_x",space="free_x")+
  scale_y_continuous(labels =comma)+
  theme(axis.text.x = element_text(size = 8))

plot + theme_bw()+
  theme(axis.text.x = element_text(siz=8, angle=270))

#align the sequences
library(DECIPHER)
seqs <- refseq(ps.top50)
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

#construct the phylogeny
library(phangorn)
#please note there are many methods to construct phylogenies - no time to go through today = https://rstudio-pubs-static.s3.amazonaws.com/345955_fba1ccbdcd8f424aa5505c15bfd75bf7.html

#first construct a neighbor-joining tree, and then fit a GTR+G+I (Generalized time-reversible with Gamma rate variation) maximum likelihood tree using the neighbor-joining tree as a starting point
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)

#this is computational bottleneck - takes about ~1 minute with 50 sequences
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

#add the phylogeny to the object
ps.top50 = merge_phyloseq(ps.top50,fitGTR$tree)
ps.top50
#phy_tree is there

#phyloseq can draw a tree - https://joey711.github.io/phyloseq/plot_tree-examples.html
plot_tree(ps.top50)
plot_tree(ps.top50, ladderize="left", color="Sample_AC",label.tips="Genus")

####There is many resources, publications and training available for advanced analysis of microbiome data in phyloseq
####Please search Dr. Susan Holmes microbiome on Youtube for lectures to learn more, or sign up for a course in microbiome analysis (we did preliminary outlier checking today)

#output a fasta file from a phyloseq object - https://rdrr.io/github/DanielSprockett/reltools/man/save_fasta.html

#install.packages("remotes")
#remotes::install_github("DanielSprockett/reltools")
library(reltools)

save_fasta(ps.top50, file = here("data/top50.fasta"), rank = "Genus")

####end of Phyloseq demo

#######################
#interactive lesson on following up - view the fasta file in a text editor and copy the >NA sequence to see if we can figure out better what type of bacteria this may be
#go to https://blast.ncbi.nlm.nih.gov/Blast.cgi and and blastn against 16S ribosomal RNA sequences (Bacteria and Arcaea)
#view the descriptions, graphic summary, alignments and taxonomy tabs
#Under 'other report' make a distance tree of results - you can make a phylogeny through genbank too!
#you can also retrieve fasta files, alignment files or other reports from Genbank


#~*~*~*end
