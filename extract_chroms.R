library(tidyverse)
library(Biostrings)

fasta_files <- list.files(path = "../", 
                          pattern = "fasta", 
                          recursive = TRUE,full.names = TRUE) %>% 
    set_names(str_extract(. , "(?<=mix_path_asm/).*(?=-pilon)"))
pilon_seqs <- fasta_files %>% map(readDNAStringSet)

seq_lengths <- pilon_seqs %>% map(width)
long_contigs <- map2(pilon_seqs, seq_lengths, ~.x[.y > 2000000])
## no contigs longer than 2 Mb for 9
long_contigs$NIST_9 <- pilon_seqs$NIST_9
long_contigs_names <- map2(.x = names(pilon_seqs), 
                           .y = long_contigs, 
                           ~paste0(.x,"|" ,names(.y))) %>% 
    set_names(names(pilon_seqs))

## Naming contigs and saving to file
for (i in names(long_contigs)) {
    names(long_contigs[[i]]) <- long_contigs_names[[i]]
    writeXStringSet(long_contigs[[i]], filepath = paste0(i,".fasta"))
}

system("cat NIST*fasta > long_contigs.fasta")
system("rm NIST*fasta")
