#------ This function elaborate the bins when reads are paired-end -----# 
# Bins are the node of the splicing graph. 


bins_paired <- function(inpf, npetype, tophat.exons, count.exons, len.exons, n.exons, readtype, gap, std){
   
   longtype <- vector(mode='list') # record each created type
   numbertype <- vector(mode='list') # record the number of each type
   compteur <- 0
   
   while(compteur<npetype){
      compteur <- compteur+1
      plouf <- scan(inpf, what=character(0), nlines=2, quiet=TRUE)
      sgone <- as.numeric(plouf[1]) # left type index
      sgtwo <- as.numeric(plouf[2]) # right type index
      nimp <- as.numeric(plouf[3]) # number of such types
      
      oneleft <- which(readtype[sgone,1:n.exons]==1) # type of left read
      oneright <- which(readtype[sgtwo,1:n.exons]==1) # type of right read
      lastleft <- oneleft[length(oneleft)]
      firstright <- oneright[1]
      
      implicated <- as.integer(unlist(strsplit(plouf[4:(4+as.numeric(plouf[3])-1)],split=':')))
      npair <- implicated[seq(2,2*nimp,2)] # number of such paired-end reads
      alldist <- implicated[seq(1,2*nimp,2)] # all insert distances
      
      # a) left and right reads overlap or stop/start into same exon or neighbour exons -- EASY case 
      if(lastleft>=(firstright-1)){
         vecbin <- unique(sort(c(oneleft,oneright)))
         newbin <- paste(vecbin, collapse='-')
         longtype[[newbin]] <- vecbin
         numbertype[[newbin]] <- c(numbertype[[newbin]],npair)
      }
      
      # b) left and right reads stop/start not into neighbours bins -- MIGHT BE AMBIGUOUS
      if(firstright>lastleft+1){
         # b1) only jump an intron (active when annotation in input)
         if(sum(count.exons[(lastleft+1):(firstright-1)])==0){ # be less stringent ????
            vecbin <- unique(sort(c(oneleft,oneright)))
            newbin <- paste(vecbin, collapse='-')
            longtype[[newbin]] <- vecbin
            numbertype[[newbin]] <- c(numbertype[[newbin]],npair)
         }
         # b2) complicated cases
         if(sum(count.exons[(lastleft+1):(firstright-1)])>0){
            djump <- alldist - (tophat.exons[firstright,1] - tophat.exons[lastleft,2]) # effective insert when jump all inside bins
            ppp <- which(count.exons[(lastleft+1):(firstright-1)]>0)
            addno <- sum(len.exons[(lastleft+1):(firstright-1)][ppp])
            dnojump <- djump + addno # effective insert when no jump at all (exept introns of course)
            seuil.nojump <- min(gap-std, gap+std-addno)
            ind.nojump <- which(djump<seuil.nojump) # indexes of no jumping reads
            if(length(ind.nojump)>0){
               vecbin <- unique(sort(c(oneleft,((lastleft+1):(firstright-1))[ppp],oneright)))
               newbin <- paste(vecbin, collapse='-')
               longtype[[newbin]] <- vecbin
               numbertype[[newbin]] <- c(numbertype[[newbin]],npair[ind.nojump])
            }
            seuil.jump <- max(gap-std, gap+std-addno)
            ind.jump <- which(djump>seuil.jump) # indexes of jumping reads
            if(length(ind.jump)>0){
               vecbin <- unique(sort(c(oneleft,oneright)))
               newbin <- paste(vecbin, collapse='-')
               longtype[[newbin]] <- vecbin
               numbertype[[newbin]] <- c(numbertype[[newbin]],npair[ind.jump])
            }
         }
      }
   }

   count.first <- sapply(numbertype, sum)

   return(list(longtype=longtype, count.first=count.first))

}
