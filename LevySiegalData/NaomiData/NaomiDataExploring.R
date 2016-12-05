# Naomi Data New

setwd('./NZ_130721_1_GR')
load('d.Rfile')

labhaploid = which(d$pID[,2] == "FY4")
strain_imno = c()
for (i in names(labhaploid)){
  id = which(i == names(d$well.list))
  strain_imno = c(strain_imno,d$well.list[[id]])
}

setwd('..')

setwd('./NZ_130722_1_GR')
load('d.Rfile')

for (i in names(labhaploid)){
  id = which(i == names(d$well.list))
  strain_imno = c(strain_imno,d$well.list[[id]])
}

print(length(strain_imno))
print(length(unique(strain_imno)))

setwd('..')
# Saving image numbers associated with that strain to file 
lapply(strain_imno, write, "FY4_ImNo.txt", append=TRUE)
