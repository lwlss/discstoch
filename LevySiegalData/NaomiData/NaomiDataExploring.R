# Naomi Data New

fpath='./NZ_130721_1_GR'
load(file.path(fpath,'/d.Rfile'))

labhaploid = which(d$pID[,2] == "FY4")
strain_impath = c()
for (i in names(labhaploid)){
  id = which(i == names(d$well.list))
  fnames=sprintf("t%02dxy%04d.jpg",replicate(length(d$well.list[[id]]),1:500),d$well.list[[id]])
  fpaths=file.path(fpath,"jpgGRimages",fnames)
  fpaths=fpaths[file.exists(fpaths)]
  strain_impath = c(strain_impath,fpaths)
}

fpath='./NZ_130722_1_GR'
load(file.path(fpath,'/d.Rfile'))

for (i in names(labhaploid)){
  id = which(i == names(d$well.list))
  fnames=sprintf("t%02dxy%04d.jpg",replicate(length(d$well.list[[id]]),1:500),d$well.list[[id]])
  fpaths=file.path(fpath,"jpgGRimages",fnames)
  fpaths=fpaths[file.exists(fpaths)]
  strain_impath = c(strain_impath,fpaths)
}

print(length(strain_impath))
print(length(unique(strain_impath)))

# Saving image numbers associated with that strain to file 
#lapply(strain_imno, write, "FY4_ImNo.txt", append=TRUE)

# Create new directory for storing links to images
outdir="~/FY4"
dir.create(outdir)
outfiles=file.path(outdir,paste(substr(strain_impath,3,16)[1],basename(strain_impath),sep="_"))
file.symlink(normalizePath(strain_impath),outfiles)

