library(metaCCA)
for (i in 1:200) {
  locusname = sprintf("inputs/Locus.%03i", i)
  S_XY_G1 = read.table(paste0(locusname, '.S_XY_G1.txt'), header = TRUE, row.names = 1)
  S_XY_G2 = read.table(paste0(locusname, '.S_XY_G2.txt'), header = TRUE, row.names = 1)
  S_XY_G3 = read.table(paste0(locusname, '.S_XY_G3.txt'), header = TRUE, row.names = 1)
  S_XY_G4 = read.table(paste0(locusname, '.S_XY_G4.txt'), header = TRUE, row.names = 1)
  S_XY_G5 = read.table(paste0(locusname, '.S_XY_G5.txt'), header = TRUE, row.names = 1)
  nsample_G1 = 2139
  nsample_G2 = 2420
  nsample_G3 = 2472
  nsample_G4 = 2084
  nsample_G5 = 3967
  ldrsid = scan(paste0(locusname, '.LD.rsid'), what="", sep="\n")
  S_XX_G1 = read.table(paste0(locusname, '.LD.G1'), header = FALSE, row.names = ldrsid, col.names = ldrsid)
  S_XX_G2 = read.table(paste0(locusname, '.LD.G2'), header = FALSE, row.names = ldrsid, col.names = ldrsid)
  S_XX_G3 = read.table(paste0(locusname, '.LD.G3'), header = FALSE, row.names = ldrsid, col.names = ldrsid)
  S_XX_G4 = read.table(paste0(locusname, '.LD.G4'), header = FALSE, row.names = ldrsid, col.names = ldrsid)
  S_XX_G5 = read.table(paste0(locusname, '.LD.G5'), header = FALSE, row.names = ldrsid, col.names = ldrsid)
  commonSNPids = Reduce(intersect, list(rownames(S_XY_G1), rownames(S_XY_G2), rownames(S_XY_G3), rownames(S_XY_G4), rownames(S_XY_G5)))
  metares = read.table(paste0(locusname, '.metares'), header = TRUE, row.names = 1)

  mySNPids = c()
  tmpSNPids = commonSNPids
  while ((length(mySNPids) < 10) & (length(tmpSNPids) > 0)) {
    tmpres = metares[tmpSNPids,]
    selectedrsid = rownames(tmpres)[which.min(tmpres$P_value)]
    mySNPids = c(mySNPids, selectedrsid)
    allowedSNPids = colnames(S_XX_G5[selectedrsid, S_XX_G5[selectedrsid,]^2 < 0.9])
    tmpSNPids = intersect(tmpSNPids, allowedSNPids)
  }

  S_YY_G1 = estimateSyy(S_XY = S_XY_G1)
  S_YY_G2 = estimateSyy(S_XY = S_XY_G2)
  S_YY_G3 = estimateSyy(S_XY = S_XY_G3)
  S_YY_G4 = estimateSyy(S_XY = S_XY_G4)
  S_YY_G5 = estimateSyy(S_XY = S_XY_G5)
  res = metaCcaGp(nr_studies = 5,
                  S_XY = list(S_XY_G1[mySNPids,], S_XY_G2[mySNPids,], S_XY_G3[mySNPids,], S_XY_G4[mySNPids,], S_XY_G5[mySNPids,]),
                  std_info = c(0, 0, 0, 0, 0),
                  S_YY = list(S_YY_G1, S_YY_G2, S_YY_G3, S_YY_G4, S_YY_G5),
                  N = c(nsample_G1, nsample_G2, nsample_G3, nsample_G4, nsample_G5),
                  analysis_type = 2,
                  SNP_id = mySNPids,
                  S_XX = list(S_XX_G1, S_XX_G2, S_XX_G3, S_XX_G4, S_XX_G5))
  print(paste(locusname, res[,2]))
  rm (locusname, mySNPids, res, metares, commonSNPids)
  rm(list = ls(pattern = "S_XX_*"))
  rm(list = ls(pattern = "S_XY_*"))
  rm(list = ls(pattern = "S_YY_*"))
  rm(list = ls(pattern = "nsample_*"))
}
