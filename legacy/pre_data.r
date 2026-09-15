args = commandArgs(trailingOnly=TRUE)
threshold=150000
ncn=30; ncd=30; maxn=0.3; mind=0.3
allmatch=c()
for(chr in 1:22){
    x1=read.table(paste("/public/group_data_2023/xiaoyn/Sprime_msea_re/output/subpops/",args[1],"/match/out.chr",chr,".mscore",sep=""),header=T)
    y1=x1[x1$SCORE>threshold,]
    x2=read.table(paste("/public/group_data_2023/xiaoyn/Sprime_msea_re/output/subpops/",args[1],"/match/out.chr",chr,".mscore",sep=""),header=T)
    y2=x2[x2$SCORE>threshold,]

    segmentids = unique(x1$SEGMENT)
    matching = sapply(segmentids,function(s){
        z1=y1[y1$SEGMENT==s,]
        z2=y2[y2$SEGMENT==s,]
        nmatch1=sum(z1$AltaiNean=="match")
        nmatch2=sum(z2$AltaiDeni=="match")
        nmis1=sum(z1$AltaiNean=="mismatch")
        nmis2=sum(z2$AltaiDeni=="mismatch")
        if(nmatch1+nmis1>=ncn & nmatch2+nmis2>=ncd){return(c(z1$CHROM[1],mean(z1$POS),nmatch1/(nmatch1+nmis1),nmatch2/(nmatch2+nmis2)))}else{return(c(z1$CHROM[1],mean(z1$POS),-1,-1))}
    })
    allmatch=cbind(allmatch,matching[,matching[3,]>=0&matching[3,]<maxn&matching[4,]>mind])
}
xvals = allmatch[4,]

# data_path = paste("/public/group_data_2023/xiaoyn/msea/output_Altai16/subpops/",args[2], "/match/match.summary.txt", sep="")
# data = read.table(data_path, header=TRUE, sep=" ")

# filtmerge = merge(data, out, by=c("chr","seg"))
# # 将结果写入指定的输出文件
out_path = paste("/public/group_data_2023/xiaoyn/Sprime_msea_re/pipeline/",args[1], "_match_rate.txt", sep="")
# write.table(filtmerge, out_path, quote=FALSE, row.names=FALSE, sep="\t")
write.table(xvals, out_path, quote=FALSE, row.names=FALSE, sep="\t")