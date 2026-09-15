#! /bin/bash
population=("CHBL" "CHDA" "CHWA" "CMBR" "CMJA" "CMKH" "CMLA" "CMME" "CMPH" "CMTO" "LALL" "LALS" "LALT" "MYBU" "MYCH" "MYNA" "MYRA" "THW2" "VICH" "VIKI" "CMKE" "CMKR" "CMKU")
# population=("CHBL" "CHDA" "CHWA" "CMBR" "CMJA" "CMKH" "CMLA" "CMME" "CMPH" "CMTO" "LALS" "LALT" "MYBU" "MYCH" "MYNA" "MYRA" "THW2" "VICH" "VIKI" "CMKE" "CMKR" "CMKU")
# population=("LALL")
run_deal(){
pop=$1
# mkdir ${pop}
# cd /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}
# 提取标题行
# bcftools view -h /public/group_data_2023/xiaoyn/Sprime_msea_re/output/subpops/${pop}/chr1.vcf.gz | \
# tail -1 | \
# awk '{for(i=1; i<=9; i++) printf $i"\t"; for(i=10; i<=NF; i++) printf $i"\t"$i"\t"; printf "\n"} '> ${pop}.head
# echo -e "SEGMENT\tALLELE\t$(head -1 /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}.head)" > ${pop}_fin.head
# awk '
# BEGIN {0
#     FS = OFS = "\t"
# }
# {
#     for (i = 1; i <= NF; i++) {
#         sample = $i
#         if (sample_count[sample]++ > 0) {
#             $i = sample ".1"
#         }
#     }
#     print
# }' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}.head > updated_${pop}.head


# 划分haplotype
# for chr in {1..22}
# do
# bcftools view -H /public/group_data_2023/xiaoyn/Sprime_msea_re/output/subpops/${pop}/chr${chr}.vcf.gz | sed 's/|/\t/g' > /public/group_data_2023/xiaoyn/Sprime_msea_re/output/subpops/${pop}/chr${chr}.haplotype
# done
# mv /public/group_data_2023/xiaoyn/Sprime_msea_re/output/subpops/${pop}/*.haplotype /public/group_dqata_2023/xiaoyn/Sprime_msea_re/ind/${pop}
 
# 合并
# for chr in {1..22}
# do
# awk 'NR==FNR{a[$2];next} ($2 in a)'  /public/group_data_2023/xiaoyn/Sprime_msea_re/output/subpops/${pop}/match/out.chr${chr}.mscore /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/chr${chr}.haplotype > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/chr${chr}.vcf
# done

# 提取allele
# for chr in {1..22}
# do
# awk 'NR >1 {print $6"\t"$7}' /public/group_data_2023/xiaoyn/Sprime_msea_re/output/subpops/${pop}/match/out.chr${chr}.mscore > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/out.chr${chr}.col67
# 合并paste
# paste /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/out.chr${chr}.col67 /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/chr${chr}.vcf | sort -k1,1n -k4,4n > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/out.chr${chr}.vcf
# cat /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_fin.head /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/out.chr${chr}.vcf > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/chr${chr}.compare
# done

# for chr in {1..22}
# do
# awk -F'\t' 'BEGIN {OFS="\t"}
# NR==1{print};NR>1{
#     for (i = 12; i <= NF; i++){
#         $i = $i -$2
#     }
#     print 
# }' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/chr${chr}.compare > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/chr${chr}.update
# done


# for chr in {1..22}
# do
# cat /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/py_work/output/${pop}/${chr}.txt >> /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/py_work/output/${pop}/${pop}.total
# done

# 提取deni片段
# # awk '{if ($6>0.4 && $5<0.3) print $0}' /public/group_data_2023/xiaoyn/Sprime_msea_re/output/subpops/${pop}/match/match.summary.txt > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_deni.txt
# # awk '{if ($6>0.3 && $5<0.3) print $0}' /public/group_data_2023/xiaoyn/Sprime_msea_re/output/subpops/${pop}/match/match.summary.txt > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_deni1.txt
# # awk '{if ($6<0.4 && $5>0.6) print $0}' /public/group_data_2023/xiaoyn/Sprime_msea_re/output/subpops/${pop}/match/match.summary.txt > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_nean.txt
# # awk 'FNR==NR {a[$1,$2];next}($1,$6) in a' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_deni.txt /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/py_work/output/${pop}/${pop}.total > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_deni
# # awk 'FNR==NR {a[$1,$2];next}($1,$6) in a' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_deni1.txt /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/py_work/output/${pop}/${pop}.total > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_deni1
# # awk 'FNR==NR {a[$1,$2];next}($1,$6) in a' /public/group_data_2023/xiaoyn/match_rate/Deni25-sprime/${pop}_deni.txt /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/py_work/output/${pop}/${pop}.total > /public/group_data_2023/xiaoyn/match_rate/Deni25-sprime/${pop}_deni
# # awk 'FNR==NR {a[$1,$2];next}($1,$6) in a' /public/group_data_2023/xiaoyn/match_rate/Deni25-sprime/${pop}_deni1.txt /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/py_work/output/${pop}/${pop}.total > /public/group_data_2023/xiaoyn/match_rate/Deni25-sprime/${pop}_deni1
# awk 'FNR==NR {a[$1,$2];next}($1,$6) in a' /public/group_data_2023/xiaoyn/match_rate/both/${pop}_deni.txt /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/py_work/output/${pop}/${pop}.total > /public/group_data_2023/xiaoyn/match_rate/both/${pop}_deni
awk 'FNR==NR {a[$1,$2];next}($1,$6) in a' /public/group_data_2023/xiaoyn/match_rate/both/${pop}_nean.txt /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/py_work/output/${pop}/${pop}.total > /public/group_data_2023/xiaoyn/match_rate/both/${pop}_nean
awk 'FNR==NR {a[$1,$2];next}($1,$6) in a' /public/group_data_2023/xiaoyn/match_rate/both/${pop}_ambiguous.txt /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/py_work/output/${pop}/${pop}.total > /public/group_data_2023/xiaoyn/match_rate/both/${pop}_ambiguous
# # awk 'FNR==NR {a[$1,$2];next}($1,$6) in a' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_nean.txt /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/py_work/output/${pop}/${pop}.total > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_nean

# # awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_deni |sed 's/\./\t/g' > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_deni_temp
# # awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_deni1 |sed 's/\./\t/g' > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_deni1_temp
# # awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' /public/group_data_2023/xiaoyn/match_rate/Deni25-sprime/${pop}_deni |sed 's/\./\t/g' >  /public/group_data_2023/xiaoyn/match_rate/Deni25-sprime/${pop}_deni_temp
# # awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' /public/group_data_2023/xiaoyn/match_rate/Deni25-sprime/${pop}_deni1 |sed 's/\./\t/g' > /public/group_data_2023/xiaoyn/match_rate/Deni25-sprime/${pop}_deni1_temp
# awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' /public/group_data_2023/xiaoyn/match_rate/both/${pop}_deni |sed 's/\./\t/g' > /public/group_data_2023/xiaoyn/match_rate/both/${pop}_deni_temp
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' /public/group_data_2023/xiaoyn/match_rate/both/${pop}_nean |sed 's/\./\t/g' > /public/group_data_2023/xiaoyn/match_rate/both/${pop}_nean_temp
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' /public/group_data_2023/xiaoyn/match_rate/both/${pop}_ambiguous |sed 's/\./\t/g' > /public/group_data_2023/xiaoyn/match_rate/both/${pop}_ambiguous_temp
# # awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5}' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_nean |sed 's/\./\t/g' > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_nean_temp

# # 去除YRI个体
# # awk 'FNR==NR {a[$1]; next}!($6 in a)' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/outgroup.txt /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_deni_temp > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_deni_new
# # awk 'FNR==NR {a[$1]; next}!($6 in a)' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/outgroup.txt /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_deni1_temp > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_deni1_new
# # awk 'FNR==NR {a[$1]; next}!($6 in a)' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/outgroup.txt /public/group_data_2023/xiaoyn/match_rate/Deni25-sprime/${pop}_deni_temp >  /public/group_data_2023/xiaoyn/match_rate/Deni25-sprime/${pop}_deni_new
# # awk 'FNR==NR {a[$1]; next}!($6 in a)' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/outgroup.txt /public/group_data_2023/xiaoyn/match_rate/Deni25-sprime/${pop}_deni1_temp > /public/group_data_2023/xiaoyn/match_rate/Deni25-sprime/${pop}_deni1_new
# awk 'FNR==NR {a[$1]; next}!($6 in a)' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/outgroup.txt /public/group_data_2023/xiaoyn/match_rate/both/${pop}_deni_temp > /public/group_data_2023/xiaoyn/match_rate/both/${pop}_deni_new
awk 'FNR==NR {a[$1]; next}!($6 in a)' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/outgroup.txt /public/group_data_2023/xiaoyn/match_rate/both/${pop}_nean_temp > /public/group_data_2023/xiaoyn/match_rate/both/${pop}_nean_new
awk 'FNR==NR {a[$1]; next}!($6 in a)' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/outgroup.txt /public/group_data_2023/xiaoyn/match_rate/both/${pop}_ambiguous_temp > /public/group_data_2023/xiaoyn/match_rate/both/${pop}_ambiguous_new
# awk 'FNR==NR {a[$1]; next}!($6 in a)' /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/outgroup.txt /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_nean_temp > /public/group_data_2023/xiaoyn/Sprime_msea_re/ind/${pop}/${pop}_nean_new

}
export -f run_deal
parallel -j 23 run_deal ::: "${population[@]}"