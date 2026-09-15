#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<zlib.h>

/*
 * Annotate every SPrime marker as match, mismatch, or notcomp for one ancient
 * diploid genome. The program loads callable-mask and ancient-genotype state
 * into arrays indexed by 1-based chromosome position, then streams the score
 * table and appends the column named by --tag.
 *
 * Inputs must contain one chromosome only. BED intervals are interpreted as
 * 0-based half-open; VCF and SPrime POS values are 1-based. Run this program
 * once per chromosome and ancient individual, as done by archaic_match.smk.
 */

void print_help(char *argv[])
{
//--deepth :add read DP\n
fprintf(stderr, "\
Usage: %s \n \
--kp/--rm/--kpall :bed file for keep/romove/no-use\n \
--sep sep :define the separator in the outpout file \n \
--mskbed  maskfile :maskfile, only one allowed as the input \n \
--vcf vcffile :archaic genome in VCF, one individuals \n \
--score scorefile : scorefile from Sprime \n \
--tag reftag :tag for the added column\n", argv[0]);
}

int main(int arg, char *argv[])
{
char mskfile[1024], vcffile[1024]="", scorefile[1024]="", cmd[10240], reftag[32], sep[32]="\\t", sp;
int i, j, k, kp, dp;
sprintf(reftag, "MATCHING");
sprintf(mskfile, "NULL");
sp='\t';
dp=0;
kp=1;
if(arg == 1){print_help(argv);exit(-1);}

i=1;
while(i<arg)
{
if((strcmp(argv[i],"--mskbed")==0 || strcmp(argv[i],"--vcf")==0 ||
    strcmp(argv[i],"--score")==0 || strcmp(argv[i],"--tag")==0 || strcmp(argv[i],"--sep")==0) && i+1>=arg){
fprintf(stderr,"Missing option value\n"); return 1;
}
if(i+1<arg && (strlen(argv[i+1])>=1024 ||
  ((strcmp(argv[i],"--tag")==0 || strcmp(argv[i],"--sep")==0) && strlen(argv[i+1])>=32))){
fprintf(stderr,"Option value too long\n"); return 1;
}
//fprintf(stderr, "%s\n", argv[i]);
if(strcmp(argv[i], "--mskbed")==0){sprintf(mskfile, "%s", argv[i+1]);i=i+2;}
else if(strcmp(argv[i], "--vcf")==0){sprintf(vcffile, "%s", argv[i+1]);i=i+2;}
else if(strcmp(argv[i], "--score")==0){sprintf(scorefile, "%s", argv[i+1]);i=i+2;}
else if(strcmp(argv[i], "--tag")==0){sprintf(reftag, "%s", argv[i+1]);i=i+2;}
else if(strcmp(argv[i], "--sep")==0){sprintf(sep, "%s", argv[i+1]);i=i+2;}
else if(strcmp(argv[i], "--kp")==0){kp=1;i=i+1;}
else if(strcmp(argv[i], "--rm")==0){kp=0;i=i+1;}
else if(strcmp(argv[i], "--kpall")==0){kp=2;i=i+1;}
else if(strcmp(argv[i], "--deepth")==0){dp=1;i=i+1;}
else {
fprintf(stderr, "wrong option: %s\n", argv[i]);
print_help(argv);
exit(-1);
}
}

if(kp!=2&&strcmp(mskfile, "NULL")==0){fprintf(stderr, "Please check mask file!\n");return -1;}

//fprintf(stderr, "%s\n", sep);
//for(i=0;i<strlen(sep);i++)fprintf(stderr, "%c*%d",sep[i], i);
if(strcmp(sep, "\\t")==0)sp='\t';
if(strcmp(sep, " ")==0)sp=' ';

FILE *ifp;
gzFile gzfp;
char buffer[10240], str[1024], tmpc;
int pos;
int start, end, L_max=0;

/* Allocate only as far as the highest SPrime marker on this chromosome. */
ifp=fopen(scorefile, "r");
if(!ifp){perror(scorefile);return 1;}
fgets(buffer, 10240, ifp);
while(fgets(buffer, 10240, ifp)!=NULL){
sscanf(buffer, "%*s %d", &pos);
if(L_max<=pos)L_max=pos;
}
fclose(ifp);
if(L_max == 0){
fprintf(stderr, "check the file %s\n", scorefile);
exit(-1);
}
fprintf(stderr, "MAX LENGTH: %d\t", L_max);

/* data[position] stores callable flag plus two ancient alleles.
 * deepth[position] == -1 means that no comparable ancient genotype exists.
 */
char **data;
int *deepth;
data=(char **)malloc((L_max+1)*sizeof(char *));
deepth=(int *)malloc((L_max+1)*sizeof(int));
if(!data || !deepth){fprintf(stderr,"Out of memory\n");return 1;}
for(i=0;i<(L_max+1);i++){
data[i]=(char *)calloc(3,sizeof(char));
if(!data[i]){fprintf(stderr,"Out of memory\n");return 1;}
if(kp==1)data[i][0]='0';
if(kp==0||kp==2)data[i][0]='1';
deepth[i]=-1;
}

if(kp==1||kp==0){
/* Expand callable BED intervals into the position-indexed lookup array. */
gzfp=gzopen(mskfile, "r");
if(gzfp == Z_NULL){
fprintf(stderr, "Error: check the file %s\n", mskfile);
exit(-1);
}
while(gzgets(gzfp, buffer, 1024)!=NULL){
if(buffer[0]=='#' || strncmp(buffer,"track",5)==0 || strncmp(buffer,"browser",7)==0)continue;
if(sscanf(buffer, "%1023s %d %d", str, &start, &end)!=3 || start<0 || end<=start){fprintf(stderr,"Invalid BED row\n");return 1;}
for(i=start+1;i<=end;i++){
if(i>L_max)break;
if(kp==1)data[i][0]='1';
if(kp==0)data[i][0]='0';
}
}
gzclose(gzfp);
}
//fprintf(stderr, "INITIATION COMPLETED!\n");

/* Read the single ancient sample and fill its two diploid alleles. */
char ref[128], alt[128], DP[12800], gt[128], format[1024], *s;
gzfp=gzopen(vcffile, "r");
if(gzfp == Z_NULL){
fprintf(stderr, "Error: check the file %s\n", vcffile);
exit(-1);
}
while(gzgets(gzfp, buffer, 10240)!=NULL){
if(buffer[0]=='#')continue;
if(sscanf(buffer, "%*s %d %*s %127s %127s %*s %*s %12799s %1023s %127s", &pos, ref, alt, DP, format, gt)!=6){
fprintf(stderr,"Invalid VCF row\n");return 1;
}
/* Extract GT by its FORMAT position instead of assuming GT is first. */
int gt_index=-1, field_index=0;
char *token=strtok(format, ":");
while(token){if(strcmp(token,"GT")==0)gt_index=field_index;field_index++;token=strtok(NULL,":");}
if(gt_index<0)continue;
s=gt;
for(j=0;j<gt_index && s;j++){s=strchr(s,':');if(s)s++;}
if(!s || strlen(s)<3)continue;
if((s[0]!='0'&&s[0]!='1') || (s[2]!='0'&&s[2]!='1') || (s[1]!='|'&&s[1]!='/'))continue;
if(s[3]!='\0'&&s[3]!=':')continue;
memmove(gt,s,3);gt[3]='\0';
//fprintf(stderr, "%d %s %s %s %s\n",pos, ref, alt, DP, gt);
if(pos>0&&pos<=L_max){
if(strlen(ref)<2&&strlen(alt)<2){
if(gt[0]=='0')data[pos][1]=ref[0];
if(gt[0]=='1')data[pos][1]=alt[0];
if(gt[2]=='0')data[pos][2]=ref[0];
if(gt[2]=='1')data[pos][2]=alt[0];
deepth[pos]=1;
/* INFO/DP is optional; matching depends on a complete callable genotype. */
s=DP;
while(s && *s){
if(strncmp(s,"DP=",3)==0){int value=atoi(s+3);if(value>0)deepth[pos]=value;break;}
s=strchr(s,';');if(s)s++;
}
}
}
}
gzclose(gzfp);
//fprintf(stderr, "DATA BASE COMPLETED!\n");
//return -1;
/* Preserve every SPrime column and append the requested ancient-reference tag. */
char SNP[2];
ifp=fopen(scorefile, "r");
if(!ifp){perror(scorefile);return 1;}
fgets(buffer, 10240, ifp);
for(i=0;i<(strlen(buffer)-1);i++)printf("%c", buffer[i]);
printf("%c%s", sp,reftag);
if(dp==1)printf("%c%s_DP",sp, reftag);
printf("\n");
while(fgets(buffer, 10240, ifp)!=NULL){
if(sscanf(buffer, "%*s %d %*s %c %c %*d %d", &pos, &SNP[0], &SNP[1], &k)!=4 || pos<1 || pos>L_max || k<0 || k>1){
fprintf(stderr,"Invalid score allele\n");return 1;
}
for(i=0;i<(strlen(buffer)-1);i++)printf("%c", buffer[i]);
if(data[pos][0]=='0'||deepth[pos]<0){
printf("%cnotcomp",sp);
if(dp==1)printf("%c%d",sp, deepth[pos]);
printf("\n");
}
else {
if(SNP[k]==data[pos][1]||SNP[k]==data[pos][2])printf("%cmatch",sp);
else printf("%cmismatch",sp);
if(dp==1)printf("%c%d", sp, deepth[pos]);
printf("\n");
}

}
fclose(ifp);


fprintf(stderr, "%s COMPLETED!\n", reftag);
return 0;
}
