# func_equiv

checking transmission / mendel errors / de novos for 19 quads and 8 trios

Experiment - parents and offspring were called in different VCFs
 - parent1_and_kids2.vcf
 - parent2_and_kids1.vcf
 - variant calling is unaffected by the familial relationship

OUTPUT file: family_transmission_denovo_summary.txt

VARIANT FILTERING: none

COLUMN description:

IID - offspring ID
FID - family number
PID - father ID
MID - mother ID
SEX - offspring sex
AFF - affection status

father_shet_count - number of AC=1 hets in father
mother_shet_count - number of AC=1 hets in mother
father_matched - number of AC=1 hets in father found in offspring
mother_matched - number of AC=1 hets in mother found in offspring
father_tr - singleton transmission rate from father to offspring
mother_tr - singleton transmission rate from mother to offspring
parental_tr - singleton transmission rate from parent to offspring
scheme - gVCF blocking scheme
N_MENDEL_ERRORS - total mendel errors (including de novo)
NSNP_MENDEL_ERRORS - total SNP mendel errors (including de novo)
denovo_count_unfiltered - denovo count (min PL = 20, min child AB = 0.2)
denovo_count_filtered - denovo count after posterier probability applied (note: posterior prob uses AF from exome sequence data)


