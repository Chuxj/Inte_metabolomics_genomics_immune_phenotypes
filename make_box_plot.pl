####set SNP and cytokine

####set Effect Allele and Alternative allele
$EA=$ARGV[0];
$AA=$ARGV[1];

$SNP=$ARGV[2];

$cytokine=$ARGV[3];

#####

$dosage=`grep -w $SNP /groups/umcg-wijmenga/tmp04/umcg-lchang/tmp02/lung_ced_T1D_Lyme_celiac_HIV/post_impu/T1D_mapping/data_used_in_mapping_T1D/snp_for_cyto_mapping.txt`;

#print "$dosage\n";

$exp=`grep -w $cytokine /groups/umcg-wijmenga/tmp04/umcg-lchang/tmp02/lung_ced_T1D_Lyme_celiac_HIV/post_impu/T1D_mapping/data_used_in_mapping_T1D/exp_for_cyto_mapping.txt`;

open OUT, ">tmp";
print OUT "$dosage";
print OUT "$exp";
close OUT;

`module load R`;

$pwd=`pwd`;

`Rscript ~/pool/boxplot_function_in_perl.R -e $EA -a $AA -p $pwd`;

`rm tmp`;
