# Calculate TSS windows
Rscript ../R/calc_wind.R \
-r ../ref/toy_gene_annot_hg19.txt \
-w ../ref/toy_gene_annot_hg19_tsswindow.txt

# Map SVs to genes
Rscript ../R/mapsv.R \
-r ../ref/toy_gene_annot_hg19_tsswindow.txt \
-i ./data/sample_ids.txt \
-b ./data/bedpe/ \
-w ./intermediate/

# Quantile normalize expression data
Rscript ../R/quantNorm.R \
-e ./data/exp.txt \
-w ./intermediate/ \
-i ./data/sample_ids.txt

# Calculate normal scores for gene expression
Rscript ../R/normal_score.R \
-e ./intermediate/DAT.exp_quant.txt \
-w ./intermediate/

# Run HYENA for 0 to 5 PCs
for i in {0..5}
  do
  Rscript ../R/HYENA.R \
  -e ./intermediate/DAT.exp.Rdata \
  -s ./intermediate/DAT.sv_mapped_filtered_numSV.txt \
  -i ./data/sample_ids.txt \
  -a ../ref/toy_gene_annot_hg19_tsswindow.txt \
  -p ./data/purity.txt \
  -c ./data/cna.txt \
  -m ./data/clindat.txt \
  -d ./results/run_PC${i}/ \
  -w ./intermediate/run_PC${i}/ \
  --pur --cn -C 10 --age --sex -f 5 --PC -n ${i}
  done

# Permute expression 5 times
for j in {1..5}
  do
  Rscript ../R/pmt_exp.R \
  -e ./intermediate/DAT.exp.Rdata \
  -w ./intermediate/exp_pmt/ \
  -x PMT${j}
  done

# Run HYENA on permuted expression
for i in {0..5}
  do
  for j in {1..5}
    do
    Rscript ../R/HYENA.R \
    -e ./intermediate/exp_pmt/PMT${j}.exp.Rdata \
    -s ./intermediate/DAT.sv_mapped_filtered_numSV.txt \
    -i ./data/sample_ids.txt \
    -a ../ref/toy_gene_annot_hg19_tsswindow.txt \
    -p ./data/purity.txt \
    -c ./data/cna.txt \
    -m ./data/clindat.txt \
    -d ./intermediate/run_PC${i}/ \
    -w ./intermediate/run_PC${i}/ \
    -x PMT${j} \
    --pur --cn -C 10 --age --sex -f 5 --PC -n ${i} --pmt
    done
  done

# Compile permuted p-values
ls -d intermediate/run_PC* | while read DIR
  do
  [ -e $DIR/pvalpmt.txt ] && rm $DIR/pvalpmt.txt
  touch $DIR/pvalpmt.txt
  tail -q -n +2 $DIR/*_posEstimate.txt >> $DIR/pvalpmt.txt
  cp $DIR/pvalpmt.txt $DIR/temp.txt
  cut -f17 $DIR/temp.txt > $DIR/pvalpmt.txt
  rm $DIR/temp.txt
  done

# Calculate empirical p-values
for i in {0..5}
  do
  Rscript ../R/empiricalp.R \
  -o ./results/run_PC${i}/*_posEstimate.txt \
  -p ./intermediate/run_PC${i}/pvalpmt.txt \
  -d ./results/run_PC${i}/
  done

# Set the PC that reaches 80% power for the final results
Rscript ../R/setpc.R \
-l ./results/run_PC0/DAT.annot_sv_pur_cn_age_sex_PC0_posEstimate_adj.txt,./results/run_PC1/DAT.annot_sv_pur_cn_age_sex_PC1_posEstimate_adj.txt,./results/run_PC2/DAT.annot_sv_pur_cn_age_sex_PC2_posEstimate_adj.txt,./results/run_PC3/DAT.annot_sv_pur_cn_age_sex_PC3_posEstimate_adj.txt,./results/run_PC4/DAT.annot_sv_pur_cn_age_sex_PC4_posEstimate_adj.txt,./results/run_PC5/DAT.annot_sv_pur_cn_age_sex_PC5_posEstimate_adj.txt \
-w ./intermediate/ \
-d ./results/ \
--emp

# Filter our eQTL driven expression changes
Rscript ../R/eqtl.R \
-Q ./data/Thyroid.signifpairs.goi.txt \
-S ./data/snv.txt \
-R ./results/*_posEstimate_adj_sig.txt \
-V ./intermediate/DAT.sv_mapped_filtered_numSV.txt \
-i ./data/sample_ids.txt \
-a ../ref/toy_gene_annot_hg19_tsswindow.txt \
-w ./intermediate/ \
-d ./results/ \
-g hg19 \
--emp

# plot expression
Rscript ../R/expr_plotter.R \
-R ./DAT.annot_sv_pur_cn_age_sex_PC0_posEstimate_adj_sig_rmeqtl.txt \
-d ./intermediate/run_PC0/matrix/ \
-w ./results/plots/
