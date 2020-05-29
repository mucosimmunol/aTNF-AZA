#qiime2 longitudinal analyses

qiime tools import \
--input-path otu_biom.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV100Format \
--output-path feature-table.qza



	qiime tools import \
	--type 'FeatureData[Taxonomy]' \
	--input-format HeaderlessTSVTaxonomyFormat \
	--input-path tax.txt \
	--output-path taxonomy.qza

qiime tools import \
--input-path rep-seqs.fna \
--type 'FeatureData[Sequence]' \
--output-path rep-seqs.qza


qiime tools export --input-path feature-table.qza --output-path exported-final
qiime tools export --input-path taxonomy.qza --output-path exported-final

cp exported-final/taxonomy.tsv biom-taxonomy.tsv
# change column names to #OTUID taxonomy confidence
/////////////////////////////////////////////////


qiime feature-table filter-samples \
  --i-table feature-table.qza \
  --m-metadata-file sample-metadata.txt \
  --p-where "DIS IN ('mc')" \
  --o-filtered-table mc-filtered-table.qza

qiime feature-table filter-samples \
  --i-table feature-table.qza \
  --m-metadata-file sample-metadata.txt \
  --p-where "DIS IN ('cu')" \
  --o-filtered-table cu-filtered-table.qza


#Calculate Diversity metricsfor mc
qiime diversity core-metrics --i-table mc-filtered-table.qza --p-sampling-depth 1000 --m-metadata-file sample-metadata.txt --output-dir core-metrics-mc
qiime metadata tabulate --m-input-file sample-metadata.txt --m-input-file core-metrics-mc/shannon_vector.qza --o-visualization mc_shannon


#Calculate Diversity metricsfor cu
qiime diversity core-metrics --i-table cu-filtered-table.qza --p-sampling-depth 1000 --m-metadata-file sample-metadata.txt --output-dir core-metrics-cu
qiime metadata tabulate --m-input-file sample-metadata.txt --m-input-file core-metrics-cu/shannon_vector.qza --o-visualization cu_shannon

#Pairwise comparisons for timepoints baseline vs early
qiime longitudinal pairwise-differences \
  --m-metadata-file sample-metadata.txt \
  --m-metadata-file core-metrics-mc/shannon_vector.qza \
  --p-metric shannon \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "baseline" \
  --p-state-2 "timepoint 1" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --o-visualization mc_shannon_pairwise-tp0_1.qzv

qiime longitudinal pairwise-differences \
  --m-metadata-file sample-metadata.txt \
  --m-metadata-file core-metrics-cu/shannon_vector.qza \
  --p-metric shannon \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "baseline" \
  --p-state-2 "timepoint 1" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --o-visualization cu_shannon_pairwise-tp0_1.qzv

qiime longitudinal pairwise-distances \
  --i-distance-matrix core-metrics-mc/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample-metadata.txt \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "baseline" \
  --p-state-2 "timepoint 1" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --o-visualization pairwise-distances_mc.qzv

qiime longitudinal pairwise-distances \
  --i-distance-matrix core-metrics-cu/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample-metadata.txt \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "baseline" \
  --p-state-2 "timepoint 1" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --p-parametric  \
  --o-visualization pairwise-distances_cu.qzv

qiime longitudinal pairwise-distances \
  --i-distance-matrix core-metrics-cu/jaccard_distance_matrix.qza \
  --m-metadata-file sample-metadata.txt \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "baseline" \
  --p-state-2 "timepoint 1" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --p-parametric  \
  --o-visualization pairwise-distances_jacc_cu.qzv

qiime longitudinal pairwise-distances \
  --i-distance-matrix core-metrics-mc/jaccard_distance_matrix.qza \
  --m-metadata-file sample-metadata.txt \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "baseline" \
  --p-state-2 "timepoint 1" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --p-parametric  \
  --o-visualization pairwise-distances_jacc_mc.qzv


#Pairwise comparisons for timepoints baseline vs late
qiime longitudinal pairwise-differences \
  --m-metadata-file sample-metadata.txt \
  --m-metadata-file core-metrics-mc/shannon_vector.qza \
  --p-metric shannon \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "baseline" \
  --p-state-2 "timepoint 2" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --o-visualization mc_shannon_pairwise-tp0_2.qzv

qiime longitudinal pairwise-differences \
  --m-metadata-file sample-metadata.txt \
  --m-metadata-file core-metrics-cu/shannon_vector.qza \
  --p-metric shannon \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "baseline" \
  --p-state-2 "timepoint 2" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --o-visualization cu_shannon_pairwise-tp0_2.qzv

qiime longitudinal pairwise-distances \
  --i-distance-matrix core-metrics-mc/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample-metadata.txt \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "baseline" \
  --p-state-2 "timepoint 2" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --o-visualization pairwise-distances_mc_0-2.qzv

qiime longitudinal pairwise-distances \
  --i-distance-matrix core-metrics-cu/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample-metadata.txt \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "baseline" \
  --p-state-2 "timepoint 2" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --p-parametric  \
  --o-visualization pairwise-distances_cu_0-2.qzv

qiime longitudinal pairwise-distances \
  --i-distance-matrix core-metrics-cu/jaccard_distance_matrix.qza \
  --m-metadata-file sample-metadata.txt \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "baseline" \
  --p-state-2 "timepoint 2" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --p-parametric  \
  --o-visualization pairwise-distances_jacc_cu_0-2.qzv

qiime longitudinal pairwise-distances \
  --i-distance-matrix core-metrics-mc/jaccard_distance_matrix.qza \
  --m-metadata-file sample-metadata.txt \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "baseline" \
  --p-state-2 "timepoint 2" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --p-parametric  \
  --o-visualization pairwise-distances_jacc_mc_0-2.qzv


#Pairwise comparisons for timepoints early vs late
qiime longitudinal pairwise-differences \
  --m-metadata-file sample-metadata.txt \
  --m-metadata-file core-metrics-mc/shannon_vector.qza \
  --p-metric shannon \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "timepoint 1" \
  --p-state-2 "timepoint 2" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --o-visualization mc_shannon_pairwise-tp1_2.qzv

qiime longitudinal pairwise-differences \
  --m-metadata-file sample-metadata.txt \
  --m-metadata-file core-metrics-cu/shannon_vector.qza \
  --p-metric shannon \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "timepoint 1" \
  --p-state-2 "timepoint 2" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --o-visualization cu_shannon_pairwise-tp1_2.qzv

qiime longitudinal pairwise-distances \
  --i-distance-matrix core-metrics-mc/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample-metadata.txt \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "timepoint 1" \
  --p-state-2 "timepoint 2" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --o-visualization pairwise-distances_mc_1-2.qzv

qiime longitudinal pairwise-distances \
  --i-distance-matrix core-metrics-cu/bray_curtis_distance_matrix.qza \
  --m-metadata-file sample-metadata.txt \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "timepoint 1" \
  --p-state-2 "timepoint 2" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --p-parametric  \
  --o-visualization pairwise-distances_cu_1-2.qzv

qiime longitudinal pairwise-distances \
  --i-distance-matrix core-metrics-cu/jaccard_distance_matrix.qza \
  --m-metadata-file sample-metadata.txt \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "timepoint 1" \
  --p-state-2 "timepoint 2" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --p-parametric  \
  --o-visualization pairwise-distances_jacc_cu_1-2.qzv

qiime longitudinal pairwise-distances \
  --i-distance-matrix core-metrics-mc/jaccard_distance_matrix.qza \
  --m-metadata-file sample-metadata.txt \
  --p-group-column TherResp \
  --p-state-column TP \
  --p-state-1 "timepoint 1" \
  --p-state-2 "timepoint 2" \
  --p-individual-id-column Case \
  --p-replicate-handling random \
  --p-parametric  \
  --o-visualization pairwise-distances_jacc_mc_1-2.qzv






