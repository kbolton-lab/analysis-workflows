bsub -oo /storage1/fs1/bolton/Active/projects/mocha/UKBB/ukbb_calls/pvcf/vep/ann.log -G compute-bolton -g /jietest -q general -M 72G -R rusage[mem=72G] -a docker(kboltonlab/jie_vep:3.0) bash -c "
    /opt/vep/src/ensembl-vep/vep \
        --format vcf \
        -i /storage1/fs1/bolton/Active/projects/mocha/UKBB/ukbb_calls/pvcf/filtered/ukb23156_c9_b1_v1_filtered.vcf.gz \
        --terms SO \
        --transcript_version \
        --offline \
        --cache \
        --symbol \
        --vcf \
        -o ukb23156_c9_b1_v1.vep.annotated.vcf \
        --fasta /storage1/fs1/bolton/Active/projects/mocha/UKBB/ukbb_calls/pvcf/vep/GRCh38_full_analysis_set_plus_decoy_hla.fa     
        --dir_cache /storage1/fs1/bolton/Active/projects/mocha/UKBB/ukbb_calls/pvcf/vep  \
        --dir_plugins /opt/vep/.vep/Plugins/ \
        --synonyms /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build50f99e75d14340ffb5b7d21b03887637/chromAlias.ensembl.txt \
        --sift p \
        --polyphen p \
        --coding_only \
        --pick \
        --plugin Frameshift \
        --plugin Wildtype \
        --plugin CADD,/storage1/fs1/bolton/Active/data/hg38/CADD/prescored/whole_genome_SNVs.tsv.gz,/storage1/fs1/bolton/Active/data/hg38/CADD/prescored/gnomad.genomes.r3.0.indel.tsv.gz \
        --plugin LoFtool,/opt/vep/.vep/Plugins/LoFtool_scores.txt \
        --everything 1 \
        --assembly GRCh38 \
        --species homo_sapiens \
        --merged \
        --cache_version 104 \
        --check_existing \
        --buffer_size 5000 \
        --custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/fixed_b38_exome.vcf.gz,gnomADe,vcf,exact,1,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS \
        --custom /storage1/fs1/bolton/Active/data/hg38/gnomad.genomes.r3.0.sites.full_trimmed_info.af_only.UKBB_intervals.vcf.bgz,gnomADg,vcf,exact,1,AF,AF_ami,AF_oth,AF_afr,AF_sas,AF_asj,AF_fin,AF_amr,AF_nfe,AF_eas \
        --custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/custom_clinvar_vcf/v20210629/clinvar.vcf.gz,clinvar,vcf,exact,1,CLINSIGN,PHENOTYPE,SCORE,RCVACC,TESTEDINGTR,PHENOTYPELIST,NUMSUBMIT,GUIDELINES \
        --force_overwrite && bgzip -f ukb23156_c9_b1_v1.vep.annotated.vcf && tabix -f ukb23156_c9_b1_v1.vep.annotated.vcf.gz
    "