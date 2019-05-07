import os
import pandas as pd

configfile: "config.yaml"

# Parse samples
sample_df = pd.read_csv(config['metadata'], index_col=1, dtype=str).dropna()


sample_dict = sample_df.to_dict()
samples = list(sample_df.index)
pdxs = list(set(sample_dict['pdx'].values()))

#  processed datasets
files10X = ['barcodes.tsv', 'genes.tsv', 'matrix.mtx']
sces = expand("data/sces/sce_{id}.rds", id=samples)
sces_qc = expand("data/sces/sce-qc_{id}.rds", id=samples)


# Let's work out what clones map to what samples

cnv_csvs = expand("data/processed_cnv/cnv_{pdx}.csv", pdx=pdxs)

# Clonealign fits
var_quantiles = [0.5]
clonealign_fits = expand("data/clonealign_fits/{id}/clonealign-{id}-var_{v}.rds",
                         id=samples, v=var_quantiles)
# Clonealign analysis
clonealign_reports = expand("reports/clonealign_analysis/{id}/clonealign-analysis-{id}-var_{v}.html",
                            id=samples, v=var_quantiles)
# de_fits = expand("data/differential_expression/{id}/limma-voom-de-{id}-var_{v}.rds",
#                  id=samples, v=var_quantiles)



rule all:
    input:
        sces_qc,
        cnv_csvs,
        clonealign_fits,
        clonealign_reports

rule read_qc_scrna:
    params:
        scrna_path=config['scrna_path'],
        curr_dir = os.getcwd(),
        max_mito = config['max_mito'],
        min_features = config['min_features'],
        input_dir = lambda wildcards: os.path.join(config['scrna_path'], wildcards.id)
    input:
        expand(os.path.join(config['scrna_path'], "{{id}}", "{f}"), f=files10X)
    output:
        sce="data/sces/sce_{id}.rds",
        report="reports/scrnaqc/scrnaqc_{id}.html"
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/read_qc_scrna.Rmd', \
        output_file='{params.curr_dir}/{output.report}', knit_root_dir='{params.curr_dir}', \
        params=list(id='{wildcards.id}', \
        max_mito={params.max_mito}, \
        min_features={params.min_features},\
        input_10X_path='{params.input_dir}', \
        output_sce_path='{output.sce}'))\" "

rule manual_qc_SA609:
    params:
        curr_dir = os.getcwd()
    input:
        [f for f in sces if "SA609" in f]
    output:
        sces=[f for f in sces_qc if "SA609" in f],
        report="reports/scrnaqc/SA609-mnn-reqc.html"
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/re-qc-SA609.Rmd', \
        output_file='{params.curr_dir}/{output.report}', knit_root_dir='{params.curr_dir}/pipeline')\" "
    

rule copy_number_to_gene:
    params:
        curr_dir = os.getcwd()
    output:
        csv="data/processed_cnv/cnv_{pdx}.csv",
        report="reports/processed_cnv/processed_cnv_{pdx}.html",
        prev_csv="data/processed_cnv/clone_prevalence_{pdx}.csv"
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/parse_cnv_data-may-retest.Rmd', \
        output_file='{params.curr_dir}/{output.report}', knit_root_dir='{params.curr_dir}', \
        params=list(pdx='{wildcards.pdx}', \
        prevalence_csv='{output.prev_csv}',\
        output_df='{output.csv}'))\" "

rule run_clonealign:
    params:
        curr_dir=os.getcwd(),
        pdx = lambda wildcards: sample_dict['pdx'][wildcards.id]
    input:
        sce="data/sces/sce-qc_{id}.rds",
        cnv=lambda wildcards: "data/processed_cnv/cnv_{}.csv".format(sample_dict['pdx'][wildcards.id]),
        prev_csv=lambda wildcards: "data/processed_cnv/clone_prevalence_{}.csv".format(sample_dict['pdx'][wildcards.id])

    output:
        fit="data/clonealign_fits/{id}/clonealign-{id}-var_{v}.rds",
        report="reports/clonealign_fits/{id}/clonealign-{id}-var_{v}.html",
        cnv_mat="data/processed_cnv/cnv_mat_clonealign_input_{id}_{v}.rds"
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/run_clonealign.Rmd',\
        output_file='{params.curr_dir}/{output.report}', knit_root_dir='{params.curr_dir}',\
        params=list(id='{wildcards.id}',\
        pdx='{params.pdx}',\
        input_sce='{input.sce}',\
        clone_prevs='{input.prev_csv}',\
        output_rds='{output.fit}',\
        output_cnv='{output.cnv_mat}',\
        gex_var_quantile={wildcards.v},\
        max_cnv_var=0.5))\" "

rule clonealign_analysis:
    params:
        curr_dir = os.getcwd()
    input:
        fit="data/clonealign_fits/{id}/clonealign-{id}-var_{v}.rds",
        cnv_mat="data/processed_cnv/cnv_mat_clonealign_input_{id}_{v}.rds",
        cnv_df=lambda wildcards: "data/processed_cnv/cnv_{}.csv".format(sample_dict['pdx'][wildcards.id]),
        sce="data/sces/sce-qc_{id}.rds",
        prev_csv=lambda wildcards: "data/processed_cnv/clone_prevalence_{}.csv".format(sample_dict['pdx'][wildcards.id])
    output:
        report="reports/clonealign_analysis/{id}/clonealign-analysis-{id}-var_{v}.html"#,
#        de_fit="data/differential_expression/{id}/limma-voom-de-{id}-var_{v}.rds"
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/clonealign_analysis.Rmd',\
        output_file='{params.curr_dir}/{output.report}', knit_root_dir='{params.curr_dir}',\
        params=list(id='{wildcards.id}',\
        var_quantile={wildcards.v},\
        input_sce='{input.sce}',\
        input_cnv_mat='{input.cnv_mat}',\
        input_cnv_df='{input.cnv_df}',\
        ca_fit='{input.fit}',\
        clone_prevs='{input.prev_csv}'))\" "



