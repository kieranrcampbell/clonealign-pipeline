import os
import pandas as pd

configfile: "config.yaml"

# Parse samples
sample_df = pd.read_csv(config['metadata'], index_col=0).dropna()
sample_dict = sample_df.to_dict()
samples = list(sample_df.index)

#  processed datasets
files10X = ['barcodes.tsv', 'genes.tsv', 'matrix.mtx']
sces = expand("data/sces/sce_{id}.rds", id=samples)
cnv_csvs = expand("data/processed_cnv/cnv_{id}.csv", id=samples)


# Clonealign fits
var_quantiles = [0.7]
clonealign_fits = expand("data/clonealign_fits/{id}/clonealign-{id}-var_{v}.rds",
                         id=samples, v=var_quantiles)
# Clonealign analysis
clonealign_reports = expand("reports/clonealign_analysis/{id}/clonealign-analysis-{id}-var_{v}.html",
                            id=samples, v=var_quantiles)
de_fits = expand("data/differential_expression/{id}/limma-voom-de-{id}-var_{v}.rds",
                 id=samples, v=var_quantiles)



rule all:
    input:
        sces,
        cnv_csvs,
        clonealign_fits,
        clonealign_reports,de_fits

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

rule copy_number_to_gene:
    params:
        curr_dir = os.getcwd()
    input:
        lambda wildcards: os.path.join(config['cnv_clone_path'],
                                       sample_dict['cnv_clone_path'][wildcards.id])
    output:
        csv="data/processed_cnv/cnv_{id}.csv",
        report="reports/processed_cnv/processed_cnv_{id}.html",
        prev_csv="data/processed_cnv/clone_prevalence_{id}.csv"
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/parse_cnv_data.Rmd', \
        output_file='{params.curr_dir}/{output.report}', knit_root_dir='{params.curr_dir}', \
        params=list(id='{wildcards.id}', \
        input_clone_path='{input}', \
        prevalence_csv='{output.prev_csv}',\
        output_df='{output.csv}'))\" "

rule run_clonealign:
    params:
        curr_dir=os.getcwd()
    input:
        sce="data/sces/sce_{id}.rds",
        cnv="data/processed_cnv/cnv_{id}.csv"
    output:
        fit="data/clonealign_fits/{id}/clonealign-{id}-var_{v}.rds",
        report="reports/clonealign_fits/{id}/clonealign-{id}-var_{v}.html",
        cnv_mat="data/processed_cnv/cnv_mat_clonealign_input_{id}_{v}.rds"
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/run_clonealign.Rmd',\
        output_file='{params.curr_dir}/{output.report}', knit_root_dir='{params.curr_dir}',\
        params=list(id='{wildcards.id}',\
        input_sce='{input.sce}',\
        input_cnv='{input.cnv}',\
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
        cnv_df="data/processed_cnv/cnv_{id}.csv",
        sce="data/sces/sce_{id}.rds",
        prev_csv="data/processed_cnv/clone_prevalence_{id}.csv"        
    output:
        report="reports/clonealign_analysis/{id}/clonealign-analysis-{id}-var_{v}.html",
        de_fit="data/differential_expression/{id}/limma-voom-de-{id}-var_{v}.rds"
    shell:
        "Rscript -e \"rmarkdown::render('pipeline/clonealign_analysis.Rmd',\
        output_file='{params.curr_dir}/{output.report}', knit_root_dir='{params.curr_dir}',\
        params=list(id='{wildcards.id}',\
        var_quantile={wildcards.v},\
        input_sce='{input.sce}',\
        input_cnv_mat='{input.cnv_mat}',\
        input_cnv_df='{input.cnv_df}',\
        ca_fit='{input.fit}',\
        cnv_prevs='{input.prev_csv}',\
        output_voom_results='{output.de_fit}'))\" "


