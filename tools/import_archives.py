"""Restore selected archived products; do not fit models or copy any sample table."""
from pathlib import Path
import hashlib, csv, shutil, sys, re
root=Path(sys.argv[1]).resolve()
out=Path(__file__).resolve().parents[1]
ch=root/'Chapter2_Lampropholis_Heredity'
r6=root/'R6_Lampropholis_candidate'
prov=[]
def sha(p):
    with p.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()
def put(src,dst,change=None,status='archived; unchanged'):
    f=out/dst
    f.parent.mkdir(exist_ok=True,parents=True)
    if change: f.write_text(change(src.read_text()))
    else: shutil.copyfile(src,f)
    prov.append([dst,str(src.relative_to(root)),sha(src),sha(f),f.stat().st_size,status])
sources=[
('geographic_context','Fig1_Geographic_context_v3','01_genomic_filtering','Fig1_Geographic_context_v3.R'),
('nuclear_structure','Fig2_Nuclear_structure','02_population_structure','Fig2_Nuclear_structure.R'),
('svdquartets','Fig3_SVDQuartets_geography_bootstrap','03_phylogenetics','Fig3_SVDQuartets_geography_v4.R'),
('contemporary_gdm','Contemporary_GDM_revised','05_contemporary_GDM','Contemporary_GDM_revised.R'),
('palaeoclimate','Figure5_palaeoclimate_final','06_palaeoclimate','Figure5_palaeoclimate_final.R'),
('treemix','TreeMix_revised','07_TreeMix','TreeMix_revised.R'),
('baseline_demography','Figure7_candidate_A_R6','08_demography','Figure7_R6_reproducible.R')]
for n,(topic,fig,folder,script) in enumerate(sources,1):
    base=r6 if n==7 else ch
    for ext in ('pdf','png'):
        put(base/'01_main_figures'/f'{fig}.{ext}',f'results/figures/Figure_{n}_{topic}.{ext}')
    put(base/'04_reproducibility/scripts'/script,f'analysis/{folder}/{script}')
scripts=ch/'04_reproducibility/scripts'
for folder,script in [
    ('02_population_structure','hierarchical_snmf.R'),
    ('03_phylogenetics','run_svdquartets_ref5_random_final.paup'),
    ('05_contemporary_GDM','gdm_definitive_spatial_pipeline.R'),
    ('06_palaeoclimate','gdm_stage6_historical_landscape.R'),
    ('08_demography','stage8b_finalize.R')]:
    x=scripts/script
    f=(lambda s:s.replace('"/projects/rstudio-export/Lampro"','Sys.getenv("LAMPRO_CHELSA_DIR","data/local/chelsa")')) if folder.startswith(('05','06')) else None
    put(x,f'analysis/{folder}/{script}',f,'only CHELSA path made configurable' if f else 'archived; unchanged')
tables={
'structure':['snmf_k_selection','snmf_cross_entropy','hierarchical_snmf_parent_counts'],
'phylogeny':['svdquartets_final_run_summary'],
'treemix':['likelihood_all_runs','likelihood_summary','curated_input_summary','population_sample_sizes','TreeMix_population_Q_summary'],
'gdm':['spatial_group_null_summary','spatial_cv_model_comparison','final_ispline_importance'],
'palaeoclimate':['corrected_historical_null_summary','historical_model_comparison'],
'moments':['baseline_moments_model_comparison','baseline_moments_convergence','SFS_haploid_projections'],
'population_genetics':['TableS_Lampropholis_population_genetic_diversity','TableS_Lampropholis_pairwise_FST','TableS_Lampropholis_pairwise_FST_matrix','TableS_Lampropholis_K4_diversity_Q070','TableS_Lampropholis_geographic_ancestry_composition','TableS_Lampropholis_candidate_exclusion_FST_sensitivity']}
for folder,names in tables.items():
    for name in names: put(ch/'03_tables'/folder/(name+'.csv'),f'results/summary_tables/{folder}/{name}.csv')
for name in ('Table_moments_demography_article_R6_checked','Figure7_R6_SFS_lineage'):
    put(r6/'03_tables/moments'/(name+'.csv'),f'results/summary_tables/moments/{name}.csv')
put(r6/'03_tables/moments/Figure7_R6_source_data.csv',
    'results/summary_tables/moments/Figure7_A_model_and_SFS_cells.csv',
    lambda x:'\n'.join(l for l in x.splitlines() if l.startswith(('"kind"','"A_model"','"A_SFS_cell"')))+'\n',
    'Figure A rows selected; original numeric values unchanged')
put(ch/'04_reproducibility/notes/chapter2_palette.csv','design/idv/palette/chapter2_palette.csv')
put(r6/'02_supplementary_figures/FigureS_demography_baseline_R6.pdf','results/figures/Figure_S_demography_baseline.pdf')
def lampro(s):
    def drop(a,b):
        nonlocal s
        s=s[:s.index(a)]+s[s.index(b):]
    drop('root <- normalizePath','set.seed')
    init='''Lout <- normalizePath(args[2],mustWork=TRUE)
input <- Sys.getenv("LAMPRO_DIVERSITY_INPUTS","data/local/diversity")
sources <- c(Lgeno="ref5.rds",Lbg="lampro_landscape_background.rds",
             Lq="snmf_selected_qmatrix_k4.csv",Lmeta="TableS1_Individual_metadata.csv",
             Lvcf="ref5_exact_demographic.vcf.gz",Lcross="TreeMix_population_label_crosswalk.csv")
sources[] <- file.path(input,sources)
stopifnot(all(file.exists(sources)))
'''
    s=s.replace('set.seed(112326)',init+'set.seed(112326)')
    drop('qa <- read("Aq")','# Decode serialized')
    s=s.replace('check_vcf(sources[["Avcf"]],A,ainfo)\n','')
    s=s.replace('stopifnot(all(A[!is.na(A)] %in% 0:2),all(L[!is.na(L)] %in% 0:2),\n          !anyDuplicated(colnames(A)),!anyDuplicated(colnames(L)))',
                'stopifnot(all(L[!is.na(L)] %in% 0:2),!anyDuplicated(colnames(L)))')
    drop('run(A,Ameta','run(L,Lmeta')
    assert not re.search(r'Ameiva|/home/|/media/',s)
    return s
put(scripts/'population_genetics/population_genetics.R','analysis/09_population_diversity/population_genetics.R',
    lampro,'isolated Lampropholis branch and paths; bootstrap draws may differ because original run advanced RNG before this branch')
def fig7(s):
    draw=s[s.index('#### Figure A:'):s.index('#### Figure B:')].rstrip()
    return '''#### Rebuild candidate A from archived summary values, not model fitting
args <- commandArgs(TRUE)
stopifnot(length(args)==1,dir.exists(args[1]))
pair <- c("K1_K3","K1_K4","K2_K4")
labels <- gsub("_","–",pair)
pal <- read.csv("design/idv/palette/chapter2_palette.csv")
idv <- setNames(pal$hex,pal$key)
K <- c(K1=idv[["violet"]],K2="#06ADF1",K3="#B8D6E2",K4=idv[["amber"]])
ink <- idv[["ink"]]
z <- read.csv("results/summary_tables/moments/Figure7_A_model_and_SFS_cells.csv")
i <- z$kind=="A_model"
mod <- data.frame(pair=z$id[i],model=z$id2[i],delta_AIC=z$value[i])
arrays <- setNames(vector("list",length(pair)),pair)
for(p in pair){
  x <- z[z$kind=="A_SFS_cell" & z$id==p,]
  n <- max(x$index1)+1L; m <- max(x$index2)+1L
  stopifnot(nrow(x)==n*m,all(x$index1==rep(0:(n-1L),m)),
    all(x$index2==rep(0:(m-1L),each=n)),x$masked[1],
    max(abs(x$observed-x$expected-x$value))<1e-8)
  arrays[[p]] <- list(res=matrix(x$value,n,m),mask=matrix(x$masked,n,m))
}
'''+draw+'''
out <- file.path(args[1],"Figure_7_baseline_demography")
cairo_pdf(paste0(out,".pdf"),width=174/25.4,height=128/25.4,family="sans",pointsize=10)
draw_A()
dev.off()
png(paste0(out,".png"),width=174,height=128,units="mm",res=600,type="cairo",pointsize=10)
draw_A()
dev.off()
'''
put(r6/'04_reproducibility/scripts/Figure7_R6_reproducible.R',
    'analysis/08_demography/plot_figure7_archived.R',fig7,'standalone plot from archived values; original draw_A retained')
if len(sys.argv)>2:
    external=Path(sys.argv[2])
    for name in ('stage8_moments_models.py','stage8b_confirmatory.py','stage8b_residual_arrays.py'):
        f=external/'analysis'/name
        if f.is_file() and not re.search(r'/home/|/media/|/projects/|Ameiva',f.read_text()):
            path=f'analysis/08_demography/{name}'
            dest=out/path; shutil.copyfile(f,dest)
            prov.append([path,'Lampropholis archive/analysis/'+name,sha(f),sha(dest),dest.stat().st_size,'archived; unchanged, requires original inputs'])
    f=external/'treemix/run_treemix.sh'
    if f.is_file():
        path='analysis/07_TreeMix/run_treemix.sh'
        dest=out/path
        dest.write_text(re.sub(r'/home/[^}\n]+/bin/treemix','treemix',f.read_text()))
        prov.append([path,'Lampropholis archive/treemix/run_treemix.sh',sha(f),sha(dest),dest.stat().st_size,'only executable default path made configurable'])
(out/'docs').mkdir(exist_ok=True)
with (out/'docs/artifact_provenance.tsv').open('w',newline='') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(('path','source_role','source_sha256','staged_sha256','bytes','notes'));w.writerows(sorted(prov))

