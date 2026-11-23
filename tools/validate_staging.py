"""Validate the local staging content; do not run any scientific fits."""
import csv,hashlib,re,subprocess
from pathlib import Path
import xml.etree.ElementTree as ET
base=Path(__file__).resolve().parents[1]
docs=list(base.rglob('*.md'))
for f in docs:
    for target in re.findall(r'\]\(([^)]+)\)',f.read_text()):
        dest=target.split('#',1)[0]
        assert not dest or dest.startswith('https://') or (f.parent/dest).exists(),(f,dest)
with (base/'docs/artifact_provenance.tsv').open() as f:
    rows=list(csv.DictReader(f,delimiter='\t'))
for row in rows:
    f=base/row['path']
    with f.open('rb') as b: assert hashlib.file_digest(b,'sha256').hexdigest()==row['staged_sha256'],f
    assert f.stat().st_size==int(row['bytes'])
for f in base.rglob('*'):
    if not f.is_file() or '.git' in f.parts: continue
    assert not f.is_symlink()
    assert f.stat().st_size < 6_000_000, f
    assert not re.search(r'(?i)\.(?:vcf|fastq|bam|rds|rdata|nex|docx)(?:\.gz)?$',f.name),f
    if f.suffix in {'.md','.R','.py','.csv','.tsv','.yml','.sh','.paup','.cff'}:
        content=f.read_text()
        if f.parent==base/'tools':
            # Known legacy server input is quoted solely as a replacement/screening pattern.
            content=content.replace('/'+'projects/rstudio-export/Lampro','')
        assert not re.search(r'/(?:home|media|projects|Users)/[A-Za-z0-9_-]+|gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|AKIA[A-Z0-9]{16}|BEGIN [A-Z ]*PRIVATE KEY',content),f
def table(name):
    with (base/'results/summary_tables'/name).open() as f:
        return list(csv.DictReader(f))
a,=table('structure/snmf_k_selection.csv')
assert (int(a['selected_k']),int(a['minimum_cross_entropy_k']))==(4,9)
a,=table('phylogeny/svdquartets_final_run_summary.csv')
assert (int(a['tree_tips']),int(float(a['random_quartets'])),int(a['bootstrap_replicates']))==(170,100000,1000)
a,=table('treemix/curated_input_summary.csv')
assert (int(a['populations']),int(a['loci']))==(14,15955)
t=table('treemix/likelihood_all_runs.csv')
assert len(t)==90 and all(sum(int(z['m'])==m for z in t)==10 for m in range(9))
t=table('gdm/spatial_group_null_summary.csv')
assert [(round(float(x['observed_mean_rmse_gain']),6),float(x['p_value'])) for x in t]==[(.006880,.141),(.000553,.019)]
a,=table('palaeoclimate/corrected_historical_null_summary.csv')
assert (round(float(a['observed_Delta_RMSE']),6),float(a['p_value']))==(.002739,.035)
t=table('moments/Table_moments_demography_article_R6_checked.csv')
assert [round(float(z['delta_AIC_IM']),2) for z in t]==[329.18,16.80,6.32] and {z['best_baseline_model'] for z in t}=={'SC'}
t=table('population_genetics/TableS_Lampropholis_population_genetic_diversity.csv')
assert len(t)==14 and sum(int(z['n']) for z in t)==163 and {int(z['n_loci']) for z in t}=={5827}
print(f'{len(docs)} Markdown files checked; {len(rows)} archive hashes and linked files verified; numeric checks passed.')
for code,folder in [('PT','pt-BR'),('EN','en')]:
    for mode,pad in [('screen',0),('print',3),('print_cropmarks',9)]:
        pdf=base/'design/pocket-guide'/folder/f'Lampropholis_pocket_guide_{code}_{mode}.pdf'
        meta=subprocess.check_output(['pdfinfo','-f','1','-l','2','-box',pdf],text=True)
        assert re.search(r'Pages:\s+2\b',meta)
        for tag,d in [('MediaBox',(0,0,297+2*pad,210+2*pad)),('TrimBox',(pad,pad,pad+297,pad+210))]:
            box=re.findall(rf'Page\s+[12]\s+{tag}:\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)',meta)
            assert len(box)==2 and all(max(abs(float(x)-y*72/25.4) for x,y in zip(b,d))<.04 for b in box),(pdf,tag)
        fonts=subprocess.check_output(['pdffonts',pdf],text=True).splitlines()[2:]
        assert fonts and all(x.split()[-5]=='yes' for x in fonts)
        if mode!='screen': continue
        d=ET.fromstring(subprocess.check_output(['pdftotext','-bbox',pdf,'-']))
        for i,page in enumerate(d.findall('.//{*}page')):
            bounds=(0,97,197,297) if i==0 else (0,100,200,297)
            for word in page.findall('.//{*}word'):
                a,b=(float(word.attrib[k])*25.4/72 for k in ('xMin','xMax'))
                assert any(a>=x+4.9 and b<=y-4.9 for x,y in zip(bounds,bounds[1:])),(code,i,word.text,a,b)
print('Six two-page PDFs: physical sizes, trim, embedded fonts and text/fold margins verified.')

