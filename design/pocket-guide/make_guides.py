"""Illustrated Lampropholis field guide: print and screen layouts. Produce English and Portuguese pocket-guide layouts from the existing scientific colour palette, artwork and geographic summaries. Keep specimen illustrations, labels and sampling geography consistent with the reviewed project materials; figure design does not change the associated genetic analyses."""
from collections import Counter
import csv,html,json,math,re,subprocess,sys,tempfile
from pathlib import Path
from PIL import Image,ImageDraw,ImageFont

base=Path(__file__).parent.resolve(); stage=base.parents[1]
ink='#251462'; violet='#6811DB'; amber='#EBB201'; wash='#F7F3FB'; rule='#C1B8CE'; white='#FFFFFF'
ks=(violet,'#06ADF1','#B8D6E2',amber)
sans='DejaVu Sans'; serif='DejaVu Serif'
font=lambda family: subprocess.check_output(('fc-match','-f','%{file}',family),text=True)
fontfile=font(sans)
def text(x,y,s,n=3.65,c=ink,ff=sans,b=False,it=False):
    return f'<text x="{x}" y="{y}" fill="{c}" font-size="{n}" font-family="{ff}" font-weight="{"bold" if b else "normal"}" font-style="{"italic" if it else "normal"}">{html.escape(s)}</text>'
def para(x,y,s,w=79,n=3.60,step=5.3):
    f=ImageFont.truetype(fontfile,round(n*20)); lines=[]
    for word in s.split():
        if lines and f.getlength(lines[-1]+' '+word)/20<=w: lines[-1]+=' '+word
        else: lines.append(word)
    return ''.join(text(x,round(y+i*step,3),l,n) for i,l in enumerate(lines)),y+len(lines)*step
def rect(x,y,w,h,col):return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" fill="{col}"/>'
def line(x,y,u,v,col=rule,width=.32):return f'<path d="M{x},{y}L{u},{v}" stroke="{col}" stroke-width="{width}" fill="none"/>'
def top(x,n,label,h):
    return text(x,13,f'{n:02d} / {label}',3,violet,b=True)+line(x,17,x+80,17)+''.join(text(x,30+i*7.4,t,6.35,ink,serif) for i,t in enumerate(h))
def svg(g,w=303,h=216,x=-3,y=-3):
    return f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}mm" height="{h}mm" viewBox="{x} {y} {w} {h}">{g}</svg>\n'
assets=base/'shared-assets'; assets.mkdir(exist_ok=True)
if len(sys.argv)>1:
    root=Path(sys.argv[1])
    with (root/'Chapter2_Lampropholis_Heredity/04_reproducibility/plot_inputs/01_sampling/Fig1_panel_sampling_audit.csv').open() as f:
        samples=list(csv.DictReader(f))
    assert len(samples)==101 and sum(int(s['n']) for s in samples)==170
    bins=Counter((math.floor(float(s['lon'])/3)*3,math.floor(float(s['lat'])/3)*3) for s in samples)
    keep=[(x,y,n) for (x,y),n in sorted(bins.items()) if n>=3]
    with (assets/'sampled_3_degree_bins.csv').open('w',newline='') as f:
        a=csv.writer(f);a.writerow(('west_deg','south_deg','sampled_localities'));a.writerows(keep)
    geodata=json.loads((root/'results/Lampropholis_Fig1_Geographic_context/data/abs_states_generalized.geojson').read_text())
    xy=lambda lon,lat:(3.4*(lon-129),2.85*(-10-lat))
    paths=[]
    for ft in geodata['features']:
        geo=ft['geometry']
        for poly in (geo['coordinates'] if geo['type']=='MultiPolygon' else [geo['coordinates']]):
            outer=poly[0]
            area=abs(sum(a[0]*b[1]-a[1]*b[0] for a,b in zip(outer,outer[1:])))/2
            if area<.075: continue
            d=''.join('M'+'L'.join(f'{x:.1f},{y:.1f}' for x,y in (xy(*v[:2]) for v in rr))+'Z' for rr in poly)
            paths.append(f'<path d="{d}"/>')
    land=''.join(paths)
    cells=''
    for west,south,n in keep:
        x,y=xy(west,south+3)
        cells+=f'<rect x="{x}" y="{y}" width="10.2" height="8.55" opacity="{.2+.68*math.sqrt(n/max(q[2] for q in keep)):.3f}" fill="{violet}" stroke="white" stroke-width=".35"/>'
    label=''.join(text(x,y,s,3.5) for x,y,s in [(51,37,'QLD'),(60,70,'NSW'),(58,85,'VIC'),(65,96,'TAS'),(21,69,'SA')])
    sub=f'<svg width="65" height="73" viewBox="0 0 91 101"><defs><clipPath id="east"><rect width="91" height="101"/></clipPath><clipPath id="land">{land}</clipPath></defs><g clip-path="url(#east)"><g fill="#F0EAF6" stroke="#9A8DA9" stroke-width=".33" fill-rule="evenodd">{land}</g><g clip-path="url(#land)">{cells}</g>{label}</g></svg>'
    (assets/'sampling_map.svg').write_text(sub.replace('<svg ','<svg xmlns="http://www.w3.org/2000/svg" ',1)+'\n')
map_svg=(assets/'sampling_map.svg').read_text().replace(' xmlns="http://www.w3.org/2000/svg"','').strip()
languages={
'PT':{
'locale':'pt-BR','badge':'GUIA CIENTÍFICO DE BOLSO','cover':['Um guia de bolso para','sua história evolutiva'],'intro':'Genômica, geografia e mudanças climáticas no leste da Austrália.','number':'19.352','ind':'indivíduos','snp':'SNPs no painel principal','private':'PROVA PRIVADA • NÃO DISTRIBUIR',
'qtag':'A PERGUNTA','qhead':['Por que estudar','este lagarto?'],'q':'De Queensland à Tasmânia e ao sul da Austrália, a amostragem atravessou paisagens contrastantes. Perguntamos quanto da diferenciação genômica se associa à distância, ao terreno e ao clima atual e passado.','map':'Mapa esquemático: tons mais escuros indicam mais localidades em células de 3°. Células com menos de três não aparecem.',
'ctag':'CRÉDITOS','chead':['Sobre o estudo'],'credit':'Material do projeto Lampropholis. O repositório original credita Sally Potter pela contribuição ao pipeline. Autoria e afiliações serão conferidas com o manuscrito aprovado.','art':'A ilustração original é assinada por Nadaline. O uso na capa aguarda confirmação dos direitos sobre a imagem de base.','cite':'Referência provisória: equipe do estudo, Genomic history. Manuscrito em preparação.','qr':'QR AGUARDA PUBLICAÇÃO',
'stag':'AMOSTRAGEM','shead':['Da amostragem','à genômica'],'sampling':'O painel nuclear combina 170 indivíduos, 101 localidades e 19.352 SNPs alinhados à referência.','flow':['Amostragem geográfica','ddRADseq','SNPs nucleares','Estrutura populacional','Paisagem e demografia'],'bg':'As análises de paisagem usaram 19.131 SNPs após retirar 221 candidatos de uma RDA recuperada.',
'ktag':'ESTRUTURA','khead':['Uma espécie, diferentes','histórias genômicas'],'k':'Quatro grupos descrevem a estrutura ampla. K = 4 foi escolhido no cotovelo; o menor valor de cross-entropy ocorreu em K = 9.','klabel':'ESTADO MAIS FREQUENTE POR GRUPO','knote':'A fração se refere a indivíduos com Q ≥ 0,70. K3 também ocorre em Victoria e Tasmânia. Os grupos não definem espécies.',
'ltag':'PAISAGEM','lhead':['O que a paisagem','nos conta?'],'results':[('Geografia','A distância geográfica deu a maior contribuição ao modelo contemporâneo.'),('Topografia','A resistência topográfica trouxe um ganho pequeno, detectável no teste condicional.'),('História climática','O clima passado acrescentou um sinal modesto. O ganho variou entre as partições espaciais.')],'last':'A diferenciação no leste australiano não se resume a um fator ambiental. Associação não demonstra seleção ou causalidade.',
'faces':['Face externa','Face interna','Aberto → fechado']},
'EN':{
'locale':'en','badge':'SCIENTIFIC POCKET GUIDE','cover':['A pocket guide to its','evolutionary history'],'intro':'Genomics, geography and climatic change across eastern Australia.','number':'19,352','ind':'individuals','snp':'SNPs in the reference panel','private':'PRIVATE PROOF • NOT FOR DISTRIBUTION',
'qtag':'THE QUESTION','qhead':['Why study','this skink?'],'q':'From Queensland to Tasmania and South Australia, sampling spanned contrasting landscapes. We asked how genomic differentiation relates to distance, terrain and present and past climates.','map':'Schematic: darker 3° cells contain more sampled localities. Cells with fewer than three are not shown.',
'ctag':'CREDITS','chead':['About the study'],'credit':'Prepared for the Lampropholis study. The original repository credits Sally Potter for her pipeline contribution. Authorship and affiliations will be checked against the approved manuscript.','art':'Nadaline signed the original illustration. Cover use awaits confirmation of rights to the underlying image.','cite':'Provisional reference: study team, Genomic history. Manuscript in preparation.','qr':'QR PENDING PUBLIC RELEASE',
'stag':'SAMPLING','shead':['From sampling','to genomics'],'sampling':'The reference-aligned nuclear panel contains 19,352 SNPs from 170 individuals across 101 localities.','flow':['Geographic sampling','ddRADseq','Nuclear SNPs','Population structure','Landscape and demography'],'bg':'Landscape analyses used 19,131 SNPs after excluding 221 candidates from a recoverable RDA.',
'ktag':'STRUCTURE','khead':['One species, different','genomic histories'],'k':'Four genetic clusters summarize broad structure. The elbow supported K = 4, while numeric cross-entropy was lowest at K = 9.','klabel':'MOST SAMPLED STATE WITHIN EACH GROUP','knote':'Fractions are for individuals with Q ≥ 0.70. K3 also occurs in Victoria and Tasmania. The clusters do not define species.',
'ltag':'LANDSCAPE','lhead':['What can landscapes','tell us?'],'results':[('Geography','Geographic distance contributed most to the contemporary model.'),('Terrain','Topographic resistance gave a small predictive gain under the conditional test.'),('Climatic history','Past climates added modest predictive information. Gains varied among spatial folds.')],'last':'Genomic differentiation across eastern Australia cannot be reduced to one factor. Association does not establish selection or causation.',
'faces':['Outside face','Inside face','Open → folded']}}
def ticks(xs):
    return ''.join(line(x,2,x,4)+line(x,206,x,208) for x in xs)
def outside(c):
    f=rect(-3,-3,303,216,white)+rect(-3,-3,100,216,wash)+rect(197,-3,103,216,'#EFE8F7')+ticks((97,197))
    f+=top(8,2,c['qtag'],c['qhead'])+para(8,48,c['q'])[0]+f'<g transform="translate(14 85)">{map_svg}</g>'+para(8,165,c['map'],76,n=3.15,step=4.5)[0]
    f+=top(107,3,c['ctag'],c['chead'])
    y=48
    for k in ('credit','art','cite'):
        g,y=para(107,y,c[k],79);f+=g;y+=7
    f+=rect(107,151,79,15,wash)+text(111,160,c['qr'],3,b=True)
    f+=rect(197,-3,103,26,ink)+text(207,14,c['badge'],3.2,white,b=True)+line(207,38,234,38,amber,.95)
    f+=text(207,55,'Lampropholis',10.7,ink,serif,it=True)+text(207,67,'delicata',10.7,ink,serif,it=True)
    f+=''.join(text(207,88+i*7,h,5.75,ink,serif) for i,h in enumerate(c['cover']))
    f+=para(207,114,c['intro'],77)[0]+line(207,146,286,146,violet)
    f+=text(207,168,'170',9.8,ink,serif)+text(244,168,c['number'],9.8,ink,serif)
    f+=text(207,177,c['ind'],3.3)+text(244,177,c['snp'],3.3)+text(207,198,c['private'],2.9,violet,b=True)
    return f
def inside(c):
    f=rect(-3,-3,303,216,white)+rect(100,-3,100,216,wash)+rect(200,-3,100,216,'#FDF9EF')+ticks((100,200))
    f+=top(8,4,c['stag'],c['shead'])+para(8,48,c['sampling'])[0]
    for i,a in enumerate(c['flow']):
        y=76+18*i; f+=rect(8,y,80,11,wash)+rect(8,y,1.4,11,violet if i<4 else amber)+text(12,y+7,a,3.45,b=True)
        if i<4:f+=line(47,y+11,47,y+17,violet,.4)
    f+=para(8,177,c['bg'],80,3.15,4.65)[0]
    f+=top(108,5,c['ktag'],c['khead'])+para(108,48,c['k'],81)[0]
    f+=text(108,79,c['klabel'],3.05,violet,b=True)+line(108,85,187,85)
    for i,(state,a,b,col) in enumerate(zip(('South Australia','Queensland','New South Wales','New South Wales'),(41,11,35,36),(43,11,56,45),ks)):
        y=97+20*i;f+=f'<circle cx="110" cy="{y-1.2}" r="1.7" fill="{col}" stroke="{ink}" stroke-width=".25"/>'
        f+=text(114,y,f'K{i+1}',3.45,b=True)+text(126,y,state,3.42)+rect(126,y+5,36,1.45,'#DCD5E4')+rect(126,y+5,36*a/b,1.45,col)+text(165,y+7,f'{a}/{b}',3.05)
    f+=para(108,181,c['knote'],79,3.13,4.6)[0]
    f+=top(208,6,c['ltag'],c['lhead'])
    for j,(head,p) in enumerate(c['results']):
        y=49+40*j;f+=line(208,y-4,286,y-4,amber,.65)+text(208,y+3,head,4.55,ink,serif)+para(208,y+12,p,76,3.55,5.1)[0]
    f+=para(208,176,c['last'],76,3.35,4.9)[0]
    return f
def boxes(path,pad):
    # Append print boxes to pdfunite's classic xref without changing the page artwork.
    b=path.read_bytes()
    old=int(re.findall(rb'startxref\s+(\d+)\s+%%EOF',b)[-1]); trailer=b[old:]
    root=re.search(rb'/Root\s+(\d+\s+\d+\s+R)',trailer)[1].decode()
    size=int(re.search(rb'/Size\s+(\d+)',trailer)[1])
    pages=re.findall(rb'(?m)^(\d+)\s+0\s+obj\s*\r?\n(<<\s*/Type\s+/Page\b[^\r\n]+)\r?\nendobj',b)
    assert len(pages)==2 and b[old:old+4]==b'xref'
    cv=lambda l:' '.join(f'{v*72/25.4:.5f}' for v in l)
    t=cv((pad,pad,pad+297,pad+210)); g=cv((pad-3,pad-3,pad+300,pad+213))
    ref=[]; b+=b'\n'
    for n,p in sorted((int(n),p) for n,p in pages):
        p=re.sub(rb'/(?:TrimBox|BleedBox)\s*\[[^\]]*\]\s*',b'',p).rstrip()
        assert p.endswith(b'>>')
        ref.append((n,len(b)));b+=f'{n} 0 obj\n'.encode()+p[:-2]+f'/TrimBox [{t}]/BleedBox [{g}] >>\nendobj\n'.encode()
    end=len(b)
    b+=b'xref\n'+b''.join(f'{n} 1\n{x:010d} 00000 n \n'.encode() for n,x in ref)
    b+=f'trailer\n<< /Size {max(size,max(n for n,_ in ref)+1)} /Root {root} /Prev {old} >>\nstartxref\n{end}\n%%EOF\n'.encode()
    path.write_bytes(b)
def review(pdf,stem,labels,tmp):
    imgs=[]
    for i in (1,2):
        f=tmp/f'preview{i}'
        subprocess.run(['pdftoppm','-f',str(i),'-l',str(i),'-singlefile','-r','95','-png',str(pdf),str(f)],check=True)
        imgs.append(Image.open(f.with_suffix('.png')).convert('RGB'))
    w,h=imgs[0].size
    canvas=Image.new('RGB',(w+60,2*h+278),wash); d=ImageDraw.Draw(canvas); title=ImageFont.truetype(font(serif),28)
    for j,im in enumerate(imgs):
        y=51+j*(h+49); d.text((30,y-38),labels[j],fill=ink,font=title);canvas.paste(im,(30,y));d.rectangle((29,y-1,30+w,y+h),outline=rule)
    y=2*h+112;d.text((30,y),labels[2],fill=ink,font=title)
    b=imgs[0].copy();b.thumbnail((230,115));canvas.paste(b,(30,y+45))
    b=imgs[0].crop((int(w*197/297),0,w,h));b.thumbnail((56,115));canvas.paste(b,(240,y+45))
    canvas.save(stem.with_name(stem.name+'_preview.png'),optimize=True)
marks=''.join(f'<path d="{x}" fill="none" stroke="#251462" stroke-width=".22"/>' for x in ('M-8,0h4','M301,0h4','M-8,210h4','M301,210h4','M0,-8v4','M297,-8v4','M0,214v4','M297,214v4'))
for lang,c in languages.items():
    folder=base/c['locale'];folder.mkdir(exist_ok=True);stem=folder/f'Lampropholis_pocket_guide_{lang}'
    faces=(outside(c),inside(c))
    for key,g in zip(('outside','inside'),faces):stem.with_name(stem.name+f'_{key}.svg').write_text(svg(g))
    with tempfile.TemporaryDirectory() as td:
        td=Path(td)
        for key,pad in [('print',3),('print_cropmarks',9),('screen',0)]:
            pdfs=[]
            for i,g in enumerate(faces):
                f=td/f'{i}.svg'; pdf=f.with_suffix('.pdf')
                f.write_text(svg(g+(marks if pad==9 else ''),297+2*pad,210+2*pad,-pad,-pad))
                subprocess.run(['inkscape',str(f),'--export-type=pdf',f'--export-filename={pdf}'],check=True)
                pdfs.append(pdf)
            out=stem.with_name(stem.name+f'_{key}.pdf')
            subprocess.run(['pdfunite',*map(str,pdfs),str(out)],check=True)
            if pad: boxes(out,pad)
            else: review(out,stem,c['faces'],td)
badges=stage/'design/idv/badges';badges.mkdir(exist_ok=True)
for n,k,v,w in [('ddrad','Markers','ddRADseq',120),('geo','Region','Eastern Australia',160),('analysis','Code','R + Python',115),('status','Stage','Private review',145)]:
    a=rect(0,0,w,24,wash)+rect(0,0,50,24,ink)+text(5,16,k,10,white)+text(57,16,v,10)+line(51,21,w-3,21,amber)
    (badges/(n+'.svg')).write_text(f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="24" viewBox="0 0 {w} 24">{a}</svg>\n')
dest=stage/'assets';dest.mkdir(exist_ok=True)
m=rect(0,0,900,174,wash)+line(30,25,120,25,amber,3)+text(30,91,'Lampropholis delicata',53,ink,serif,it=True)+text(31,129,'Genomics, geography and history in eastern Australia',20)+text(31,152,'RESEARCH REPOSITORY / PRIVATE STAGING',10)+f'<g transform="translate(739 9) scale(2.05)">{map_svg}</g>'
(dest/'masthead.svg').write_text('<svg xmlns="http://www.w3.org/2000/svg" width="900" height="174" viewBox="0 0 900 174">'+m+'</svg>\n')

