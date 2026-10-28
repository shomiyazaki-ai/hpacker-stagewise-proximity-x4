from pathlib import Path
import json,re,textwrap
from xml.sax.saxutils import escape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,Image,KeepTogether,PageBreak,Preformatted
R=Path(__file__).resolve().parents[1];D=R/'data/authoritative';T=json.loads((D/'tables.json').read_text())
fontdir=Path('/usr/share/fonts/truetype/dejavu')
for name,file in [('DV','DejaVuSans.ttf'),('DV-B','DejaVuSans-Bold.ttf'),('DVM','DejaVuSansMono.ttf')]:pdfmetrics.registerFont(TTFont(name,str(fontdir/file)))
pdfmetrics.registerFontFamily('DV',normal='DV',bold='DV-B',italic='DV',boldItalic='DV-B')
W=A4[0]-104
styles={
 'body':ParagraphStyle('body',fontName='DV',fontSize=9.2,leading=13.15,spaceAfter=6.2,splitLongWords=True),
 'title':ParagraphStyle('title',fontName='DV-B',fontSize=20,leading=25,spaceAfter=12,textColor=colors.HexColor('#153b54')),
 'h1':ParagraphStyle('h1',fontName='DV-B',fontSize=12.3,leading=16,spaceBefore=11,spaceAfter=6,keepWithNext=True,textColor=colors.HexColor('#153b54')),
 'h2':ParagraphStyle('h2',fontName='DV-B',fontSize=10,leading=13.6,spaceBefore=8,spaceAfter=5,keepWithNext=True),
 'cell':ParagraphStyle('cell',fontName='DV',fontSize=7.8,leading=10.2,splitLongWords=True),
 'head':ParagraphStyle('head',fontName='DV-B',fontSize=7.8,leading=10.2,textColor=colors.white,splitLongWords=True),
 'caption':ParagraphStyle('caption',fontName='DV',fontSize=8.1,leading=11.1,spaceBefore=5,spaceAfter=7),
 'tabletitle':ParagraphStyle('tabletitle',fontName='DV-B',fontSize=8.6,leading=11.5,spaceBefore=8,spaceAfter=5,keepWithNext=True),
 'foot':ParagraphStyle('foot',fontName='DV',fontSize=7.6,leading=10.4,spaceBefore=4,spaceAfter=8),
 'code':ParagraphStyle('code',fontName='DVM',fontSize=7.3,leading=10,spaceBefore=4,spaceAfter=8,backColor=colors.HexColor('#f1f5f7'),borderPadding=7)
}
widths={'Table_1':[110,300,81],'Table_2':[285,95,111],'Table_3':[43,43,124,127,154],'Table_4':[50,145,95,105,96],'Table_5':[48,101,105,146,91],'Table_S1':[140,166,185],'Table_S2A_1':[50,50,78,66,73,78,96],'Table_S2A_2':[45,70,93,283],'Table_S2B':[38,62,75,80,72,88,76],'Table_S3':[95,53,46,54,55,97,91],'Table_S4A':[116,125,125,125],'Table_S4B':[50,141,91,110,99],'Table_S5':[90,130,160,111]}
def txt(s):return escape(str(s)).replace('\n','<br/>')
def footer(c,doc):
 c.saveState();c.setStrokeColor(colors.HexColor('#d3dde3'));c.line(52,39,A4[0]-52,39);c.setFont('DV',7);c.setFillColor(colors.HexColor('#526572'));c.drawString(52,27,'X4 | 1.1.0 | 28 Oct 2026');c.drawRightString(A4[0]-52,27,str(doc.page));c.restoreState()
def make(kind,out):
 elems=json.loads((R/'source'/(kind+'.json')).read_text());story=[]
 figure_heading=[]
 for index,e in enumerate(elems):
  typ=e['type']
  if typ=='title':story.append(Paragraph(txt(e['text']),styles['title']))
  elif typ=='heading':
   hp=Paragraph(txt(e['text']),styles['h1' if e['level']==1 else 'h2'])
   if index+1<len(elems) and elems[index+1]['type']=='figure':figure_heading=[hp]
   else:story.append(hp)
  elif typ=='paragraph':story.append(Paragraph(txt(e['text']),styles['body']))
  elif typ=='pagebreak':story.append(PageBreak())
  elif typ=='code':
   ls=[]
   for l in e['text'].splitlines():ls.extend(textwrap.wrap(l,108,break_long_words=True,break_on_hyphens=False) or [''])
   story.append(Preformatted('\n'.join(ls),styles['code']))
  elif typ=='figure':
   from PIL import Image as PI
   f=R/'figures'/e['file'];ww,hh=PI.open(f).size;im=Image(str(f),width=W,height=W*hh/ww)
   story.append(KeepTogether(figure_heading+[Spacer(1,4),im,Paragraph(txt(e['caption']),styles['caption'])]));figure_heading=[]
  elif typ=='table':
   d=T[e['id']];cw=[x*W/sum(widths[e['id']]) for x in widths[e['id']]]
   rows=[[Paragraph(txt(h),styles['head']) for h in d['headers']]]+[[Paragraph(txt(v),styles['cell']) for v in row] for row in d['rows']]
   t=Table(rows,colWidths=cw,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#244e65')),('VALIGN',(0,0),(-1,-1),'TOP'),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('ROWBACKGROUNDS',(0,1),(-1,-1),[colors.HexColor('#f0f4f6'),colors.white]),('LINEBELOW',(0,-1),(-1,-1),.5,colors.HexColor('#b4c4ce'))]))
   block=[Paragraph(txt(d['caption']),styles['tabletitle']),t,Paragraph(txt(d['footnote']),styles['foot'])]
   story.append(KeepTogether(block))
 doc=SimpleDocTemplate(str(out),pagesize=A4,leftMargin=52,rightMargin=52,topMargin=45,bottomMargin=51,title=kind.replace('_',' ')+' - X4',author='Sho Miyazaki',pageCompression=1)
 doc.build(story,onFirstPage=footer,onLaterPages=footer)
make('Technical_Note',R/'manuscript/Technical_Note_v1.1.0.pdf')
make('Supplementary_Information',R/'supplement/Supplementary_Information_v1.1.0.pdf')
print('Built two final-candidate PDFs')
