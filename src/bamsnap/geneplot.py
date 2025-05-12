from PIL import ImageFont, Image, ImageDraw
from .conf import COLOR, GENE_ANNOT_FILE, REFER_SEQ_VERSION
from .util import getTemplatePath, getDataPath, comma, gzopen, decodeb, convert_int_list, getrgb, get_scale
import tabix

LINETYPE = []
LINETYPE.append('exon')

class GeneAnnot():
    def __init__(self, rec, header):
        self.data = {}
        for i in range(len(header)):
            self.data[header[i]] = rec[i]
        self.chrom = rec[header.index('CHROM')]
        self.spos = int(rec[header.index('SPOS')])
        self.epos = int(rec[header.index('EPOS')])
        self.strand = rec[header.index('strand')]
        self.is_negative = True
        if self.strand == "+":
            self.is_negative = False
        self.gene_id = rec[header.index('gene_id')]
        self.gene_name = rec[header.index('gene_name')]
        self.gene_biotype = rec[header.index('gene_biotype')]
        self.transcripts = []
        self.visible_transcripts = []
        self.load_transcript()

    def __str__(self):
        return self.gene_name + "("+ self.gene_id +")"

    def load_transcript(self):
        transcript_id_list = [self.data['gene_id']]
        transcript_biotype_list = self.data['gene_biotype']
        transcript_spos_list = self.data['SPOS']
        transcript_epos_list = self.data['EPOS']
        ltype_list = {}
        for ltype in LINETYPE:
            ltype_list['SPOS'] = self.data['SPOS']
            ltype_list['EPOS'] = self.data['EPOS']
        for i, tid in enumerate(transcript_id_list):
            t1 = TranscriptAnnot(tid, self.gene_biotype, self.spos, self.epos)
            for ltype in LINETYPE:
                t1.subregion['SPOS'] = ltype_list['SPOS']
                t1.subregion['EPOS'] = ltype_list['EPOS']
            t1.set_subregion()
            self.transcripts.append(t1)

    def set_visible_transcript(self, spos, epos):
        for t1 in self.transcripts:
            if t1.epos >= spos and t1.spos <= epos:
                self.visible_transcripts.append(t1)
        

class TranscriptAnnot():
    def __init__(self, tid, biotype, spos, epos):
        self.transcript_id = tid
        self.biotype = biotype
        self.subregion = {}
        self.spos = int(spos)
        self.epos = int(epos)

    def __str__(self):
        return self.transcript_id

    def set_subregion(self):
        for ltype in LINETYPE:
            self.subregion['SPOS'] = convert_int_list(self.subregion['SPOS'])
            self.subregion['EPOS'] = convert_int_list(self.subregion['EPOS'])


class GenePlot():
    def __init__(self, chrom, spos, epos, xscale, w, refversion="hg38", show_transcript = True):
        self.chrom = chrom
        self.nchrom = chrom.replace('chr', '')
        self.spos = spos
        self.epos = epos
        self.g_len = self.epos - self.spos + 1
        self.font = None
        self.gene_annot_file = getDataPath(GENE_ANNOT_FILE.replace("#REFSEQVERSION#", REFER_SEQ_VERSION[refversion]))
        self.gene_annot_tb = tabix.open(self.gene_annot_file)
        self.gene_annot_header = []
        self.gene_annot = []
        self.show_transcript = show_transcript
        self.w = w
        self.h = 0
        self.im = None
        self.bgcolor = "FFFFFF"
        self.noline = 0
        self.margin = 5
        self.lineheight = 10
        self.gene_pos_color = "ffac9c"
        self.gene_neg_color = "A19Cff"
        self.xscale = xscale
        self.load_gene_structure()

    def load_gene_structure(self):
        if len(self.gene_annot_header) == 0:
            for line in gzopen(self.gene_annot_file):
                line = decodeb(line)
                if line[0] == "#":
                    self.gene_annot_header = line[1:].split('\t')
                    self.gene_annot_header[-1] = self.gene_annot_header[-1].strip()
                    break
        pos_str = self.nchrom + ':' + str(self.spos) + '-' + str(self.epos)

        self.gene_annot = []
        for rec in self.gene_annot_tb.querys(pos_str):
            ga = GeneAnnot(rec, self.gene_annot_header)
            ga.set_visible_transcript(self.spos, self.epos)
            self.gene_annot.append(ga)
            self.noline += len(ga.visible_transcripts)
            # self.noline += len(ga.transccripts)
        

    def draw(self, dr):
        y = self.h - 1
        x1 = 0
        x2 = self.w

        yi = 0
        fontsize = self.font.getsize('C')
        for ga in self.gene_annot:
            # x1 = int((ga.spos - self.spos) * self.scale_x)
            # x2 = int((ga.epos - self.spos) * self.scale_x)
            # if x1 < 0:
            #     x1 = 0
            # if x2 > panel_xy[1][0]:
            #     x2 = panel_xy[1][0]
            # col1 = COLOR['GENE_NEG'] if ga.is_negative else COLOR['GENE_POS']
            # dr.line([(x1, yi), (x2, yi)], fill=col1, width=2)
            # x = int((min(ga.epos, self.epos) - max(ga.spos, self.spos)) / 2) * self.scale_x
            # dr.text(( x - 50 , yi-15), ga.gene_name, font=self.font, fill=COLOR['COORDINATE'])

            # for t1 in ga.transcripts:
            for t1 in ga.visible_transcripts:
                # yi += 30
                # yi += margin + fontsize[1] + lineheight
                yi += self.margin

                # x = int((min(ga.epos, self.epos) - max(ga.spos, self.spos)) / 2) * self.xscale.scale_x
                x = self.xscale.get_x( (min(t1.epos, self.epos) + max(t1.spos, self.spos))/2 )['cpos']
                
                txt = ga.gene_name + " ("  + t1.transcript_id + ")"
                x1 = min(max (x - int((len(txt) * fontsize[0])/2) , 0), self.w-len(txt) * fontsize[0])
                dr.text( (x1, yi), txt, font=self.font, fill=COLOR['COORDINATE'])

                yi += fontsize[1] + int(self.lineheight/2)
                col1 = getrgb(self.gene_neg_color, whitening=50) if ga.is_negative else getrgb(self.gene_pos_color, whitening=50)
                for i, s1 in enumerate(t1.subregion['SPOS']):
                    x1 = max(self.xscale.get_x(t1.subregion['SPOS'][i])['spos'], 0)
                    x2 = max(min(self.xscale.get_x(t1.subregion['EPOS'][i])['epos'], self.w), 0)
                    if x1 > 0 or x2 > 0:
                        dr.line([(x1, yi), (x2, yi)], fill=col1, width=self.lineheight)

                x1 = self.xscale.get_x(t1.spos)['spos']
                x2 = self.xscale.get_x(t1.epos)['epos']
                if x1 < 0:
                    x1 = 0
                if x2 > self.w:
                    x2 = self.w

                col1 = getrgb(self.gene_neg_color) if ga.is_negative else getrgb(self.gene_pos_color)
                dr.line([(x1, yi), (x2, yi)], fill=col1, width=2)
                
                xi = 10
                d = 3
                for i in range(100):
                    xi = i * 60 + 20
                    if xi > x2:
                        break
                    if xi >= x1:
                        if ga.is_negative:
                            dr.polygon([(xi, yi), (xi + d, yi + d), (xi + d, yi - d)], fill=col1)
                        else:
                            dr.polygon([(xi, yi), (xi - d, yi + d), (xi - d, yi - d)], fill=col1)
                        
                yi += int(self.lineheight/2)
    
    def get_image(self):
        if self.im is None:
            # if self.font is None:
            #     self.set_font()
            fontsize = self.font.getsize('C')
            self.h = self.noline * (self.margin + fontsize[1] + self.lineheight)

            self.im = Image.new('RGBA', (self.w, self.h), getrgb(self.bgcolor))
            dr = ImageDraw.Draw(self.im)
            self.draw(dr)
            
        return self.im        
            
