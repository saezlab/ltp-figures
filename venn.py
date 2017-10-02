#!/usr/bin/env/python

# Dénes Türei EMBL 2017
# turei.denes@gmail.com

import matplotlib as mpl
import matplotlib.backends.backend_pdf
import matplotlib_venn as venn
import collections

class Venn(object):
    
    def __init__(self,
                 infile,
                 outfile,
                 labels = None,
                 figsize = (7, 7),
                 colors = ('#6EA945', '#FCCC06', '#007B7F'),
                 **kwargs):
        
        
        self.infile = infile
        self.outfile = outfile
        self.figsize = figsize
        self.labels = labels
        self.colors = colors
        self.param = kwargs
        self.read_sets()
        self.plot_venn()
    
    def read_sets(self):
        
        self.sets = collections.defaultdict(lambda: set([]))
        
        with open(self.infile, 'r') as fp:
            
            _ = fp.readline()
            
            for l in fp:
                
                l = l.strip().split('\t')
                
                for sc in l[2]:
                    
                    self.sets[sc].add((l[0], l[1]))
        
        self.sets = dict(self.sets)

    def plot_venn(self):
        
        mpl.rc('font', family = 'DINPro', size = 48)
        self.pdf = mpl.backends.backend_pdf.PdfPages(self.outfile)
        self.fig = mpl.figure.Figure(figsize=self.figsize)
        self.cvs = mpl.backends.backend_pdf.FigureCanvasPdf(self.fig)
        self.ax  = self.fig.add_axes((0.05, 0.05, 0.90, 0.90))
        
        self.labels = (
            list(zip(sorted(self.sets.keys()), sorted(self.sets.keys())))
            if self.labels is None else self.labels
        )
        
        self.lsets = [self.sets[l[0]] for l in self.labels]
        self.method = 'venn2' if len(self.lsets) == 2 else 'venn3'
        
        self.venn = getattr(venn, self.method)(
            subsets = self.lsets,
            set_labels = [l[1] for l in self.labels],
            set_colors = self.colors,
            ax = self.ax,
            **self.param
        )
        
        # self.fig.tight_layout()
        self.cvs.draw()
        self.cvs.print_figure(self.pdf)
        self.pdf.close()
        self.fig.clf()


if __name__ == '__main__':
    
    # import venn
    v1 = venn.Venn('screens_venn.tsv', 'screens_venn.pdf',
                  labels = [('A', 'In vivo MS'), ('E', 'In vitro MS'), ('L', 'Literature')])
    v2 = venn.Venn('screens_venn.tsv', 'screens_hptlc_venn.pdf',
                  labels = [('A', 'In vivo MS'), ('T', 'In vivo HPTLC'), ('L', 'Literature')],
                  colors = ('#6EA945', '#DA0025', '#007B7F')
                  )
