""" Virana plot tool for plotting. Part of the Virana package.

    (c) 2013, Michael Zeidler, MPI for Informatics.
"""

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import colorsys
from numpy import arange

from matplotlib.patches import Rectangle

from matplotlib import rcParams
try:
    from plumbum import cli
except ImportError:
    message = 'This script requires the plumbum python package\n'
    sys.stderr.write(message)
    sys.exit(1)




class Family:
    
    def __init__(self,name,typ,samples=None):
        self.name=name
        if samples!=None:
            self.samples=samples
        else:
            self.samples={}
        self.typ=typ
    
    def add(self,sample,region,state,reads,bp,length):
        if sample not in self.samples:
            s = Sample(sample)
            self.samples[sample]=s

        return self.samples[sample].add(region,state,reads,bp,length)

        
    def __str__(self):
        return "{'type':'%s', 'samples':%s}"%(self.typ,str(self.samples))
    
    def __repr__(self):
        return "{'type':'%s', 'samples':%s}"%(self.typ,str(self.samples))
              
class Sample:

    def __init__(self,name,regions=None):
        self.name=name
        if regions!=None:
            self.regions=regions
        else:
            self.regions={}
    def add(self,region,state,reads,bp,length):
        if region not in self.regions:
            r = Region(region,state,reads,bp,length)
            self.regions[region]=r
            return True
        else:
            return False
    
    def __str__(self):
        return str(self.regions)
    
    def __repr__(self):
        return str(self.regions)            
                

class Region:
    
    def __init__(self,id,state,reads,bp,length):
        self.id=id
        self.state=state
        self.reads=reads
        self.bp=bp
        self.length=length
    
    def __str__(self):
        return "{'state':'%s', 'reads':%s, 'bp':%s, 'length':%s}" %(self.state,self.reads,self.bp,self.length)
    
    def __repr__(self):
        return "{'state':'%s', 'reads':%s, 'bp':%s, 'length':%s}" %(self.state,self.reads,self.bp,self.length)





class CLI(cli.Application):
    """Plot."""
    PROGNAME = "vref"
    VERSION = "1.0.0"
    DESCRIPTION = \
"""DESCRIPTION: Virana vplot - plot.

The Virana plot utility ('vplot') generates plots for the statistic file output of vhom. 

https://github.com/schelhorn/virana

Schelhorn S-E, Fischer M, Tolosi L, Altmueller J, Nuernberg P, et al. (2013)
Sensitive Detection of Viral Transcripts in Human Tumor Transcriptomes.
PLoS Comput Biol 9(10): e1003228. doi:10.1371/journal.pcbi.1003228"""

    USAGE = """USAGE: The program has one mode that can be accessed by
       [vplot | python vplot.py] plot
       """

    def main(self, *args):

        if args:
            print self.DESCRIPTION
            print
            print self.USAGE
            print("ERROR: Unknown command %r" % (args[0]))
            return 1

        if not self.nested_command:
            print self.DESCRIPTION
            print
            print self.USAGE
            print("ERROR : No command given")
            return 1

@CLI.subcommand("plot")
class Plotter(cli.Application):
    """ Generate plot for vhom statistic file."""

    stat_file = cli.SwitchAttr(['-s', '--stats'],str,
                                   mandatory=True,
                                   help="Path to statistic file generated with vhom.")

    plot_path  = cli.SwitchAttr(['-o', '--output_file'], str, mandatory=True,
                                 help="Sets the pdf output file. Note that all the plots are stored within a single pdf file.")

    debug       = cli.Flag(["-d", "--debug"], help="Enable debug messages")

    rcParams['figure.dpi'] = 100
    rcParams['axes.facecolor'] = 'white'
    rcParams['font.size'] = 10
    rcParams['patch.edgecolor'] = 'white'
    rcParams['patch.linewidth']=0.5
    rcParams['patch.facecolor'] = 'black'
    rcParams['font.family'] = 'StixGeneral'

    def means(self,start,stop):
    	list = arange(start,stop)
    	return (float(sum(list))/float(len(list)))  

    def remove_border(self,axes=None, top=False, right=False, left=True, bottom=True):
	    """
	    Minimize chartjunk by stripping out unnecesasry plot borders and axis ticks
	    
	    The top/right/left/bottom keywords toggle whether the corresponding plot border is drawn
	    """
	    ax = axes or plt.gca()
	    ax.spines['top'].set_visible(top)
	    ax.spines['right'].set_visible(right)
	    ax.spines['left'].set_visible(left)
	    ax.spines['bottom'].set_visible(bottom)
	    
	    #turn off all ticks
	    ax.yaxis.set_ticks_position('none')
	    ax.xaxis.set_ticks_position('none')
	    
	    #now re-enable visibles
	    if top:
	        ax.xaxis.tick_top()
	    if bottom:
	        ax.xaxis.tick_bottom()
	    if left:
	        ax.yaxis.tick_left()
	    if right:
	        ax.yaxis.tick_right()


    def customeLegend(self,color_map):
	    legend_map=[[],[]]
	    for sample in color_map:
	        legend_map[0].append(Rectangle((0,0),0,0,visible=False))
	        legend_map[1].append(sample)
	        for label  in color_map[sample]:
	            box=Rectangle((0,0),1,1,color=color_map[sample][label],fill=True)
	            legend_map[0].append(box)
	            legend_map[1].append(label)
	    plt.figlegend(legend_map[0],legend_map[1],loc=9,ncol=len(color_map),prop={'size':8})  

    def rgb_color_variants(self,stat_dict):
    
	    n_samples=0
	    sample_names=[]
	    
	    for fam,obj in stat_dict.iteritems():
	        if len(obj.samples)>n_samples:
	            n_samples=len(obj.samples)
	        for sample in obj.samples.iterkeys():
	            if sample not in sample_names:
	                sample_names.append(sample)

	    
	    assert n_samples == len(sample_names) #iterate over samples -> color_map
	    rgbs={}
	    i=1
	    for sample in sample_names:
	        rgbs[sample]={}
	        h=1.*(i/n_samples)
	        h=float(h)/2.
	        l=[0.1,0.4,0.7]
	        rgbs[sample]["pathogen"]=colorsys.hls_to_rgb(h,l[0],1.0)
	        rgbs[sample]["ambiguous"]=colorsys.hls_to_rgb(h,l[1],1.0)
	        rgbs[sample]["human"]=colorsys.hls_to_rgb(h,l[2],1.0)

	        i+=1
	    return rgbs

    def parseStatFile(self,filename):

		with open(filename,'r') as file:     
		    values=[ map(str,line.split('\t')) for line in file ]

		    family_dict={}
		    for line in values[1:]:
		        fam_type,fam,region,sample,state,reads,length,bp,coverage=line
		        if fam not in family_dict:
		            family= Family(fam,fam_type)
		            family_dict[fam]=family
		        family_dict[fam].add(sample,region,state,reads,bp,length)        
		    
		return family_dict

    def generatePyPlot(self,stat_dict,output_file):
	    pp = PdfPages(output_file)
	       
	    color_map=self.rgb_color_variants(stat_dict)
	    n_samples=len(color_map)
	    for plot in ["reads","bp"]:
	        fig=plt.figure()        
	        i=0    
	        for f,fam in stat_dict.iteritems():
	            j=0
	            for s,sample in fam.samples.iteritems():
	                values = []
	                color = []
	                for r,region in sample.regions.iteritems():
	                    values.append(int(getattr(region,plot)))
	                    color.append(color_map[s][region.state])  
	                if(len(values)>1):
	                	b = values[:-1]
	                	b.insert(0,0);
	                	for u in range(1,len(b)):
	                		b[u]+=b[u-1] 	
	                	plt.bar([i+j]*len(values),values,bottom=b,width=0.8,color=color)
	                else:
	                	plt.bar([i+j]*len(values),values,width=0.8,color=color)
	                
	                j+=1
	            i+= 1+n_samples
	        pos=[self.means(x,x+n_samples)+0.4 for x in arange(0,i+n_samples,1+n_samples)]    
	        plt.xticks(pos, stat_dict.keys(), rotation='vertical')
	        if(plot=="reads"):
	            plt.ylabel("Cumulative reads assigned to family")
	        if(plot=="bp"):
	            plt.ylabel("Cumulative basepairs assigned to family")
	        plt.xlabel("")
	        fig.subplots_adjust(bottom=0.25,wspace=0.5,top=0.8)
	        self.remove_border()
	        
	        self.customeLegend(color_map)   
	        pp.savefig(fig)
	    pp.close()

    def main(self):
        """plots"""

        if self.debug:
            logging.getLogger().setLevel(logging.DEBUG)


        stat_dict=self.parseStatFile(self.stat_file)
        self.generatePyPlot(stat_dict,self.plot_path)

if __name__ == "__main__":
    CLI.run()
