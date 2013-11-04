import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from matplotlib import rcParams
from matplotlib import colors
from matplotlib import cm

from matplotlib.patches import Rectangle

from vhom import RegionRunner

from argparse import ArgumentParser
from matplotlib.backends.backend_pdf import PdfPages

import os
import shutil

import sys

import colorsys

from numpy import arange

import copy




rcParams['figure.figsize'] = (6, 6)
rcParams['figure.dpi'] = 100
rcParams['axes.facecolor'] = 'white'
rcParams['font.size'] = 10
rcParams['patch.edgecolor'] = 'white'
rcParams['patch.linewidth']=0.5
rcParams['patch.facecolor'] = 'black'
rcParams['font.family'] = 'StixGeneral'

    
def rgb_color_variants(stat_dict):
    
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
def convert_color_map(color_map):
    deep = copy.deepcopy(color_map)
    for sample in deep:
        for state in deep[sample]:
            deep[sample][state] = colors.rgb2hex(deep[sample][state])
    return deep

def means(start,stop):
    list = arange(start,stop)
    return (float(sum(list))/float(len(list)))                
     

def remove_border(axes=None, top=False, right=False, left=True, bottom=True):
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
        
        
def generateHTML(stat_dict,file,directory):
    rval=["<html><head><title>Homologous groups and regions derived from hit files </title></head><body><p/>\n"]
    #rval.append('<link rel="stylesheet" type="text/css" href="/static/style/base.css">')
    rval.append("<a href='plot.pdf'>Download Plot</a>")
    rval.append("<div id='graph' class='graph'></div>")
    curr_dir = os.path.dirname(os.path.realpath(__file__))
    shutil.copy(os.path.join(curr_dir,"jquery-1.10.2.min.js"),directory)
    shutil.copy(os.path.join(curr_dir,"jqBarGraph.2.1.js"),directory)
    rval.append("<script type='text/javascript' src='jquery-1.10.2.min.js'></script>")
    rval.append("<script type='text/javascript' src='jqBarGraph.2.1.js'></script>")
    rval.append('<script>')	
    rval.append("$('#graph').jqBarGraph({")
    rval.append("data: %s,"%stat_dict)
    
    color_map=convert_color_map(rgb_color_variants(stat_dict))

    rval.append("colors: %s," %color_map)
    rval.append("legend: true,")
    rval.append("tab: 'reads',")
    rval.append("title: '<H3>Visualisation of Data</H3>',")
    rval.append("width: %d"%((len(stat_dict)*50)+150))
    rval.append("});")

    rval.append("</script>")
    
    rval.append( '</body></html>' )
    with open(file,'w') as file:
        file.write("\n".join( rval ))
                    
def customeLegend(color_map):
    legend_map=[[],[]]
    for sample in color_map:
        legend_map[0].append(Rectangle((0,0),0,0,visible=False))
        legend_map[1].append(sample)
        for label  in color_map[sample]:
            box=Rectangle((0,0),1,1,color=color_map[sample][label],fill=True)
            legend_map[0].append(box)
            legend_map[1].append(label)
    plt.figlegend(legend_map[0],legend_map[1],loc=9,ncol=len(color_map),prop={'size':8})        
    

    
    
def generatePyPlot(stat_dict,output_file):
    pp = PdfPages(output_file)
       
    color_map=rgb_color_variants(stat_dict)
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
        pos=[means(x,x+n_samples)+0.4 for x in arange(0,i+n_samples,1+n_samples)]    
        plt.xticks(pos, stat_dict.keys(), rotation='vertical')
        if(plot=="reads"):
            plt.ylabel("Cumulative reads assigned to family")
        if(plot=="bp"):
            plt.ylabel("Cumulative basepairs assigned to family")
        plt.xlabel("")
        fig.subplots_adjust(bottom=0.25,wspace=0.5,top=0.8)
        remove_border()
        
        customeLegend(color_map)   
        pp.savefig(fig)
    pp.close()
    


def parseStatFile(filename, directory):

    with open(filename,'r') as file:     
        values=[ map(str,line.split('\t')) for line in file ]

        family_dict={}
        for line in values[1:]:
            fam_type,fam,region,sample,state,reads,length,bp,coverage=line
            if fam not in family_dict:
                family= Family(fam,fam_type,directory=directory)
                family_dict[fam]=family
            family_dict[fam].add(sample,region,state,reads,bp,length)        
        
    return family_dict

class Family:
    
    def __init__(self,name,typ,directory=None,samples=None):
        self.name=name
        if samples!=None:
            self.samples=samples
        else:
            self.samples={}
        self.directory=directory  
        self.typ=typ
        if self.directory:
            self.files=self.__getFiles(directory)
        else:
            self.files=None
    
    def add(self,sample,region,state,reads,bp,length):
        if sample not in self.samples:
            s = Sample(sample)
            self.samples[sample]=s

        return self.samples[sample].add(region,state,reads,bp,self.files[region],length)

        
    def __getFiles(self,directory):
        dlist = os.listdir(os.path.join(directory,self.typ+"_"+self.name))
        dict={}
        for file in dlist:
            try:
                wc,region,suf = file.split("_")
            except ValueError:
                continue
            if region not in dict:
                dict[region]={}
            dict[region][suf]=file #bit redundant?! just save TRUE or FALSE?
        return dict
        
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
    def add(self,region,state,reads,bp,files,length):
        if region not in self.regions:
            r = Region(region,state,reads,bp,files,length)
            self.regions[region]=r
            return True
        else:
            return False
    
    def __str__(self):
        return str(self.regions)
    
    def __repr__(self):
        return str(self.regions)            
                

class Region:
    
    def __init__(self,id,state,reads,bp,files,length):
        self.id=id
        self.state=state
        self.reads=reads
        self.bp=bp
        self.files=files
        self.length=length
    
    def __str__(self):
        return "{'state':'%s', 'reads':%s, 'bp':%s, 'files':%s, 'length':%s}" %(self.state,self.reads,self.bp,self.files,self.length)
    
    def __repr__(self):
        return "{'state':'%s', 'reads':%s, 'bp':%s, 'files':%s, 'length':%s}" %(self.state,self.reads,self.bp,self.files,self.length)
      

if __name__ == "__main__":
    
    
    parser = ArgumentParser()
    
    a = parser.add_argument
    a("--html_file",dest="html_file")
    a("--directory",dest="directory")
    a("--stats",dest="stat_file")
      
    (options,args)= parser.parse_known_args()
    
    args.insert(0,"dummy")
    try:
        sys.exit(1)#RegionRunner.run(argv=args)
    except SystemExit:
        stat_dict=parseStatFile(options.stat_file,options.directory)

        generatePyPlot(stat_dict,os.path.join(options.directory,"plot.pdf"))
        generateHTML(stat_dict,options.html_file,options.directory)
    
        
                        
 
 