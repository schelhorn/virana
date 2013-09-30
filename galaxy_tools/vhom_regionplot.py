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




rcParams['figure.figsize'] = (6, 6)
rcParams['figure.dpi'] = 100
rcParams['axes.facecolor'] = 'white'
rcParams['font.size'] = 10
rcParams['patch.edgecolor'] = 'white'
rcParams['patch.linewidth']=0.5
rcParams['patch.facecolor'] = 'black'
rcParams['font.family'] = 'StixGeneral'


def color_variants(r,g,b):
    h,l2,s= colorsys.rgb_to_hls(r,g,b)
    l=[l2-0.4,l2-0.1,l2+0.2]
    hexs=[]
    for i in l:
        hexs.append(colors.rgb2hex(colorsys.hls_to_rgb(h,i,s)))
    return hexs
    
def rgb_color_variants(r,g,b):
    h,l2,s= colorsys.rgb_to_hls(r,g,b)
    l=[l2-0.4,l2-0.1,l2+0.2]
    rgbs={}
    rgbs["pathogen"]=colorsys.hls_to_rgb(h,l[0],s)
    rgbs["ambiguous"]=colorsys.hls_to_rgb(h,l[1],s)
    rgbs["human"]=colorsys.hls_to_rgb(h,l[2],s)
    return rgbs

                
     

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
    rval.append('<link rel="stylesheet" type="text/css" href="/static/style/base.css">')
    rval.append("<a href='plot.pdf'>Download Plot</a>")
    #sven create <div> here for some text
    n_samples=0
    for sample in stat_dict.iterkeys():
        n_samples+=1
        rval.append("<div id=%s_graph></div><p/>"%(sample))
    curr_dir = os.path.dirname(os.path.realpath(__file__))
    shutil.copy(os.path.join(curr_dir,"jquery-1.10.2.min.js"),directory)
    shutil.copy(os.path.join(curr_dir,"jqBarGraph.2.1.js"),directory)
    rval.append("<script type='text/javascript' src='jquery-1.10.2.min.js'></script>")
    rval.append("<script type='text/javascript' src='jqBarGraph.2.1.js'></script>")
    rval.append('<script>')
    i=0
    for sample in stat_dict.iterkeys():
        i+=1
        rval.append("%s_array = new Array(" %sample)
        for family in stat_dict[sample].iterkeys():
            values=[]
            for region,data in stat_dict[sample][family].iteritems():
                  values.append([int(data[0][1]),int(data[0][2]), data[0][0], data[0][3],int(data[0][4]),int(region)])
            rval.append(str([values,family])+",")
        rval[-1]=rval[-1][:-1]
        rval.append(");")
        rval.append("%s_files=%s;"%(sample,str(getFiles(directory))))	
        rval.append("$('#%s_graph').jqBarGraph({"%sample)
        rval.append("sample: '%s',"%sample)
        rval.append("data: %s_array,"%sample)
        rval.append("files: %s_files,"%sample)
        r,g,b,a=cm.hsv(1.*(1-(i/n_samples)))
        rval.append("colors: %s," %color_variants(r,g,b))
        rval.append("legend: true,")
        rval.append("tab: 'reads',")
        rval.append("title: '<H3>Visualisation of Sample %s</H3>',"%sample)
        rval.append("width: %d"%((len(stat_dict[sample])*50)+150))
        rval.append("});")

    rval.append("</script>")
    
    rval.append( '</body></html>' )
    with open(file,'w') as file:
        file.write("\n".join( rval ))
def transformDict(stat_dict):
    dict={}
    position={}
    dict_label={}
    for sample in stat_dict.iterkeys():
        dict[sample]=[{},{}]
        position[sample]={}
        dict_label[sample]={}
        j=0
        for family in stat_dict[sample].iterkeys():
         
            i=0
            for region,values in stat_dict[sample][family].iteritems():
                
                if i not in dict[sample][0]:
                    dict[sample][0][i]=[]
                    dict[sample][1][i]=[]
                    position[sample][i]=[]
                    dict_label[sample][i]=[]
                dict[sample][0][i].append(int(values[0][1]))
                dict[sample][1][i].append(int(values[0][2]))
                position[sample][i].append(j)
                dict_label[sample][i].append(values[0][0])
                i+=1
            j+=1
    return dict,position,dict_label
                    
def customeLegend(color_map):
    legend_map=[[],[]]
    for label  in color_map:
        box=Rectangle((0,0),1,1,color=color_map[label],fill=True)
        legend_map[0].append(box)
        legend_map[1].append(label)
    plt.figlegend(legend_map[0],legend_map[1],loc=8,ncol=len(color_map))        
    

    
    
def generatePyPlot(stat_dict,output_file):
    pp = PdfPages(output_file)
    fig = plt.figure()
    
    n_samples=len(stat_dict)
    dict,position,dict_label=transformDict(stat_dict)
                
    i=0    
    for sample in dict.iterkeys():
        fig = plt.figure()
        
        fig.subplots_adjust(bottom=0.25,hspace=0.5)
        
        r,g,b,a=cm.hsv(1.*(1-(i/n_samples)))
        color_map=rgb_color_variants(r,g,b)
        
        for sub in [0,1]:
          plt.subplot(1,2,sub+1)
          bottom={}
          for level,values in dict[sample][sub].iteritems():
              colors=[color_map[r] for r in dict_label[sample][level]]

              if level==0:

                  plt.bar(position[sample][level],values,bottom=0, width=0.8 ,color=colors)
              else:
                  lastvalues=[]
                  for oldpos in range(len(values)):
                      lastvalues.append(bottom[position[sample][level][oldpos]])         

                  plt.bar(position[sample][level],values, width=0.8,bottom=lastvalues ,color=colors)
              for pos in range(len(values)):
                  if position[sample][level][pos] not in bottom:
                        bottom[position[sample][level][pos]]=values[pos]
                  else:
                        bottom[position[sample][level][pos]]+=values[pos]
          pos=[x+0.4 for x in range(len(stat_dict[sample]))]    
          plt.xticks(pos, stat_dict[sample].keys(), rotation='vertical')
          if(sub==0):
            plt.ylabel("Cumulative reads assigned to family")
          else:
            plt.ylabel("Cumulative basepairs assigned to family")
          plt.xlabel("")
          fig.subplots_adjust(bottom=0.25,wspace=0.5)
          remove_border()
        
        customeLegend(color_map)   
        plt.suptitle(sample)
        pp.savefig(fig)
        i+=1
    pp.close()
    


def parseStatFile(filename):

    with open(filename,'r') as file:     
        values=[ map(str,line.split('\t')) for line in file ]

        dict={}
        for line in values[1:]:
            if line[3] not in dict:
                dict[line[3]]={}
                dict[line[3]][line[1]]={}
                dict[line[3]][line[1]][line[2]]=[[line[4],line[5],line[7],line[0],line[6]]]
            else:
                if line[1] not in dict[line[3]]:
                    dict[line[3]][line[1]]={}
                    dict[line[3]][line[1]][line[2]]=[[line[4],line[5],line[7],line[0],line[6]]]
                else:
                    if line[2] not in dict[line[3]][line[1]]:
                        dict[line[3]][line[1]][line[2]]=[[line[4],line[5],line[7],line[0],line[6]]]
                    else:
                         dict[line[3]][line[1]][line[2]].append([line[4],line[5],line[7],line[0],line[6]])        
#        for key in dict.iterkeys():
#            for family_key,values in dict[key].iteritems():
#                for line in values:
#                    print key,family_key,line
        return dict
def getFiles(directory):
    rval={}
    dlist = os.listdir(directory)
    for dire in dlist:
        if(os.path.isdir(os.path.join(directory,dire))):
             rval[dire]={}
             flist = os.listdir(os.path.join(directory,dire))
             for file in flist:
                split=file.split("_")
                region=split[1]
                if region not in rval[dire]:
                    rval[dire][region]=[file]
                else:
                    rval[dire][region].append(file)
                         
    return rval

        

if __name__ == "__main__":
    
    
    parser = ArgumentParser()
    
    a = parser.add_argument
    a("--html_file",dest="html_file")
    a("--directory",dest="directory")
    a("--stats",dest="stat_file")
      
    (options,args)= parser.parse_known_args()
    
    args.insert(0,"dummy")
    try:
        RegionRunner.run(argv=args)
    except SystemExit:
        stat_dict=parseStatFile(options.stat_file)
        generatePyPlot(stat_dict,os.path.join(options.directory,"plot.pdf"))
        generateHTML(stat_dict,options.html_file,options.directory)
    
        
                        
 
 