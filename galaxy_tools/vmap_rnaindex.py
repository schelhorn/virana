from argparse import ArgumentParser
import os
from vmap import RNAIndex

if __name__ == "__main__":
    
    
    parser = ArgumentParser()
    
    a = parser.add_argument
    a("-o","--html_file",dest="html_file")
    a("-d","--dir",dest="directory")
      
    (options,args)= parser.parse_known_args()
    
    args.insert(0,"dummy")
    try:
        RNAIndex.run(argv=args)
    except SystemExit:
        f = open(options.html_file,'w')
        rval = ["<html><head><title>STAR Index Galaxy Composite Dataset </title></head><body><p/>\n"]
        rval.append('<div>This composite dataset is composed of the following files:<p/><ul>')
        flist = os.listdir(options.directory)
        for file in flist:
           rval.append( '<li><a href="%s">%s</a></li>' % ( file, file) )
        rval.append( '</ul></body></html>' )
    
        f.write("\n".join( rval ))
        f.close() 
