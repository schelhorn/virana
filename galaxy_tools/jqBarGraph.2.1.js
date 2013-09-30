/**
 * jqBarGraph - jQuery plugin
 * @version: 1.1 (2011/04/03)
 * @requires jQuery v1.2.2 or later 
 * @author Ivan Lazarevic
 * Examples and documentation at: http://www.workshop.rs/jqbargraph/
 * 
 * Dual licensed under the MIT and GPL licenses:
 *   http://www.opensource.org/licenses/mit-license.php
 *   http://www.gnu.org/licenses/gpl.html
 * 
 * @param data: arrayOfData // array of data for your graph
 * @param title: false // title of your graph, accept HTML
 * @param barSpace: 10 // this is default space between bars in pixels
 * @param width: 400 // default width of your graph ghhjgjhg
 * @param height: 200 //default height of your graph
 * @param color: '#000000' // if you don't send colors for your data this will be default bars color
 * @param colors: false // array of colors that will be used for your bars and legends
 * @param lbl: '' // if there is no label in your array
 * @param sort: false // sort your data before displaying graph, you can sort as 'asc' or 'desc'
 * @param position: 'bottom' // position of your bars, can be 'bottom' or 'top'. 'top' doesn't work for multi type
 * @param prefix: '' // text that will be shown before every label
 * @param postfix: '' // text that will be shown after every label
 * @param animate: true // if you don't need animated appereance change to false
 * @param speed: 2 // speed of animation in seconds
 * @param legendWidth: 100 // width of your legend box
 * @param legend: false // if you want legend change to true
 * @param legends: false // array for legend. for simple graph type legend will be extracted from labels if you don't set this
 * @param type: false // for multi array data default graph type is stacked, you can change to 'multi' for multi bar type
 * @param showValues: true // you can use this for multi and stacked type and it will show values of every bar part
 * @param showValuesColor: '#fff' // color of font for values 

 * @example  $('#divForGraph').jqBarGraph({ data: arrayOfData });  
  
**/

(function($) {
	var opts = new Array;
	var level = new Array;
	
	$.fn.jqBarGraph = $.fn.jqbargraph = function(options){
	
	init = function(el){

		opts[el.id] = $.extend({}, $.fn.jqBarGraph.defaults, options);
		$(el).css({ 'width': opts[el.id].width, 'height': opts[el.id].height, 'position':'relative', 'text-align':'center' });
		

		doGraph(el);

	};
	
	// sum of array elements
	sum = function(ar,index){
		total = 0;
		for(val in ar){
			total += ar[val][index];
		}
		return total.toFixed(2);
	};
	
	// count max value of array
	max = function(ar,index){
		maxvalue = 0;
		for(var val in ar){
			value = ar[val][0];
			if(value instanceof Array) value = sum(value,index);	
			if (parseFloat(value) > parseFloat(maxvalue)) maxvalue=value;
		}	
		return maxvalue;	
	};

	// max value of multi array
	maxMulti = function(ar,index){
		maxvalue = 0;
		maxvalue2 = 0;
		
		for(var val in ar){
			ar2 = ar[val][0];
			
			for(var val2 in ar2){
				if(ar2[val2][index]>maxvalue2) maxvalue2 = ar2[val2][index];
			}

			if (maxvalue2>maxvalue) maxvalue=maxvalue2;
		}
	
		return maxvalue;		
	};
	
		
	doGraph = function(el){
		
		arr = opts[el.id];
		data = arr.data;
		if(arr.tab=='reads'){
		  val_index=0;
		}
		else{
		  val_index=1;
		}

		//check if array is bad or empty
		if(data == undefined) {
			$(el).html('There is not enought data for graph');
			return;
		}
		//sorting ascending or descending
		if(arr.tab == 'reads'){
		  if(arr.sort == 'asc'){
		      data.sort(sortReadsAsc);    
		  }
		  if(arr.sort == 'desc'){
		      data.sort(sortReadsDesc);        
		  }
		}
		if(arr.tab == 'basepairs'){
		  if(arr.sort == 'asc'){
		      data.sort(sortBasesAsc);
		   }
		   if(arr.sort == 'desc'){
		      data.sort(sortBasesDesc);
		   }
		}
		if(arr.sortBar == 'asc') sortBars(data, barAsc);
		if(arr.sortBar == 'desc')sortBars(data, barDesc);
    
		
		
		legend = '';
		prefix = arr.prefix;
		postfix = arr.postfix;
		space = arr.barSpace; //space between bars
		legendWidth = arr.legend ? arr.legendWidth : 0; //width of legend box
		fieldWidth = ($(el).width()-legendWidth)/data.length; //width of bar
		totalHeight =  $(el).height(); //total height of graph box
		var leg = new Array(); //legends array
		
		//max value in data, I use this to calculate height of bar
		max = max(data,val_index);
	
		color_map={};
		color_map["pathogen"]=arr.colors[0];
		color_map["ambiguous"]=arr.colors[1];
		color_map["human"]=arr.colors[2];
		
		if(arr.tab=='reads'){
		  val_index=0;
		}
		else{
		  val_index=1;

		}

 		for(var val in data){
 			
 			valueData = data[val][0];
 			if (valueData instanceof Array) 
 				value = sum(valueData,val_index);
 			else
 				value = valueData;
 			
 			lbl = data[val][1];
			unique = val+el.id; //unique identifier
			divid=""+el.id;
 			
 			if(arr.type == 'multi') color = 'none';
 				
 			if (lbl == undefined) lbl = arr.lbl;
 			
 			margin_top=14/(data.length);

 		
 			out  = "<div class='graphField"+el.id+"' id='graphField"+unique+"' style='position: absolute'>";
 			out += "<div class='graphValue"+el.id+"' id='graphValue"+unique+"'>"+parseInt(value)+"</div>";
 			
 			out += "<div class='graphBar"+el.id+"' id='graphFieldBar"+unique+"' style='background-color:#fff;position: relative; overflow: hidden;'></div>";
            //console.log(color)
			// if there is no legend or exist legends display lbl at the bottom
 			
 				out += "<a class='graphLabelLink' href=#files"+el.id+" onclick=fillDiv('"+el.id+"','"+valueData[0][3]+"_"+lbl+"')><div class='graphLabel"+el.id+"' id='graphLabel"+unique+"'style='margin-top:"+margin_top+"px;'>"+lbl+"</div></a>";
 			out += "</div>";
 			
 			
 			 			
			$(el).append(out);
			$(".graphLabel"+el.id).css({'-webkit-transform': 'rotate(30deg)', '-moz-transform': 'rotate(30deg)','-o-transform': 'rotate(30deg)', '-ms-transform': 'rotate(30deg)','transform': 'rotate(30deg)', 'height':'100'});
 			$('a.graphLabel').css({
 			    'text-decoration': 'none',
 			    'color':'black'    
 			});
 			
 			//$('#graphLabel'+el.id).rotateLeft();
 			
 			//size of bar
 			totalHeightBar = totalHeight - $('.graphLabel'+el.id).height() - $('.graphValue'+el.id).height()-margin_top; 
 			fieldHeight = (totalHeightBar*value)/max;	
 			$('#graphField'+unique).css({ 
 				'left': (fieldWidth)*val, 
 				'width': fieldWidth-space, 
 				'margin-left': space});
 	
 			// multi array
 			if(valueData instanceof Array){
 				
				if(arr.type=="multi"){

					maxe = maxMulti(data,val_index);
					totalHeightBar = fieldHeight = totalHeight - $('.graphLabel'+el.id).height()-margin_top;
					$('.graphValue'+el.id).remove();
				} else {

					maxe = max;
				}
				
 				for (i in valueData){
 					heig = (totalHeightBar*valueData[i][val_index]/maxe);
                    if(navigator.userAgent.match('Safari') && !navigator.userAgent.match('Chrom')){
                        
                        heig = heig + (arr.borderSize/1.85);
                    } 
                    
 					wid = parseInt((fieldWidth-space)/valueData.length);

 					sv = ''; // show values
 					fs = 0; // font size
 					if (arr.showValues){
 						sv = arr.prefix+valueData[i][0]+arr.postfix;
 						fs = 12; // font-size is 0 if showValues = false
 					}
 					o = "<a class='subBarLink' href=#files"+el.id+" onclick=fillDiv('"+el.id+"','"+valueData[i][3]+"_"+lbl+"','"+valueData[i][5]+"')><div class='subBars"+el.id+"' style=' box-sizing:border-box; -moz-box-sizing:border-box; -ms-box-sizing:border-box; -webkit-box-sizing:border-box; height:"+heig+"px;  border-top:"+arr.borderSize+"px solid; border-color: "+arr.showValuesColor+"; background-color: "+color_map[valueData[i][2]]+"; left:"+wid*i+"px; color:"+arr.showValuesColor+"; font-size:"+fs+"px' >"+sv+"</div></a>";
 					$('#graphFieldBar'+unique).prepend(o);
 				}
 			}
 			
 			if(arr.type=='multi')
 				$('.subBars'+el.id).css({ 'width': wid, 'position': 'absolute', 'bottom': 0 });
 
 			//position of bars
 			if(arr.position == 'bottom') $('.graphField'+el.id).css('bottom',0);
            
            
			
 			
 			
 			// animated apearing
 			if(arr.animate){
 				$('#graphFieldBar'+unique).css({ 'height' : 0 });
 				$('#graphFieldBar'+unique).animate({'height': fieldHeight},arr.speed*1000);
 			} else {
 				$('#graphFieldBar'+unique).css({'height': fieldHeight});
 			}
 			
 		}
 			
 		
 		
 		createLegend(color_map,el.id); // create legend from array
 		createLinks(el.id);
 		//position of legend
 		
			$(el).append("<div id='legendHolder"+unique+"'></div>");
	 		$('#legendHolder'+unique).css({ 'width': legendWidth, 'float': 'right', 'margin-left':'100px','text-align' : 'left'});
	 		$('#legendHolder'+unique).append(legend);
	 		$('#legendHolder'+unique).append(links);
	 		$('.legendBar'+el.id).css({ 'float':'left', 'margin': 3, 'margin-left':'10','height': 12, 'width': 20, 'font-size': 0});
	 		$('.linkBar'+el.id).css({'margin-left':'10'});
	 		
	 		
	 		$("#sortAsc"+el.id).click(function(){
	           if(opts[el.id].sort!='asc'){
	               opts[el.id].sort='asc';
	               $('#graphHolder'+el.id).html('');
	               $('#graphHolder'+el.id).jqbargraph(opts[el.id]);
	           }
    
	        });
	        $("#sortDesc"+el.id).click(function(){
	           if(opts[el.id].sort!='desc'){
	               opts[el.id].sort='desc';
	               $('#graphHolder'+el.id).html('');
	               $('#graphHolder'+el.id).jqbargraph(opts[el.id]);
	           }
    
	        });
	        $("#sortBarAsc"+el.id).click(function(){
	           if(opts[el.id].sortBar!='asc'){
	               opts[el.id].sortBar='asc';
	               $('#graphHolder'+el.id).html('');
	               $('#graphHolder'+el.id).jqbargraph(opts[el.id]);
	           }
    
	        });
	 		$("#sortBarDesc"+el.id).click(function(){
	           if(opts[el.id].sortBar!='desc'){
	               opts[el.id].sortBar='desc';
	               $('#graphHolder'+el.id).html('');
	               $('#graphHolder'+el.id).jqbargraph(opts[el.id]);
	           }
    
	        });
	        
	        $("#showBasepairs"+el.id).click(function(){
	               opts[el.id].tab='basepairs';
	               $('#graphHolder'+el.id).html('');
	               $('#graphHolder'+el.id).jqbargraph(opts[el.id]);
	               
	      
    
	        });
	        $("#showReads"+el.id).click(function(){
	               opts[el.id].tab='reads';
	               $('#graphHolder'+el.id).html('');
	               $('#graphHolder'+el.id).jqbargraph(opts[el.id]);
	               
	           
    
	        });
 		    
 		    

 		//position of title
 		if(arr.title){
 			$(el).wrap("<div id='graphHolder"+el.id+"'></div>");
 			if(arr.tab=='reads'){
 			    $('#graphHolder'+el.id).prepend("<div id='label"+el.id+"'>Cumulative reads assigned to family</div>").css({ 'width' : arr.width+'px', 'text-align' : 'center' });
 			}else{
   			    $('#graphHolder'+el.id).prepend("<div id='label"+el.id+"'>Cumulative basepairs assigned to family</div>").css({ 'width' : arr.width+'px', 'text-align' : 'center' });
 			}
 			$('#graphHolder'+el.id).prepend("<a href=#files"+el.id+" onclick=fillDiv('"+el.id+"')>"+arr.title+"</a>").css({ 'width' : arr.width+'px', 'text-align' : 'center' });
 		    
  		}
 		    $("#graphHolder"+el.id).append("<div class='files"+el.id+"' id='files"+el.id+"' ></div><p/>");
 		    $('.files'+el.id).css({'width' : arr.width+'px','background-color':'silver', 'border':'2px solid gray', 'display':'none'})
 		    
 		     $("#graphHolder"+el.id).append("<div class='image"+el.id+"' id='image"+el.id+"' ></div>");
 		    $('.image'+el.id).css({'width' : arr.width+'px','background-color':'silver', 'border':'2px solid gray', 'display':'none'})
 		
	};


	//creating legend from array
	createLegend = function(color_map,id){
		legend = '';
		for(var val in color_map){
	 			legend += "<div id='legend"+val+"' style='overflow: hidden; zoom: 1;'>";
	 			legend += "<div class='legendBar"+id+"' id='legendColor"+val+"' style='background-color:"+color_map[val]+"'></div>";
	 			legend += "<div class='legendLabel"+id+"' id='graphLabel"+unique+"'>"+val+"</div>";
	 			legend += "</div>";			
		}
	};
	
	createLinks = function(id){
	   links="<div class='linkBar"+id+"'>";
	   links+="<div ><a href='javascript:void(0);' id='sortAsc"+id+"' >sort  ascending</a></div>"
	   links+="<div ><a href='javascript:void(0);' id='sortDesc"+id+"'>sort  descending</a></div>"
	   links+="<div ><a href='javascript:void(0);' id='sortBarAsc"+id+"'>sort bars ascending</a></div>"
	   links+="<div ><a href='javascript:void(0);' id='sortBarDesc"+id+"'>sort bars descending</a></div>"
	   if(opts[id].tab=='reads'){
	       links+="<div ><a href='javascript:void(0);' id='showBasepairs"+id+"'>show basepair chart </a></div>"    
	   }else{
	       links+="<div ><a href='javascript:void(0);' id='showReads"+id+"'>show read chart </a></div>"    
	   }
	   links+="</div>"
	   
	     
	};
	
	fillDiv = function (elid, dir, region) {
	   
	        
	   
	        arr=opts[elid];
	        
	        dict=arr.files;
	        
	        generateDownloadLink = function(family,region,filetype){
	           out=" no file ";

	           if($.inArray("region_"+region+"_"+filetype,dict[family][region])!=-1){
	               out="<a href="+family+"/region_"+region+"_"+filetype+">download</a>";
	           }
	           return out;
	           
	        };
	        generateShowLink = function(family,region){
	           out=" no image ";

	           if($.inArray("region_"+region+"_consensus.png",dict[family][region])!=-1){
	               out="<a href=#image"+elid+" onclick=showImage('"+elid+"','"+region+"','"+ family +"/region_"+region+"_consensus.png')>show</a>";
	           }
	           return out;
	           
	        };

	        
	        data = arr.data
	        for(var val in data){
	           for(var element in data[val][0]){     
	                   dict[data[val][0][0][3]+"_"+data[val][1]][data[val][0][element][5]].push([data[val][0][element][0],data[val][0][element][1],data[val][0][element][4]]);
	               
	           }      
	        }
	        
            var div = document.getElementById("files" + elid);
            div.innerHTML="";
            
            out = "<div><H3>Files for Sample "+arr.sample+"</H3></div><p/>";
            
            out += "<div class='" + elid + "_files' id='" + elid + "_files'><table id='"+elid+"table' class='"+elid+"table' border=1 cellpadding=3 cellspacing=0 style=' border: 1pt solid #000000; border-Collapse: collapse; font-size: 12px;'>";
            out += "<tr><th>family</th><th>region</th><th>#reads</th><th>#basepairs</th><th>region length</th><th>unaligned fasta</th><th>bam alignment</th><th>consensus fasta</th><th>consensus diagram</th></tr>"
            if (!dir) {
                for (var directory in dict) {
                    //out += "<b>"+directory+"</b>";
                    //out += "<div class='dirlist' id='" + directory + "_files'";
                    for (var region in dict[directory]) {
                        last = dict[directory][region][dict[directory][region].length-1];
                        if(last instanceof Array && last[0]){
                            out += "<tr><td>"+directory+"</td><td>"+region+"</td><td>"+last[0]+"</td><td>"+last[1]+"</td><td>"+last[2]+"</td><td>"+generateDownloadLink(directory,region,"unaligned.fa.bzip2")+"</td><td>"+generateDownloadLink(directory,region,"alignment.bam")+"</td><td>"+generateDownloadLink(directory,region,"consensus.fa")+"</td><td>"+generateShowLink(directory,region)+" "+generateDownloadLink(directory,region,"consensus.png")+"</td></tr>";
                        }
                    }
                    out += "</div>";
                }
            } else {
                out += "<b>"+ dir+"</b>" ;
                
                if (!region) {
                    for (var region in dict[dir]) {
                        last = dict[dir][region][dict[dir][region].length-1];
                        if(last instanceof Array && last[0]){
                            out += "<tr><td>"+dir+"</td><td>"+region+"</td><td>"+last[0]+"</td><td>"+last[1]+"</td><td>"+last[2]+"</td><td>"+generateDownloadLink(dir,region,"unaligned.fa.bzip2")+"</td><td>"+generateDownloadLink(dir,region,"alignment.bam")+"</td><td>"+generateDownloadLink(dir,region,"consensus.fa")+"</td><td>"+generateShowLink(dir,region)+" "+generateDownloadLink(dir,region,"consensus.png")+"</td></tr>";
                        }
                    }
                } else {

                    last = dict[dir][region][dict[dir][region].length-1];
                    if(last instanceof Array && last[0]){
                        out += "<tr><td>"+dir+"</td><td>"+region+"</td><td>"+last[0]+"</td><td>"+last[1]+"</td><td>"+last[2]+"</td><td>"+generateDownloadLink(dir,region,"unaligned.fa.bzip2")+"</td><td>"+generateDownloadLink(dir,region,"alignment.bam")+"</td><td>"+generateDownloadLink(dir,region,"consensus.fa")+"</td><td>"+generateShowLink(dir,region)+"     "+generateDownloadLink(dir,region,"consensus.png")+"</td></tr>";
                    }
                }
            }
            out += "</div>";
            
            $(div).append(out);
            $(div).append("<div class='close_files"+elid+"'style='cursor:pointer;'>x</div>");
            $("a").css({'color':'black'});
            $("table").css({'border':'collapse'});
            
            $(".close_files"+elid).click(function(){
                $('.files'+elid).animate({'height':'0px'},1000,'',function(){
                    $('.files'+elid).css({'display':'none'});
                    $('.files'+elid).removeClass("selected");
                });
                
                
            });
            $("div.close_files"+elid).css({'color':'white','position':'absolute','top':'2','left':'2','background-color':arr.colors[1],'border':'1px solid white','border-radius':'20px','height':'20px','width':'20px','text-align': 'center'});
            $("a.close_files"+elid).css({'color':'white', 'text-decoration': 'none'});
            $('.files'+elid).css({'width': 'auto','display':'none'});
            wi=$('.files'+elid).width()+15;
            

           
            if(!$('.files'+elid).hasClass("selected")){
                $('.files'+elid).css({'height': 'auto','display':'none'});
                hei=$('.files'+elid).height();
                hei=parseInt(hei)+10;
                hei=hei>400? hei=400 : hei;
                $('.files'+elid).css({'width':wi+'px','position':'relative','height':'0px', 'overflow' : 'auto','background-color':'rgb(223, 229, 249)', 'border':'2px solid silver', 'display':'block' });
                $('.files'+elid).animate({'height':hei+'px'}, 1000,'',function(){
                   
                    $('.files'+elid).addClass("selected");
                });
            }
            else{
                
                curr_hei=$('.files'+elid).height();
                $('.files'+elid).css({'height': 'auto','display':'none'});
                hei=$('.files'+elid).height();
                hei=parseInt(hei)+10;
                hei=hei>400? hei=400 : hei;
                $('.files'+elid).css({'height': curr_hei+'px', 'width':wi+'px','position':'relative', 'overflow' : 'auto','background-color':'rgb(223, 229, 249)', 'border':'2px solid silver', 'display':'block' }).animate({'height':hei+'px'}, 1000);	           
            }
	};

	 
	 showImage = function (elid,region, consensus_image,wi) {
            var div = document.getElementById("image" + elid);
            div.innerHTML="";
            out="<div><H3>Consensus Image for Region "+region+"</H3></div>";
            newImg= new Image();
            newImg=consensus_image
            wid=newImg.naturalWidth
            hig=newImg.naturalHeight
            out+="<img id='image_file"+elid+"' src="+consensus_image+" alt="+consensus_image+" style='width:"+wid+";height:"+hig+";'>";
            $(div).append(out)
     
     
            $(div).append("<div class='close_image"+elid+"'style='cursor:pointer;'>x</div>");
            $(".close_image"+elid).click(function(){
                
                    $('.image'+elid).css({'display':'none'});
                    $('.image'+elid).removeClass("selected");
                
                
            });
            $("div.close_image"+elid).css({'color':'white','position':'absolute','top':'2','left':'2','background-color':arr.colors[1],'border':'1px solid white','border-radius':'20px','height':'20px','width':'20px','text-align': 'center'});
            $("a.close_image"+elid).css({'color':'white', 'text-decoration': 'none'});
            
            wi=$('.files'+elid).width();
            
            
            if(!$('.image'+elid).hasClass("selected")){
                $('.image'+elid).css({'width':wi+'px','position':'relative','height':'auto', 'overflow' : 'auto','background-color':'rgb(223, 229, 249)', 'border':'2px solid silver', 'display':'block' });
                
            }
       
	 };


	this.each (
		function()
		{ init(this); }
	)
	
};

	// default values
	$.fn.jqBarGraph.defaults = {	
		sample: 'no_sample_id',
		barSpace: 10,
		width: 400,
		height: 600,
		color: '#000000',
		colors: false,
		lbl: '',
		sort: false, // 'asc' or 'desc'
		sortBar: false,
		tab: 'reads',
		position: 'bottom', // or 'top' doesn't work for multi type
		prefix: '',
		postfix: '',
		animate: true,
		speed: 3.0,
		legendWidth: 150,
		legend: false,
		type: false, // or 'multi'
		showValues: false,
		borderSize: 1,
		showValuesColor: '#fff',
		title: false
	};
	
	
	//sorting functions
	function sortReadsAsc(a,b){
	    sum_a=0
	    for(var values in a){
	       if(a[values] instanceof Array){
	           for(var val in a[values]){
	               sum_a+=a[values][val][0];
	           }
	       }
	    }
	    sum_b=0
	    for(var values in b){
	       if(b[values] instanceof Array){
	           for(var val in b[values]){
	               sum_b+=b[values][val][0];
	           }
	       }
	    }
	    
		if (sum_a<sum_b) return -1;
		if (sum_a>sum_b) return 1;
		return 0;
	}
	function sortBasesAsc(a,b){
	    sum_a=0
	    for(var values in a){
	       if(a[values] instanceof Array){
	           for(var val in a[values]){
	               sum_a+=a[values][val][1];
	           }
	       }
	    }
	    sum_b=0
	    for(var values in b){
	       if(b[values] instanceof Array){
	           for(var val in b[values]){
	               sum_b+=b[values][val][1];
	           }
	       }
	    }
	    
		if (sum_a<sum_b) return -1;
		if (sum_a>sum_b) return 1;
		return 0;
	}
	
	function sortReadsDesc(a,b){
		sum_a=0
	    for(var values in a){
	       if(a[values] instanceof Array){
	           for(var val in a[values]){
	               sum_a+=a[values][val][0];
	           }
	       }
	    }
	    sum_b=0
	    for(var values in b){
	       if(b[values] instanceof Array){
	           for(var val in b[values]){
	               sum_b+=b[values][val][0];
	           }
	       }
	    }
	    
		if (sum_a<sum_b) return 1;
		if (sum_a>sum_b) return -1;
		return 0;
	}
	function sortBasesDesc(a,b){
		sum_a=0
	    for(var values in a){
	       if(a[values] instanceof Array){
	           for(var val in a[values]){
	               sum_a+=a[values][val][1];
	           }
	       }
	    }
	    sum_b=0
	    for(var values in b){
	       if(b[values] instanceof Array){
	           for(var val in b[values]){
	               sum_b+=b[values][val][1];
	           }
	       }
	    }
	    
		if (sum_a<sum_b) return 1;
		if (sum_a>sum_b) return -1;
		return 0;
	}
	
    function sortBars(data,fun){
        for(var values in data){
                last = data[values].pop();
                for(var val in data[values]){
                    data[values][val].sort(fun);
                }
                data[values].push(last);       
        }           
    }
	
    function sortBars(data,fun){
        for(var values in data){
                last = data[values].pop();
                for(var val in data[values]){
                    data[values][val].sort(fun);
                }
                data[values].push(last);       
        }           
    }


    
    function barAsc(a,b){
        if(a[0]<b[0]) return -1;
        if(a[0]>b[0]) return 1;
        return 0;   
    }
	
	function barDesc(a,b){
	    if(a[0]<b[0]) return 1;
        if(a[0]>b[0]) return -1;
	    return 0;  
	}

})(jQuery);