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
 			$(el).css({'height': opts[el.id].height, 'position':'relative', 'text-align':'center' });


 			doGraph(el);

 		};

	// sum of array elements
	sum = function(ar,tab){
		total = 0;
		for(var val in ar){
			total += ar[val][tab];
		}
		return total.toFixed(2);
	};
	
	// count max value of array
	max = function(ar,tab){
		maxvalue = 0;
		for(var fam in ar){
			for(var sample in ar[fam]["samples"]){

				value = sum(ar[fam]["samples"][sample],tab);
				if (parseFloat(value) > parseFloat(maxvalue)) maxvalue=value;
			}
		}	
		return maxvalue;	
	};
	

	doGraph = function(el){
		
		arr = opts[el.id];
		data = arr.data;
		data_keys=Object.keys(data);
		color_map=arr.colors;
		
		n_bars=0;

		for(fam in data){
			for(sample in data[fam]){
				n_bars+=1;
			}
		}


		//check if array is bad or empty
		if(data == undefined) {
			$(el).html('There is not enought data for graph');
			return;
		}
		//sorting ascending or descending
		
		if(arr.sort == 'asc'){
			data_keys=sortGraph(arr,sortAsc);    
		}
		if(arr.sort == 'desc'){
			data_keys=sortGraph(arr,sortDesc);        
		}
	

		max1 = max(data,arr.tab);

		m = Math.max(max(data,'reads'),max(data,'bp'));

		w = (m+"").length * 10;
		legend = '';
		prefix = arr.prefix;
		postfix = arr.postfix;
		space = arr.barSpace; //space between bars
		legendWidth = arr.legend ? arr.legendWidth : 0; //width of legend box
		fieldWidth = w //width of bar
		
		totalHeight =  $(el).height(); //total height of graph box
		var leg = new Array(); //legends array
		
		//max value in data, I use this to calculate height of bar
		

		var i=0;
		var fields=0
		for(var key in data_keys){
			space=arr.barSpace;
			fam = data_keys[key];

			for(var sample in data[fam]["samples"]){ 

				valueData = data[fam]["samples"][sample];
				value = sum(valueData,arr.tab);
				sample_keys=Object.keys(valueData);

				lbl = fam
				unique = fam+sample+el.id; //unique identifier
				fieldWidth = w+space;
				fields = fields+w+arr.barSpace;

				if (lbl == undefined) lbl = arr.lbl;

				margin_top=4;

				unique_fam = fam+el.id;
				out  = "<div class='graphField"+el.id+"' id='graphField"+unique+"' style='position: absolute'>";

				out += "<a class='graphLabelLink' href=#files"+el.id+" onclick=fillDiv('"+el.id+"','"+fam+"','"+sample+"')><div class='graphValue"+el.id+"' id='graphValue"+unique+"'>"+parseInt(value)+"</div></a>";

				out += "<div class='graphBar"+el.id+"' id='graphFieldBar"+unique+"' style='background-color:#fff;position: relative; overflow: hidden;'></div>";
	            //console.log(color)
				// if there is no legend or exist legends display lbl at the bottom

				out += "<a class='graphLabelLink' href=#files"+el.id+" onclick=fillDiv('"+el.id+"','"+fam+"')><div class='graphLabel"+el.id+"' id='graphLabel"+unique+"'style='margin-top:"+margin_top+"px;'>"+lbl+"</div></a>";
				out += "</div>";



				$(el).append(out);
				$(".graphLabel"+el.id).css({'-webkit-transform': 'rotate(30deg)', '-moz-transform': 'rotate(30deg)','-o-transform': 'rotate(30deg)', '-ms-transform': 'rotate(30deg)','transform': 'rotate(30deg)', 'height':'100'});
				$('a.graphLabel').css({
					'text-decoration': 'none',
					'color':'black'    
				});

				totalHeightBar = totalHeight - $('.graphLabel'+el.id).height() - $('.graphValue'+el.id).height()-margin_top; 
				fieldHeight = (totalHeightBar*value)/max1;

				$('#graphField'+unique).css({ 
					'left': (w+arr.barSpace)*i, 
					'width': w, 
					'margin-left': space});


				maxe = max1;
				k=0;
				for(var r in valueData){
					k+=1;
				}


				if(arr.sortBar=='asc'){
					sample_keys=sortBar(valueData,arr.tab, barAsc);
				}
				if(arr.sortBar=='desc'){
					sample_keys=sortBar(valueData,arr.tab, barDesc);
				}


				for(var s_key in sample_keys){
					region = sample_keys[s_key];

					heig = (totalHeightBar*valueData[region][arr.tab]/maxe);

					if(navigator.userAgent.match('Safari') && !navigator.userAgent.match('Chrom')){

						heig = heig + (arr.borderSize/1.85);
					} 

					wid = parseInt((fieldWidth-space)/k);
	 				sv = ''; // show values
	 				fs = 0; // font size
	 				if (arr.showValues){
	 					sv = arr.prefix+valueData[region][arr.tab]+arr.postfix;
	 					fs = 12; // font-size is 0 if showValues = false
	 				}

	 				o = "<a class='subBarLink' href=#files"+el.id+" onclick=fillDiv('"+el.id+"','"+fam+"','"+sample+"','"+region+"')><div class='subBars"+el.id+"' style=' box-sizing:border-box; -moz-box-sizing:border-box; -ms-box-sizing:border-box; -webkit-box-sizing:border-box; height:"+heig+"px;  border-top:"+arr.borderSize+"px solid; border-color: "+arr.showValuesColor+"; background-color: "+color_map[sample][valueData[region]["state"]]+"; left:"+wid*parseInt(region)+"px; color:"+arr.showValuesColor+"; font-size:"+fs+"px' >"+sv+"</div></a>";
	 				$('#graphFieldBar'+unique).prepend(o);
	 			}

	 			//position of bars
	 			if(arr.position == 'bottom') $('.graphField'+el.id).css('bottom',0);



	 			
	 			
	 			// animated apearing
	 			if(arr.animate){
	 				$('#graphFieldBar'+unique).css({ 'height' : 0 });
	 				$('#graphFieldBar'+unique).animate({'height': fieldHeight},arr.speed*1000);
	 			} else {
	 				$('#graphFieldBar'+unique).css({'height': fieldHeight});
	 			}
	 			i=i+1;
	 			space=0;
	 		}
	 		
	 	}
	 	ws = fields + legendWidth + 40;
		$(el).css({'width': ws})
		arr.width=ws;

 		createLegend(color_map,el.id); // create legend from array
 		createLinks(el.id);
 		//position of legend
 		
 		$(el).append("<div id='legendHolder"+el.id+"'></div>");
 		$('#legendHolder'+el.id).css({ 'width': 200, 'float': 'right','text-align' : 'left'});
 		$('#legendHolder'+el.id).append(legend);
 		$('#legendHolder'+el.id).append(links);
 		$('.legendBar'+el.id).css({ 'float':'left', 'margin': 3, 'margin-left':'40','height': 12, 'width': 20, 'font-size': 0});
 		$('.linkBar'+el.id).css({'margin-left':'40'});


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
 			opts[el.id].tab='bp';
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
		for(var sample in color_map){
			legend += "<div id='legend"+sample+"' style='overflow: hidden; zoom: 1; margin-left: 20;'><b>"+sample+"</b></div>";
			for(var state in color_map[sample]){
				legend += "<div id='legend"+sample+state+"' style='overflow: hidden; zoom: 1;'>";
				legend += "<div class='legendBar"+id+"' id='legendColor"+sample+state+"' style='background-color:"+color_map[sample][state]+"'></div>";
				legend += "<div class='legendLabel"+id+"' id='graphLabel"+unique+"'>"+state+"</div>";
				legend += "</div>";	
			}			
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
	
	fillDiv = function (elid, family,sample, region) {



		arr=opts[elid].data;

		generateDownloadLink = function(family,sample,region,filetype){
			out=" no file ";

			if(arr[family]["samples"][sample][region]["files"][filetype]){
				out="<a href="+arr[family]["type"]+"_"+family+"/region_"+region+"_"+filetype+">download</a>";
			}
			return out;

		};
		generateShowLink = function(family,sample,region){
			out=" no image ";
			if(arr[family]["samples"][sample][region]["files"]["consensus.png"]){
				out="<a href=#image"+elid+" onclick=showImage('"+elid+"','"+region+"','"+arr[family]["type"]+"_"+family +"/region_"+region+"_consensus.png')>show</a>";
			}
			return out;

		};
		generateLine=function(family,sample,region){

			return "<tr><td>"+family+"</td><td>"+sample+"</td><td>"+region+"</td><td>"+arr[family]["samples"][sample][region]["reads"]+"</td><td>"+arr[family]["samples"][sample][region]["bp"]+"</td><td>"+arr[family]["samples"][sample][region]["length"]+"</td><td>"+generateDownloadLink(family,sample,region,"unaligned.fa.bzip2")+"</td><td>"+generateDownloadLink(family,sample,region,"alignment.bam")+"</td><td>"+generateDownloadLink(family,sample,region,"consensus.fa")+"</td><td>"+generateShowLink(family,sample,region)+" "+generateDownloadLink(family,sample,region,"consensus.png")+"</td></tr>";

		};


		var div = document.getElementById("files" + elid);
		div.innerHTML="";

		out = "<div><H3>Files for selected data "+sample+"</H3></div><p/>";

		out += "<div class='" + elid + "_files' id='" + elid + "_files'><table id='"+elid+"table' class='"+elid+"table' border=1 cellpadding=3 cellspacing=0 style=' border: 1pt solid #000000; border-Collapse: collapse; font-size: 12px;'>";
		out += "<tr><th>family</th><th>sample</th><th>region</th><th>#reads</th><th>#basepairs</th><th>region length</th><th>unaligned fasta</th><th>bam alignment</th><th>consensus fasta</th><th>consensus diagram</th></tr>"
		if (!family) {
			if(sample){
				for (var fam in arr) {
					for (var region in arr[fam]["samples"][sample]) {
						out+= generateLine(fam,sample,region);
					}

					out += "</div>";
				}

			}else{
				for (var fam in arr) {
					for(var sample in arr[fam]["samples"]){
						for (var region in arr[fam]["samples"][sample]) {
							out+= generateLine(fam,sample,region);
						}

						out += "</div>";
					}
				}
			}
		} else {
			out += "<b>"+ family+"</b>" ;
			if (!sample) {
				for (var sample in arr[family]["samples"]) {
					for(var region in arr[family]["samples"][sample]){

						out+= generateLine(family,sample,region);            			
					}
				}	
			} else {
				if(!region){
					for(var region in arr[family]["samples"][sample]){
						out+= generateLine(family,sample,region);
					}
				} else {

					out+= generateLine(family,sample,region);
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
		$("div.close_files"+elid).css({'color':'white','position':'absolute','top':'2','left':'2','background-color':'red','border':'1px solid white','border-radius':'20px','height':'20px','width':'20px','text-align': 'center'});
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
		newImg.src=consensus_image;
		wid=newImg.width;
		hig=newImg.height;
		out+="<img id='image_file"+elid+"' src="+consensus_image+" alt="+consensus_image+" style='width:"+wid+";height:"+hig+";'>";
		$(div).append(out)


		$(div).append("<div class='close_image"+elid+"'style='cursor:pointer;'>x</div>");
		$(".close_image"+elid).click(function(){

			$('.image'+elid).css({'display':'none'});
			$('.image'+elid).removeClass("selected");


		});
		$("div.close_image"+elid).css({'color':'white','position':'absolute','top':'2','left':'2','background-color':'red','border':'1px solid white','border-radius':'20px','height':'20px','width':'20px','text-align': 'center'});
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
	
	function sortGraph(arr,fun){
		
		sum = function(ar,tab){
			total = 0;
			for(var val in ar){
				total += ar[val][tab];
			}
			return total;
		};


		data = arr.data;
		keys=Object.keys(data);
		tab = arr.tab;
		keys2 = new Array(keys.length);
		
		for(var key in keys){
			fam = keys[key];
			ss=0;
			for(var sample in data[fam]['samples']){

				s = sum(data[fam]["samples"][sample],tab);
			}
			ss = ss+s;
			keys2[key]={'key':keys[key], 'sum':ss};
		}
		keys2.sort(fun);

		for (key in keys){

			keys[key] = keys2[key]['key'];
		}
		return keys;
	};
	
	//sorting functions
	function sortAsc(a,b){

		if (a["sum"]<b["sum"]) return -1;
		if (a["sum"]>b["sum"]) return 1;
		return 0;
	}
	function sortDesc(a,b){

		if (a["sum"]<b["sum"]) return 1;
		if (a["sum"]>b["sum"]) return -1;
		return 0;
	}
	
	
	function sortBar(sample,tab,fun){
		keys=Object.keys(sample);
		keys2 = new Array(data_keys.length);
		for(var key in keys){
			keys2[key]={ 'key':keys[key],'value':sample[keys[key]][tab]};

		}

		keys2.sort(fun);
		for (key in keys){
			keys[key]= keys2[key]['key'];
		}

		return keys;
	}
	



	function barAsc(a,b){
		if(a['value']<b['value']) return -1;
		if(a['value']>b['value']) return 1;
		return 0;   
	}
	
	function barDesc(a,b){
		if(a['value']<b['value']) return 1;
		if(a['value']>b['value']) return -1;
		return 0;  
	}

})(jQuery);