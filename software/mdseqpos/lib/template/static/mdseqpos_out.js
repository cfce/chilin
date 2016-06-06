/*
 * mdseqpos_out.js - script to control mdseqpos_out.html UI
 * 
 * Copyright (c) 2012 Len Taing
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * Last Modified: 2012-06-12
 *
 * Description:
 * A set of fns which I wished was included in prototype.js
 */

/**
 * Class: ModelEvent
 * Description: This class does two important things: 1. register listeners
 * for a model event, and 2. when the model event is "invoked" (via the 
 * notify method), all of the registered listeners are informed.
 *
 * @param: modelRef -- the object that is passed as the message of notify
 * NOTE: it is customary for modelRef to be MODEL that uses the Model Event
 */
function ModelEvent(modelRef) {
    this.modelRef = modelRef;
    this.listeners = [];
    var outer = this;
    
    this.register = function(listenerFn) {
	var i = outer.listeners.indexOf(listenerFn);
	//make sure there are no duplicates
	if (i == -1) { outer.listeners.push(listenerFn);}
    }
    
    //WARNING: i think this is dependent on protoype and shouldn't be used!
    this.unregister = function(listenerFn) {
	var i = outer.listeners.indexOf(listenerFn);
	if (i != -1) { outer.listeners.splice(i, 1);}
    }
    
    this.notify = function() {
	for (var i = 0; i < outer.listeners.length; i++) {
	    outer.listeners[i](outer.modelRef);
	}
	//outer.listeners.each(function(lstnr) {lstnr(outer.modelRef);});
    }
}

// ****Models***

/**
 * Class: Motif
 * this function contains all of the information for any individual motif
 * fields: id, factors, consensus sequence, pssm - motif matrix, logoImg,
 * hit score, cutoff score, zscore, pvalue, and gene group
 * 
 * note: we don't have gene group information yet
 * 
 * The input should be a hashtable of the fields and their associated value
 */
function Motif(paramObjs) {	
    this.fields = ["id", "factors", "dbd", "entrezs", "refseqs", "species", 
		   "pssm", "numhits", "cutoff", "zscore", "pvalue",
		   "meanposition", 'seqpos_results'];
    //var fieldsToCopy = ["numhits", "cutoff", "zscore","pvalue","meanposition"];
    var outer = this;

    //1. set the values
    //e.g. this is short-hand for: this.id = paramObjs.id
    //2. create "get" functions, e.g. getid()
    for (var i = 0; i < outer.fields.length; i++) {
	/*
	if (fieldsToCopy.indexOf(outer.fields[i]) != -1) {
	    //EXCEPTION: seqpos_results
	    this[outer.fields[i]] = paramObjs['seqpos_results'][outer.fields[i]];
	} else {
	    this[outer.fields[i]] = paramObjs[outer.fields[i]];
	}
	*/
	this[outer.fields[i]] = paramObjs[outer.fields[i]];
	this["get"+outer.fields[i]] = function(field) {
	    return function(){return outer[field]; }
	}(outer.fields[i]);
    }

    this.toString = function() {
	tmp = "";
	for (var i = 0; i < outer.fields.length; i++) {
	    tmp += outer.fields[i]+":"+paramObjs[outer.fields[i]]+"\n";
	}
	return tmp;
    }
}

/**
 * Class: MotifModel
 * Description: model for the mdseqpos_out page
 * 
 * Input:
 * motifList: List of Motif objects (see Motif class above)
 * motifDists: Hashtable of pairwise motif distances
 */
function MotifModel(motifList, motifDists) {
    //BEGIN READ-ONLY fields
    //SET sort motifList by zscore--and create a helper hastable to map 
    //id->motif
    this.motifMap = {}
    this.motifList = motifList.sort(function(m1,m2){return m1.zscore-m2.zscore;});
    this.motifDists = motifDists;

    this.defaultCutoff = -11;
    this.defaultGroups = 5;
    this.defaultSpecies = "Any";
    //END READ-ONLY fields

    this.currCutoff = 0;//this.defaultCutoff;
    this.currGroups = this.defaultGroups;
    this.currSpecies = this.defaultSpecies;

    this.currMotifList = motifList;
    this.currClusters = null;
    
    //VIOLATING MVC: need to register the groups slider b/c groups may 
    //sometimes auto-set it to 1
    this.groupSldr = null;

    var outer = this;
    
    //init the motifMap
    for (var i = 0; i < this.motifList.length; i++) {
	outer.motifMap[outer.motifList[i].id] = outer.motifList[i];
    }

    var getsetters = ['currCutoff', 'currGroups', 'currSpecies', 
		      'currMotifList', 'currClusters'];
    
    //set the incoming params: e.g. this.entryId = paramObj.entryId         
    //and create get functions                                              
    for (var i = 0; i < getsetters.length; i++) {
	var name = upperCase1stLtr(getsetters[i]);
	//we need a closure for help                                  
	this["get"+name] = function(field) {
	    return function(){return outer[field]; }
	}(getsetters[i]);
	
	//create the event and set functions                                
	this[getsetters[i]+"Event"] = new ModelEvent(this);
	this["set"+name] = function(field) {
	    return function(param) {
		if (param != outer[field]) {                              
		    outer[field] = param;
		    outer[field+"Event"].notify();
		}
	    }
	}(getsetters[i])
    }

    //Listeners
    this.currCutoffEventLstnr = function() {
	//Filter out by the z-score, and set currMotifList
	var tmp = [];
	for (var i = 0; i < outer.motifList.length; i++) {
	    if (outer.motifList[i].zscore < outer.currCutoff) {
		tmp.push(outer.motifList[i]);
	    }
	}
	outer.setCurrMotifList(tmp);
    }
    this.currCutoffEvent.register(this.currCutoffEventLstnr);

    this.currGroupsEventLstnr = function() {
	//recluster using the currMotifList
	var tmp = outer.cluster(outer.currMotifList, outer.currGroups);
	outer.setCurrClusters(tmp);
    }
    this.currGroupsEvent.register(this.currGroupsEventLstnr);

    this.currSpeciesEventLstnr = function() {
	//filter the motifList by species
	//special case: 'Any' = no filtering
	var tmp = [];
	if (outer.currSpecies == "Any") {
	    tmp = outer.motifList;
	} else {
	    var list = outer.motifList;
	    for (var i = 0; i < list.length; i++) {
		if (checkStrInList(outer.currSpecies, list[i].species)) {
		    tmp.push(list[i]);
		}
	    }
	}
	outer.setCurrMotifList(tmp);
    }
    this.currSpeciesEvent.register(this.currSpeciesEventLstnr);

    this.currMotifListLstnr = function() {
	//might need to re-cluster
	var tmp = outer.cluster(outer.currMotifList, outer.currGroups);
	outer.setCurrClusters(tmp);
    }
    this.currMotifListEvent.register(this.currMotifListLstnr);

    //k means clustering based on motif distances
    this.cluster = function(mlist, k) {
	if (mlist.length < k) {
	    //IF the list is too small--set k to mlist.length
	    if (mlist.length > 0) {
		k = mlist.length;
		outer.setCurrGroups(k);
		if (outer.groupSldr) {
		    outer.groupSldr.value = k;
		    outer.groupSldr.notify(); //FAKING an event!
		}
	    } else {
		return [];
	    }
	}

	//INIT: ret to array of size k of empty arrays
	var ret = [];
	for (var i = 0; i < k; i++) {
	    ret.push([]);
	}

	//INIT: randomly throw motifs into clusters--throw first K into
	// the first k bins
	for (var i = 0; i < k; i++) {
	    ret[i].push(mlist[i].id);
	}
	for (var i = k; i < mlist.length; i++) {
	    var bin = Math.floor(Math.random()*k);
	    ret[bin].push(mlist[i].id);
	}
	 
	var loopCount = 0;
	do {
	    //calculate cluster means
	    var means = [];
	    for (var i = 0; i < ret.length; i++) {
		var tmp = outer.clusterMediod(ret[i]);
		//clusterMediod returns [dist, motif id]
		means.push(tmp[1]);
	    }

	    //INIT: a new temporary array
	    var tmp = [];
	    for (var i = 0; i < k; i++) {
		tmp.push([]);
	    }
	    //assignment step: assign it to the closest mediod
	    var newAssignment = false;
	    for (var i = 0; i < ret.length; i++) {
		for (var j = 0; j < ret[i].length; j++) {
		    //if mediod, then assign it to the curr cluster:
		    if (ret[i][j] == means[i]) {
			tmp[i].push(ret[i][j]);
		    } else {
			//assign it to the closest mediod
			var min = outer.getDist(ret[i][j], means[0]);
			var clss = 0;
			for (var l = 1; l < means.length; l++) {
			    var d = outer.getDist(ret[i][j], means[l]);
			    if (d < min) {
				min = d;
				clss = l;
			    }
			}
			if (clss != i) {
			    newAssignment = true;
			}
			tmp[clss].push(ret[i][j]);
		    }
		}
	    }
	    ret = tmp;
	    loopCount +=1;
	} while (newAssignment);
	
	//NOTE: returning motif IDS only
	return ret;
    }
    
    //given two motif ids, return the distance between the two motifs
    this.getDist = function(m1, m2) {
	var d = outer.motifDists[m1+":"+m2] ? outer.motifDists[m1+":"+m2] : outer.motifDists[m2+":"+m1];
	return d[0];
    }
    
    //given a list of motif ids, calculate the mediod--the motif whose distance
    //to all other motifs is minimal
    this.clusterMediod = function(lst) {
	var ret = [];
	for (var i = 0; i < lst.length; i++) {
	    var sum = 0;
	    for (var j = 0; j < lst.length; j++) {
		if (i != j) {
		    var d = outer.getDist(lst[i], lst[j]);
		    sum += d;
		}
	    }
	    ret.push([sum/lst.length, lst[i]]);
	}
	return ret.sort()[0];
    }

    //initial Cluster!
    this.setCurrCutoff(this.defaultCutoff);
    this.currClusters = this.cluster(this.currMotifList, this.currGroups);

}
// ****VIEW***
function MotifTableView(motifModel, container) {
    this.motifModel = motifModel;
    this.container = container;

    var outer = this;
    var _logoHeight = 90;
    var _logoWidth = 288;

    //save the TRs so that we can display them or not
    this.trByGroup;
    this.fields = ["group", "id", "factors", "dbd", "logo", "plot", "zscore", "pvalue"];

    this.showState;

    this.mClusters;

    this.makeHTML = function() {
	//CLEAR the container
	//check for empty list:
	if (outer.motifModel.currMotifList.length <= 0) {
	    outer.container.innerHTML = "";
	    var tmp = $D('h1',{innerHTML:'No Motifs to show'});
	    outer.container.appendChild(tmp);
	    return;
	}

	outer.container.innerHTML = "";
	//reset the showState: initially all of them are not shown
	outer.showState = [];
	for (var i = 0; i < outer.motifModel.getCurrClusters().length; i++) {
	    outer.showState.push(false);
	}

	var clusters = outer.motifModel.getCurrClusters();
	outer.mClusters = [];
	outer.trByGroup = [];
	//Convert list of motifIDs to list of motifs (sorted by zscore)
	for (var i = 0; i < clusters.length; i++) {
	    outer.mClusters.push(outer.mapMotifs(clusters[i]));
	}
	//sort groups by zscore--take the top hit i.e. first elm
	outer.mClusters.sort(function(a,b){return a[0].zscore - b[0].zscore;});

	//DRAW!
	for (var c = 0; c < clusters.length; c++) {
	    outer.trByGroup[c] = [];
	    var group = outer.mClusters[c];
	    
	    //create the table
	    var tbl = document.createElement('table');
	    tbl.className = "datatable";
	    var header = document.createElement('tr');
	    tbl.appendChild(header);
	    
	    //create the table header
	    for (var i = 0; i < outer.fields.length; i++) {
		var field = outer.fields[i];
		var tmp = document.createElement('th');
		var span = document.createElement('span');
		
		if (field == "dbd") {
		    span.innerHTML = "DNA BindDom";
		} else if (field == "group") {
		    span.innerHTML = "group"+c+" ("+group.length+")";
		} else {
		    span.innerHTML = field;
		}
		tmp.appendChild(span);
		header.appendChild(tmp);
	    }   

	    //add the motifs!
	    for (var i = 0; i < group.length; i++) {
		var newTr = document.createElement('tr');
		//add it to the TR list
		outer.trByGroup[c].push(newTr);

		newTr.id = group[i].id;
		//Alternating rows effect
		newTr.className = (i % 2 == 0) ? "altrow":"";
		tbl.appendChild(newTr);
		for (var j = 0; j < outer.fields.length; j++) {
		    var tmp = document.createElement('td');
		    
		    if (outer.fields[j] == "group") {//SPECIAL CASES
			tmp.innerHTML = "group"+c;

			tmp.onmouseover = function() {
			    this.style.cursor = "pointer";
			}
			tmp.onmouseout = function() {
			    this.style.cursor = "pointer";
			}
			tmp.onclick = function(groupN) {
			    return function() { outer.toggle(groupN)}
			} (c);

		    } else if (outer.fields[j] == 'factors') { 
			var factors = group[i][outer.fields[j]];
			if (factors) {
			    //FACTORS is a list
			    var fact = "";
			    for (var k = 0; k < factors.length; k++) {
				fact += factors[k]+ "\n";
			    }
			    tmp.innerHTML = fact;
			}
		    } else if (outer.fields[j] == 'logo') {
			var canvas = makeDNALogo(group[i].pssm, _logoHeight, 
						 _logoWidth);
			var img = $D("img");
			//img.src = canvas.toDataURL('image/png');
			tmp.appendChild(canvas);
			
			/*
			tmp.appendChild($D('br'));
			var btn = $D('input', {type:'button', value:'show pssm'});
			tmp.appendChild(btn);
			btn.style.display="none";
			*/

			canvas.ondblclick = function(motif) {
			    return function () {
				//SHOW PSSM window
				//NOTE:this call works in safari but not in chrome!
				popUpPssm(motif);
			    }
			}(group[i]);
			tmp.onmouseover=function() {this.style.cursor="pointer";};
			tmp.onmouseout=function() {this.style.cursor="auto";};
		    } else if (outer.fields[j] == 'plot') {
			//PLOT HERE
			//x-labels
			var xCat = [];
			for (var z = 0; z< group[i].seqpos_results.plot.bin_avg.length; z++) {
			    xCat.push(z*group[i].seqpos_results.plot.chunk);
			}
                        var chart = new Highcharts.Chart({
                            chart: {
                                renderTo: tmp,
                                backgroundColor:'rgba(255, 255, 255, 0.1)',
                                type: 'line',
                                width: _logoWidth,
                                height: _logoHeight,
                            },
                            title: {text: ''},
                            legend: {enabled:false},
			    xAxis: {categories:xCat, tickInterval:10},
                            yAxis: {
                                title: {text: '',},
                                min:0, max:200,
                                tickInterval:50,
                            },
                            plotOptions: {
                                line: {
                                    dataLabels: {enabled: false},
                                    enableMouseTracking: false,
                                    animation:false,
                                    color:'#A60000',
                                },
				series: { marker: {radius:0}},
                            },
                            series: [{
				name: '',
                                data: group[i].seqpos_results.plot.bin_avg
                            }],
                            credits:{enabled:false}
                        });
		    } else if (outer.fields[j] == 'pvalue') {
			var val = group[i][outer.fields[j]];
			if (val < 0.0001) {
			    //convert 0.00000257564144033 --> 2.5e-6
			    var pow = 0;
			    while(val * Math.pow(10, pow) < 1) {
				pow++;
			    }
			    tmp.innerHTML = Math.floor(val*Math.pow(10,pow+1))/10 +"e-"+pow;
			} else {
			    tmp.innerHTML = val;
			}
		    } else {
			tmp.innerHTML = group[i][outer.fields[j]];
		    }
		    newTr.appendChild(tmp);
		}
	    }
	    
	    outer.container.appendChild(tbl);
	    //add a break
	    outer.container.appendChild($D('br'));
	 
	    //show/hide the group
	    outer.showGroup(c);
	}

    }

    this.init = function() {
	if (outer.motifModel.currMotifList.length > 0) {
            outer.makeHTML();
	} else { //No motifs found, remove the reset btn
            var reset_btn = document.getElementById("reset_btn");
            reset_btn.parentNode.removeChild(reset_btn);
	}
    }

    this.mapMotifs = function(motifIDs) {
	var ret = [];
	for (var i = 0; i < motifIDs.length; i++) {
	    ret.push(outer.motifModel.motifMap[motifIDs[i]]);
	}
	//sort by zscore--the smallest, i.e. the most negative/sig. is first!
	ret.sort(function(a, b) { return a.zscore - b.zscore;});
	return ret;
    }

    this.showGroup = function(groupN) {
	//Given a group, checks the showState--if true, then the other TRs
	//are shown, otherwise they are hidden
	var show = outer.showState[groupN];
	for (var i = 1; i < outer.trByGroup[groupN].length; i++) {
	    outer.trByGroup[groupN][i].style.display = (show)?"table-row":"none";
	}
    }

    this.toggle = function(groupN) {
	outer.showState[groupN] = !outer.showState[groupN];
	outer.showGroup(groupN);
    }

}


/*
var foo = [{name: 'a', zscore:-71}, {name:'b', zscore: -50}, {name:'c', zscore:-33}, {name:'d', zscore:-84}, {name:'e', zscore:-83}];

var tmp = foo.sort(function(a,b) {return a.zscore - b.zscore;});


for (var i = 0; i < foo.length; i++) {
    alert(foo[i].name);
}
*/

var motifModel;
var motifTblView;

function initPage() {

    //Convert the json to Motif Objects
    var motifList = [];
    var fieldsToCopy = ["numhits", "cutoff", "zscore","pvalue","meanposition"];
    for (var i = 0; i < motifList_json.length; i++) {
	//NOTE: zscore, pvalue, cutoff, numhits, meanposition are in 
	//motifList_json[i][seqpos_results]--copy it over
	for (var j = 0; j < fieldsToCopy.length; j++) {
	    motifList_json[i][fieldsToCopy[j]] = motifList_json[i].seqpos_results[fieldsToCopy[j]];
	}
	//alert(motifList_json[i]['zscore']);
	motifList.push(new Motif(motifList_json[i]));
    }

    motifModel = new MotifModel(motifList, motifDists);
    
    motifTblView = new MotifTableView(motifModel, 
				      document.getElementById('motif_table'));
    motifTblView.init();

    var redraw = function() { 
	motifTblView.makeHTML();
    }
    motifModel.currGroupsEvent.register(redraw);
    motifModel.currMotifListEvent.register(redraw);

    //init the sliders
    var zscoreSldr = document.getElementById('zscoreSlider');
    zscoreSldr.max = 0;
    document.getElementById('zMax').innerHTML = 0;
    //using the fact that the motifModel.motifList is sorted
    var zMin = Math.floor(motifModel.motifList[0].zscore);
    zscoreSldr.min = zMin;
    document.getElementById('zMin').innerHTML = zMin;
    zscoreSldr.value = motifModel.defaultCutoff;
    var zVal = document.getElementById('zscoreSldrVal');
    zVal.innerHTML = motifModel.defaultCutoff;

    zscoreSldr.onmouseup = function() { 
	zVal.innerHTML = this.value;
	motifModel.setCurrCutoff(this.value);
    }

    //set the number of motifs
    var numMotifSpan = document.getElementById('numMotifs');
    numMotifSpan.innerHTML = motifModel.currMotifList.length;
    motifModel.currMotifListEvent.register(function() {
	    numMotifSpan.innerHTML = motifModel.currMotifList.length;
	});

    var groupsSldr = document.getElementById('groupsSlider');
    groupsSldr.min = 1;
    document.getElementById('grpMin').innerHTML = 1;
    groupsSldr.max = 11;
    document.getElementById('grpMax').innerHTML = 11;
    groupsSldr.value = motifModel.defaultGroups;
    var groupsVal = document.getElementById('grpSldrVal');
    groupsVal.innerHTML = motifModel.defaultGroups;

    groupsSldr.onmouseup = function() {
	groupsVal.innerHTML = this.value;
	motifModel.setCurrGroups(this.value);
    }
    groupsSldr.notify = function() {
	groupsVal.innerHTML = this.value;
    }

    //register the groupSldr w/ motifModel
    motifModel.groupSldr = groupsSldr;

    //activate the dump btn
    var dumpTablesBtn = document.getElementById('txt_btn');
    dumpTablesBtn.onclick = function () { dumpTable();}

    //activate the reset btn
    var resetBtn = document.getElementById('reset_btn');
    resetBtn.onclick = function() {
	motifModel.setCurrCutoff(motifModel.defaultCutoff);
	zscoreSldr.value = motifModel.defaultCutoff;
	zVal.innerHTML = motifModel.defaultCutoff;

	motifModel.setCurrGroups(motifModel.defaultGroups);
	groupsSldr.value = motifModel.defaultGroups;
	groupsVal.innerHTML = motifModel.defaultGroups;
    }

}

function dumpTable() {
    //based on popupWindow and popupPSSM in focus_frame.js
    var win = window.open("focus_frame_2.html", '',"width=900, height=480,"+
                          "resizable=yes, scrollbars=yes,toolbar=no,"+
                          "location=no, menubar=no, status=yes");

    //NEED to create the motif elements AFTER the page load
    win.onload = function(event) {
        //BUILD the new window's DOM
        win.document.title = "Seqpos Output";
        var tmp = win.document.createElement('div');
        win.document.body.appendChild(tmp);

        var pre = win.document.createElement('pre');
        //pre.innerHTML = "Hello World";
        pre.style.fontSize="10pt";
        tmp.appendChild(pre);

        //build the table:
        var txtTbl = "";
		
	if (motifModel.motifList.length > 0) {
	    //copy the fields, and drop pssm, b/c they're not needed
	    var fields = motifModel.motifList[0].fields.slice(0);
	    fields.splice(fields.indexOf('pssm'), 1);

	    //dump the data -- use the clusters in the motifTblView
	    var mClusters = motifTblView.mClusters;
	    for (var c = 0; c < mClusters.length; c++) {
		var lst = mClusters[c];
		//build table header 
		for (var i=0; i < fields.length; i++) {
		    txtTbl += fields[i]+"\t";
		}
		txtTbl += "\n";

		for (var i = 0; i < lst.length; i++) {
		    for (var j = 0; j < fields.length; j++) {
			var val = lst[i][fields[j]];
			if (val) {
			    txtTbl += val+"\t";
			} else {
			    txtTbl += "\t";
			}
		    }
		    txtTbl += "\n";
		}
		txtTbl += "\n";
	    }
	    pre.innerHTML = txtTbl;
	}
    }
}

function popUpPssm(motif) {
    var win = window.open("focus_frame_2.html", '',"width=300, height=600,"+
			  "resizable=yes, scrollbars=yes, toolbar=no,"+
			  "location=no, menubar=no, status=yes");

    //NEED to create the PSSM matrix AFTER the page is loaded
    win.onload = function(event) {
	win.document.title = motif.id + " PSSM";
	//load the pssm_popup.css style
	var css = win.document.createElement('link');
	css.href = 'pssm_popup.css';
	css.rel = "stylesheet";
	css.type = "text/css";
	win.document.body.appendChild(css);
	//BUILD the table
	var tbl = win.document.createElement('table');
	tbl.className = "datatable";
	win.document.body.appendChild(tbl);

	var header = win.document.createElement('tr');
	tbl.appendChild(header);
	//BUILD the header
	var order = ['', 'A', 'C', 'G', 'T'];
	for (var i = 0; i < order.length; i++) {
	    var th = win.document.createElement('th');
	    th.innerHTML = order[i];
	    header.appendChild(th);
	}
	var pssm = "["; //represent the pssm as a matrix
	for (var i = 0; i < motif.pssm.length; i++) {
	    var tr = win.document.createElement('tr');
	    tbl.appendChild(tr);
	    //The nucleotide index
	    var tmp = win.document.createElement('td');
	    tmp.innerHTML = i+1;
	    if (i % 2 == 0) { tr.className="altrow"; }
	    tr.appendChild(tmp);
	    var row_pssm = "["; //represent the pssm row
	    for (var j = 0; j < motif.pssm[i].length; j++) {
		var td = win.document.createElement('td');
		tr.appendChild(td);
		td.innerHTML = motif.pssm[i][j];
		
		//build row_pssm
		row_pssm += (j != 0) ? ", " : "";
		row_pssm += motif.pssm[i][j];
	    }
	    row_pssm += "]"; //close row
	    pssm += (i != 0) ? ",\n" : "";
	    pssm += row_pssm;
	}
	pssm += "]"; //close the matrix
	//var build the txt dump
	var div = win.document.createElement('div');
	win.document.body.appendChild(div);
	var p = win.document.createElement('p');
	p.innerHTML="If you would like to input this motif into our \'Screen all motifs\' tool: <br/>OPTION A: <br/>1. copy the text below, and <br/>2. paste it into \"PSSM Raw Text\" in the Motif:Screen All Motifs tool<br/>--OR--<br/></br>OPTION B:<br/>1. copy and paste the txt below into a text file <br/>2. upload the file onto cistrome, and <br/>3. run Motif:Screen All Motifs selecting the text file as the PSSM file";
	p.style.fontSize="10pt";
	div.appendChild(p);
	var pre = win.document.createElement('pre');
	pre.innerHTML = pssm;
	pre.style.fontSize="10pt";
	div.appendChild(pre);	
    }
}

