/*
 * goLogo.js - a script to dynamically generate sequence logos
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
 * Function to return the heights of each base-pair given a pssm
 * ref: http://en.wikipedia.org/wiki/Sequence_logo
 * height(a,i) - where a is the base, and i is the pos = f(a,i) * Ri
 * f(a,i) - frequency of base a in pos i, i.e. pssm[i][a]
 * Ri = for nucleic acids, Ri = 2 − (Hi + en)
 * Hi = - SUM(f(a,i) * log2(f(a,i)))
 * en = (s - 1) / (2 * ln(2) * n) - where s = 4 base pairs, n = length of pssm0
 *
 * heights works for nucleic OR amino acids!
 */
function heights(pssm) {
    //s = 4 for nucleic acids, or 20 for amino acids--this depends on how 
    //many cols are in pssm
    var s = pssm[0].length;
    var en = (s - 1.0)/(2 * Math.log(2) * pssm.length);
    var R = [];
    for (var i = 0; i < pssm.length; i++) {
	var row = pssm[i];
	var sum = 0;
	for (var j = 0; j < row.length; j++) {
	    //check for 0 (ln of 0 is Nan!--and cast as float
	    var val = row[j] + 10e-20;
	    sum += val * Math.log(val, 2);
	}
	//Ri  = log2(s) - (Hi + en), but for some reason, Math.log(s, 2)
	//is incorrect, so we use the log division rule	
	R.push((Math.log(s)/Math.log(2)) - (-1.0*sum + en));
    }

    //calc heights
    var tmp = [];
    for (var i = 0; i < pssm.length; i++) {
	var row = pssm[i];
	var h = [];
	for (var j = 0; j < row.length; j++) {
	    h.push(row[j] * R[i]);
	}
	tmp.push(h);
    }
    return tmp;
}

function drawYAxis(ctx, numBits, height, margin) {
    //draw y-axis
    ctx.strokeStyle = 'black';
    ctx.beginPath();
    ctx.moveTo(margin, 0);
    ctx.lineTo(margin, height - margin);
    ctx.stroke();
    //Enumerate y-axis
    ctx.fillStyle="black";
    var dh = Math.round((height - margin) / numBits);
    for (var i = 0; i <= numBits; i++) {
	if (i != numBits) {
	    ctx.fillText(i, margin - 10, (height - margin) - dh*i);
	} else {
	    ctx.fillText(i, margin - 10, 10);
	}
    }
}

function loadImgs(subdir, imgNames) {
    var imgs = [];
    for (var i = 0; i < imgNames.length; i++) {
	var tmp = document.createElement("img");
	tmp.src = subdir+"/"+imgNames[i]+".png";
	imgs.push(tmp);
    }
    return imgs;
}
/**
 * Function to generate a sequence logo 
 * Generalize: numBits = 2.0 or 10.0 (for amino)???
 * -Can't hardcode y-axis
 * -images hardcoded to be ACGT
 * returns the canvas that was made that can be attached to the DOM or 
 * further manipulated
 */
function makeDNALogo(pssm, height, width) {
    //Constants
    var imgNames = ["a","c","g","t"];
    var imgs = loadImgs("dna", imgNames);

    var numBits = 2.0;
    var canvas = document.createElement('canvas');
    canvas.height = height;
    canvas.width = width;
    
    var ctx = canvas.getContext('2d');
    var margin = (height > width) ? (0.05 * height) : (0.05 * width);
    margin = Math.round(margin);

    drawYAxis(ctx, numBits, height, margin);

    //ctx.textAlign='center';

    //x-axis
    ctx.beginPath();
    ctx.moveTo(margin, height - margin);
    ctx.lineTo(width, height - margin);
    ctx.stroke();

    var w = Math.round((width - margin) / pssm.length);
    var h = heights(pssm);

    for (var i = 0; i < h.length; i++) {
	var tmp = [];
	for (var j = 0; j < h[i].length; j++) {
	    tmp.push({'img':imgs[j], 'height':h[i][j]});
	}
	//sort ascending
	tmp.sort(function(a,b) { return a.height - b.height;});
	var start = height - margin;
	for (var j = 0; j < tmp.length; j++) {
	    var pxH = Math.round((height - margin) * (tmp[j].height/numBits));
	    ctx.drawImage(tmp[j].img, w*i + margin, (start - pxH), w, pxH);
	    ctx.fillText(i + 1, w*i + margin + w/2.0, (height - margin + 10));
	    start -= pxH;
	}
    }

    return canvas;
}

function makeAminoLogo(container, pssm, height, width) {
    var imgNames=["a","r","n","d","c","e","q","g","h","i",
    		  "l","k","m","f","p","s","t","w","y","v"];
    var imgs = loadImgs("amino",imgNames);
    var numBits = 4.0;
    var canvas = document.createElement('canvas');
    canvas.height = height;
    canvas.width = width;
    
    var ctx = canvas.getContext('2d');
    var margin = (height > width) ? (0.05 * height) : (0.05 * width);
    margin = Math.round(margin);

    drawYAxis(ctx, numBits, height, margin);

    //x-axis
    ctx.beginPath();
    ctx.moveTo(margin, height - margin);
    ctx.lineTo(width, height - margin);
    ctx.stroke();

    var w = Math.round((width - margin) / pssm.length);
    var h = heights(pssm);

    for (var i = 0; i < h.length; i++) {
	var tmp = [];
	for (var j = 0; j < h[i].length; j++) {
	    tmp.push({'img':imgs[j], 'height':h[i][j]});
	}
	//sort ascending
	tmp.sort(function(a,b) { return a.height - b.height;});
	var start = height - margin;
	for (var j = 0; j < tmp.length; j++) {
	    var pxH = Math.round((height - margin) * (tmp[j].height/numBits));
	    ctx.drawImage(tmp[j].img, w*i + margin, (start - pxH), w, pxH);
	    ctx.fillText(i + 1, w*i + margin + w/2.0, (height - margin + 10));
	    start -= pxH;
	}
    }

    container.appendChild(canvas);
}

