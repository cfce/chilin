/*
 * Utils.js
 * 
 * Copyright (c) 2011 Len Taing
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
 * Function: $D 
 * Description: dom element creator shortcut/convenience
 * var tmp = $D('span', {'class':'foo'});
 * RATHER than:
 * var tmp = document.createElement('span');
 * tmp.class = "foo";
 */
function $D(elmType, opts) {
    var tmp = document.createElement(elmType);
    if (opts) {
	//set optional fields
	for (var i in opts) {
	    tmp[i] = opts[i];
	}
    }

    return tmp;
}

/**
 * Function: getattr - my homage to the python fn
 * Description: dereferencing fn
 * @param: obj - javascript object
 * @param: field - string of the field to get
 * KEY note: this fn can follow dots., e.g. obj = {'foo':{'bar':5}}
 * getattr(obj, "foo.bar") --> 5
 */
function getattr(obj, field, debug) {
    var fldLst = field.split(".");
    var curr = obj;

    //just keep deref until we exhaust the list, or we break the link
    for (var i = 0; i < fldLst.length; i++) {
	//whenever we reach null, return
	if (curr == null) {
	    return null;
	}

	//otherwise, if we can't find the key
	var keys = Object.keys(curr);
	if (keys.indexOf(fldLst[i]) == -1) {
	    return null;
	} else {
	    curr = curr[fldLst[i]];
	}
    }
    return curr;
}

/** 
 * Function: inspect - inspects jscript objects
 */
function inspect(obj) {
    var flds = "";
    for (x in obj) {
	flds += x + "\n";
    }
    return flds;
}

/**                                                                            
 * Fn: takes a string and returns the string w/ the first letter uppercased
 * NOTE: maybe I should move the following to a file called StringUtils. 
 */
function upperCase1stLtr(s) {
    if (s == null) { return null;}
    if (s.length > 0) {
        var ltr = s.substring(0,1);
        var rest = s.substring(1,s.length);
        return ltr.toUpperCase() + rest;
    } else {//empty string
        return "";
    }
}

// returns true if string s is in list (of strings)
function checkStrInList(s, list) {
    if (list != null) {
	for (var i = 0; i < list.length; i++) {
	    if (s == list[i]) {
		return true;
	    }
	}
    }
    return false;
}
