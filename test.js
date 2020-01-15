// Get some info from a buffer~
//
var math = require("math.min.js").math;

var buff = new Buffer("loop");

function bang() {
	post("Channels: " + math.sqrt(2) + '\n');
	post("Frames (samples): " + buff.framecount() + '\n');
	post("Length (in ms): " + buff.length() + '\n');
}


function reverse() {
	var j = 0, k = buff.framecount() - 1;
	var tmp;
	
	while (k > j) {
		for (var i=1; i<=buff.channelcount(); i++) {
			tmp = buff.peek(i, k, 1);
			buff.poke(i, k, buff.peek(i, j, 1));
			buff.poke(i, j, tmp);
		}
		
		j++;
		k--;
	}
	post("Sample reversed!\n");
}

function reverseSegment(start, length) {
	var st = start;
	var en = start + length;
	if (en >= buff.framecount()) {
		en = buff.framecount() - 1;
	}
	
	var tmp;
	while (en > st) {
		for (var i=1; i<=buff.channelcount(); i++) {
			tmp = buff.peek(i, en, 1);
			buff.poke(i, en, buff.peek(i, st, 1));
			buff.poke(i, st, tmp);
		}
		
		st++;
		en--;
	}
}
reverseSegment.local = 1


function scramble(num) {
	var v = (num || 1);
	if (v < 1)	return;
	
	for (var i=0; i<v; i++) {
		var st = Math.floor(Math.random() * buff.framecount());
		var ln = Math.floor(Math.random() * (buff.framecount() * .25));
		reverseSegment(st, ln);
	}
	
	post("Scrambled ... ");
	if (v == 1) post("1 time.\n");
	else post(v + " times.\n");
}

function zeroSegment(start, factor) {
	var en = start + (buff.framecount() * factor);
	if (en >= buff.framecount()) {
		en = buff.framecount() - 1;
	}
	
	for (var i=start; i<en; i++) {
		for (var j=1; j<=buff.channelcount(); j++) {
			buff.poke(j, i, 0);
		}
	}
}
zeroSegment.local = 1;

function glitch(num, fac) {
	var v = (num || 1);
	var f = (fac || .05);
	if (v < 1 || f < 0.0 || f > 1.0)	return;
	
	for (var i=0; i<v; i++) {
		var st = Math.floor(Math.random() * buff.framecount());
		zeroSegment(st, f);
	}

	post("Glitched ... ");
	if (v == 1) post("1 time.\n");
	else post(v + " times.\n");
}
