#!/usr/bin/env k8

var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

var _hic_re_cigar = /(\d+)([MIDSH])/g;

function _hic_resolve_frag(opt, a)
{
	if (a.length < 2) return;

	// test single- or paired-end
	var is_pe = false, is_se = false, qname = a[0][0];
	for (var i = 0; i < a.length; ++i) {
		if (a[i][1] & 0x1) is_pe = true;
		else is_se = true;
	}
	if (is_pe && is_se)
		throw Error(qname + " shouldn't be both PE and SE");

	var b = [], read_nums = 0;
	for (var i = 0; i < a.length; ++i) {
		var t = a[i];
		var mapq = parseInt(t[4]);
		if (mapq < opt.min_mapq) continue; // skip hits with low mapping quality
		var m, clip = [0, 0], x = 0, y = 0;
		var flag = t[1];
		var rev = flag&16? true : false;
		while ((m = _hic_re_cigar.exec(t[5])) != null) {
			var op = m[2], len = parseInt(m[1]);
			if (op == 'M') {
				x += len, y += len;
			} else if (op == 'I') {
				y += len;
			} else if (op == 'D') {
				x += len;
			} else if (op == 'S' || op == 'H') {
				clip[y == 0? 0 : 1] = len;
				y += len;
			}
		}
		var rs = parseInt(t[3]), re = rs + x;
		var qs, qe;
		if (!rev) qs = clip[0], qe = y - clip[1];
		else qs = clip[1], qe = y - clip[0];
		var read_num = flag >> 6 & 0x3;
		if (read_num == 3)
			throw Error(t[0] + ": incorrect read number flags");
		read_nums |= 1 << read_num;
		if (read_num == 2) { // NB: qlen1 could be zero here
			var qs1 = y - qe;
			qe = y - qs, qs = qs1;
			rev = !rev;
		}
		b.push([read_num, rev, y, qs, qe, t[2], rs, re, mapq]);
	}
	if (b.length < 2) return;

	b.sort(function(x,y) { return x[0] != y[0]? x[0] - y[0] : x[3] - y[3] });

	// identify segments
	var segs = [[b[0][5], b[0][6], b[0][7], b[0][1]? '-' : '+', 0, b[0][8], 1]];
	for (var i = 1; i < b.length; ++i) {
		var p = b[i - 1], q = b[i];
		var new_seg = true;
		if (p[1] == q[1] && p[5] == q[5]) { // same strand and chr
			var dist = !p[1]? q[6] - p[7] : q[7] - p[6];
			if (dist < opt.min_dist && dist >= -opt.min_dist) new_seg = false;
		}
		if (new_seg) {
			segs.push([q[5], q[6], q[7], q[1]? '-' : '+', 0, q[8], 1]);
		} else {
			var last = segs[segs.length - 1];
			if (p[1]) last[1] = q[6];
			else last[2] = q[7];
			last[5] = last[5] > q[8]? last[5] : q[8];
			++last[6];
		}
	}
	if (segs.length < 2) return;

	function ovlp_ratio(s1, s2) {
		var min_st = s1[1] < s2[1]? s1[1] : s2[1];
		var max_st = s1[1] > s2[1]? s1[1] : s2[1];
		var min_en = s1[2] < s2[2]? s1[2] : s2[2];
		var max_en = s1[2] > s2[2]? s1[2] : s2[2];
		return max_st >= min_en? 0 : (min_en - max_st) / (s1[2] - s1[1] < s2[2] - s2[1]? s1[2] - s1[1] : s2[2] - s2[1]);
	}

	// fixed unmerged mates, a corner case
	if (segs.length == 4 && segs[0][0] == segs[2][0] && segs[1][0] == segs[3][0] && segs[0][3] == segs[2][3] && segs[1][3] == segs[3][3]) {
		if (ovlp_ratio(segs[0], segs[2]) > 0.9 && ovlp_ratio(segs[1], segs[3]) > 0.9)
			segs.length = 2;
	}

	// debugging
	if (opt.verbose >= 4) {
		for (var i = 0; i < b.length; ++i)
			print(qname, b[i].join("\t"));
	}

	if (opt.fmt_pairs) {
		for (var i = 0; i < segs.length - 1; ++i)
			print(qname, segs[i][0], segs[i][2], segs[i+1][0], segs[i+1][1], segs[i][3], segs[i+1][3]);
	} else {
		var out = [];
		for (var i = 0; i < segs.length; ++i)
			out.push(segs[i].join(":"));
		print(qname, out.join("\t"));
	}
}

function hic_sam2seg(args)
{
	var c, opt = { min_mapq:20, min_dist:500, fmt_pairs:false, verbose:3 };
	while ((c = getopt(args, "q:v:d:p")) != null) {
		if (c == 'q') opt.min_mapq = parseInt(getopt.arg);
		else if (c == 'v') opt.verbose = parseInt(getopt.arg);
		else if (c == 'd') opt.min_dist = parseInt(getopt.arg);
		else if (c == 'p') opt.fmt_pairs = true;
	}

	if (args.length - getopt.ind == 0) {
		print("Usage: hickit.js sam2seg [options] <in.sam>");
		print("Options:");
		print("  -q INT     min mapping quality [" + opt.min_mapq + "]");
		print("  -d INT     min distance between segments [" + opt.min_dist + "]");
		print("  -p         output .pairs format (segments by default)");
		return 1;
	}

	var buf = new Bytes();
	var file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
	var a = [];
	if (opt.fmt_pairs)
		print("## pairs format v1.0");
	var hdr_done = false;
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (t[0][0] == '@') {
			if (t[0] == '@SQ') {
				var sn = null, ln = null;
				for (var i = 1; i < t.length; ++i) {
					var m;
					if ((m = /(LN|SN):(\S+)/.exec(t[i])) != null) {
						if (m[1] == 'SN') sn = m[2];
						else if (m[1] == 'LN') ln = parseInt(m[2]);
					}
				}
				if (sn == null || ln == null)
					throw Error("missing SN or LN at an @SQ line");
				if (opt.fmt_pairs)
					print("#chromosome: " + sn + " " + ln);
			}
			continue;
		}
		if (!hdr_done) {
			if (opt.fmt_pairs)
				print("#columns: readID chr1 pos1 chr2 pos2 strand1 strand2");
			hdr_done = true;
		}
		t[1] = parseInt(t[1]);
		if (a.length > 0 && a[0][0] != t[0]) {
			_hic_resolve_frag(opt, a);
			a.length = 0;
		}
		if ((t[1] & 0x4) == 0)
			a.push(t);
	}
	_hic_resolve_frag(opt, a);
	file.close();
	buf.destroy();
	return 0;
}

function main(args)
{
	if (args.length == 0) {
		print("Usage: hickit.js <command> [arguments]");
		print("Commands:");
		print("  sam2seg        convert SAM to segmentations/pairs");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'sam2seg') hic_sam2seg(args);
	else throw Error("unrecognized command: " + cmd);
}
main(arguments);

