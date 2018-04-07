#!/usr/bin/env k8

/*****************************
 ***** Library functions *****
 *****************************/

/*******************************
 * Command line option parsing *
 *******************************/

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

/***********************
 * Interval operations *
 ***********************/

Interval = {};

Interval.sort = function(a)
{
	if (typeof a[0] == 'number')
		a.sort(function(x, y) { return x - y });
	else a.sort(function(x, y) { return x[0] != y[0]? x[0] - y[0] : x[1] - y[1] });
}

Interval.merge = function(a, sorted)
{
	if (typeof sorted == 'undefined') sorted = true;
	if (!sorted) Interval.sort(a);
	var k = 0;
	for (var i = 1; i < a.length; ++i) {
		if (a[k][1] >= a[i][0])
			a[k][1] = a[k][1] > a[i][1]? a[k][1] : a[i][1];
		else a[++k] = a[i].slice(0);
	}
	a.length = k + 1;
}

Interval.index_end = function(a, sorted)
{
	if (a.length == 0) return;
	if (typeof sorted == 'undefined') sorted = true;
	if (!sorted) Interval.sort(a);
	a[0].push(0);
	var k = 0, k_en = a[0][1];
	for (var i = 1; i < a.length; ++i) {
		if (k_en <= a[i][0]) {
			for (++k; k < i; ++k)
				if (a[k][1] > a[i][0])
					break;
			k_en = a[k][1];
		}
		a[i].push(k);
	}
}

Interval.find_intv = function(a, x)
{
	var left = -1, right = a.length;
	if (typeof a[0] == 'number') {
		while (right - left > 1) {
			var mid = left + ((right - left) >> 1);
			if (a[mid] > x) right = mid;
			else if (a[mid] < x) left = mid;
			else return mid;
		}
	} else {
		while (right - left > 1) {
			var mid = left + ((right - left) >> 1);
			if (a[mid][0] > x) right = mid;
			else if (a[mid][0] < x) left = mid;
			else return mid;
		}
	}
	return left;
}

Interval.find_ovlp = function(a, st, en)
{
	if (a.length == 0 || st >= en) return [];
	var l = Interval.find_intv(a, st);
	var k = l < 0? 0 : a[l][a[l].length - 1];
	var b = [];
	for (var i = k; i < a.length; ++i) {
		if (a[i][0] >= en) break;
		else if (st < a[i][1])
			b.push(a[i]);
	}
	return b;
}

/************************
 *** hickit functions ***
 ************************/

var _hic_re_cigar = /(\d+)([MIDSH])/g;

function hic_vcf2tsv(args)
{
	var c, vcf_no_chr = false, inc_indel = false;
	while ((c = getopt(args, 'cg')) != null) {
		if (c == 'c') vcf_no_chr = true;
		else if (c == 'g') inc_indel = true;
	}
	if (args.length - getopt.ind == 0) {
		print("Usage: hickit.js vcf2tsv [options] <in.vcf>");
		print("Options:");
		print("  -c     convert 'chr1' in VCF to '1'");
		print("  -g     include INDELs");
		exit(1);
	}
	var file = new File(args[getopt.ind]);
	var buf = new Bytes();
	while (file.readline(buf) >= 0) {
		var m, t = buf.toString().split("\t");
		if (t[0][0] == '#') continue;
		if (t[6] != '.' && t[6] != 'PASS') continue;
		var s = t[4].split(",");
		s.unshift(t[3]);
		var max = 0;
		for (var i = 0; i < s.length; ++i)
			max = max > s[i].length? max : s[i].length;
		if (!inc_indel && max > 1) continue;
		if ((m = /^(\d+)\|(\d+)/.exec(t[9])) == null) continue;
		var a1 = parseInt(m[1]), a2 = parseInt(m[2]);
		if (a1 == a2) continue;
		if (a1 >= s.length || a2 >= s.length) throw Error("incorrect VCF");
		var chr = vcf_no_chr? t[0].replace(/^chr/, "") : t[0];
		var pos = parseInt(t[1]) - 1;
		print(chr, pos + 1, s[a1], s[a2]);
	}
	buf.destroy();
	file.close();
}

function _hic_resolve_frag(opt, a)
{
	if (a.length < 2) return;

	// test single- or paired-end
	var is_pe = false, is_se = false, qname = opt.no_qname? '.' : a[0][0];
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
		var rs = parseInt(t[3]) - 1, re = rs + x;
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
	var segs = [[b[0][5], b[0][6], b[0][7], b[0][1]? '-' : '+', '.', b[0][8], 1]];
	for (var i = 1; i < b.length; ++i) {
		var p = b[i - 1], q = b[i];
		var new_seg = true;
		if (p[1] == q[1] && p[5] == q[5]) { // same strand and chr
			var dist = !p[1]? q[6] - p[7] : q[7] - p[6];
			if (dist < opt.min_dist && dist >= -opt.min_dist) new_seg = false;
		}
		if (new_seg) {
			segs.push([q[5], q[6], q[7], q[1]? '-' : '+', '.', q[8], 1]);
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
		if (ovlp_ratio(segs[0], segs[2]) > 0.5 && ovlp_ratio(segs[1], segs[3]) > 0.5)
			segs.length = 2;
	}

	// debugging
	if (opt.verbose >= 4) {
		for (var i = 0; i < b.length; ++i)
			print(qname, b[i].join("\t"));
	}

	// print out
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
	var c, opt = { min_mapq:20, min_dist:500, fmt_pairs:false, no_qname:false, verbose:3, fn_var:null };
	while ((c = getopt(args, "q:V:d:pNv:")) != null) {
		if (c == 'q') opt.min_mapq = parseInt(getopt.arg);
		else if (c == 'V') opt.verbose = parseInt(getopt.arg);
		else if (c == 'd') opt.min_dist = parseInt(getopt.arg);
		else if (c == 'p') opt.fmt_pairs = true;
		else if (c == 'N') opt.no_qname = true;
		else if (c == 'v') opt.fn_var = getopt.arg;
	}

	if (args.length - getopt.ind == 0) {
		print("Usage: hickit.js sam2seg [options] <in.sam>");
		print("Options:");
		print("  -q INT     min mapping quality [" + opt.min_mapq + "]");
		print("  -d INT     min distance between segments [" + opt.min_dist + "]");
		print("  -p         output .pairs format (segments by default)");
		print("  -N         don't print fragment name");
		return 1;
	}

	var buf = new Bytes();

	if (opt.fn_var != null) {
		opt.snp = {};
		var file = new File(opt.fn_var);
		while (file.readline(buf) >= 0) {
			var t = buf.toString().split("\t");
			if (t[2].length != 1 || t[3].length != 1) continue;
			if (opt.snp[t[0]] == null) opt.snp[t[0]] = [];
			var pos = parseInt(t[1]);
			opt.snp[t[0]].push([pos - 1, pos, t[2], t[3]]);
		}
		file.close();
		for (var c in opt.snp)
			Interval.index_end(opt.snp[c], true);
	}

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
		print("  vcf2tsv        convert phased VCF to simple TSV (chr, pos1, al1, al2)");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'sam2seg') hic_sam2seg(args);
	else if (cmd == 'vcf2tsv') hic_vcf2tsv(args);
	else throw Error("unrecognized command: " + cmd);
}
main(arguments);

