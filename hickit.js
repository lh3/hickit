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

var hic_sub_delim = '!';
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
		var qs, qe, qlen = y;
		if (!rev) qs = clip[0], qe = qlen - clip[1];
		else qs = clip[1], qe = qlen - clip[0];
		var read_num = flag >> 6 & 0x3;
		if (read_num == 3)
			throw Error(t[0] + ": incorrect read number flags");
		read_nums |= 1 << read_num;
		if (read_num == 2) { // NB: qlen1 could be zero here
			var qs1 = qlen - qe;
			qe = qlen - qs, qs = qs1;
			rev = !rev;
		}
		var phase = '.';
		if (opt.snp != null && opt.snp[t[2]] != null) { // phasing
			var v = Interval.find_ovlp(opt.snp[t[2]], rs, re);
			if (v.length > 0) {
				var x = rs, y = 0;
				var phases = [];
				while ((m = _hic_re_cigar.exec(t[5])) != null) {
					var op = m[2], len = parseInt(m[1]);
					if (op == 'M') {
						for (var j = 0; j < v.length; ++j) {
							var p = v[j][0];
							if (x <= p && p < x + len) {
								p = p - x + y;
								if (p < 0 || p >= t[9].length)
									throw Error("CIGAR parsing error");
								var c = t[9][p];
								var q = t[10].length == t[9].length? t[10].charCodeAt(p) - 33 : opt.min_baseq;
								if (q >= opt.min_baseq) {
									if (c == v[j][2]) phases.push(0);
									else if (c == v[j][3]) phases.push(1);
									else if (opt.verbose >= 2)
										warn('WARNING: a new allele ' + c + ' on read ' + qname + ' at position ' + t[2] + ':' + v[j][1] + ' (not ' + v[j][2] + '/' + v[j][3] + ')');
								}
							}
						}
						x += len, y += len;
					} else if (op == 'I') y += len;
					else if (op == 'D') x += len;
					else if (op == 'S') y += len;
				}
				var n = [0, 0];
				for (var k = 0; k < phases.length; ++k)
					++n[phases[k]];
				if (n[0] > 0 && n[1] == 0) phase = 0;
				else if (n[0] == 0 && n[1] > 0) phase = 1;
				else if (n[0] > 0 && n[1] > 0 && opt.verbose >= 2)
					warn('WARNING: conflicting phase at a segment of read ' + qname);
			}
		}
		b.push([read_num, rev, qlen, qs, qe, t[2], rs, re, mapq, phase]);
	}
	if (b.length < 2) return;

	b.sort(function(x,y) { return x[0] != y[0]? x[0] - y[0] : x[3] - y[3] });

	// identify segments
	var segs = [[b[0][5], b[0][6], b[0][7], b[0][1]? '-' : '+', b[0][9], b[0][8], 1]];
	for (var i = 1; i < b.length; ++i) {
		var p = b[i - 1], q = b[i];
		var new_seg = true;
		if (p[1] == q[1] && p[5] == q[5]) { // same strand and chr
			var dist = !p[1]? q[6] - p[7] : q[7] - p[6];
			if (dist < opt.min_dist && dist >= -opt.min_dist) new_seg = false;
		}
		if (new_seg) {
			segs.push([q[5], q[6], q[7], q[1]? '-' : '+', q[9], q[8], 1]);
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
			out.push(segs[i].join(hic_sub_delim));
		print(qname, out.join("\t"));
	}
}

function hic_sam2seg(args)
{
	var c, opt = { min_mapq:20, min_baseq:20, min_dist:500, fmt_pairs:false, no_qname:false, verbose:3, fn_var:null };
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
		print("  -v FILE    phased SNPs (typically vcf2tsv output)");
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

function hic_chronly(args)
{
	var c, re = new RegExp("^(chr)?([0-9]+|[XY])$");
	while ((c = getopt(args, "r:y")) != null) {
		if (c == 'r') re = new RegExp(getopt.arg);
		else if (c == 'y') re = new RegExp("^(chr)?([0-9]+|X)$");
	}
	if (getopt.ind == args.length) {
		print("Usage: hicket.js chronly [options] <in.pairs>|<in.seg>");
		print("Options:");
		print("  -r STR     regexp to keep [^(chr)?([0-9]+|[XY])]");
		print("  -y         filter out Y chromosome");
		exit(1);
	}
	var buf = new Bytes();
	var file = args[getopt.ind] == '-'? new File() : new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var m, line = buf.toString();
		if ((m = /^#chromosome:\s+(\S+)/.exec(line)) != null) {
			if (re.test(m[1]))
				print(line);
		} else if (line[0] == '#') {
			print(line);
		} else if ((m = /^\S+\t(\S+)\t\d+\t(\S+)\t\d+/.exec(line)) != null) { // .pairs
			if (re.test(m[1]) && re.test(m[2]))
				print(line);
		} else { // .seg
			var t = line.split("\t"), n = 0;
			for (var i = 1; i < t.length; ++i)
				if ((m = /^([^\s!]+)/.exec(t[i])) != null && re.test(m[1]))
					++n;
			if (n == t.length - 1)
				print(line);
		}
	}
	file.close();
	buf.destroy();
}

function hic_bedflt(args)
{
	var c, min_ov_len = 30, min_ov_ratio = 0.5;
	while ((c = getopt(args, "l:r:")) != null) {
		if (c == 'l') min_ov_len = parseInt(getopt.arg);
		else if (c == 'r') min_ov_ratio = parseInt(getopt.arg);
	}
	if (args.length - getopt.ind < 1) {
		print("Usage: hickit.js bedflt [options] <flt.bed> <in.seg>");
		print("Options:");
		print("  -l INT     min overlap length [" + min_ov_len + "]");
		print("  -r FLOAT   min overlap ratio [" + min_ov_ratio + "]");
		exit(1);
	}

	var file, buf = new Bytes();

	var bed = {};
	file = new File(args[getopt.ind]);
	while (file.readline(buf) >= 0) {
		var t = buf.toString().split("\t");
		if (bed[t[0]] == null) bed[t[0]] = [];
		bed[t[0]].push([parseInt(t[1]), parseInt(t[2])]);
	}
	file.close();
	for (var chr in bed) {
		Interval.sort(bed[chr]);
		Interval.merge(bed[chr]);
		Interval.index_end(bed[chr]);
	}

	var re = /\t([^\s!]+)!(\d+)!(\d+)/g;
	file = getopt.ind + 1 < args.length && args[getopt.ind+1] != '-'? new File(args[getopt.ind+1]) : new File();
	while (file.readline(buf) >= 0) {
		var m, flt = false, line = buf.toString();
		if (line[0] == '#') {
			print(line);
			continue;
		}
		while ((m = re.exec(line)) != null) {
			if (bed[m[1]] == null) continue;
			var st = parseInt(m[2]), en = parseInt(m[3]);
			var ov = Interval.find_ovlp(bed[m[1]], st, en);
			var ov_len = 0;
			for (var i = 0; i < ov.length; ++i) {
				var max_st = st > ov[i][0]? st : ov[i][0];
				var min_en = en < ov[i][1]? en : ov[i][1];
				if (min_en < max_st) throw Error("Bug!");
				ov_len += min_en - max_st;
			}
			//print(m[1], m[2], m[3], ov.length, ov_len);
			if (ov_len >= min_ov_len && ov_len >= (en - st) * min_ov_ratio)
				flt = true;
		}
		if (!flt) print(line);
	}
	file.close();

	buf.destroy();
}

var seq_nt4_table = [
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
];

function hic_gfeat(args)
{
	var c, fn_ref = null, conf = { action:1, ref:null, win:100000 };
	while ((c = getopt(args, "c:r:w:")) != null) {
		if (c == 'r') fn_ref = getopt.arg;
		else if (c == 'w') conf.win = atoi(optarg);
	}
	if (getopt.ind == args.length) {
		print("Usage: hickit.js gfeat [options] <in.3dg>");
		print("Options:");
		print("  -r FILE    FASTA of the reference genome []");
		print("  -w INT     window size [" + conf.win + "]");
		exit(1);
	}
	if (conf.action == 1 && fn_ref == null)
		throw Error("option '-r' required to compute CpG density");

	var file, buf = new Bytes();

	if (fn_ref != null) {
		var gt = '>'.charCodeAt(0);
		conf.ref = {};
		file = new File(fn_ref);
		var name = null, seq = null;
		warn("Reading the reference genome...");
		while (file.readline(buf) >= 0) {
			if (buf[0] == gt) {
				var m, s = buf.toString();
				if ((m = />(\S+)/.exec(s)) == null)
					throw Error("malformated FASTA");
				if (name != null) conf.ref[name] = seq;
				name = m[1], seq = new Bytes();
			} else seq.set(buf);
		}
		if (name != null) conf.ref[name] = seq;
		file.close();
		warn("Converting the base encoding...");
		for (var name in conf.ref) {
			var seq = conf.ref[name];
			warn(name, seq.length);
			for (var i = 0; i < seq.length; ++i)
				seq[i] = seq_nt4_table[seq[i]];
		}
	}

	function process_interval(chr, st0, en0, conf) {
		var ref = null;
		if (conf.ref[chr] == null) {
			var alt = chr.replace(/[ab]$/, "").replace(/\(\S+\)$/, "");
			if (conf.ref[alt] == null)
				throw Error("unable to determine the chromosome name '" + alt + "'");
			ref = conf.ref[alt];
		} else ref = conf.ref[chr];
		var st = st0, en = en0 != null? en0 : ref.length;
		if (en - st < conf.win) {
			var x = (conf.win - (en - st)) >> 1;
			st = st - x > 0? st - x : 0;
			en = en + x < ref.length? en + x : ref.length;
		}
		var n_cpg = 0, n_gc = 0, n_tot = 0;
		for (var i = st; i < en; ++i) {
			var c = ref[i];
			if (c < 4) {
				++n_tot;
				if (c == 1 && i < en - 1 && ref[i+1] == 2)
					++n_cpg;
				if (c == 1 || c == 2)
					++n_gc;
			}
		}
		return n_tot? (n_cpg / n_tot).toFixed(6) : 0;
	}

	warn("Processing 3dg...");
	file = new File(args[getopt.ind]);
	var line0 = null, chr0 = null, pos0 = null;
	while (file.readline(buf) >= 0) {
		var line = buf.toString();
		if (line[0] == '#') {
			print(line);
			continue;
		}
		var t = line.split("\t");
		if (t.length < 5) continue;
		t[1] = parseInt(t[1]);
		if (chr0 != null)
			print(line0, process_interval(chr0, pos0, chr0 == t[0]? t[1] : null, conf));
		chr0 = t[0], pos0 = t[1], line0 = line;
	}
	if (line0 != null)
		print(line0, process_interval(chr0, pos0, null, conf));
	file.close();

	if (conf.ref)
		for (var name in conf.ref)
			conf.ref[name].destroy();
	buf.destroy();
}

function main(args)
{
	if (args.length == 0) {
		print("Usage: hickit.js <command> [arguments]");
		print("Commands:");
		print("  sam2seg        convert SAM to segments/pairs");
		print("  vcf2tsv        convert phased VCF to simple TSV (chr, pos1, al1, al2)");
		print("  chronly        filter out non-chromosomal segments/pairs");
		print("  bedflt         filter out segments overlapping a BED");
		print("  gfeat          compute genomic features");
		exit(1);
	}

	var cmd = args.shift();
	if (cmd == 'sam2seg') hic_sam2seg(args);
	else if (cmd == 'vcf2tsv') hic_vcf2tsv(args);
	else if (cmd == 'chronly') hic_chronly(args);
	else if (cmd == 'bedflt') hic_bedflt(args);
	else if (cmd == 'gfeat') hic_gfeat(args);
	else throw Error("unrecognized command: " + cmd);
}

main(arguments);
