#include <getopt.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "hickit.h"

#define HICKIT_VERSION "r263"

#include <sys/resource.h>
#include <sys/time.h>

double cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

static struct option long_options[] = {
	{ "min-leg-dist",   required_argument, 0, 0 },   // 0
	{ "max-seg",        required_argument, 0, 0 },   // 1
	{ "min-mapq",       required_argument, 0, 0 },   // 2
	{ "out-pairs",      required_argument, 0, 'o' }, // 3
	{ "out-seg",        required_argument, 0, 0 },   // 4
	{ "highlight",      required_argument, 0, 0 },   // 5
	{ "tads",           optional_argument, 0, 0 },   // 6
	{ "bead-radius",    required_argument, 0, 0 },   // 7
	{ "line-width",     required_argument, 0, 0 },   // 8
	{ "dup-dist",       required_argument, 0, 0 },   // 9
	{ "imput-nei",      required_argument, 0, 0 },   // 10
	{ "val-frac",       required_argument, 0, 0 },   // 11
	{ "version",        no_argument,       0, 0 },   // 12
	{ "out-val",        required_argument, 0, 0 },   // 13
	{ "out-png",        required_argument, 0, 0 },   // 14
	{ "png-no-dim",     no_argument,       0, 0 },   // 15
	{ "view",           no_argument,       0, 0 },   // 16
	{ "loops",          optional_argument, 0, 0 },   // 17
	{ "loop-inner",     required_argument, 0, 0 },   // 18
	{ "loop-mid",       required_argument, 0, 0 },   // 19
	{ "loop-outer",     required_argument, 0, 0 },   // 20
	{ "loop-p",         required_argument, 0, 0 },   // 21
	{ 0, 0, 0, 0}
};

static inline int64_t hk_parse_num(const char *str)
{
	double x;
	char *p;
	x = strtod(optarg, &p);
	if (*p == 'G' || *p == 'g') x *= 1e9;
	else if (*p == 'M' || *p == 'm') x *= 1e6;
	else if (*p == 'K' || *p == 'k') x *= 1e3;
	return (int64_t)(x + .499);
}

int main(int argc, char *argv[])
{
	int c, long_idx, has_options = 0;
	FILE *fp;
	struct hk_map *m = 0;
	struct hk_bmap *d3 = 0;
	int n_tads = 0, n_tads_prev = 0, n_loops = 0;
	struct hk_pair *tads = 0, *tads_prev = 0, *loops = 0;
	krng_t rng;

	// general parameters
	int ploidy = 2, seed = 1, max_iter = 1000, radius = 10000000, width = 780;
	float phase_thres = 0.75f;
	// immediate input filters
	int max_seg = 3, min_mapq = 20, min_leg_dist = 1000, dup_dist = 100;
	// TAD calling parameters
	float tad_area_weight = 15.0f, tad_min_cnt_weight = 0.1f;
	// loop calling parameters
	int loop_radius[3] = { 5000, 25000, 55000 };
	float loop_pv = 1e-5;
	// imputation parameters
	int imput_max_nei = 50, imput_min_radius = 50000;
	float imput_val_frac = 0.1f, imput_pseudo_cnt = 0.4f;
	// PNG generation
	int png_no_dim = 0;
	// 3D modeling
	struct hk_fdg_conf fdg_opt;
	// 3D viewing
	struct hk_v3d_opt v3d_opt;
	char *v3d_hl = 0;

	kr_srand_r(&rng, seed);
	hk_fdg_conf_init(&fdg_opt);
	hk_v3d_opt_init(&v3d_opt);

	while ((c = getopt_long(argc, argv, "i:o:r:c:T:P:n:w:p:b:e:k:R:a:s:I:O:D:Suz:L:", long_options, &long_idx)) >= 0) {
		has_options = 1;
		if (c == 'i') {
			if (m) hk_map_destroy(m);
			m = hk_map_read(optarg);
			assert(m);
			if (ploidy == 2) hk_map_phase_male_XY(m);
			if (m->pairs == 0 && m->segs)
				m->pairs = hk_seg2pair(m->n_segs, m->segs, min_leg_dist, max_seg, min_mapq, &m->n_pairs);
			else if (m->pairs && min_leg_dist > 0)
				m->n_pairs = hk_pair_filter_close_legs(m->n_pairs, m->pairs, min_leg_dist);
			if (dup_dist)
				m->n_pairs = hk_pair_dedup(m->n_pairs, m->pairs, dup_dist);
		} else if (c == 'o') {
			fp = strcmp(optarg, "-") == 0? stdout : fopen(optarg, "w");
			assert(fp);
			hk_print_pair(fp, m->cols, m->d, m->n_pairs, m->pairs);
			if (fp != stdout) fclose(fp);
		} else if (c == 'r') {
			radius = hk_parse_num(optarg);
			assert(radius > 0);
		} else if (c == 'n') {
			max_iter = fdg_opt.n_iter = atoi(optarg);
			assert(max_iter > 0);
		} else if (c == 'p') {
			phase_thres = atof(optarg);
			assert(phase_thres > 0.0f && phase_thres <= 1.0f);
		} else if (c == 'P') {
			ploidy = atoi(optarg);
			assert(ploidy == 1 || ploidy == 2);
		} else if (c == 'w') {
			width = atoi(optarg);
			assert(width > 0);
		} else if (c == 'u') {
			hk_impute(m->n_pairs, m->pairs, radius, imput_min_radius, imput_max_nei, max_iter, imput_pseudo_cnt, 1);
			m->cols |= 0x3c;
		} else if (c == 'z') {
			tad_min_cnt_weight = atof(optarg);
		} else if (c == 'T') {
			assert(m && m->pairs);
			if (tads_prev) free(tads_prev);
			if (tads) tads_prev = tads, n_tads_prev = n_tads;
			tads = hk_pair2tad(m->d, m->n_pairs, m->pairs, tad_min_cnt_weight, tad_area_weight, &n_tads);
			m->cols |= 1<<7;
			fp = strcmp(optarg, "-") == 0? stdout : fopen(optarg, "w");
			hk_print_pair(fp, 1<<7, m->d, n_tads, tads);
			if (fp != stdout) fclose(fp);
		} else if (c == 'L') {
			assert(m && m->pairs);
			if (loops) free(loops);
			loops = hk_pair2loop(m->d, m->n_pairs, m->pairs, loop_radius, loop_pv, &n_loops);
			m->cols |= 1<<6 | 1<<9 | 1<<10;
			fp = strcmp(optarg, "-") == 0? stdout : fopen(optarg, "w");
			hk_print_pair(fp, m->cols, m->d, n_loops, loops);
			if (fp != stdout) fclose(fp);
		} else if (c == 'c') {
			int min_flt_cnt;
			min_flt_cnt = atoi(optarg);
			if (min_flt_cnt > 0)
				m->n_pairs = hk_pair_filter_isolated(m->n_pairs, m->pairs, radius, min_flt_cnt, 0.0f);
			else
				hk_pair_count_nei(m->n_pairs, m->pairs, radius, radius);
			m->cols |= 1<<6 | 1<<9;
		} else if (c == 'a') {
			tad_area_weight = atof(optarg);
			assert(tad_area_weight > 0.0f);
		} else if (c == 'I') {
			if (d3) hk_bmap_destroy(d3);
			d3 = hk_3dg_read(optarg);
			assert(d3);
		} else if (c == 'O') {
			fp = strcmp(optarg, "-") == 0? stdout : fopen(optarg, "w");
			assert(fp);
			hk_print_3dg(fp, d3);
			if (fp != stdout) fclose(fp);
		} else if (c == 'D') {
			if (d3) {
				float max_dist;
				max_dist = atof(optarg);
				assert(max_dist > 1.0f);
				m->n_pairs = hk_pair_flt_3d(d3, m->n_pairs, m->pairs, max_dist);
			}
		} else if (c == 'S') {
			if (ploidy == 2) {
				struct hk_map *tmp;
				assert((m->cols & 0x3c) == 0x3c);
				tmp = hk_pair_split_phase(m, phase_thres);
				hk_map_destroy(m);
				m = tmp;
			}
		} else if (c == 'b') {
			int bin_size;
			struct hk_bmap *b;
			bin_size = hk_parse_num(optarg);
			assert(bin_size > 0);
			if (!(m->cols & 1<<6)) hk_pair_count_nei(m->n_pairs, m->pairs, radius, radius);
			b = hk_bmap_gen(m->d, m->n_pairs, m->pairs, bin_size);
			hk_fdg(&fdg_opt, b, d3, &rng);
			if (d3) hk_bmap_destroy(d3);
			d3 = hk_bmap_bead_dup(b);
			hk_bmap_destroy(b);
		} else if (c == 'e') {
			fdg_opt.step = atof(optarg);
			assert(fdg_opt.step > 0.0f && fdg_opt.step < 1.0f);
		} else if (c == 'k') {
			fdg_opt.k_rel_rep = atof(optarg);
			assert(fdg_opt.k_rel_rep > 0.0f);
		} else if (c == 'R') {
			fdg_opt.d_r = atof(optarg);
			assert(fdg_opt.d_r > 0.0f);
		} else if (c == 's') {
			seed = atol(optarg);
			kr_srand_r(&rng, seed);
		} else if (c == 0) {
			if (long_idx == 0) min_leg_dist = hk_parse_num(optarg); // --min-leg-dist
			else if (long_idx ==  1) max_seg = atoi(optarg); // --max-seg
			else if (long_idx ==  2) min_mapq = atoi(optarg); // --min-mapq
			else if (long_idx ==  5) v3d_hl = optarg; // --highlight
			else if (long_idx ==  7) v3d_opt.bead_radius = atof(optarg); // --bead-radius
			else if (long_idx ==  8) v3d_opt.line_width = atof(optarg); // --line-width
			else if (long_idx == 9)  dup_dist = hk_parse_num(optarg); // --dup-dist
			else if (long_idx == 10) imput_max_nei = atoi(optarg); // --imput-nei
			else if (long_idx == 11) imput_val_frac = atof(optarg); // --val-frac
			else if (long_idx == 15) png_no_dim = 1; // --png-no-dim
			else if (long_idx == 18) loop_radius[0] = hk_parse_num(optarg); // --loop-inner
			else if (long_idx == 19) loop_radius[1] = hk_parse_num(optarg); // --loop-mid
			else if (long_idx == 20) loop_radius[2] = hk_parse_num(optarg); // --loop-outer
			else if (long_idx == 21) loop_pv = atof(optarg); // --loop-p
			else if (long_idx ==  4) { // --out-seg
				assert(m && m->segs);
				fp = strcmp(optarg, "-") == 0? stdout : fopen(optarg, "w");
				assert(fp);
				hk_print_seg(fp, m->d, m->n_segs, m->segs);
				if (fp != stdout) fclose(fp);
			} else if (long_idx == 6) { // --tads
				if (tads_prev) free(tads_prev);
				if (tads) tads_prev = tads, n_tads_prev = n_tads;
				if (!optarg) {
					assert(m && m->pairs);
					tads = hk_pair2tad(m->d, m->n_pairs, m->pairs, tad_min_cnt_weight, tad_area_weight, &n_tads);
					m->cols |= 1<<7;
				} else {
					struct hk_map *t;
					t = hk_map_read(optarg);
					assert(t);
					n_tads = t->n_pairs;
					tads = t->pairs;
					t->pairs = 0;
					hk_map_destroy(t);
				}
			} else if (long_idx == 17) { // --loops
				if (loops) free(loops);
				if (!optarg) {
					assert(m && m->pairs);
					loops = hk_pair2loop(m->d, m->n_pairs, m->pairs, loop_radius, loop_pv, &n_loops);
					m->cols |= 1<<6 | 1<<9 | 1<<10;
				}
			} else if (long_idx == 12) { // --version
				printf("%s\n", HICKIT_VERSION);
				exit(0);
			} else if (long_idx == 13) { // --out-val
				fp = strcmp(optarg, "-") == 0? stdout : fopen(optarg, "w");
				assert(fp);
				if (tads == 0) tads = hk_pair2tad(m->d, m->n_pairs, m->pairs, tad_min_cnt_weight, tad_area_weight, &n_tads);
				m->cols |= 1<<7;
				hk_mark_by_tad(n_tads, tads, m->n_pairs, m->pairs);
				hk_validate_holdback(&rng, imput_val_frac, m->n_pairs, m->pairs);
				hk_impute(m->n_pairs, m->pairs, radius, imput_min_radius, imput_max_nei, max_iter, imput_pseudo_cnt, 1);
				hk_validate_roc(fp, m->n_pairs, m->pairs);
				if (fp != stdout) fclose(fp);
			} else if (long_idx == 14) { // --out-png
				hk_pair_image(m->d, m->n_pairs, m->pairs, width, (m->cols&0x3c) == 0x3c? phase_thres : -1.0, png_no_dim,
							  n_tads, tads, n_tads_prev, tads_prev, optarg);
#ifdef HAVE_GL
			} else if (long_idx == 16) { // --view
				int fake_argc = 1;
				char *fake_argv[1];
				fake_argv[0] = "hickit";
				hk_v3d_prep(&fake_argc, fake_argv);
				hk_v3d_view(d3, &v3d_opt, seed, v3d_hl);
#endif
			}
		}
	}

	if (!has_options) {
		FILE *fp = stderr;
		fprintf(fp, "Usage: hickit [options]\n");
		fprintf(fp, "Options:\n");
		fprintf(fp, "  Actions:\n");
		fprintf(fp, "    -i FILE             read .pairs or .seg FILE and apply input filters []\n");
		fprintf(fp, "    -o FILE             write the pairs to FILE []\n");
		fprintf(fp, "    -I FILE             read .3dg FILE []\n");
		fprintf(fp, "    -O FILE             write the 3D model to FILE []\n");
		fprintf(fp, "    -c INT              filter pairs if within -r, #neighbors<INT []\n");
		fprintf(fp, "    -T FILE             call TADs and write to FILE []\n");
		fprintf(fp, "    -L FILE             call loops and write to FILE []\n");
		fprintf(fp, "    -u                  impute missing phases\n");
		fprintf(fp, "    --out-val=FILE      save validation to FILE []\n");
		fprintf(fp, "    --out-png=FILE      write 2D contact map to FILE in PNG []\n");
		fprintf(fp, "    -S                  separate homologous chromosomes\n");
		fprintf(fp, "    -b NUM              3D modeling at NUM resolution []\n");
		fprintf(fp, "    -D FLOAT            filter contacts with large 3D distance []\n");
#ifdef HAVE_GL
		fprintf(fp, "    --view              3D view\n");
#endif
		fprintf(fp, "  General settings:\n");
		fprintf(fp, "    -s INT              random number seed [%d]\n", seed);
		fprintf(fp, "    -P INT              ploidy: 1 or 2 [%d]\n", ploidy);
		fprintf(fp, "    -r NUM              radius [10m]\n");
		fprintf(fp, "    -n INT              max iterations [%d]\n", max_iter);
		fprintf(fp, "    -w INT              image or viewer width [%d]\n", width);
		fprintf(fp, "    -p FLOAT            prob. threshold for a contact considered phased [%g]\n", phase_thres);
		fprintf(fp, "  Input filters:\n");
		fprintf(fp, "    --max-seg=NUM       ignore fragments with >INT segments (.seg only) [%d]\n", max_seg);
		fprintf(fp, "    --min-mapq=NUM      min mapping quality (.seg only) [%d]\n", min_mapq);
		fprintf(fp, "    --min-leg-dist=NUM  min base-pair distance between the two legs [%d]\n", min_leg_dist);
		fprintf(fp, "    --dup-dist=NUM      remove contacts within NUM-bp (dedup; 0 to disable) [%d]\n", dup_dist);
		fprintf(fp, "  TAD calling:\n");
		fprintf(fp, "    -a FLOAT            area weight (larger for smaller TADs) [%g]\n", tad_area_weight);
		fprintf(fp, "    -z INT              min TAD count weight [%g]\n", tad_min_cnt_weight);
		fprintf(fp, "  Loop calling:\n");
		fprintf(fp, "    --loop-inner=NUM    peak radius [%d]\n", loop_radius[0]);
		fprintf(fp, "    --loop-mid=NUM      gap radius [%d]\n", loop_radius[1]);
		fprintf(fp, "    --loop-outer=NUM    outer radius [%d]\n", loop_radius[2]);
		fprintf(fp, "    --loop-p=FLOAT      pseudo P-value threshold [%g]\n", loop_pv);
		fprintf(fp, "  Imputation:\n");
		fprintf(fp, "    --imput-nei=INT     max neighbors [%d]\n", imput_max_nei);
		fprintf(fp, "    --val-frac=FLOAT    fraction to hold out for validation [%g]\n", imput_val_frac);
		fprintf(fp, "  3D modeling:\n");
		fprintf(fp, "    -e FLOAT            step size [%g]\n", fdg_opt.step);
		fprintf(fp, "    -k FLOAT            relative repulsive stiffness [%g]\n", fdg_opt.k_rel_rep);
		fprintf(fp, "    -R FLOAT            relative repulsive radius [%g]\n", fdg_opt.d_r);
#ifdef HAVE_GL
		fprintf(fp, "  3D viewing:\n");
		fprintf(fp, "    --line-width=FLOAT  line width [%g]\n", v3d_opt.line_width);
		fprintf(fp, "    --bead-radius=FLOAT bead radius [%g]\n", v3d_opt.bead_radius);
		fprintf(fp, "    --highlight=STR     comma-delimited list of chromosome names to highlight []\n");
#endif
		fprintf(fp, "\n");
		fprintf(fp, "Examples:\n");
		fprintf(fp, "  hickit -i raw.pairs.gz -U1k -u -o imput.pairs\n");
		fprintf(fp, "  hickit -i imput.pairs -Sr1m -c1 -r10m -c5 -b4m -b1m -b200k -D5 -b50k -D5 -b20k -O out.3dg\n");
		return 1;
	}

	if (tads) free(tads);
	if (d3) hk_bmap_destroy(d3);
	if (m) hk_map_destroy(m);

	if (hk_verbose >= 3) {
		int i;
		fprintf(stderr, "[M::%s] Version: %s\n", __func__, HICKIT_VERSION);
		fprintf(stderr, "[M::%s] CMD:", __func__);
		for (i = 0; i < argc; ++i)
			fprintf(stderr, " %s", argv[i]);
		fputc('\n', stderr);
		fprintf(stderr, "[M::%s] CPU time: %.3f sec\n", __func__, cputime());
	}
	return 0;
}
