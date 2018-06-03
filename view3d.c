#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hkpriv.h"

#ifdef HAVE_GL
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <gl/glut.h>
#endif

struct {
	struct hk_bmap *m;
	fvec3_t *color;
	float *alpha;
	krng_t rng;
	int line_width;
	int n_hl;
	int is_white;
	int use_ball;
	char **hl;
	float feat_lo, feat_hi;
	int feat_color;
} global;

static void prep_color(const struct hk_sdict *d);

static inline void set_feat_color(float x, float alpha)
{
	if (x < global.feat_lo) x = global.feat_lo;
	if (x > global.feat_hi) x = global.feat_hi;
	x = (x - global.feat_lo) / (global.feat_hi - global.feat_lo);
	glColor4f(1.0f - x, x, 0, alpha);
}

static void cb_draw(void)
{
	static const double X = .525731112119133606, Z = .850650808352039932;
	static GLfloat vdata[12][3] = {
		{-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},
		{0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},
		{Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}
	};
	static GLuint tindices[20][3] = {
		{1,4,0}, {4,9,0}, {4,5,9}, {8,5,4}, {1,8,4},
		{1,10,8}, {10,3,8}, {8,3,5}, {3,2,5}, {3,7,2},
		{3,10,7}, {10,6,7}, {6,11,7}, {6,0,11}, {6,1,0},
		{10,1,6}, {11,0,9}, {2,11,9}, {5,2,9}, {11,2,7}
	};

	struct hk_bmap *m = global.m;
	int32_t i;

	if (global.is_white) glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	else glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLineWidth(global.line_width);
	for (i = 0; i < m->d->n; ++i) {
		int32_t j, cnt = (int32_t)m->offcnt[i];
		int32_t off = m->offcnt[i] >> 32;
		if (global.alpha[i] < 0.9) continue;
		glColor4f(global.color[i][0], global.color[i][1], global.color[i][2], global.alpha[i]);
		glBegin(GL_LINE_STRIP);
		for (j = 0; j < cnt; ++j) {
			if (m->feat && global.feat_color) set_feat_color(m->feat[off+j], global.alpha[i]);
			glVertex3f(m->x[off+j][0], m->x[off+j][1], m->x[off+j][2]);
		}
		glEnd();
		if (global.use_ball) {
			for (j = 0; j < cnt; ++j) {
				int k;
				if (m->feat && global.feat_color) set_feat_color(m->feat[off+j], global.alpha[i]);
				glPushMatrix();
				glTranslatef(m->x[off+j][0], m->x[off+j][1], m->x[off+j][2]);
				glScalef(0.01, 0.01, 0.01);
#if 1
				glBegin(GL_TRIANGLES);
				for (k = 0; k < 20; ++k) {
					glVertex3fv(&vdata[tindices[k][0]][0]);
					glVertex3fv(&vdata[tindices[k][1]][0]);
					glVertex3fv(&vdata[tindices[k][2]][0]);
				}
				glEnd();
#else
				glutSolidSphere(1, 5, 4);
#endif
				glPopMatrix();
			}
		}
	}
	for (i = 0; i < m->d->n; ++i) {
		int32_t j, cnt = (int32_t)m->offcnt[i];
		int32_t off = m->offcnt[i] >> 32;
		if (global.alpha[i] > 0.9) continue;
		glColor4f(global.color[i][0], global.color[i][1], global.color[i][2], global.alpha[i]);
		glBegin(GL_LINE_STRIP);
		for (j = 0; j < cnt; ++j) {
			if (m->feat && global.feat_color) set_feat_color(m->feat[off+j], global.alpha[i]);
			glVertex3f(m->x[off+j][0], m->x[off+j][1], m->x[off+j][2]);
		}
		glEnd();
	}
	glFlush();
	glutSwapBuffers();
}

static void rotate(int axis, float x)
{
	const float step = 5.0f;
	double mat_modview[16], mat_proj[16];
	double xx, yy, zz, s;
	int view_point[4];
	glGetDoublev(GL_MODELVIEW_MATRIX, mat_modview);
	glGetDoublev(GL_PROJECTION_MATRIX, mat_proj);
	glGetIntegerv(GL_VIEWPORT, view_point);
	gluProject(0.0, 0.0, 0.0, mat_modview, mat_proj, view_point, &xx, &yy, &zz);
	s = step * x;
	if (axis == 0)      gluUnProject(xx + 10.0, yy, zz, mat_modview, mat_proj, view_point, &xx, &yy, &zz);
	else if (axis == 1) gluUnProject(xx, yy + 10.0, zz, mat_modview, mat_proj, view_point, &xx, &yy, &zz);
	else if (axis == 2) gluUnProject(xx, yy, zz + 10.0, mat_modview, mat_proj, view_point, &xx, &yy, &zz);
	glMatrixMode(GL_MODELVIEW);
	glRotatef(s, xx, yy, zz);
	cb_draw();
}

static void cb_key(unsigned char key, int x, int y)
{
	if (key == 27 || key == 'q' || key == 'Q') {
		hk_bmap_destroy(global.m);
		free(global.color);
		exit(0);
	} else if (key == 'c' || key == 'C') {
		prep_color(global.m->d);
		cb_draw();
	} else if (key == 'f' || key == 'F') {
		global.feat_color = !global.feat_color;
		cb_draw();
	} else if (key == '>' || key == '.') {
		glMatrixMode(GL_PROJECTION);
		gluLookAt(0, 0, 0.05,  0, 0, -1.0,  0, 1, 0);
		cb_draw();
	} else if (key == '<' || key == ',') {
		glMatrixMode(GL_PROJECTION);
		gluLookAt(0, 0, -0.05,  0, 0, -1.0,  0, 1, 0);
		cb_draw();
	} else if (key == '0') {
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		cb_draw();
	} else if (key == '[') {
		rotate(2, 1.0f);
	} else if (key == ']') {
		rotate(2, -1.0f);
	} else if (key == 'b' || key == 'B') {
		global.use_ball = !global.use_ball;
		cb_draw();
	} else if (key == 'z' || key == 'Z') {
		global.is_white = !global.is_white;
		cb_draw();
	}
}

static void cb_special_key(int key, int x, int y)
{
	if (key == GLUT_KEY_UP)         rotate(0, 1.0f);
	else if (key == GLUT_KEY_DOWN)  rotate(0, -1.0f);
	else if (key == GLUT_KEY_LEFT)  rotate(1, 1.0f);
	else if (key == GLUT_KEY_RIGHT) rotate(1, -1.0f);
}

void hk_v3d_prep(int *argc, char *argv[])
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_RGB|GLUT_DEPTH|GLUT_ALPHA|GLUT_DOUBLE);
}

static void prep_color(const struct hk_sdict *d)
{
	int i, j;
	global.color = CALLOC(fvec3_t, d->n);
	global.alpha = CALLOC(float, d->n);
	for (i = 0; i < d->n; ++i) {
		global.color[i][0] = kr_drand_r(&global.rng);
		global.color[i][1] = kr_drand_r(&global.rng);
		global.color[i][2] = kr_drand_r(&global.rng);
		if (global.n_hl > 0) {
			global.alpha[i] = 0.15f;
			for (j = 0; j < global.n_hl; ++j)
				if (strcmp(d->name[i], global.hl[j]) == 0)
					global.alpha[i] = 1.0f;
		} else global.alpha[i] = 1.0f;
	}
}

static void split_hl(const char *hl)
{
	const char *p, *q;
	int k;
	for (q = hl, k = 0; *q; ++q)
		if (*q == ',') ++k;
	global.n_hl = k + 1;
	global.hl = CALLOC(char*, global.n_hl);
	for (p = q = hl, k = 0;; ++q) {
		if (*q == 0 || *q == ',') {
			global.hl[k] = CALLOC(char, q - p + 1);
			memcpy(global.hl[k], p, q - p);
			++k;
			p = q + 1;
			if (*q == 0) break;
		}
	}
}

void hk_v3d_view(struct hk_bmap *m, int width, int line_width, int color_seed, const char *hl)
{
	memset(&global, 0, sizeof(global));
	global.m = m;
	global.line_width = line_width;
	global.is_white = 1;
	kr_srand_r(&global.rng, color_seed);
	hk_fdg_normalize(m);
	if (hl) split_hl(hl);
	prep_color(m->d);
	if (m->feat) {
		float *a;
		int i;
		a = CALLOC(float, m->n_beads);
		for (i = 0; i < m->n_beads; ++i)
			a[i] = m->feat[i];
		global.feat_lo = ks_ksmall_float(m->n_beads, a, (int)(m->n_beads * 0.05f));
		global.feat_hi = ks_ksmall_float(m->n_beads, a, (int)(m->n_beads * 0.95f));
		if (hk_verbose >= 3)
			fprintf(stderr, "[M::%s] low = %f, high = %f\n", __func__, global.feat_lo, global.feat_hi);
		global.feat_color = 1;
	}

	glutInitWindowSize(width, width);
	glutCreateWindow("View3D");
	glutDisplayFunc(cb_draw);
	glutKeyboardFunc(cb_key);
	glutSpecialFunc(cb_special_key);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);

	glutMainLoop();
}
#endif
