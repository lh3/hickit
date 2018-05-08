#include <stdlib.h>
#include "hkpriv.h"

#ifdef HAVE_GL
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <gl/glut.h>
#endif

struct {
	struct hk_bmap *m;
} global = { 0 };

static void cb_draw(void)
{
	struct hk_bmap *m = global.m;
	int32_t i;

	srand48(1);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLineWidth(2);
	for (i = 0; i < m->d->n; ++i) {
		int32_t j, cnt = (int32_t)m->offcnt[i];
		int32_t off = m->offcnt[i] >> 32;
		glColor3f(drand48(), drand48(), drand48());
		glBegin(GL_LINE_STRIP);
		for (j = 0; j < cnt; ++j)
			glVertex3f(m->x[off+j][0], m->x[off+j][1], m->x[off+j][2]);
		glEnd();
	}
	glFlush();
}

static void cb_key(unsigned char key, int x, int y)
{
	if (key == 27 || key == 'q' || key == 'Q') {
		hk_bmap_destroy(global.m);
		exit(0);
	}
}

static void cb_special_key(int key, int x, int y)
{
	const float step = 5.0f;
	if (key == GLUT_KEY_UP || key == GLUT_KEY_DOWN || key == GLUT_KEY_RIGHT || key == GLUT_KEY_LEFT) {
		double mat_modview[16], mat_proj[16];
		double xx, yy, zz, s;
		int view_point[4];
		glGetDoublev(GL_MODELVIEW_MATRIX, mat_modview);
		glGetDoublev(GL_PROJECTION_MATRIX, mat_proj);
		glGetIntegerv(GL_VIEWPORT, view_point);
		gluProject(0.0, 0.0, 0.0, mat_modview, mat_proj, view_point, &xx, &yy, &zz);
		s = key == GLUT_KEY_LEFT || key == GLUT_KEY_UP? step : -step;
		if (key == GLUT_KEY_RIGHT || key == GLUT_KEY_LEFT)
			gluUnProject(xx, yy + 10.0, zz, mat_modview, mat_proj, view_point, &xx, &yy, &zz);
		else
			gluUnProject(xx + 10.0, yy, zz, mat_modview, mat_proj, view_point, &xx, &yy, &zz);
		glRotatef(s, xx, yy, zz);
		cb_draw();
	}
}

void hk_v3d_prep(int *argc, char *argv[])
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_RGB|GLUT_DEPTH);
}

void hk_v3d_view(struct hk_bmap *m, int width)
{
	global.m = m;
	hk_fdg_normalize(m);
	glutInitWindowSize(width, width);
	glutCreateWindow("View3D");
	glutDisplayFunc(cb_draw);
	glutKeyboardFunc(cb_key);
	glutSpecialFunc(cb_special_key);

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);

	glutMainLoop();
}
#else
void hk_v3d_prep(int *argc, char *argv[]) {}
void hk_v3d_view(struct hk_bmap *m, int width) {}
#endif
