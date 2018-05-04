#include <stdlib.h>
#include "hkpriv.h"

#ifdef HAVE_GL
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <gl/glut.h>
#endif

static struct hk_bmap *g_m;

static void cb_draw(void)
{
	struct hk_bmap *m = g_m;
	int32_t i;
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);

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

	glEnd();
	glFlush();
}

void hk_v3d_prep(int *argc, char *argv[])
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_RGB|GLUT_DEPTH);
}

void hk_v3d_view(struct hk_bmap *m, int width)
{
	g_m = m;
	hk_fdg_normalize(m);
	glutInitWindowSize(width, width);
	glutCreateWindow("View3D");
	glutDisplayFunc(cb_draw);
	glutMainLoop();
}
#else
void hk_v3d_prep(int *argc, char *argv[]) {}
void hk_v3d_view(struct hk_bmap *m, int width) {}
#endif
