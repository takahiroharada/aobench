#include<metal_stdlib>
using namespace metal;

#define WIDTH        256
#define HEIGHT       256
#define NSUBSAMPLES  2
#define NAO_SAMPLES  8

#define M_PI 3.141526f


#define R_ARGS rSeeds, rStep
#define R_ARG_LIST float2 rSeeds, thread int* rStep

float draw(float2 seed, thread int* ss)
{
	float c = 0.01f;
	seed.x*=cos(ss[0]*c);
	seed.y*=sin(ss[0]*c);
	ss[0]++;
	return fract(sin(seed.x*12.9898 + seed.y*78.233) * 43758.5453);
}

#define drand48() draw(R_ARGS)

typedef struct
{
	float x;
	float y;
	float z;
} vec3;


typedef struct
{
	float t;
	vec3    p;
	vec3    n;
	int    hit;
} Isect;

typedef struct
{
	vec3    center;
	float radius;
	
} Sphere;

typedef struct
{
	vec3    p;
	vec3    n;
	
} Plane;

typedef struct
{
	vec3    org;
	vec3    dir;
} Ray;


float vdot(vec3 v0, vec3 v1)
{
	return v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
}

void vcross(thread vec3 *c, vec3 v0, vec3 v1)
{
	
	c->x = v0.y * v1.z - v0.z * v1.y;
	c->y = v0.z * v1.x - v0.x * v1.z;
	c->z = v0.x * v1.y - v0.y * v1.x;
}

void vnormalize(thread vec3 *c)
{
	float length = sqrt(vdot((*c), (*c)));
	
	if (fabs(length) > 1.0e-17) {
		c->x /= length;
		c->y /= length;
		c->z /= length;
	}
}

void
ray_sphere_intersect(thread Isect *isect, const thread Ray *ray, const thread Sphere *sphere)
{
	vec3 rs;
	
	rs.x = ray->org.x - sphere->center.x;
	rs.y = ray->org.y - sphere->center.y;
	rs.z = ray->org.z - sphere->center.z;
	
	float B = vdot(rs, ray->dir);
	float C = vdot(rs, rs) - sphere->radius * sphere->radius;
	float D = B * B - C;
	
	if (D > 0.0) {
		float t = -B - sqrt(D);
		
		if ((t > 0.0) && (t < isect->t)) {
			isect->t = t;
			isect->hit = 1;
			
			isect->p.x = ray->org.x + ray->dir.x * t;
			isect->p.y = ray->org.y + ray->dir.y * t;
			isect->p.z = ray->org.z + ray->dir.z * t;
			
			isect->n.x = isect->p.x - sphere->center.x;
			isect->n.y = isect->p.y - sphere->center.y;
			isect->n.z = isect->p.z - sphere->center.z;
			
			vnormalize(&(isect->n));
		}
	}
}

void
ray_plane_intersect(thread Isect *isect, const thread Ray *ray, const thread Plane *plane)
{
	float d = -vdot(plane->p, plane->n);
	float v = vdot(ray->dir, plane->n);
	
	if (fabs(v) < 1.0e-17) return;
	
	float t = -(vdot(ray->org, plane->n) + d) / v;
	
	if ((t > 0.0) && (t < isect->t)) {
		isect->t = t;
		isect->hit = 1;
		
		isect->p.x = ray->org.x + ray->dir.x * t;
		isect->p.y = ray->org.y + ray->dir.y * t;
		isect->p.z = ray->org.z + ray->dir.z * t;
		
		isect->n = plane->n;
	}
}

void
orthoBasis(thread vec3 *basis, vec3 n)
{
	basis[2] = n;
	basis[1].x = 0.0; basis[1].y = 0.0; basis[1].z = 0.0;
	
	if ((n.x < 0.6) && (n.x > -0.6)) {
		basis[1].x = 1.0;
	} else if ((n.y < 0.6) && (n.y > -0.6)) {
		basis[1].y = 1.0;
	} else if ((n.z < 0.6) && (n.z > -0.6)) {
		basis[1].z = 1.0;
	} else {
		basis[1].x = 1.0;
	}
	
	vcross(&basis[0], basis[1], basis[2]);
	vnormalize(&basis[0]);
	
	vcross(&basis[1], basis[2], basis[0]);
	vnormalize(&basis[1]);
}


void ambient_occlusion(thread Sphere* spheres, Plane plane, thread vec3 *col, const thread Isect *isect, R_ARG_LIST)
{
	int    i, j;
	int    ntheta = NAO_SAMPLES;
	int    nphi   = NAO_SAMPLES;
	float eps = 0.0001;
	
	vec3 p;
	
	p.x = isect->p.x + eps * isect->n.x;
	p.y = isect->p.y + eps * isect->n.y;
	p.z = isect->p.z + eps * isect->n.z;
	
	vec3 basis[3];
	orthoBasis(basis, isect->n);
	
	float occlusion = 0.0;
	
	for (j = 0; j < ntheta; j++) {
		for (i = 0; i < nphi; i++) {
			float theta = sqrt(drand48());
			float phi   = 2.0 * M_PI * drand48();
			
			float x = cos(phi) * theta;
			float y = sin(phi) * theta;
			float z = sqrt(1.0 - theta * theta);
			
			// local -> global
			float rx = x * basis[0].x + y * basis[1].x + z * basis[2].x;
			float ry = x * basis[0].y + y * basis[1].y + z * basis[2].y;
			float rz = x * basis[0].z + y * basis[1].z + z * basis[2].z;
			
			Ray ray;
			
			ray.org = p;
			ray.dir.x = rx;
			ray.dir.y = ry;
			ray.dir.z = rz;
			
			Isect occIsect;
			occIsect.t   = 1.0e+17;
			occIsect.hit = 0;
			
			ray_sphere_intersect(&occIsect, &ray, &spheres[0]);
			ray_sphere_intersect(&occIsect, &ray, &spheres[1]);
			ray_sphere_intersect(&occIsect, &ray, &spheres[2]);
			ray_plane_intersect (&occIsect, &ray, &plane);
			
			if (occIsect.hit) occlusion += 1.0;
		}
	}
	
	occlusion = (ntheta * nphi - occlusion) / (float)(ntheta * nphi);
	
	col->x = occlusion;
	col->y = occlusion;
	col->z = occlusion;
}

kernel
void AoKernel(device float* Out[[buffer(0)]],
			  uint2 tid [[thread_position_in_grid]])
{
	device float* fimg = Out;
	const int h = HEIGHT;
	const int w = WIDTH;
	const int nsubsamples = NSUBSAMPLES;

	const int x = tid.x;
	const int y = tid.y;

	Sphere spheres[3];
	Plane  plane;
	
	spheres[0].center.x = -2.0;
	spheres[0].center.y =  0.0;
	spheres[0].center.z = -3.5;
	spheres[0].radius = 0.5;
	
	spheres[1].center.x = -0.5;
	spheres[1].center.y =  0.0;
	spheres[1].center.z = -3.0;
	spheres[1].radius = 0.5;
	
	spheres[2].center.x =  1.0;
	spheres[2].center.y =  0.0;
	spheres[2].center.z = -2.2;
	spheres[2].radius = 0.5;
	
	plane.p.x = 0.0;
	plane.p.y = -0.5;
	plane.p.z = 0.0;
	
	plane.n.x = 0.0;
	plane.n.y = 1.0;
	plane.n.z = 0.0;
	

	float2 rSeed = (float2)(x/(float)w, y/(float)h);
	int rStep = 0;
	
	for (int v = 0; v < nsubsamples; v++) {
		for (int u = 0; u < nsubsamples; u++) {
			float px = (x + (u / (float)nsubsamples) - (w / 2.0)) / (w / 2.0);
			float py = -(y + (v / (float)nsubsamples) - (h / 2.0)) / (h / 2.0);
			
			Ray ray;
			
			ray.org.x = 0.0;
			ray.org.y = 0.0;
			ray.org.z = 0.0;
			
			ray.dir.x = px;
			ray.dir.y = py;
			ray.dir.z = -1.0;
			vnormalize(&(ray.dir));
			
			Isect isect;
			isect.t   = 1.0e+17;
			isect.hit = 0;
			
			ray_sphere_intersect(&isect, &ray, &spheres[0]);
			ray_sphere_intersect(&isect, &ray, &spheres[1]);
			ray_sphere_intersect(&isect, &ray, &spheres[2]);
			ray_plane_intersect (&isect, &ray, &plane);
			
			if (isect.hit)
			{
				vec3 col;
				ambient_occlusion(spheres, plane, &col, &isect, rSeed, &rStep);
				
				fimg[3 * (y * w + x) + 0] += col.x;
				fimg[3 * (y * w + x) + 1] += col.y;
				fimg[3 * (y * w + x) + 2] += col.z;
			}
		}
	}
	
	fimg[3 * (y * w + x) + 0] /= (float)(nsubsamples * nsubsamples);
	fimg[3 * (y * w + x) + 1] /= (float)(nsubsamples * nsubsamples);
	fimg[3 * (y * w + x) + 2] /= (float)(nsubsamples * nsubsamples);
	
}




