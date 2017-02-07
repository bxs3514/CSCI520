/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

  Xiangshun Bei
*/

#include "jello.h"
#include "physics.h"

#define STRUCTURAL_SPRING_LENGTH	0.142857
#define SHEER_SPRING_LENGTH_SHORT	0.202030
#define SHEER_SPRING_LENGTH_LONG	0.247435
#define BEND_SPRING_LENGTH			0.285714

#define FORCEFIELD_INDEX(i, j, k) ((i * jello->resolution * jello->resolution) + (j * jello->resolution) + (k))

typedef enum SpringType
{
	SPRING_STRUCTURAL		=		0x00000000,
	SPRING_SHEER_SHORT		=		0x00000001,
	SPRING_SHEER_LONG		=		0x00000002,
	SPRING_BEND				=		0x00000004,
	SPRING_COLLISION		=		0x00000008,
}SpringType;

typedef enum CollisionType
{
	COLLISION_NONE				=		0x00000000,
	COLLISION_X_POS_BOUND		=		0x00000001,
	COLLISION_X_NAG_BOUND		=		0x00000002,
	COLLISION_Y_POS_BOUND		=		0x00000004,
	COLLISION_Y_NAG_BOUND		=		0x00000008,
	COLLISION_Z_POS_BOUND		=		0x00000010,
	COLLISION_Z_NAG_BOUND		=		0x00000020,
	COLLISION_PLANE_BOUND		=		0x00000040,
}CollisionType;

struct PosVelocityPtrKey
{
	Vector3* pos;
	Vector3* vel;
};

struct CollisionPoint
{
	Vector3 pos;
	bool collision = false;
}worldCollisionPoint[8][8][8];

const Vector3 surfaceNormals[6] =
{
	{ -1, 0, 0 },//x+
	{ 1, 0, 0 }, //x-
	{ 0, -1, 0 },//y+
	{ 0, 1, 0 }, //y-
	{ 0, 0, -1 },//z+
	{ 0, 0, 1 }, //z-
};

bool checkRange(double i, double j, double k, double range)
{
	return (i <= range) && (i >= -range) &&
		(j <= range) && (j >= -range) &&
		(k <= range) && (k >= -range);
}

void clamp(int* i, int range)
{
	*i = *i < 0 ? *i = 0 : *i;
	*i = *i > range - 1 ? range - 1 : *i;
}

void setPosVelPtrKey(
	struct PosVelocityPtrKey* pSurround,
	struct Vector3 p[8][8][8],
	struct Vector3 v[8][8][8],
	int i, int j, int k,
	int di, int dj, int dk
)
{
	int it = i + di, jt = j + dj, kt = k + dk;
	if ((it >= 0 && it <= 7) &&
		(jt >= 0 && jt <= 7) &&
		(kt >= 0 && kt <= 7))
	{
		pSurround->pos = &p[it][jt][kt];
		pSurround->vel = &v[it][jt][kt];
	}
	else
	{
		pSurround->pos = pSurround->vel = nullptr;
	}
}

void copySurroundPointPtrs(
	SpringType type,
	struct PosVelocityPtrKey surround[],
	struct Vector3 p[8][8][8],
	struct Vector3 v[8][8][8],
	int i, int j, int k)
{
	int d = 1 + (type == SPRING_BEND);
	if (type == SPRING_STRUCTURAL || type == SPRING_BEND)
	{
		setPosVelPtrKey(&surround[0], p, v, i, j, k, -d, 0, 0);
		setPosVelPtrKey(&surround[1], p, v, i, j, k, d, 0, 0);
		setPosVelPtrKey(&surround[2], p, v, i, j, k, 0, -d, 0);
		setPosVelPtrKey(&surround[3], p, v, i, j, k, 0, d, 0);
		setPosVelPtrKey(&surround[4], p, v, i, j, k, 0, 0, -d);
		setPosVelPtrKey(&surround[5], p, v, i, j, k, 0, 0, d);
	}
	else if (type == SPRING_SHEER_SHORT)
	{
		//Short sheer springs
		setPosVelPtrKey(&surround[0], p, v, i, j, k, -d, -d, 0);
		setPosVelPtrKey(&surround[1], p, v, i, j, k, d, d, 0);
		setPosVelPtrKey(&surround[2], p, v, i, j, k, -d, d, 0);
		setPosVelPtrKey(&surround[3], p, v, i, j, k, d, -d, 0);
		setPosVelPtrKey(&surround[4], p, v, i, j, k, 0, -d, -d);
		setPosVelPtrKey(&surround[5], p, v, i, j, k, 0, d, d);
		setPosVelPtrKey(&surround[6], p, v, i, j, k, 0, d, -d);
		setPosVelPtrKey(&surround[7], p, v, i, j, k, 0, -d, d);
		setPosVelPtrKey(&surround[8], p, v, i, j, k, -d, 0, -d);
		setPosVelPtrKey(&surround[9], p, v, i, j, k, d, 0, d);
		setPosVelPtrKey(&surround[10], p, v, i, j, k, -d, 0, d);
		setPosVelPtrKey(&surround[11], p, v, i, j, k, d, 0, -d);
	}
	else if (type == SPRING_SHEER_LONG)
	{
		//Long sheer springs
		setPosVelPtrKey(&surround[0], p, v, i, j, k, -d, -d, -d);
		setPosVelPtrKey(&surround[1], p, v, i, j, k, d, d, d);
		setPosVelPtrKey(&surround[2], p, v, i, j, k, d, -d, -d);
		setPosVelPtrKey(&surround[3], p, v, i, j, k, -d, d, d);
		setPosVelPtrKey(&surround[4], p, v, i, j, k, d, d, -d);
		setPosVelPtrKey(&surround[5], p, v, i, j, k, -d, -d, d);
		setPosVelPtrKey(&surround[6], p, v, i, j, k, -d, d, -d);
		setPosVelPtrKey(&surround[7], p, v, i, j, k, d, -d, d);
	}
	else
	{
		//Exception
	}
}

int checkCollision(struct Vector3* p, int a, int b, int c, int d)
{
	double x = p->x, y = p->y, z = p->z;
	int side = a * x + b * y + c * z + d >= 0;

	int collisionBit = COLLISION_NONE;

	if (x > boundingBoxRange) collisionBit |= COLLISION_X_POS_BOUND;
	if (x < -boundingBoxRange) collisionBit |= COLLISION_X_NAG_BOUND;

	if (y > boundingBoxRange) collisionBit |= COLLISION_Y_POS_BOUND;
	if (y < -boundingBoxRange) collisionBit |= COLLISION_Y_NAG_BOUND;

	if (z > boundingBoxRange) collisionBit |= COLLISION_Z_POS_BOUND;
	if (z < -boundingBoxRange) collisionBit |= COLLISION_Z_NAG_BOUND;

	if (side != originalSide) collisionBit |= COLLISION_PLANE_BOUND;

	return collisionBit;
}

Vector3 calPlasticForce(
	SpringType type,
	struct world* jello,
	struct Vector3* A, struct Vector3* B,
	struct Vector3* Va, struct Vector3* Vb,
	int forceFieldIndex
)
{
	Vector3 L = { 0, 0, 0 };
	Vector3 LN = { 0, 0, 0 };
	double length = 0;
	Vector3 hookForce = { 0 ,0 ,0 };

	Vector3 V = { 0, 0, 0 };
	Vector3 dampingForce = { 0, 0, 0 };

	Vector3 force = { 0, 0, 0 };

	if (B && Vb)
	{
		//Calculate hook force
		pDIFFERENCE(*A, *B, L);
		LENGTH(L, length);
		LN = L;
		pNORMALIZE(LN);

		switch (type)
		{
		case SPRING_STRUCTURAL:
			pMULTIPLY(LN, -jello->kElastic * (length - STRUCTURAL_SPRING_LENGTH), hookForce);
			break;
		case SPRING_SHEER_SHORT:
			pMULTIPLY(LN, -jello->kElastic * (length - SHEER_SPRING_LENGTH_SHORT), hookForce);
			break;
		case SPRING_SHEER_LONG:
			pMULTIPLY(LN, -jello->kElastic * (length - SHEER_SPRING_LENGTH_LONG), hookForce);
			break;
		case SPRING_BEND:
			pMULTIPLY(LN, -jello->kElastic * (length - BEND_SPRING_LENGTH), hookForce);
			break;
		case SPRING_COLLISION:
			pMULTIPLY(LN, -jello->kCollision * length, hookForce);
			break;
		}

		//Calculate damping force
		pDIFFERENCE(*Va, *Vb, V);
		double dot = 0.0;
		DOTPRODUCT(V, LN, dot);

		if (type != SPRING_COLLISION)
		{
			pMULTIPLY(LN, -jello->dElastic * dot, dampingForce);
		}
		else
		{
			pMULTIPLY(LN, -jello->dCollision * dot, dampingForce);
		}

		//Combine Hook, Damping and Force field into one force
		pSUM(hookForce, dampingForce, force);
		if (forceFieldIndex >= 0)
		{
			pSUM(force, jello->forceField[forceFieldIndex], force);
		}
	}
	else
	{
		force = { 0, 0, 0 };
	}
	return force;
}

void calCollisionPoint(const Vector3* P, const Vector3* V, const Vector3* N, double D, Vector3* originalCollisionPoint)
{
	Vector3 collisionVector = { 0, 0, 0 };
	Vector3 collisionPoint = { 0, 0, 0 };
	Vector3 VN = *V;
	double length = 0;
	pNORMALIZE(VN);
	pMULTIPLY(VN, -1, VN);
	double VoN = 0.0;

	DOTPRODUCT(VN, (*N), VoN);
	D /= VoN;
	pMULTIPLY(VN, D, collisionVector);
	pSUM((*P), collisionVector, collisionPoint);

	if (checkRange(collisionPoint.x, collisionPoint.y, collisionPoint.z, boundingBoxRange))
	{
		memcpy(originalCollisionPoint, &collisionPoint, sizeof(Vector3));
	}
}


Vector3 calSingleMassCollisionForce(struct world* jello, int i, int j, int k)
{
	int collisionBit = checkCollision(&(jello->p[i][j][k]), jello->a, jello->b, jello->c, jello->d);

	Vector3 collisionForce = { 0, 0, 0 };

	//Vector3 collisionPoint = { 0, 0, 0 };
	//Vector3 collisionVector = { 0, 0, 0 };
	Vector3 tempForce = { 0, 0, 0 };
	Vector3* V = &jello->v[i][j][k];
	Vector3* P = &jello->p[i][j][k];
	double D = 0.0;

	if (collisionBit | COLLISION_NONE)
	{
		//Vector3 tempForce = { 0, 0, 0 };
		if (collisionBit & COLLISION_X_POS_BOUND)
		{
			D = jello->p[i][j][k].x - boundingBoxRange;
			if (!worldCollisionPoint[i][j][k].collision)
			{
				//calCollisionPoint(P, V, &surfaceNormals[0], D, &collisionPoint);
				calCollisionPoint(P, V, &surfaceNormals[0], D, &worldCollisionPoint[i][j][k].pos);
				worldCollisionPoint[i][j][k].collision = true;
			}
			/*D = boundingBoxRange - jello->p[i][j][k].x;
			DOTPRODUCT(VN, surfaceNormals[0], VoN);
			D /= VoN;
			pMULTIPLY(VN, D, collisionVector);
			pSUM(jello->p[i][j][k], collisionVector, collisionPoint);*/


			tempForce = { boundingBoxRange - jello->p[i][j][k].x, 0, 0 };
			pMULTIPLY(tempForce, jello->kCollision, tempForce);
			/*pSUM(collisionForce, tempForce, collisionForce);*/
		}
		if (collisionBit & COLLISION_X_NAG_BOUND)
		{
			D = -jello->p[i][j][k].x - boundingBoxRange;
			if (!worldCollisionPoint[i][j][k].collision)
			{
				calCollisionPoint(P, V, &surfaceNormals[1], D, &worldCollisionPoint[i][j][k].pos);
				worldCollisionPoint[i][j][k].collision = true;
			}
			/*DOTPRODUCT(VN, surfaceNormals[1], VoN);
			D /= VoN;
			pMULTIPLY(VN, D, collisionVector);
			pSUM(jello->p[i][j][k], collisionVector, collisionPoint);*/
			tempForce = { -jello->p[i][j][k].x - boundingBoxRange, 0, 0 };
			pMULTIPLY(tempForce, jello->kCollision, tempForce);
			/*pSUM(collisionForce, tempForce, collisionForce);*/
		}
		if (collisionBit & COLLISION_Y_POS_BOUND)
		{
			D = jello->p[i][j][k].y - boundingBoxRange;
			if (!worldCollisionPoint[i][j][k].collision)
			{
				calCollisionPoint(P, V, &surfaceNormals[2], D, &worldCollisionPoint[i][j][k].pos);
				worldCollisionPoint[i][j][k].collision = true;
			}
			/*DOTPRODUCT(VN, surfaceNormals[2], VoN);
			D /= VoN;
			pMULTIPLY(VN, D, collisionVector);
			pSUM(jello->p[i][j][k], collisionVector, collisionPoint);*/
			tempForce = { 0, boundingBoxRange - jello->p[i][j][k].y, 0 };
			pMULTIPLY(tempForce, jello->kCollision, tempForce);
			/*pSUM(collisionForce, tempForce, collisionForce);*/
		}
		if (collisionBit & COLLISION_Y_NAG_BOUND)
		{
			D = -jello->p[i][j][k].y - boundingBoxRange;
			if (!worldCollisionPoint[i][j][k].collision)
			{
				calCollisionPoint(P, V, &surfaceNormals[3], D, &worldCollisionPoint[i][j][k].pos);
				worldCollisionPoint[i][j][k].collision = true;
			}
			/*DOTPRODUCT(VN, surfaceNormals[3], VoN);
			D /= VoN;
			pMULTIPLY(VN, D, collisionVector);
			pSUM(jello->p[i][j][k], collisionVector, collisionPoint);*/
			tempForce = { 0, -jello->p[i][j][k].y - boundingBoxRange, 0 };
			pMULTIPLY(tempForce, jello->kCollision, tempForce);
			/*pSUM(collisionForce, tempForce, collisionForce);*/
		}
		if (collisionBit & COLLISION_Z_POS_BOUND)
		{
			D = jello->p[i][j][k].z - boundingBoxRange;

			if (!worldCollisionPoint[i][j][k].collision)
			{
				calCollisionPoint(P, V, &surfaceNormals[4], D, &worldCollisionPoint[i][j][k].pos);
				worldCollisionPoint[i][j][k].collision = true;
			}
			/*DOTPRODUCT(VN, surfaceNormals[4], VoN);
			D /= VoN;
			pMULTIPLY(VN, D, collisionVector);
			pSUM(jello->p[i][j][k], collisionVector, collisionPoint);*/
			tempForce = { 0, 0, boundingBoxRange - jello->p[i][j][k].z };
			pMULTIPLY(tempForce, jello->kCollision, tempForce);
			/*pSUM(collisionForce, tempForce, collisionForce);*/
		}
		if (collisionBit & COLLISION_Z_NAG_BOUND)
		{
			D = -jello->p[i][j][k].z - boundingBoxRange;
			if (!worldCollisionPoint[i][j][k].collision)
			{
				calCollisionPoint(P, V, &surfaceNormals[5], D, &worldCollisionPoint[i][j][k].pos);
				worldCollisionPoint[i][j][k].collision = true;
			}
			/*DOTPRODUCT(VN, surfaceNormals[5], VoN);
			D /= VoN;
			pMULTIPLY(VN, D, collisionVector);
			pSUM(jello->p[i][j][k], collisionVector, collisionPoint);*/
			tempForce = { 0, 0, -jello->p[i][j][k].z - boundingBoxRange };
			pMULTIPLY(tempForce, jello->kCollision, tempForce);
			/*pSUM(collisionForce, tempForce, collisionForce);*/
		}
		if (collisionBit & COLLISION_PLANE_BOUND)
		{
			D = jello->a * jello->p[i][j][k].x + jello->b * jello->p[i][j][k].y + jello->c * jello->p[i][j][k].z + jello->d;
			D *= preComputeInclinePara;
			if (D < 0) D = -D;

			const Vector3 N = { jello->a, jello->b, jello->c };
			if (!worldCollisionPoint[i][j][k].collision)
			{
				calCollisionPoint(P, V, &N, D, &worldCollisionPoint[i][j][k].pos);
				worldCollisionPoint[i][j][k].collision = true;
			}
			/*DOTPRODUCT(VN, N, VoN);
			D /= VoN;
			pMULTIPLY(VN, D, collisionVector);
			pSUM(jello->p[i][j][k], collisionVector, collisionPoint);*/

			/*tempForce = { jello->a, jello->b, jello->c };
			pMULTIPLY(tempForce, D, tempForce);
			pMULTIPLY(tempForce, jello->kCollision, tempForce);

			DOTPRODUCT(jello->v[i][j][k], tempForce, VoN);
			if (VoN > 0) pMULTIPLY(tempForce, -1, tempForce);

			/*pSUM(collisionForce, tempForce, collisionForce);*/
		}

		Vector3 Vb = { 0, 0, 0 };
		collisionForce = calPlasticForce(SPRING_COLLISION, jello, &jello->p[i][j][k], &worldCollisionPoint[i][j][k].pos, &jello->v[i][j][k], &Vb, -1);
	}
	else
	{
		worldCollisionPoint[i][j][k].collision = false;
	}

	return collisionForce;
}

/* Calculate the field force by interpolating */
Vector3 calForceFieldForce(const Vector3& p)
{
	Vector3 forceIndex = { 0, 0, 0 };
	forceIndex.x = (p.x + boundingBoxRange) * invForceFieldGridLength;
	forceIndex.y = (p.y + boundingBoxRange) * invForceFieldGridLength;
	forceIndex.z = (p.z + boundingBoxRange) * invForceFieldGridLength;

	return forceIndex;
}

Vector3 calSingleMassForce(
	SpringType type,
	struct world* jello,
	struct PosVelocityPtrKey surround[],
	const int surroundCount,
	int i, int j, int k
)
{
	Vector3 forces[12];

	for (int m = 0; m < surroundCount; ++m)
	{
		Vector3 forceFieldIndex = calForceFieldForce(jello->p[i][j][k]);
		int x = (int)forceFieldIndex.x;
		int y = (int)forceFieldIndex.y;
		int z = (int)forceFieldIndex.z;
		clamp(&x, jello->resolution);
		clamp(&y, jello->resolution);
		clamp(&z, jello->resolution);

		forces[m] = calPlasticForce(
			type, jello,
			&(jello->p[i][j][k]), surround[m].pos,
			&(jello->v[i][j][k]), surround[m].vel,
			FORCEFIELD_INDEX(x, y, z)
		);
	}

	//Forces Integration
	Vector3 force = { 0, 0, 0 };
	for (int i = 0; i < surroundCount; ++i)
	{
		pSUM(force, forces[i], force);
	}

	return force;
}

/* Computes acceleration to every control point of the jello cube,
   which is in state given by 'jello'.
   Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct Vector3 a[8][8][8])
{
	//Calculate drag force by mouse
	Vector3 mouseDragForce = { 0, dragForce[0], -dragForce[1] };

	/* for you to implement ... */
	for (int i = 0; i < 8; ++i)
	{
		for (int j = 0; j < 8; ++j)
		{
			for (int k = 0; k < 8; ++k)
			{
				a[i][j][k].x = a[i][j][k].y = a[i][j][k].z = 0;//init

				struct Vector3 F = { 0, 0, 0 };
				struct Vector3 fSturctural = { 0, 0, 0 };
				struct Vector3 fSheerShort = { 0, 0, 0 };
				struct Vector3 fSheerLong = { 0, 0, 0 };
				struct Vector3 fBend = { 0, 0, 0 };

				struct PosVelocityPtrKey structuralSurround[6];
				struct PosVelocityPtrKey sheerShortSurround[12];
				struct PosVelocityPtrKey sheerLongSurround[8];
				struct PosVelocityPtrKey bendSurround[6];

				copySurroundPointPtrs(SPRING_STRUCTURAL, structuralSurround, jello->p, jello->v, i, j, k);
				fSturctural = calSingleMassForce(SPRING_STRUCTURAL, jello, structuralSurround, 6, i, j, k);

				copySurroundPointPtrs(SPRING_SHEER_SHORT, sheerShortSurround, jello->p, jello->v, i, j, k);
				fSheerShort = calSingleMassForce(SPRING_SHEER_SHORT, jello, sheerShortSurround, 12, i, j, k);

				copySurroundPointPtrs(SPRING_SHEER_LONG, sheerLongSurround, jello->p, jello->v, i, j, k);
				fSheerLong = calSingleMassForce(SPRING_SHEER_LONG, jello, sheerLongSurround, 8, i, j, k);

				copySurroundPointPtrs(SPRING_BEND, bendSurround,
					jello->p, jello->v, i, j, k);
				fBend = calSingleMassForce(SPRING_BEND, jello, bendSurround, 6, i, j, k);

				//Calculate collision force
				Vector3 collisionForce = calSingleMassCollisionForce(jello, i, j, k);

				pSUM(fSturctural, fSheerShort, F);
				pSUM(F, fSheerLong, F);
				pSUM(F, fBend, F);
				pSUM(F, collisionForce, F);
				pSUM(F, mouseDragForce, F);
				P_DIVIDE(F, jello->mass, a[i][j][k]);
			}
		}
	}
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
	int i, j, k;
	Vector3 a[8][8][8];

	computeAcceleration(jello, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
				jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
				jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
				jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
				jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
				jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

			}
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
	Vector3 F1p[8][8][8], F1v[8][8][8],
		F2p[8][8][8], F2v[8][8][8],
		F3p[8][8][8], F3v[8][8][8],
		F4p[8][8][8], F4v[8][8][8];

	Vector3 a[8][8][8];


	struct world buffer;

	int i, j, k;

	buffer = *jello; // make a copy of jello

	computeAcceleration(jello, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				pMULTIPLY(jello->v[i][j][k], jello->dt, F1p[i][j][k]);
				pMULTIPLY(a[i][j][k], jello->dt, F1v[i][j][k]);
				pMULTIPLY(F1p[i][j][k], 0.5, buffer.p[i][j][k]);
				pMULTIPLY(F1v[i][j][k], 0.5, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	computeAcceleration(&buffer, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// F2p = dt * buffer.v;
				pMULTIPLY(buffer.v[i][j][k], jello->dt, F2p[i][j][k]);
				// F2v = dt * a(buffer.p,buffer.v);     
				pMULTIPLY(a[i][j][k], jello->dt, F2v[i][j][k]);
				pMULTIPLY(F2p[i][j][k], 0.5, buffer.p[i][j][k]);
				pMULTIPLY(F2v[i][j][k], 0.5, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	computeAcceleration(&buffer, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// F3p = dt * buffer.v;
				pMULTIPLY(buffer.v[i][j][k], jello->dt, F3p[i][j][k]);
				// F3v = dt * a(buffer.p,buffer.v);     
				pMULTIPLY(a[i][j][k], jello->dt, F3v[i][j][k]);
				pMULTIPLY(F3p[i][j][k], 0.5, buffer.p[i][j][k]);
				pMULTIPLY(F3v[i][j][k], 0.5, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	computeAcceleration(&buffer, a);


	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// F3p = dt * buffer.v;
				pMULTIPLY(buffer.v[i][j][k], jello->dt, F4p[i][j][k]);
				// F3v = dt * a(buffer.p,buffer.v);     
				pMULTIPLY(a[i][j][k], jello->dt, F4v[i][j][k]);

				pMULTIPLY(F2p[i][j][k], 2, buffer.p[i][j][k]);
				pMULTIPLY(F3p[i][j][k], 2, buffer.v[i][j][k]);
				pSUM(buffer.p[i][j][k], buffer.v[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F1p[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F4p[i][j][k], buffer.p[i][j][k]);
				pMULTIPLY(buffer.p[i][j][k], 1.0 / 6, buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], jello->p[i][j][k], jello->p[i][j][k]);

				pMULTIPLY(F2v[i][j][k], 2, buffer.p[i][j][k]);
				pMULTIPLY(F3v[i][j][k], 2, buffer.v[i][j][k]);
				pSUM(buffer.p[i][j][k], buffer.v[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F1v[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F4v[i][j][k], buffer.p[i][j][k]);
				pMULTIPLY(buffer.p[i][j][k], 1.0 / 6, buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], jello->v[i][j][k], jello->v[i][j][k]);
			}

	return;
}
